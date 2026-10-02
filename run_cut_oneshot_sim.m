% "K 를 정하고 고유벡터 K 개로 한 번에 자르기" 비교 + 계산 시간 (합성 타깃 12종 x 점군 4개 = 48개 전체)
% 가중치: C8a 와 같은 법선 가중(sigma 0.2) + alpha=1
%
% K 정하기
%   K1   : 전체 스펙트럼 최대 비율 간격 (영고유값 2개 이상이면 K = 연결성분 수)
%   K1s  : 거의 0 인 고유값은 각각 별도 군집으로 확정, 그 뒤 구조적 고유값에서만 최대 비율.
%          최대 비율이 2 미만이면 더 안 나눔 (K = 1 도 가능)      selectKRatio 'skipnear'
%   K1sc : 진짜 연결성분(끊어진 덩어리)마다 따로 라플라시안을 풀고 K1s 를 적용, 자르기도 성분별로.
%          끊어진 덩어리들의 고유값이 한 줄로 섞이면 비율 간격이 흐려지는 문제를 피한다
%   orc  : 정답 부품 수 (참고: K 를 맞혔을 때의 상한)
%
% 고유벡터 K 개로 한 번에 자르기
%   pcca  : PCCA+ — simplex 꼭짓점 K 개로 회전한 뒤 소속도 argmax       (cutPCCA)
%   km    : k-means, 20 회 반복 시작 (NJW: L_sym 고유벡터 행 정규화)
%   ncut  : 임베딩 MST 를 후보로, 원 그래프 NCut 증가가 최소인 간선부터   (cutTreeGreedy 'ncut')
%   len   : 임베딩 MST 최장 간선부터                                     (cutTreeGreedy 'length')
%   sign  : 회전 없이 고유벡터 2..K 의 부호 패턴으로 바로 나누고 큰 K 개만 남김
%           (0 근처 점은 부호가 흔들려 애매해짐 — 회전(PCCA+)이 왜 필요한지 보기 위한 비교)
%   기준  : 재귀 이분할 + 병목 판정 (cutRecursive, K 불필요)
%
% 시간: 그래프 구성, 전역 고유분해(20개), 방식별 자르기 시간을 따로 잰다.
%       K1sc 는 성분별 고유분해 시간을 자르기 시간에 포함한다.
% 출력: out_cutting_simcut/<yyyyMMdd_HHmmss>_oneshot/  summary.csv, results.mat, labels.mat, summary.png
rng(1);

%% 0) 설정
tuneSet = {'cubesat6u', 'cubesat6u_mid', 'hexsat'};
if ~exist('targets','var') || isempty(targets)
    targets = {'cubesat6u','cubesat6u_mid','hexsat','cubesat3u_petal','spin_drum','geo_comsat', ...
               'eo_sat','telescope','soyuz','capsule','station_t','upper_stage'};
end
clouds  = {'complete', 's1_mrbar_h3', 's1_oblique_h3', 's1_pvbar_h3'};
kEig    = 20; tau = 0.25;
kRules  = {'K1','K1s','K1sc','orc'};
cutters = {'pcca','km','ncut','len','sign'};
saveLab = {'K1sc_pcca','K1sc_ncut','K1s_pcca','orc_ncut','rec'};          % 그림용으로 라벨 저장

doMerge = exist('mergeDirs','var') && ~isempty(mergeDirs);   % 나눠 돌린 결과 폴더들을 합쳐 요약만
outDir = fullfile('out_cutting_simcut', [char(datetime('now','Format','yyyyMMdd_HHmmss')) '_oneshot']);
if ~exist(outDir, 'dir'), mkdir(outDir); end
fprintf('출력: %s\n', outDir);

%% 1) 실행
rows = struct([]); figdata = struct([]);
if doMerge, targets = {}; end
for it = 1:numel(targets)
    d = dir(fullfile('out_cutting_sim', ['*_' targets{it}]));
    if isempty(d), continue; end
    runDir = fullfile(d(end).folder, d(end).name);
    for ic = 1:numel(clouds)
        C = loadSimCloud(fullfile(runDir, [clouds{ic} '.mat']));
        t = tic;
        [W, gi] = buildGraphConvex(C.P, C.sensorPos(C.view_id,:), struct('mode','unsigned','minCompSize',50));
        tGraph = toc(t);
        P = C.P(gi.idxKeep,:); gt = C.label(gi.idxKeep);
        N = size(W,1); ms = max(60, round(0.01*N));
        comp = conncomp(graph(W, 'upper')).'; n0 = max(comp);
        Kgt = numel(unique(gt));
        q  = full(sum(W,2)); Dq = spdiags(1./q, 0, N, N);
        Wt = Dq * W * Dq; Wt = (Wt + Wt.')/2;

        t = tic;
        [Ln, Dis] = normalizeGraph(W, 1);
        kE = min(kEig, N-2);
        [U, lam] = embedSpectral(Ln, kE);
        Psi = Dis * U;
        tEig = toc(t);

        Ks = struct('K1',  selectKRatio(lam, struct('rule','max')), ...
                    'K1s', selectKRatio(lam, struct('rule','skipnear')), 'orc', Kgt);
        ctx = struct('target', targets{it}, 'cloud', clouds{ic}, 'isTune', double(any(strcmp(targets{it}, tuneSet))), ...
                     'N', N, 'n0', n0, 'Kgt', Kgt, 'tGraph', tGraph, 'tEig', tEig);
        fd = struct('target', targets{it}, 'cloud', clouds{ic}, 'P', single(P), 'gt', uint8(gt), ...
                    'Kgt', Kgt, 'K1', Ks.K1, 'K1s', Ks.K1s, 'labels', struct());

        for ik = 1:numel(kRules)
            kr = kRules{ik};
            for icut = 1:numel(cutters)
                name = [kr '_' cutters{icut}];
                t = tic;
                try
                    if strcmp(kr, 'K1sc')
                        [lab, Kuse] = cutPerComponent(W, Wt, comp, cutters{icut}, kE, ms);
                    else
                        Kuse = min(Ks.(kr), kE-1);
                        lab  = cutOneShot(cutters{icut}, Psi, U, Wt, Kuse, ms);
                    end
                catch err
                    warning('%s 실패: %s', name, err.message); lab = ones(N,1); Kuse = NaN;
                end
                tm = toc(t);
                [rows, lab] = addRow(rows, ctx, name, kr, cutters{icut}, Kuse, tm, lab, gt, Wt, lam);
                if any(strcmp(name, saveLab)), fd.labels.(name) = uint8(lab); end
            end
        end
        t = tic; lab = cutRecursive(W, struct('tau', tau, 'minSize', ms)); tm = toc(t);
        [rows, lab] = addRow(rows, ctx, 'rec', '-', 'rec', NaN, tm, lab, gt, Wt, lam);
        fd.labels.rec = uint8(lab);
        figdata = [figdata; fd]; %#ok<AGROW>

        h = rows(strcmp({rows.target}, targets{it}) & strcmp({rows.cloud}, clouds{ic}));
        g = @(m) h(strcmp({h.method}, m));
        fprintf('%-16s %-14s N=%5d 성분=%d Kgt=%2d | K1 %2d K1s %2d K1sc %2d | pcca %.2f/%.2f/%.2f  ncut %.2f/%.2f/%.2f  rec %.2f | 고유분해 %.1fs\n', ...
            targets{it}, clouds{ic}, N, n0, Kgt, Ks.K1, Ks.K1s, g('K1sc_pcca').Kused, ...
            g('K1_pcca').ARI, g('K1s_pcca').ARI, g('K1sc_pcca').ARI, g('K1_ncut').ARI, g('K1s_ncut').ARI, g('K1sc_ncut').ARI, ...
            g('rec').ARI, tEig);
    end
    writetable(struct2table(rows), fullfile(outDir, 'summary.csv'));
end
if doMerge
    T = table();
    for m = 1:numel(mergeDirs)
        R = load(fullfile(mergeDirs{m}, 'results.mat')); Fm = load(fullfile(mergeDirs{m}, 'labels.mat'));
        T = [T; R.T]; figdata = [figdata; Fm.figdata]; %#ok<AGROW>
    end
    targets = unique(T.target, 'stable').';
else
    T = struct2table(rows);
end
writetable(T, fullfile(outDir, 'summary.csv'));
shown = saveLab;
save(fullfile(outDir, 'results.mat'), 'T', 'tuneSet', 'targets', 'clouds', 'kRules', 'cutters');
save(fullfile(outDir, 'labels.mat'), 'figdata', 'shown');

%% 2) 요약 (48개 전체)
fprintf('\n[K 정확도] 48개 점군 (K1sc 는 성분별 K 의 합)\n');
for kr = {'K1','K1s','K1sc'}
    q = strcmp(T.method, [kr{1} '_pcca']);
    Ku = T.Kused(q); Kg = T.Kgt(q);
    fprintf('  %-5s 정확 %2d/48  과대 평균 %.2f (최대 %d)  과소 평균 %.2f (최대 %d)\n', kr{1}, nnz(Ku==Kg), ...
        mean(max(Ku-Kg,0)), max(Ku-Kg), mean(max(Kg-Ku,0)), max(Kg-Ku));
end
mlist = unique(T.method, 'stable');
fprintf('\n%-11s | 48개: ARI   과분할 미분할 조각수 | 튜닝 ARI 테스트 ARI | 자르기 시간 중앙값 [s] (최대)\n', '방식');
Sm = zeros(numel(mlist), 6);
for im = 1:numel(mlist)
    q = strcmp(T.method, mlist{im});
    Sm(im,:) = [mean(T.ARI(q)), mean(max(T.K(q)-T.Kgt(q),0)), mean(max(T.Kgt(q)-T.K(q),0)), mean(T.frag(q)), ...
                median(T.timeSec(q)), max(T.timeSec(q))];
    fprintf('%-11s |      %.3f  %5.2f  %5.2f  %5.2f |   %.3f    %.3f   |   %.3f  (%.2f)\n', mlist{im}, Sm(im,1:4), ...
        mean(T.ARI(q & T.isTune==1)), mean(T.ARI(q & T.isTune==0)), Sm(im,5:6));
end
u = unique(T(:, {'target','cloud','N','tGraph','tEig'}));
fprintf('\n공통 단계 시간 중앙값: 그래프 구성 %.2fs, 전역 고유분해(20개) %.2fs  (N 중앙값 %d)\n', ...
    median(u.tGraph), median(u.tEig), round(median(u.N)));

f = figure('Visible','off', 'Position',[30 30 1900 1100]);
tl = tiledlayout(f, 2, 2, 'TileSpacing','compact', 'Padding','compact');
title(tl, '48개 점군 전체: K 정하기 x 한 번에 자르기', 'FontSize', 13);
A = nan(numel(kRules), numel(cutters)); O = A; Ud = A; Tm = A;
for ik = 1:numel(kRules)
    for icut = 1:numel(cutters)
        q = strcmp(T.method, [kRules{ik} '_' cutters{icut}]);
        A(ik,icut) = mean(T.ARI(q)); O(ik,icut) = mean(max(T.K(q)-T.Kgt(q),0));
        Ud(ik,icut) = mean(max(T.Kgt(q)-T.K(q),0)); Tm(ik,icut) = median(T.timeSec(q));
    end
end
cl = {'PCCA+','k-means','MST+NCut','MST 최장간선','부호 패턴'};
qr = strcmp(T.method, 'rec');
ax = nexttile(tl); bar(ax, A.'); ylim(ax, [0 1]); grid(ax,'on'); ylabel(ax, '평균 ARI');
yline(ax, mean(T.ARI(qr)), '--k', sprintf('재귀+병목 %.2f', mean(T.ARI(qr))));
set(ax, 'XTickLabel', cl); legend(ax, kRules, 'Location','northeastoutside'); title(ax, '정확도');
ax = nexttile(tl); bar(ax, O.'); grid(ax,'on'); ylabel(ax, '평균 과분할 (K 초과)');
yline(ax, mean(max(T.K(qr)-T.Kgt(qr),0)), '--k', '재귀+병목');
set(ax, 'XTickLabel', cl); legend(ax, kRules, 'Location','northeastoutside'); title(ax, '과분할');
ax = nexttile(tl); bar(ax, Ud.'); grid(ax,'on'); ylabel(ax, '평균 미분할 (K 부족)');
yline(ax, mean(max(T.Kgt(qr)-T.K(qr),0)), '--k', '재귀+병목');
set(ax, 'XTickLabel', cl); legend(ax, kRules, 'Location','northeastoutside'); title(ax, '미분할');
ax = nexttile(tl); bar(ax, Tm.'); set(ax, 'YScale', 'log'); grid(ax,'on'); ylabel(ax, '자르기 시간 중앙값 [s] (로그)');
yline(ax, median(T.timeSec(qr)), '--k', '재귀+병목'); yline(ax, median(u.tEig), ':r', '전역 고유분해');
set(ax, 'XTickLabel', cl); legend(ax, kRules, 'Location','northeastoutside'); title(ax, '시간 (그래프·전역 고유분해 제외)');
exportgraphics(f, fullfile(outDir, 'summary.png'), 'Resolution', 110); close(f);
fprintf('\n저장 완료: %s\n', outDir);

%% ===================== local functions =====================
function lab = cutOneShot(name, Psi, U, Wt, K, ms)
% 고유벡터 K 개로 한 번에 K 조각
N = size(Psi,1);
if K <= 1, lab = ones(N,1); return; end
switch name
    case 'pcca', lab = cutPCCA(Psi(:,1:K));
    case 'km'
        X = U(:,1:K); X = X ./ max(vecnorm(X,2,2), eps);
        lab = kmeans(X, K, 'Replicates', 20, 'MaxIter', 300);
    case 'ncut', lab = cutTreeGreedy(Psi(:,1:K), Wt, K, struct('score','ncut'));
    case 'len',  lab = cutTreeGreedy(Psi(:,1:K), Wt, K, struct('score','length','minSize',ms));
    case 'sign', lab = cutSign(Psi(:,2:K), K, Wt);
end
end

function lab = cutSign(Y, K, Wt)
% 부호 패턴으로 바로 나눔: 2^(K-1) 칸 중 점이 많은 K 칸만 남기고 나머지는 이웃 칸에 붙임
[~, ~, code] = unique(Y > 0, 'rows');
sz = accumarray(code, 1);
[~, o] = sort(sz, 'descend');
keepC = o(1:min(K, numel(o)));
lab = zeros(size(code));
for k = 1:numel(keepC), lab(code == keepC(k)) = k; end
for it = 1:100                                          % 나머지를 가장 강하게 이어진 칸으로
    un = find(lab == 0);
    if isempty(un), break; end
    [ii, jj, ww] = find(Wt(un,:)); ii = ii(:); jj = jj(:); ww = ww(:);
    ok = lab(jj) > 0; if ~any(ok), break; end
    A = accumarray([ii(ok), lab(jj(ok))], ww(ok), [numel(un) numel(keepC)]);
    [mx, b] = max(A, [], 2); if ~any(mx > 0), break; end
    lab(un(mx > 0)) = b(mx > 0);
end
lab(lab == 0) = numel(keepC) + 1;
end

function [lab, Ksum] = cutPerComponent(W, Wt, comp, name, kE, ms)
% 끊어진 덩어리마다 따로: 그 덩어리의 라플라시안 → K1s → 한 번에 자르기
N = size(W,1); lab = zeros(N,1); off = 0; Ksum = 0;
for c = 1:max(comp)
    idx = find(comp == c); n = numel(idx);
    if n < 3*ms
        lc = ones(n,1); Kc = 1;
    else
        [Lc, Dc] = normalizeGraph(W(idx,idx), 1);
        kc = min(kE, n-2);
        [Uc, lamc] = embedSpectral(Lc, kc);
        Kc = min(selectKRatio(lamc, struct('rule','skipnear')), kc-1);
        lc = cutOneShot(name, Dc*Uc, Uc, Wt(idx,idx), Kc, ms);
    end
    [~, ~, lc] = unique(lc);
    lab(idx) = off + lc; off = off + max(lc); Ksum = Ksum + Kc;
end
end

function [rows, lab] = addRow(rows, ctx, name, kRule, cutter, Kused, tm, lab, gt, Wt, lam)
e = evalCutResult(lab, gt, Wt, lam);
r = ctx;
r.method = name; r.kRule = kRule; r.cutter = cutter; r.Kused = Kused;
r.K = e.K; r.ARI = e.ARI; r.purity = e.purity; r.frag = e.frag; r.ncutK = e.ncutK;
r.upsilon = e.upsilon; r.timeSec = tm;
rows = [rows; r];
lab = e.labels;
end
