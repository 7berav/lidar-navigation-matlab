% 합성 타깃 전체에서 "가중치 x 자르는 방식" 비교
%
% 가중치 (둘 다 mutual kNN k=30 + self-tuning, alpha=1)
%   u : 법선이 다르면 무조건 약화 (sigma 0.2)                 = C8a 의 가중치
%   c : 오목 접합만 약화, 볼록 모서리는 유지 (buildGraphConvex) = 제안
%
% 자르는 방식
%   mst : 임베딩(15차원) MST 에서 phi 작은 컷부터 K-1 개      (segmentEmbedding 'phiorder')
%   km  : 앞쪽 고유벡터 K 개에서 k-means                        (표준 spectral clustering)
%   rec : 원 그래프 재귀 이분할, 병목(rho <= tau)일 때만 자름   (cutRecursive, K 불필요)
%
% K 정하기 (mst, km)
%   gap : 가장 큰 고유값 간격  K = argmax_k (lambda_{k+1} - lambda_k),  k >= 연결성분 수
%   orc : 정답 부품 수 (보이는 것)를 알려줌 — 방식 자체의 상한을 보기 위한 참고용
%
% tau = 0.25 는 튜닝 타깃(cubesat6u, cubesat6u_mid, hexsat)에서만 정했다 (경계 rho 0~0.09,
% 부품 내부 rho >= 0.33). 나머지 9종은 테스트.
%
% 지표: ARI / K / 순도(부품을 덜 합칠수록 ↑) / 조각수(부품의 90% 를 덮는 데 필요한 군집 수, 1 이 정답)
% 출력: out_cutting_simcut/<yyyyMMdd_HHmmss>_methods/
%   summary.csv (타깃 x 점군 x 방식), results.mat, overview_<target>.png, scores.png
rng(1);

%% 0) 설정
tuneSet = {'cubesat6u', 'cubesat6u_mid', 'hexsat'};
if ~exist('targets','var') || isempty(targets)
    targets = {'cubesat6u','cubesat6u_mid','hexsat','cubesat3u_petal','spin_drum','geo_comsat', ...
               'eo_sat','telescope','soyuz','capsule','station_t','upper_stage'};
end
clouds = {'complete', 's1_mrbar_h3', 's1_oblique_h3', 's1_pvbar_h3'};
tau    = 0.25;
kEig   = 20;
wModes = struct('tag', {'u','c'}, 'mode', {'unsigned','convex'});
shown  = {'u_mst_gap','u_km_gap','u_rec','c_mst_gap','c_km_gap','c_rec'};   % overview 그림에 넣을 방식

outDir = fullfile('out_cutting_simcut', [char(datetime('now','Format','yyyyMMdd_HHmmss')) '_methods']);
if ~exist(outDir, 'dir'), mkdir(outDir); end
segHash = fileMD5('segmentEmbedding.m');
fprintf('출력: %s\nsegmentEmbedding.m md5 = %s\n', outDir, segHash);

%% 1) 타깃 x 점군 x 방식
rows = struct([]);
for it = 1:numel(targets)
    d = dir(fullfile('out_cutting_sim', ['*_' targets{it}]));
    if isempty(d), fprintf('[건너뜀] %s 점군 없음\n', targets{it}); continue; end
    runDir = fullfile(d(end).folder, d(end).name);
    fig = figure('Visible','off', 'Position',[20 20 2600 1500]);
    tl  = tiledlayout(fig, numel(clouds), numel(shown)+1, 'TileSpacing','compact', 'Padding','compact');
    title(tl, sprintf('%s  (tau=%.2f)', targets{it}, tau), 'Interpreter','none');

    for ic = 1:numel(clouds)
        C = loadSimCloud(fullfile(runDir, [clouds{ic} '.mat']));
        S = C.sensorPos(C.view_id, :);
        L = struct();                                   % 방식 이름 → 라벨
        for iw = 1:numel(wModes)
            tg = wModes(iw).tag;
            [W, gi] = buildGraphConvex(C.P, S, struct('mode', wModes(iw).mode, 'minCompSize', 50));
            P  = C.P(gi.idxKeep, :); gt = C.label(gi.idxKeep);
            N  = size(W,1); ms = max(60, round(0.01*N));
            nComp = max(conncomp(graph(W, 'upper')));
            Kgt   = numel(unique(gt));

            [Ln, Dis] = normalizeGraph(W, 1);
            kE = min(kEig, N-2);
            [U, lam] = embedSpectral(Ln, kE);
            Y = Dis * U;
            Kgap = gapK(lam, nComp, min(15, kE-1));

            Ks = struct('gap', Kgap, 'orc', Kgt);
            for kn = {'gap','orc'}
                K = Ks.(kn{1});
                t = tic; lab = cutMST(Y, W, nComp, K, ms);      tm = toc(t);
                rows = addRow(rows, targets{it}, clouds{ic}, [tg '_mst_' kn{1}], lab, gt, N, nComp, Kgt, Kgap, tm, tuneSet);
                L.([tg '_mst_' kn{1}]) = lab;
                t = tic; lab = cutKmeans(Y, K);                 tm = toc(t);
                rows = addRow(rows, targets{it}, clouds{ic}, [tg '_km_' kn{1}], lab, gt, N, nComp, Kgt, Kgap, tm, tuneSet);
                L.([tg '_km_' kn{1}]) = lab;
            end
            t = tic; lab = cutRecursive(W, struct('tau', tau, 'minSize', ms)); tm = toc(t);
            rows = addRow(rows, targets{it}, clouds{ic}, [tg '_rec'], lab, gt, N, nComp, Kgt, Kgap, tm, tuneSet);
            L.([tg '_rec']) = lab;
            if iw == 1, Pu = P; gtu = gt; else, Pc = P; end
        end

        drawLab(nexttile(tl), Pu, gtu);
        title(gca, sprintf('%s: GT (K=%d)', clouds{ic}, numel(unique(gtu))), 'Interpreter','none');
        for is = 1:numel(shown)
            r = rows(strcmp({rows.target}, targets{it}) & strcmp({rows.cloud}, clouds{ic}) & strcmp({rows.method}, shown{is}));
            if shown{is}(1) == 'u', Pp = Pu; else, Pp = Pc; end
            drawLab(nexttile(tl), Pp, L.(shown{is}));
            title(gca, sprintf('%s: K=%d, ARI=%.2f', shown{is}, r.K, r.ARI), 'Interpreter','none');
        end
        rr = rows(strcmp({rows.target}, targets{it}) & strcmp({rows.cloud}, clouds{ic}));
        fprintf('%-16s %-14s Kgt=%2d |', targets{it}, clouds{ic}, rr(1).Kgt);
        for is = 1:numel(shown)
            r = rr(strcmp({rr.method}, shown{is}));
            fprintf(' %s %2d/%.2f', shown{is}, r.K, r.ARI);
        end
        fprintf('\n');
    end
    exportgraphics(fig, fullfile(outDir, ['overview_' targets{it} '.png']), 'Resolution', 90); close(fig);
    T = struct2table(rows); writetable(T, fullfile(outDir, 'summary.csv'));     % 중간 저장
end
T = struct2table(rows);
writetable(T, fullfile(outDir, 'summary.csv'));
save(fullfile(outDir, 'results.mat'), 'T', 'tau', 'tuneSet', 'targets', 'clouds', 'segHash');

%% 2) 방식별 요약 (튜닝 / 테스트 분리)
methods = unique(T.method, 'stable');
fprintf('\n%-12s | %-38s | %-38s\n', '방식', '튜닝 3종: ARI  과분할  미분할  조각수', '테스트 9종: ARI  과분할  미분할  조각수');
Sm = zeros(numel(methods), 8);
for im = 1:numel(methods)
    for s = 0:1
        q = strcmp(T.method, methods{im}) & T.isTune == (1-s);
        Sm(im, s*4 + (1:4)) = [mean(T.ARI(q)), mean(max(T.K(q)-T.Kgt(q),0)), mean(max(T.Kgt(q)-T.K(q),0)), mean(T.frag(q))];
    end
    fprintf('%-12s |     %.3f   %5.2f   %5.2f   %5.2f     |     %.3f   %5.2f   %5.2f   %5.2f\n', methods{im}, Sm(im,:));
end

f = figure('Visible','off', 'Position',[30 30 1700 520]);
tl = tiledlayout(f, 1, 3, 'TileSpacing','compact', 'Padding','compact');
ax = nexttile(tl); bar(ax, Sm(:,[1 5])); ylabel(ax,'평균 ARI'); legend(ax, {'튜닝','테스트'}, 'Location','northwest');
ax.XTick = 1:numel(methods); ax.XTickLabel = methods; ax.TickLabelInterpreter = 'none'; ax.XTickLabelRotation = 40; grid(ax,'on');
ax = nexttile(tl); bar(ax, Sm(:,[6 7]), 'stacked'); ylabel(ax,'테스트: 평균 |K - 정답|'); legend(ax, {'과분할 (K 초과)','미분할 (K 부족)'}, 'Location','northwest');
ax.XTick = 1:numel(methods); ax.XTickLabel = methods; ax.TickLabelInterpreter = 'none'; ax.XTickLabelRotation = 40; grid(ax,'on');
ax = nexttile(tl); bar(ax, Sm(:,8)); ylabel(ax,'테스트: 평균 조각수 (1 이 정답)'); yline(ax, 1, '--');
ax.XTick = 1:numel(methods); ax.XTickLabel = methods; ax.TickLabelInterpreter = 'none'; ax.XTickLabelRotation = 40; grid(ax,'on');
exportgraphics(f, fullfile(outDir, 'scores.png'), 'Resolution', 110); close(f);
fprintf('\n저장 완료: %s\n', outDir);

%% ===================== local functions =====================
function K = gapK(lam, nComp, kMax)
% 가장 큰 고유값 간격. k 는 연결성분 수(영고유값 개수) 이상
kMin = max(nComp, 1);
g = diff(lam(:));
[~, i] = max(g(kMin:kMax));
K = kMin + i - 1;
end

function lab = cutMST(Y, W, nComp, K, ms)
% 임베딩 MST 에서 phi 작은 컷부터 (K - 임베딩 성분 수)개. 실패하면 0 라벨
N = size(W,1);
dim  = min(15, size(Y,2) - nComp);
Yuse = Y(:, nComp+1 : nComp+dim);
ky   = min(40, N-2);
m    = embedComponents(Yuse, ky);
o = struct('ky', ky, 'minSize', ms, 'cutMethod', 'phiorder', 'phiMax', 1, ...
           'candidates', 'all', 'maxCuts', max(K - m, 0));
try
    lab = segmentEmbedding(Yuse, W, o);
    lab = attachZero(lab, W);
catch err
    warning('segmentEmbedding 실패: %s', err.message);
    lab = zeros(N, 1, 'uint32');
end
end

function lab = cutKmeans(Y, K)
N = size(Y,1);
if K <= 1, lab = ones(N, 1, 'uint32'); return; end
X = Y(:, 1:K);
X = X ./ max(vecnorm(X, 2, 2), eps);
lab = uint32(kmeans(X, K, 'Replicates', 5, 'MaxIter', 300));
end

function m = embedComponents(Y, ky)
N = size(Y,1);
[idxNN, distNN] = knnsearch(Y, Y, 'K', ky+1);
Sg = sparse(repmat((1:N)', ky, 1), reshape(idxNN(:,2:end),[],1), reshape(distNN(:,2:end),[],1), N, N);
m = max(conncomp(graph(max(Sg, Sg.'), 'upper')));
end

function lab = attachZero(lab, W)
% 라벨 0 을 W 로 가장 강하게 이어진 군집에 붙임 (모든 방식을 coverage 1 로 맞춰 비교)
lab = double(lab); K = max(lab);
if K == 0, lab = ones(size(lab), 'uint32'); return; end
for it = 1:100
    un = find(lab == 0);
    if isempty(un), break; end
    [ii, jj, ww] = find(W(un, :)); ii = ii(:); jj = jj(:); ww = ww(:);
    ok = lab(jj) > 0;
    if ~any(ok), break; end
    A = accumarray([ii(ok), lab(jj(ok))], ww(ok), [numel(un) K]);
    [mx, b] = max(A, [], 2);
    if ~any(mx > 0), break; end
    lab(un(mx > 0)) = b(mx > 0);
end
lab(lab == 0) = K + 1;
lab = uint32(lab);
end

function rows = addRow(rows, target, cloud, method, lab, gt, N, nComp, Kgt, Kgap, tm, tuneSet)
e = evalSegmentation(lab, gt);
[~, ~, li] = unique(double(lab)); [gu, ~, gi] = unique(double(gt));
Ct = accumarray([li gi], 1);
purity = sum(max(Ct, [], 2)) / sum(Ct(:));
fr = zeros(numel(gu), 1);
for q = 1:numel(gu)
    col = sort(Ct(:,q), 'descend'); fr(q) = find(cumsum(col) >= 0.9*sum(col), 1);
end
major = sum(Ct, 1).' >= 30;
r = struct('target', target, 'cloud', cloud, 'method', method, 'isTune', double(any(strcmp(target, tuneSet))), ...
           'N', N, 'nComp', nComp, 'Kgt', Kgt, 'Kgap', Kgap, 'K', numel(unique(lab)), ...
           'ARI', e.ARI, 'purity', purity, 'frag', mean(fr(major)), 'timeSec', tm);
rows = [rows; r];
end

function h = fileMD5(f)
fid = fopen(f, 'r'); b = fread(fid, Inf, '*uint8'); fclose(fid);
md = java.security.MessageDigest.getInstance('MD5'); md.update(b);
h = sprintf('%02x', typecast(md.digest, 'uint8'));
end

function drawLab(ax, P, lab)
[~, ~, li] = unique(double(lab)); K = max(li); cm = lines(max(K,1));
hold(ax,'on');
for k = 1:K
    m = li == k;
    scatter3(ax, P(m,1), P(m,2), P(m,3), 3, cm(k,:), 'filled');
end
hold(ax,'off'); axis(ax,'equal','tight','off'); view(ax, [1 -1.2 0.8]);
end
