% ISS 실데이터 스펙트럴 분할 파이프라인 검증 스크립트 (피팅 없음)
% 전처리는 ISS_normcut_4_modified.m 1~57행과 동일:
%   ISS_stationary.xyz 크롭(P1~P4) + DATA4.mat 세그먼트 5~10 보강 → 반경 thinning
%
% legacy(C0) → improved(C7) 사이를 한 축씩 바꾼 설정들을 전부 돌려 비교하고 저장
%   C0 legacy      : alpha=0, u_4..u_18 (레거시 cluster=3 하드코딩), MST 최장 간선 (80-m)개 절단
%                    ← ISS_normcut_4_modified.m 294~387행과 같은 동작 (비트 단위 일치 확인됨)
%   C1 skipfix     : 건너뛸 고유벡터 수를 실제 연결성분 수로
%   C2/C3 alpha    : Coifman–Lafon alpha = 0.5 / 1
%   C4 cond_t80    : conductance 판정, 컷 수 예산은 레거시와 동일(80-m)
%   C5 cond_free   : conductance 판정만으로 컷 수 결정 (상한 없음)
%   C6 cond_strict : C5 + phiMax 0.02
%   C7 improved    : eigengap K* 로 임베딩 차원/컷 수 결정
%   C9L phiorder   : phi 가 작은 컷부터 (후보는 긴 간선 상위 10%), phiMax = eta*sqrt(2*lambda_2)
%   C8a normal     : 법선 일치도로 W 재가중 (buildGraph 의 normalSigma)
%   C8b/c pref4/6  : 4차/6차 다항식 모델 선호 커널로 W 재가중 (modelPreferenceWeights)
%   C12 final      : 법선 sigma 0.3 + 국소화 모드 제외 + phi <= 0.005 인 간선을
%                    phi/길이 가 작은 순으로 절단, 최소 군집 크기 100 (시드 평균으로 고른 설정)
%
% 정답 라벨이 없으므로 진단 지표로 비교:
%   - 군집 수 / coverage(라벨>0 비율) / 군집 크기
%   - 군집별 conductance phi_k = cut(S_k, S_k~)/min(vol S_k, vol S_k~)  (작을수록 잘 분리)
%   - ARI(설정 vs C0 legacy)
%
% 출력: out_cutting_ISS/<yyyyMMdd_HHmmss>/
%   summary.csv    설정별 파라미터 + 지표 표
%   results.mat    T(표), cfg, labels(N x 설정 수), Pg(점 좌표), spec(alpha별 고유값/K*)
%   spectrum.png, overview.png, seg_<설정>.png
rng(1);

%% 0) 데이터 경로
xyzFile   = 'ISS_stationary.xyz';
data4File = fullfile('resources', 'DATA4.mat');   % resources/는 MATLAB 예약 폴더명이라 addpath 불가 → 경로로 로드

%% 1) 전처리 (ISS_normcut_4_modified.m 4~57행)
P  = readmatrix(xyzFile, 'FileType', 'text');
P1 = P(P(:,1) >= 13 & P(:,1) <= 17 & P(:,2) >= -5 & P(:,2) <= 15 & P(:,3) <= 3 & P(:,3) >= -21, :);
P2 = P(P(:,1) >= -17 & P(:,1) <= -13 & P(:,2) >= -5 & P(:,2) <= 15 & P(:,3) <= 3 & P(:,3) >= -21, :);
P3 = P(P(:,1) >= -20 & P(:,1) <= 25 & P(:,2) >= 2.5 & P(:,2) <= 9 & P(:,3) >= 3 & P(:,3) <= 8, :);
P4 = P(P(:,1) >= -12 & P(:,1) <= 12 & P(:,2) <= 2.5 & P(:,3) >= 2.5 & P(:,3) <= 13, :);

S4 = load(data4File);
v4 = fieldnames(S4);
D4 = S4.(v4{1});
Pi = vertcat(D4(5:10).P);

P_use = [P1(1:7:end,:); P2(1:6:end,:); P3(1:11:end,:); P4(1:9:end,:); Pi(1:10:end,:)];

% 반경 thinning: 무작위 순서로 훑으며 반경 r 내 이웃 억제 (레거시와 동일 규칙, 이웃 목록만 일괄 계산)
r_thin = 0.15;
N0   = size(P_use,1);
nbrs = rangesearch(P_use, P_use, r_thin);
ord  = randperm(N0);
keep = false(N0,1); blocked = false(N0,1);
for t = 1:N0
    i = ord(t);
    if blocked(i), continue; end
    keep(i) = true;
    blocked(nbrs{i}) = true;
end
P_use = P_use(keep,:);
fprintf('[thin] kept %d / %d (%.1f%%), r=%.3f\n', size(P_use,1), N0, 100*size(P_use,1)/N0, r_thin);

%% 2) 그래프 (레거시 62~134행: mutual kNN k=30, self-tuning gamma=1, tau=3번째 이웃, 성분<150 제거)
gOpts = struct('k',30, 'gamma',1.0, 'tauIdx',3, 'minCompSize',150);
[W, gInfo] = buildGraph(P_use, gOpts);
Pg    = P_use(gInfo.idxKeep, :);
nComp = max(conncomp(graph(W, 'upper')));        % 영고유값 개수 = 연결 성분 수
fprintf('그래프: N=%d, 간선=%d, 연결성분=%d\n', size(W,1), nnz(triu(W,1)), nComp);

kEig = 45;

%% 3) 비교 설정 (legacy → improved 사이를 한 축씩 바꿔가며)
%   skip   : 건너뛸 고유벡터 수. 'legacy'=3 하드코딩, 'ncomp'=연결성분 수
%   dim    : 임베딩 차원. 숫자 또는 'eigengap' (= K* - skip)
%   weight : 가중치 W. 'none' | 'normal'(법선 일치도) | 'pref4'/'pref6'(다항식 선호 커널)
%   wpar   : normal → normalSigma, pref → sigmaD
%   cuts   : 컷 수 상한. 'target80'=레거시(80-m), 'eigengap'=K*-1, 'none'=상한 없음
%   eta    : phiorder 전용. phiMax = eta * sqrt(2*lambda_2)
specRows = {
%  name                alpha  skip      dim         weight   wpar  method         qThr  phiMax  cuts        eta
   'C0_legacy',        0,     'legacy', 15,         'none',  NaN,  'length',      0,    NaN,    'target80', NaN
   'C1_skipfix',       0,     'ncomp',  15,         'none',  NaN,  'length',      0,    NaN,    'target80', NaN
   'C2_alpha05',       0.5,   'ncomp',  15,         'none',  NaN,  'length',      0,    NaN,    'target80', NaN
   'C3_alpha1',        1,     'ncomp',  15,         'none',  NaN,  'length',      0,    NaN,    'target80', NaN
   'C4_cond_t80',      1,     'ncomp',  15,         'none',  NaN,  'conductance', 0.90, 0.05,   'target80', NaN
   'C5_cond_free',     1,     'ncomp',  15,         'none',  NaN,  'conductance', 0.90, 0.05,   'none',     NaN
   'C6_cond_strict',   1,     'ncomp',  15,         'none',  NaN,  'conductance', 0.90, 0.02,   'none',     NaN
   'C7_improved',      1,     'ncomp',  'eigengap', 'none',  NaN,  'conductance', 0.90, 0.05,   'eigengap', NaN
   'C9L_phiord_e0.5',  1,     'ncomp',  15,         'none',  NaN,  'phiorder',    0.90, NaN,    'none',     0.5
   'C9L_phiord_e1.0',  1,     'ncomp',  15,         'none',  NaN,  'phiorder',    0.90, NaN,    'none',     1
   'C8a_normal_s0.2',  1,     'ncomp',  15,         'normal', 0.2, 'phiorder',    0.90, NaN,    'none',     1
   'C8b_pref4_s0.3',   1,     'ncomp',  15,         'pref4',  0.3, 'phiorder',    0.90, NaN,    'none',     1
   'C8c_pref6_s0.3',   1,     'ncomp',  15,         'pref6',  0.3, 'phiorder',    0.90, NaN,    'none',     1
};
cfg = cell2struct(specRows, {'name','alpha','skip','dim','weight','wpar', ...
                             'method','qThr','phiMax','cuts','eta'}, 2);
[cfg.extra] = deal(struct());

% C12 최종 후보 — 시드 6개 평균으로 고른 설정 (run_seed_check_ISS.m 참고)
%   법선 sigma 0.3, 국소화된 고유벡터 제외(minEff), 후보 = 전체 간선,
%   phi <= 0.005 (절대 기준) 인 간선 중 phi/(길이)^1 이 작은 것부터, 최소 군집 크기 100
cfg(end+1) = struct('name','C12_final', 'alpha',1, 'skip','ncomp', 'dim',15, ...
    'weight','normal', 'wpar',0.3, 'method','phiorder', 'qThr',0, 'phiMax',0.005, ...
    'cuts','none', 'eta',NaN, ...
    'extra', struct('minEff',60, 'candidates','all', 'lenPower',1, 'minSize',100));
nC = numel(cfg);

%% 4) 가중치별 W + (가중치, alpha)별 스펙트럼
wkey = arrayfun(@(c) sprintf('%s_%g', c.weight, c.wpar), cfg, 'UniformOutput', false);
[uw, ~, wIdx] = unique(wkey, 'stable');
Wall = cell(numel(uw), 1);
for q = 1:numel(uw)
    c = cfg(find(wIdx == q, 1));
    t = tic;
    switch c.weight
        case 'none'
            Wall{q} = W;
        case 'normal'
            gO = gOpts; gO.normalSigma = c.wpar;
            [Wn, giN] = buildGraph(P_use, gO);
            assert(isequal(giN.idxKeep, gInfo.idxKeep), '법선 가중 후 점 집합이 달라졌습니다');
            Wall{q} = Wn;
        case {'pref4', 'pref6'}
            ordp = str2double(c.weight(end));
            Wall{q} = modelPreferenceWeights(Pg, W, struct('order', ordp, 'sigmaD', c.wpar));
        otherwise
            error('알 수 없는 weight "%s"', c.weight);
    end
    if ~strcmp(c.weight, 'none')
        fprintf('가중치 %s 준비 (%.1fs)\n', uw{q}, toc(t));
    end
end

skey = arrayfun(@(i) sprintf('%s|%g', wkey{i}, cfg(i).alpha), (1:nC)', 'UniformOutput', false);
[us, ~, sIdx] = unique(skey, 'stable');
spec = struct('key', {}, 'alpha', {}, 'lambda', {}, 'Y', {}, 'Kgap', {}, 'gaps', {}, 'nComp', {}, 'Neff', {});
for q = 1:numel(us)
    i = find(sIdx == q, 1); Wc = Wall{wIdx(i)};
    [L, Dis] = normalizeGraph(Wc, cfg(i).alpha);
    [U, lam] = embedSpectral(L, kEig);
    [Kg, gp] = selectEigenGap(lam, struct('kMin',2, 'kMax',40));
    spec(end+1) = struct('key', us{q}, 'alpha', cfg(i).alpha, 'lambda', lam, 'Y', Dis*U, ...
                         'Kgap', Kg, 'gaps', gp, 'nComp', max(conncomp(graph(Wc,'upper'))), ...
                         'Neff', 1 ./ sum(U.^4, 1));   %#ok<SAGROW> 고유벡터별 유효 점 수
    fprintf('%-14s : lambda 1~4 = %s, K* = %d, sqrt(2*l2) = %.4f\n', ...
        us{q}, mat2str(lam(1:4)',3), Kg, sqrt(2*lam(2)));
end

%% 5) 설정별 분할 + 지표
N = size(W,1);
labels = zeros(N, nC, 'uint32');
gtFile = fullfile('out_cutting_ISS','gt','gt_labels.mat');    % GLB 부품 라벨 (있으면 채점)
useGT = false;
if exist(gtFile, 'file')
    GT = load(gtFile);
    if isequal(GT.Pg, Pg)
        [gu, ~, gGrp] = unique(GT.groupId);
        gMajor = find(accumarray(gGrp,1) >= 60);
        useGT = true;
    end
end
rows = struct([]);
for c = 1:nC
    s  = spec(sIdx(c));
    Wc = Wall{wIdx(c)};
    ex = cfg(c).extra;
    if strcmp(cfg(c).skip, 'legacy'), skip = 3; else, skip = s.nComp; end
    if ischar(cfg(c).dim), dim = max(s.Kgap - skip, 2); else, dim = cfg(c).dim; end
    modes = (skip+1):size(s.Y,2);
    modes = modes(s.Neff(modes) >= getx(ex, 'minEff', 0));    % 국소화된 모드 제외 (기본: 제외 안 함)
    Yuse  = s.Y(:, modes(1:dim));

    switch cfg(c).cuts
        case 'target80', maxCuts = max(80 - embedComponents(Yuse, 40), 0);
        case 'eigengap', maxCuts = max(s.Kgap - 1, 1);
        case 'none',     maxCuts = Inf;
    end
    sOpts = struct('ky',40, 'cutMethod',cfg(c).method, 'qThr',cfg(c).qThr, ...
                   'minSize',getx(ex,'minSize',60), 'maxCuts',maxCuts);
    phiUsed = cfg(c).phiMax;
    if strcmp(cfg(c).method, 'phiorder')
        if ~isnan(cfg(c).eta), phiUsed = cfg(c).eta * sqrt(2*s.lambda(2)); end   % eta 없으면 절대 기준
        sOpts.lenQuantile = cfg(c).qThr;
        sOpts.candidates  = getx(ex, 'candidates', 'length');
        sOpts.lenPower    = getx(ex, 'lenPower', 0);
    end
    if ~isnan(phiUsed), sOpts.phiMax = phiUsed; end

    t = tic;
    [lab, dg] = segmentEmbedding(Yuse, Wc, sOpts);
    tc = toc(t);
    labels(:,c) = lab;

    phi = clusterConductance(Wc, lab);
    sz  = accumarray(double(lab(lab>0)), 1);
    ari = evalSegmentation(lab, labels(:,1));      % C0 legacy 대비 (예측 0라벨 제외)

    purity = NaN; fragAvg = NaN; ariGT = NaN;
    if useGT
        [purity, fragAvg, ariGT] = scoreAgainstGT(lab, gGrp, numel(gu), gMajor);
    end

    r = struct('name', cfg(c).name, 'alpha', cfg(c).alpha, 'skip', skip, 'dim', dim, ...
               'weight', cfg(c).weight, 'wpar', cfg(c).wpar, ...
               'method', cfg(c).method, 'qThr', cfg(c).qThr, 'phiMax', phiUsed, ...
               'cutRule', cfg(c).cuts, 'maxCuts', maxCuts, 'nCuts', size(dg.cutEdges,1), ...
               'K', double(max(lab)), 'coverage', mean(lab>0), ...
               'phiMed', median(phi), 'phiWorst', max(phi), ...
               'sizeMin', min(sz), 'sizeMed', median(sz), 'sizeMax', max(sz), ...
               'purity', purity, 'fragAvg', fragAvg, 'ARIgt', ariGT, ...
               'ARIvsLegacy', ari.ARI, 'timeSec', tc);
    rows = [rows; r]; %#ok<AGROW>
    fprintf(['%-17s K=%3d cov=%.3f cuts=%4d dim=%2d phi med=%.4f ' ...
             'purity=%.3f frag=%.2f ARIgt=%.3f ARI(C0)=%.3f (%.1fs)\n'], ...
        r.name, r.K, r.coverage, r.nCuts, r.dim, r.phiMed, ...
        r.purity, r.fragAvg, r.ARIgt, r.ARIvsLegacy, r.timeSec);
end
T = struct2table(rows);

%% 6) 저장
outDir = fullfile('out_cutting_ISS', char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')));
if ~exist(outDir, 'dir'), mkdir(outDir); end
writetable(T, fullfile(outDir, 'summary.csv'));
spec = rmfield(spec, 'Y');                          % 임베딩 좌표는 용량 때문에 제외
save(fullfile(outDir, 'results.mat'), 'T', 'cfg', 'labels', 'Pg', 'spec', ...
     'nComp', 'gOpts', 'r_thin', 'kEig');

% 스펙트럼
figure(20); clf;
subplot(1,2,1); hold on;
for i = 1:numel(spec), plot(1:kEig, spec(i).lambda, 'o-', 'LineWidth', 1.2); end
hold off; grid on; xlabel('index k'); ylabel('\lambda_k'); title('L spectrum');
legend({spec.key}, 'Location', 'northwest', 'Interpreter', 'none', 'FontSize', 8);
subplot(1,2,2); hold on;
for i = 1:numel(spec), plot(1:kEig-1, spec(i).gaps, 'o-', 'LineWidth', 1.2); end
hold off; grid on; xlabel('k'); ylabel('\lambda_{k+1}-\lambda_k');
title(sprintf('eigengap (K* = %s)', mat2str([spec.Kgap])));
exportgraphics(figure(20), fullfile(outDir, 'spectrum.png'), 'Resolution', 120);

% 설정별 개별 그림
for c = 1:nC
    figure(30); clf;
    plotSeg(gca, Pg, labels(:,c), segTitle(T(c,:)));
    exportgraphics(figure(30), fullfile(outDir, ['seg_' cfg(c).name '.png']), 'Resolution', 150);
end

% 전체 한눈에
figure(21); clf;
set(gcf, 'Position', [50 50 1600 800]);
for c = 1:nC
    plotSeg(subplot(2, ceil(nC/2), c), Pg, labels(:,c), segTitle(T(c,:)));
end
exportgraphics(figure(21), fullfile(outDir, 'overview.png'), 'Resolution', 120);

fprintf('\n저장 완료: %s\n', outDir);

%% 로컬 함수
function v = getx(s, f, d)
% 설정별 추가 옵션(extra) 읽기: 없으면 기본값
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function m = embedComponents(Y, ky)
% segmentEmbedding 내부와 같은 임베딩 kNN 그래프의 연결 성분 수 (레거시 need_cut 계산용)
N = size(Y,1);
[idxNN, distNN] = knnsearch(Y, Y, 'K', ky+1);
S = sparse(repmat((1:N)', ky, 1), reshape(idxNN(:,2:end),[],1), ...
           reshape(distNN(:,2:end),[],1), N, N);
m = max(conncomp(graph(max(S, S.'), 'upper')));
end

function [purity, fragAvg, ariGT] = scoreAgainstGT(lab, gGrp, nGrp, gMajor)
% GLB 부품 라벨 기준 채점
%   purity  : 군집에서 가장 많은 부품이 차지하는 비율 (서로 다른 부품을 합칠수록 낮음)
%   fragAvg : 부품 점의 90%를 덮는 데 필요한 군집 수의 평균 (쪼갤수록 높음)
lab = double(lab); v = lab > 0;
Ct = accumarray([lab(v) gGrp(v)], 1, [max(lab) nGrp]);
purity = sum(max(Ct,[],2)) / sum(Ct(:));
fr = zeros(numel(gMajor),1);
for q = 1:numel(gMajor)
    col = sort(Ct(:,gMajor(q)), 'descend');
    fr(q) = find(cumsum(col) >= 0.9*sum(col), 1);
end
fragAvg = mean(fr);
e = evalSegmentation(uint32(lab), gGrp);
ariGT = e.ARI;
end

function phi = clusterConductance(W, labels)
% 군집별 conductance (라벨 0 제외)
d = full(sum(W,2));
volTot = sum(d);
K = double(max(labels));
phi = zeros(K,1);
for k = 1:K
    S = (labels == k);
    volS = sum(d(S));
    cutS = volS - full(sum(sum(W(S,S))));
    phi(k) = cutS / max(min(volS, volTot - volS), eps);
end
end

function str = segTitle(r)
str = sprintf('%s: K=%d, cov %.2f, phi_med %.3f', r.name{1}, r.K, r.coverage, r.phiMed);
end

function plotSeg(ax, P, lab, ttl)
axes(ax); hold on;
idx0 = (lab == 0);
if any(idx0)
    scatter3(P(idx0,1), P(idx0,2), P(idx0,3), 1, [0.6 0.5 0.5], ...
             'filled', 'MarkerFaceAlpha', 0.4);
end
Ks = double(max(lab));
cmap = lines(max(Ks,1));
for kk = 1:Ks
    m = (lab == kk);
    scatter3(P(m,1), P(m,2), P(m,3), 3, cmap(kk,:), 'filled', 'MarkerFaceAlpha', 0.7);
end
hold off; axis equal tight; grid on; view([1 1 1]);
title(ttl, 'Interpreter', 'none'); xlabel('X'); ylabel('Y'); zlabel('Z');
end
