% 시뮬레이터 점군(sim/raycast_target.py → out_cutting_sim)에 대한 스펙트럴 커팅 진단
%
% 설정 (그래프/임베딩/컷)
%   C7_eigengap : 법선 가중 없음, 임베딩 차원 = eigengap K* - skip,
%                 컷 = conductance 판정(qThr 0.90, phiMax 0.05), 컷 수 K*-1
%                 ← test_occlusion_spectral.m / run_pipeline_ISS.m 의 C7 규칙
%   C8a_K2      : PCA 법선 가중 sigma=0.2, 임베딩 15차원,
%                 컷 = phi 작은 것부터(후보 = 긴 간선 상위 10%), phi > sqrt(2*lambda) 에서 정지
%                 ← run_cut_compare_ISS.m 의 K2_phi_len (ISS 기준 현재 최적)
%   둘 다 alpha = 1 (Coifman–Lafon 밀도 보정), mutual kNN k=30, self-tuning
%
% 연결성분 처리 (같은 설정을 두 방식으로)
%   global  : 현행. 영고유벡터 nComp 개를 건너뛰고 전체를 한 번에 임베딩/분할
%   percomp : W 의 연결성분마다 따로 라플라시안 → 임베딩 → 분할 (성분은 절대 섞이지 않음)
%
% 출력: out_cutting_simcut/<yyyyMMdd_HHmmss>_<target>/
%   summary.csv, results.mat, diag_<cloud>_<cfg>.png
%     그림 8칸: GT / 스펙트럼+간격 / 고유벡터 2개 / 임베딩 공간의 포레스트와 컷 /
%               3D 에 그린 포레스트와 컷 / global 분할 / percomp 분할
rng(1);

%% 0) 입력
if ~exist('target','var') || isempty(target), target = 'cubesat6u'; end   % 밖에서 target 을 주면 그 타깃
clouds = {'complete', 's1_mrbar_h3', 's1_oblique_h3', 's1_pvbar_h3'};
kEig   = 20;
gBase  = struct('k',30, 'gamma',1.0, 'tauIdx',3, 'minCompSize',50);
cfg = struct( ...
    'name',        {'C7_eigengap', 'C8a_K2'}, ...
    'normalSigma', {0, 0.2}, ...
    'dim',         {'eigengap', 15}, ...
    'cut',         {'cond', 'phiorder'}, ...
    'minSize',     {60, 60});

d = dir(fullfile('out_cutting_sim', ['*_' target]));
assert(~isempty(d), 'out_cutting_sim 에 %s 점군이 없습니다 (sim/raycast_target.py 먼저 실행)', target);
runDir = fullfile(d(end).folder, d(end).name);
outDir = fullfile('out_cutting_simcut', [char(datetime('now','Format','yyyyMMdd_HHmmss')) '_' target]);
if ~exist(outDir, 'dir'), mkdir(outDir); end
fprintf('입력: %s\n출력: %s\n', runDir, outDir);

%% 1) 점군 x 설정
rows = struct([]); res = struct([]);
for f = 1:numel(clouds)
    C = loadSimCloud(fullfile(runDir, [clouds{f} '.mat']));
    for c = 1:numel(cfg)
        gOpts = gBase; gOpts.normalSigma = cfg(c).normalSigma;
        [W, gi] = buildGraph(C.P, gOpts);
        P  = C.P(gi.idxKeep, :);
        gt = C.label(gi.idxKeep);
        comp  = conncomp(graph(W, 'upper')).';
        nComp = max(comp);

        tG = tic; [labG, dG] = segOne(W, cfg(c), nComp, kEig); tG = toc(tG);
        tP = tic; labP = segPerComp(W, comp, cfg(c), kEig);    tP = toc(tP);

        mG = evalSegmentation(labG, gt);
        mP = evalSegmentation(labP, gt);
        mC = evalSegmentation(uint32(comp), gt);
        r = struct('cloud', clouds{f}, 'cfg', cfg(c).name, 'N', size(W,1), 'nComp', nComp, ...
                   'Kgt', numel(unique(gt)), 'Kgap', dG.Kgap, 'dim', dG.dim, 'phiStop', dG.phiC, ...
                   'K_global', double(max(labG)), 'ARI_global', mG.ARI, 'cov_global', mG.coverage, ...
                   'K_percomp', double(max(labP)), 'ARI_percomp', mP.ARI, 'cov_percomp', mP.coverage, ...
                   'ARI_compOnly', mC.ARI, 'tGlobal', tG, 'tPercomp', tP);
        rows = [rows; r]; %#ok<AGROW>
        res(end+1).cloud = clouds{f}; res(end).cfg = cfg(c).name; %#ok<SAGROW>
        res(end).P = P; res(end).gt = gt; res(end).labG = labG; res(end).labP = labP;
        res(end).lambda = dG.lambda; res(end).cutEdges = dG.cutEdges;
        fprintf('%-14s %-12s N=%5d comp=%d Kgt=%d K*=%2d | global K=%2d ARI=%.3f | percomp K=%2d ARI=%.3f | 성분만 ARI=%.3f\n', ...
            r.cloud, r.cfg, r.N, r.nComp, r.Kgt, r.Kgap, r.K_global, r.ARI_global, ...
            r.K_percomp, r.ARI_percomp, r.ARI_compOnly);

        drawDiag(fullfile(outDir, sprintf('diag_%s_%s.png', clouds{f}, cfg(c).name)), ...
                 P, gt, C.partNames, comp, dG, labG, labP, r);
    end
end
T = struct2table(rows);
writetable(T, fullfile(outDir, 'summary.csv'));
save(fullfile(outDir, 'results.mat'), 'T', 'res', 'cfg', 'gBase', 'kEig', 'runDir');
fprintf('\n저장 완료: %s\n', outDir);

%% ===================== local functions =====================
function [lab, d] = segOne(W, cfg, skip, kEig)
% 연결성분 skip 개짜리 그래프를 한 번에 임베딩/분할 (현행 방식)
N  = size(W,1);
kE = min(kEig, N-2);
[L, Dis] = normalizeGraph(W, 1);
[U, lam] = embedSpectral(L, kE);
Y    = Dis * U;
Kgap = selectEigenGap(lam, struct('kMin', max(2,skip), 'kMax', min(15, kE-1)));
if ischar(cfg.dim), dim = max(Kgap - skip, 2); else, dim = cfg.dim; end
dim  = min(dim, kE - skip);
Yuse = Y(:, skip+1 : skip+dim);
phiC = sqrt(2 * max(lam(skip+1), 0));           % 첫 비영 고유값 기준 Cheeger 규모

o = struct('ky', min(40, N-2), 'minSize', cfg.minSize);
switch cfg.cut
    case 'cond'
        o.cutMethod = 'conductance'; o.qThr = 0.90; o.phiMax = 0.05; o.maxCuts = max(Kgap-1, 1);
    case 'phiorder'
        o.cutMethod = 'phiorder'; o.phiMax = phiC; o.candidates = 'length'; o.lenQuantile = 0.90;
end
[lab, dg] = segmentEmbedding(Yuse, W, o);
d = struct('lambda', lam, 'Y', Y, 'Yuse', Yuse, 'Kgap', Kgap, 'dim', dim, 'skip', skip, ...
           'phiC', phiC, 'cutEdges', dg.cutEdges, 'forest', dg.T_all.Edges.EndNodes);
end

function lab = segPerComp(W, comp, cfg, kEig)
% 연결성분마다 따로 분할. 작은 성분(3*minSize 미만)은 자르지 않고 군집 하나로 둔다
lab = zeros(size(W,1), 1, 'uint32'); off = 0;
for c = 1:max(comp)
    idx = find(comp == c);
    if numel(idx) < 3*cfg.minSize
        lc = ones(numel(idx), 1, 'uint32');
    else
        lc = segOne(W(idx,idx), cfg, 1, kEig);
        if ~any(lc), lc(:) = 1; end
    end
    lab(idx(lc > 0)) = off + lc(lc > 0);
    off = off + double(max(lc));
end
end

function drawDiag(file, P, gt, partNames, comp, dG, labG, labP, r)
f = figure('Visible','off', 'Position',[30 30 2200 1000]);
tl = tiledlayout(f, 2, 4, 'TileSpacing','compact', 'Padding','compact');
title(tl, sprintf('%s  |  %s  |  N=%d, 연결성분 %d, 보이는 부품 %d', ...
    r.cloud, r.cfg, r.N, r.nComp, r.Kgt), 'Interpreter','none');

ax = nexttile(tl); drawLab(ax, P, gt);
title(ax, ['GT: ' strjoin(partNames, ' / ')], 'Interpreter','none');

ax = nexttile(tl); lam = dG.lambda; n = numel(lam);
yyaxis(ax,'left');  plot(ax, 1:n, lam, 'o-', 'LineWidth',1.2); ylabel(ax,'\lambda_k');
yyaxis(ax,'right'); bar(ax, (1:n-1)+0.5, diff(lam), 0.5, 'FaceAlpha',0.35, 'EdgeColor','none');
ylabel(ax,'\lambda_{k+1}-\lambda_k');
xline(ax, r.Kgap+0.5, '--r', sprintf('eigengap K*=%d', r.Kgap), 'LabelVerticalAlignment','bottom');
xline(ax, r.Kgt+0.5, ':', sprintf('정답 K=%d', r.Kgt), 'Color',[0 0.5 0], 'LabelVerticalAlignment','middle');
grid(ax,'on'); xlabel(ax,'k'); title(ax, sprintf('스펙트럼 (영고유값 %d개)', r.nComp));

for j = 1:2
    ax = nexttile(tl); k = dG.skip + j;
    scatter3(ax, P(:,1), P(:,2), P(:,3), 4, dG.Y(:,k), 'filled');
    axis(ax,'equal','tight'); view(ax, [1 -1.2 0.8]); colormap(ax, turbo); colorbar(ax);
    title(ax, sprintf('고유벡터 y_%d (\\lambda=%.2e)', k, lam(k)));
end

% 임베딩 공간: 포레스트(회색) + 실제로 자른 간선(빨강)
ax = nexttile(tl); Yu = dG.Yuse; if size(Yu,2) < 3, Yu(:,3) = 0; end
hold(ax,'on');
drawEdges(ax, Yu, dG.forest, [0.6 0.6 0.6], 0.4);
cm = lines(max(double(max(gt)),1));
scatter3(ax, Yu(:,1), Yu(:,2), Yu(:,3), 5, cm(double(gt),:), 'filled');
if ~isempty(dG.cutEdges), drawEdges(ax, Yu, dG.cutEdges(:,1:2), [1 0 0], 2.5); end
hold(ax,'off'); grid(ax,'on'); view(ax, 3);
xlabel(ax, sprintf('y_%d', dG.skip+1)); ylabel(ax, sprintf('y_%d', dG.skip+2)); zlabel(ax, sprintf('y_%d', dG.skip+3));
title(ax, sprintf('임베딩 공간 (앞 3축 / 총 %d차원), 색 = GT, 빨강 = 컷 %d개', r.dim, size(dG.cutEdges,1)));

% 3D 공간에 그린 같은 포레스트
ax = nexttile(tl); hold(ax,'on');
drawEdges(ax, P, dG.forest, [0.45 0.45 0.45], 0.4);
if ~isempty(dG.cutEdges), drawEdges(ax, P, dG.cutEdges(:,1:2), [1 0 0], 2.5); end
hold(ax,'off'); axis(ax,'equal','tight'); grid(ax,'on'); view(ax, [1 -1.2 0.8]);
title(ax, '임베딩 MST 를 3D 에 그린 것 (빨강 = 컷)');

ax = nexttile(tl); drawLab(ax, P, labG);
title(ax, sprintf('global(현행): K=%d, ARI=%.3f', r.K_global, r.ARI_global));
ax = nexttile(tl); drawLab(ax, P, labP);
title(ax, sprintf('percomp: K=%d, ARI=%.3f  (성분만: %.3f)', r.K_percomp, r.ARI_percomp, r.ARI_compOnly));

exportgraphics(f, file, 'Resolution', 100); close(f);
end

function drawLab(ax, P, lab)
lab = double(lab); K = max(lab); cm = lines(max(K,1));
hold(ax,'on');
m0 = lab == 0;
if any(m0), scatter3(ax, P(m0,1), P(m0,2), P(m0,3), 3, [0.7 0.7 0.7], 'filled'); end
for k = 1:K
    m = lab == k;
    if any(m), scatter3(ax, P(m,1), P(m,2), P(m,3), 4, cm(k,:), 'filled'); end
end
hold(ax,'off'); axis(ax,'equal','tight'); grid(ax,'on'); view(ax, [1 -1.2 0.8]);
xlabel(ax,'x R-bar'); ylabel(ax,'y V-bar'); zlabel(ax,'z H-bar');
end

function drawEdges(ax, X, E, col, lw)
if isempty(E), return; end
n = size(E,1);
xs = [X(E(:,1),1), X(E(:,2),1), nan(n,1)].';
ys = [X(E(:,1),2), X(E(:,2),2), nan(n,1)].';
zs = [X(E(:,1),3), X(E(:,2),3), nan(n,1)].';
plot3(ax, xs(:), ys(:), zs(:), '-', 'Color', col, 'LineWidth', lw);
end
