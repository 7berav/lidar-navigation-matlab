% 가림(occlusion)이 스펙트럴 임베딩/분할에 주는 영향 테스트 (샘플 1개)
%
% 씬: ISS 유사 구조 — 원통 모듈 A, 연결 원통(connector), 원통 모듈 B, 모듈 B에 붙은 태양전지판(얇은 판)
% 가림: Hidden Point Removal (Katz et al. 2007) 로 시점별 가시 점만 남김
%
% 비교 (전체 그래프 vs 가려진 그래프, 공통 점 기준)
%   1) 고유값 lambda_1..15
%   2) 개별 고유벡터 상관 |corr(y_full_k, y_occ_k)|  ← "값"이 유지되는가
%   3) 앞쪽 K개 고유벡터 부분공간의 principal angle  ← "부분공간"이 유지되는가
%   4) 정답 지시벡터가 앞쪽 K개 부분공간에 담기는 비율 (cluster capture)
%   5) 파이프라인(eigengap K* + conductance 컷) 분할 ARI
%        vs 원 GT / vs 가시 연결성분으로 쪼갠 GT / vs 전체 그래프 분할
%
% 출력: out_cutting_occlusion/<yyyyMMdd_HHmmss>/ overview.png, subspace.png, results.mat
rng(1);

%% 0) 파라미터
rho    = 100;     % 표면 점 밀도 [pts / unit^2]
noise  = 0.01;    % 가우시안 노이즈 sigma
alpha  = 1;       % Coifman–Lafon alpha
kEig   = 20;      % 계산할 고유쌍 수
hprExp = 2.5;     % HPR 반경 지수 (클수록 더 많이 보임)
gOpts  = struct('k',30, 'gamma',1.0, 'tauIdx',3, 'minCompSize',50);

%% 1) 합성 씬 (x축 방향으로 A - connector - B, B 옆에 판)
Rx = [0 0 1; 0 1 0; -1 0 0];            % 로컬 z축 -> 전역 x축
cyl = struct('name',{'moduleA','connector','moduleB'}, ...
             'r',   {1.0, 0.4, 0.7}, ...
             'h',   {3.0, 0.5, 2.5}, ...          % 반길이
             't',   {[-3.5 0 0], [0 0 0], [3.0 0 0]}, ...
             'caps',{true, false, true});
P = zeros(0,3); gt = zeros(0,1,'uint32');
for i = 1:numel(cyl)
    c = cyl(i);
    area = 2*pi*c.r*2*c.h + c.caps*2*pi*c.r^2;
    Pi = sampleCylinder(round(rho*area), c.r, c.h, c.caps) * Rx.' + c.t;
    P  = [P; Pi]; gt = [gt; repmat(uint32(i), size(Pi,1), 1)]; %#ok<AGROW>
end
panelHalf = [1.2 2.5 0.02];
panelArea = 8*(panelHalf(1)*panelHalf(2) + panelHalf(2)*panelHalf(3) + panelHalf(1)*panelHalf(3));
panel = struct('type','box', 'N',round(rho*panelArea), 'size',panelHalf, ...
               'R',eye(3), 't',[3.0, 0.7+panelHalf(2), 0]);   % 판 가장자리가 모듈 B 옆면에 닿음
Pp = synthPointCloud(panel, 0);
P  = [P; Pp]; gt = [gt; repmat(uint32(4), size(Pp,1), 1)];
P  = P + noise*randn(size(P));
partNames = {'moduleA','connector','moduleB','panel'};
Kgt = numel(partNames);
fprintf('씬: %d점, 파트 %d개 (%s)\n', size(P,1), Kgt, strjoin(partNames, ', '));

%% 2) 전체(가림 없음) 그래프 스펙트럼 + 분할
full = analyzeCloud(P, (1:size(P,1))', gt, gOpts, alpha, kEig);
fprintf('\n[full]  N=%d  comp=%d  K*=%d  ARI(GT)=%.3f  K_pred=%d\n', ...
    numel(full.idx), full.nComp, full.Kgap, full.ariGT, full.Kpred);

%% 3) 시점별 가림
views = struct('name', {'top(+z)', 'side(-y)', 'oblique', 'end(+x)'}, ...
               'C',    {[0 0 25], [0 -25 0], [15 15 10], [25 0 2]});
nV = numel(views);
res = struct([]);
for v = 1:nV
    vis    = hprVisible(P, views(v).C, hprExp);
    visIdx = find(vis);
    occ    = analyzeCloud(P(visIdx,:), visIdx, gt(visIdx), gOpts, alpha, kEig);

    % 공통 점 (두 그래프 모두 살아남은 원본 인덱스)
    [~, iaF, ibO] = intersect(full.idx, occ.idx);
    Yf = full.Y(iaF, :);
    Yo = occ.Y(ibO, :);
    gtC = gt(full.idx(iaF));

    % (2) 개별 고유벡터 상관: 같은 번호끼리 / 가장 잘 맞는 번호
    kk = 2:min(6, kEig);
    Cm = abs(corr(Yf(:, kk), Yo(:, 2:10)));
    corrSame = diag(Cm(:, kk-1)).';
    corrBest = max(Cm, [], 2).';

    % (3) 부분공간 principal angle: 앞쪽 K개끼리, 그리고 full K개 ⊂ occ K+2개?
    angK    = principalAngles(Yf(:, 1:Kgt), Yo(:, 1:Kgt));
    angKsub = principalAngles(Yf(:, 1:Kgt), Yo(:, 1:Kgt+2));

    % (4) 정답 지시벡터 capture (공통 점 기준, 같은 GT)
    capF = clusterCapture(Yf(:, 1:Kgt), gtC);
    capO = clusterCapture(Yo(:, 1:Kgt), gtC);

    % (5) 분할 비교
    ariFullVsOcc = evalSegmentation(occ.labels(ibO), full.labels(iaF));

    r = struct('view', views(v).name, 'C', views(v).C, 'Nvis', numel(occ.idx), ...
               'visFrac', numel(occ.idx)/numel(full.idx), 'nComp', occ.nComp, ...
               'partVisFrac', accumarray(double(gt(occ.idx)), 1, [Kgt 1]).' ./ ...
                              accumarray(double(gt(full.idx)), 1, [Kgt 1]).', ...
               'lambda', occ.lambda.', 'corrSame', corrSame, 'corrBest', corrBest, ...
               'angK', angK.', 'angKsub', angKsub.', 'capFull', capF, 'capOcc', capO, ...
               'Kgap', occ.Kgap, 'Kpred', occ.Kpred, 'KgtRefined', occ.KgtRef, ...
               'ariGT', occ.ariGT, 'ariGTref', occ.ariGTref, 'ariVsFull', ariFullVsOcc.ARI);
    res = [res; r]; %#ok<AGROW>
    occAll{v} = occ; %#ok<SAGROW>

    fprintf('\n[%s]  가시 %d점 (%.0f%%), 파트별 가시율 %s, comp=%d\n', r.view, r.Nvis, ...
        100*r.visFrac, mat2str(round(100*r.partVisFrac)), r.nComp);
    fprintf('   lambda2..5   full %s  occ %s\n', mat2str(full.lambda(2:5).',3), mat2str(occ.lambda(2:5).',3));
    fprintf('   |corr| 같은번호 u2..u6 %s   best-match %s\n', mat2str(corrSame,2), mat2str(corrBest,2));
    fprintf('   principal angle K=%d  [deg] %s   (full K ⊂ occ K+2: %s)\n', Kgt, ...
        mat2str(round(angK.',1)), mat2str(round(angKsub.',1)));
    fprintf('   cluster capture  full %.3f  occ %.3f\n', capF, capO);
    fprintf('   분할 K*=%d K_pred=%d | ARI GT=%.3f  GT(가시성분 %d개)=%.3f  vs full분할=%.3f\n', ...
        r.Kgap, r.Kpred, r.ariGT, r.KgtRefined, r.ariGTref, r.ariVsFull);
end

%% 4) 저장 + 그림
outDir = fullfile('out_cutting_occlusion', char(datetime('now','Format','yyyyMMdd_HHmmss')));
if ~exist(outDir,'dir'), mkdir(outDir); end
save(fullfile(outDir,'results.mat'), 'P', 'gt', 'full', 'res', 'views', 'occAll');

f1 = figure('Visible','off', 'Position',[50 50 1600 330*(nV+1)]);
tl = tiledlayout(nV+1, 4, 'TileSpacing','compact', 'Padding','compact');
sets = [{full}, occAll];
rowNames = [{'full'}, {views.name}];
for r = 1:nV+1
    s  = sets{r}; Ps = P(s.idx,:);
    nexttile; plotLabels(Ps, gt(s.idx)); title([rowNames{r} ' : GT']);
    if r > 1, hold on; plotCam(views(r-1).C); hold off; end
    nexttile; scatter3(Ps(:,1),Ps(:,2),Ps(:,3), 2, s.Y(:,2), 'filled');
    axis equal; view([1 -1.5 1]); colormap(gca, turbo); title('y_2 (Fiedler)');
    nexttile; plotLabels(Ps, s.labels);
    title(sprintf('분할 K_{pred}=%d, ARI=%.2f', s.Kpred, s.ariGT));
    nexttile; plot(1:15, full.lambda(1:15), 'ko-'); hold on;
    if r > 1, plot(1:15, s.lambda(1:15), 'rs-'); legend('full','occ','Location','northwest'); end
    hold off; grid on; xlabel('k'); ylabel('\lambda_k'); title('spectrum');
end
exportgraphics(f1, fullfile(outDir,'overview.png'), 'Resolution', 110);

f2 = figure('Visible','off', 'Position',[50 50 1200 400]);
tiledlayout(1,3, 'TileSpacing','compact');
nexttile; bar(vertcat(res.angK)); set(gca,'XTickLabel',{res.view}); ylabel('deg');
title(sprintf('principal angles (앞쪽 %d개)', Kgt)); grid on;
nexttile; bar([vertcat(res.corrSame), nan(nV,1), vertcat(res.corrBest)]);
set(gca,'XTickLabel',{res.view}); ylim([0 1]); grid on;
title('|corr| u_2..u_6: 같은번호 | best-match');
nexttile; bar([[res.capFull].', [res.capOcc].']); set(gca,'XTickLabel',{res.view});
ylim([0 1]); legend('full','occ','Location','south'); grid on; title('cluster capture');
exportgraphics(f2, fullfile(outDir,'subspace.png'), 'Resolution', 110);
fprintf('\n저장: %s\n', outDir);

%% ===================== local functions =====================
function out = analyzeCloud(Pc, idxOrig, gtc, gOpts, alpha, kEig)
% 그래프 → 스펙트럼 → eigengap + conductance 분할 → GT 평가
[W, gi] = buildGraph(Pc, gOpts);
nComp   = max(conncomp(graph(W,'upper')));
[L, Dis] = normalizeGraph(W, alpha);
[U, lam] = embedSpectral(L, kEig);
Y = Dis * U;

Kgap = selectEigenGap(lam, struct('kMin', max(2,nComp), 'kMax', 15));
skip = nComp;                                   % run_pipeline_ISS C7 규칙
dim  = max(Kgap - skip, 2);
sOpts = struct('ky',40, 'cutMethod','conductance', 'qThr',0.90, 'phiMax',0.05, ...
               'minSize',60, 'maxCuts', max(Kgap-1,1));
labels = segmentEmbedding(Y(:, skip+1:skip+dim), W, sOpts);

g    = gtc(gi.idxKeep);
gRef = refineGT(W, g);
mA = evalSegmentation(labels, g);
mB = evalSegmentation(labels, gRef);
out = struct('idx', idxOrig(gi.idxKeep), 'W', W, 'Y', Y, 'lambda', lam, 'nComp', nComp, ...
             'Kgap', Kgap, 'labels', labels, 'Kpred', double(max(labels)), ...
             'ariGT', mA.ARI, 'ariGTref', mB.ARI, 'KgtRef', double(max(gRef)));
end

function g2 = refineGT(W, g)
% GT 라벨을 가시 그래프의 연결성분으로 다시 쪼갬 (가림으로 끊긴 파트는 별개 군집으로 인정)
g2 = zeros(size(g), 'uint32'); off = 0;
for l = unique(g).'
    ii = find(g == l);
    c  = conncomp(graph(W(ii,ii)));
    g2(ii) = uint32(off + c(:));
    off = off + max(c);
end
end

function ang = principalAngles(A, B)
Qa = orth(A - 0); Qb = orth(B);
s  = svd(Qa.' * Qb);
ang = acosd(min(s(1:min(size(Qa,2), size(Qb,2))), 1));
end

function cap = clusterCapture(Yk, g)
% 정규화한 GT 지시벡터들이 span(Yk)에 담기는 에너지 비율 (1 = 완전히 piecewise-constant)
Q = orth(Yk);
labs = unique(g);
H = zeros(numel(g), numel(labs));
for j = 1:numel(labs), h = double(g == labs(j)); H(:,j) = h / norm(h); end
cap = norm(Q.'*H, 'fro')^2 / norm(H, 'fro')^2;
end

function vis = hprVisible(P, C, gexp)
% Hidden Point Removal (Katz, Tal, Basri 2007): 구면 반전 후 볼록껍질 꼭짓점 = 가시 점
p = P - C;
n = vecnorm(p, 2, 2);
R = max(n) * 10^gexp;
f = p + 2*(R - n) .* p ./ n;
K = convhulln([f; 0 0 0]);
v = unique(K(:)); v(v > size(P,1)) = [];
vis = false(size(P,1), 1); vis(v) = true;
end

function Q = sampleCylinder(N, r, h, caps)
% 로컬 z축 원통 표면 균일 샘플링 (면적 가중)
aSide = 2*pi*r*2*h; aCap = caps*2*pi*r^2;
isCap = rand(N,1) < aCap/(aSide + aCap);
th = 2*pi*rand(N,1);
rr = r*ones(N,1); rr(isCap) = r*sqrt(rand(nnz(isCap),1));
z  = h*(2*rand(N,1) - 1);
z(isCap) = h*sign(rand(nnz(isCap),1) - 0.5);
Q = [rr.*cos(th), rr.*sin(th), z];
end

function plotLabels(Ps, lab)
lab = double(lab); K = max(lab); cm = lines(max(K,1));
hold on;
m0 = lab == 0;
if any(m0), scatter3(Ps(m0,1),Ps(m0,2),Ps(m0,3), 1, [0.7 0.7 0.7], 'filled'); end
for k = 1:K
    m = lab == k;
    if any(m), scatter3(Ps(m,1),Ps(m,2),Ps(m,3), 2, cm(k,:), 'filled'); end
end
hold off; axis equal; view([1 -1.5 1]); grid on;
end

function plotCam(C)
d = -C / norm(C) * 3; c0 = -d*2.2;
quiver3(c0(1),c0(2),c0(3), d(1),d(2),d(3), 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 1);
end
