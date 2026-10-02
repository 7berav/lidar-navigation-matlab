% ISS 컷 방식 비교 — 가중치는 C8a 로 고정 (α=1, skip=연결성분 수, 법선 σ=0.2, 임베딩 15차원)
%
% 모두 같은 임베딩 MST 위에서 "어떻게 자르고 언제 멈추는가"만 바꾼다.
%
% 컷 규칙
%   K1 length79      : 가장 먼 간선부터 (80 - 성분수)개                  ← C3 규칙
%   K1n              : K1 + 버려진 점(라벨 0)을 W로 이웃 군집에 붙임 (coverage 1 로 공정 비교)
%   K2 phi_len       : φ 작은 컷부터, 후보 = 긴 간선 상위 10%             ← C8a (현재 최적)
%   K3 phi_all       : φ 작은 컷부터, 후보 = 크기 조건만
%   K4 phi_junction  : φ 작은 컷부터, 후보 = 유의미한 분기점(가지 3개 이상)에 붙은 간선
%   K5 hdb           : HDBSCAN — 계층 전체에서 안정도가 최대인 군집 조합을 한 번에 선택
%   K6 hdb_mr15      : 〃, 상호도달거리(k=15)로 만든 MST (사슬 억제)
%   K7 hdb_mr15_m120 : 〃, 최소 군집 크기 120
%
% 정지 기준 (컷 규칙은 K2 와 같게 두고 멈추는 곳만 바꿈)
%   K2  : φ > √(2λ₂) 이면 정지 (Cheeger 상한 규모)
%   K8  : φ 순서로 끝까지 자른 뒤 log φ 가 가장 크게 뛰는 곳에서 정지 (φ-gap)
%   K9  : eigengap K* = argmax(λ_{k+1}-λ_k) 개로 정지
%   K10 : log-eigengap K* = argmax log(λ_{k+1}/λ_k) 개로 정지
%
% 평가는 GLB 부품 라벨(out_cutting_ISS/gt/gt_labels.mat)로만 하며, 방법 안에는 쓰지 않는다.
% 출력: out_cutting_ISS/cutcmp_<yyyyMMdd_HHmmss>/  summary.csv, labels.mat, overview.png, gap.png

%% 1) 데이터 + 그래프 + 임베딩 (C8a)
P_use = prepareISS();
gOpts = struct('k',30, 'gamma',1.0, 'tauIdx',3, 'minCompSize',150, 'normalSigma',0.2);
[W, gInfo] = buildGraph(P_use, gOpts);
Pg    = P_use(gInfo.idxKeep, :);
nComp = max(conncomp(graph(W, 'upper')));
kEig  = 45;
[L, Dis] = normalizeGraph(W, 1);
[U, lam] = embedSpectral(L, kEig);
Y    = Dis * U;
Yuse = Y(:, nComp+1 : nComp+15);
phiC = sqrt(2*lam(2));
fprintf('N=%d, 성분=%d, lambda2=%.3e, sqrt(2*lambda2)=%.4f\n', size(W,1), nComp, lam(2), phiC);

%% 2) 정지 기준 후보 계산
base = struct('ky',40, 'minSize',60);

% φ-gap: 끝까지 자르며 φ 수열을 얻고, log φ 가 가장 크게 뛰는 곳에서 멈춤
oFull = base; oFull.cutMethod = 'phiorder'; oFull.phiMax = 1; oFull.lenQuantile = 0.90; oFull.maxCuts = 200;
[~, dgFull] = segmentEmbedding(Yuse, W, oFull);
phiSeq = dgFull.cutEdges(:,4);
[~, tg] = max(diff(log(phiSeq)));
phiGap  = sqrt(phiSeq(tg) * phiSeq(tg+1));
fprintf('phi 수열 %d개, 가장 큰 log 점프: %d번째 컷 뒤 (%.4f -> %.4f)\n', ...
    numel(phiSeq), tg, phiSeq(tg), phiSeq(tg+1));

% eigengap / log-eigengap
Kgap = selectEigenGap(lam, struct('kMin',2, 'kMax',40));
kk   = (nComp+1):40;
[~, ig] = max(log(lam(kk+1) ./ lam(kk)));
KlogGap = kk(ig);
fprintf('eigengap K*=%d, log-eigengap K*=%d\n', Kgap, KlogGap);

%% 3) 설정
mk = @(varargin) mergeOpts(base, struct(varargin{:}));
cfg = {
    'K1_length79',        mk('cutMethod','length', 'qThr',0, 'maxCuts',max(80 - embedComponents(Yuse,40), 0))
    'K1n_length79_assign',mk('cutMethod','length', 'qThr',0, 'maxCuts',max(80 - embedComponents(Yuse,40), 0), 'assignNoise',true)
    'K2_phi_len',         mk('cutMethod','phiorder', 'phiMax',phiC, 'candidates','length', 'lenQuantile',0.90)
    'K3_phi_all',         mk('cutMethod','phiorder', 'phiMax',phiC, 'candidates','all')
    'K4_phi_junction',    mk('cutMethod','phiorder', 'phiMax',phiC, 'candidates','junction')
    'K5_hdb',             mk('cutMethod','hdbscan')
    'K6_hdb_mr15',        mk('cutMethod','hdbscan', 'mreachK',15)
    'K7_hdb_mr15_m120',   mk('cutMethod','hdbscan', 'mreachK',15, 'minSize',120)
    'K8_phi_gap',         mk('cutMethod','phiorder', 'phiMax',phiGap, 'lenQuantile',0.90)
    'K9_eigengap',        mk('cutMethod','phiorder', 'phiMax',1, 'lenQuantile',0.90, 'maxCuts',Kgap-1)
    'K10_logeigengap',    mk('cutMethod','phiorder', 'phiMax',1, 'lenQuantile',0.90, 'maxCuts',KlogGap-1)
};
nC = size(cfg,1);

%% 4) 실행 + 채점
GT = load(fullfile('out_cutting_ISS','gt','gt_labels.mat'));
assert(isequal(GT.Pg, Pg), 'GT 점 집합이 다릅니다');
[gu, ~, gGrp] = unique(GT.groupId);
gn     = GT.groupNames(gu);
gMajor = find(accumarray(gGrp,1) >= 60);
iH = find(contains(gn,'Harmony')); iU = find(contains(gn,'Unity'));

N = size(W,1);
labels = zeros(N, nC, 'uint32');
rows = struct([]); dgs = cell(nC,1);
for c = 1:nC
    t = tic;
    [lab, dg] = segmentEmbedding(Yuse, W, cfg{c,2});
    tc = toc(t);
    labels(:,c) = lab; dgs{c} = dg;
    covRaw = mean(lab > 0);
    if isfield(dg, 'coverageRaw'), covRaw = dg.coverageRaw; end
    s = scoreAll(double(lab), gGrp, numel(gu), gMajor, iH, iU, W);
    r = struct('name', cfg{c,1}, 'K', double(max(lab)), 'coverage', mean(lab>0), ...
               'coverageRaw', covRaw, 'nCuts', size(dg.cutEdges,1), ...
               'purity', s.purity, 'fragAvg', s.frag, 'ARIgt', s.ari, ...
               'Harmony', s.fH, 'Unity', s.fU, 'worstMixedN', s.worstN, ...
               'phiMed', s.phiMed, 'timeSec', tc);
    rows = [rows; r]; %#ok<AGROW>
    fprintf('%-17s K=%3d cov=%.3f(raw %.3f) cuts=%3d purity=%.3f frag=%.2f ARI=%.3f Har=%d Uni=%d worst=%4d (%.1fs)\n', ...
        r.name, r.K, r.coverage, r.coverageRaw, r.nCuts, r.purity, r.fragAvg, r.ARIgt, ...
        r.Harmony, r.Unity, r.worstMixedN, r.timeSec);
end
T = struct2table(rows);

%% 5) 저장 + 그림
outDir = fullfile('out_cutting_ISS', ['cutcmp_' char(datetime('now','Format','yyyyMMdd_HHmmss'))]);
if ~exist(outDir,'dir'), mkdir(outDir); end
writetable(T, fullfile(outDir,'summary.csv'));
names = cfg(:,1);
save(fullfile(outDir,'labels.mat'), 'T','names','labels','Pg','lam','phiSeq','phiC','phiGap','Kgap','KlogGap');

pal = [31 119 180; 255 127 14; 44 160 44; 214 39 40; 148 103 189; 140 86 75; 227 119 194; ...
       127 127 127; 188 189 34; 23 190 207; 174 199 232; 255 187 120; 152 223 138; 255 152 150; ...
       197 176 213; 196 156 148; 247 182 210; 199 199 199; 219 219 141; 158 218 229; ...
       0 0 128; 128 128 0; 0 128 128; 128 0 64; 64 64 64]/255;
f = figure('Position',[30 30 2300 950]);
tl = tiledlayout(f, 2, ceil(nC/2), 'TileSpacing','compact', 'Padding','compact');
for c = 1:nC
    drawSeg(nexttile(tl), Pg, double(labels(:,c)), pal, ...
        sprintf('%s\nK=%d ARI=%.3f 순도=%.3f', names{c}, T.K(c), T.ARIgt(c), T.purity(c)));
end
exportgraphics(f, fullfile(outDir,'overview.png'), 'Resolution', 110); close(f);

f = figure('Position',[30 30 1800 560]);
tl = tiledlayout(f, 1, 3, 'TileSpacing','compact', 'Padding','compact');
ax = nexttile(tl);
yyaxis(ax,'left');  plot(ax, 2:40, diff(lam(2:41)), 'o-', 'LineWidth',1.1); ylabel(ax,'\lambda_{k+1}-\lambda_k');
yyaxis(ax,'right'); plot(ax, kk, log(lam(kk+1)./lam(kk)), 's-', 'LineWidth',1.1); ylabel(ax,'log(\lambda_{k+1}/\lambda_k)');
xline(ax, Kgap, '--', sprintf('eigengap K*=%d', Kgap)); xline(ax, KlogGap, ':', sprintf('log K*=%d', KlogGap));
grid(ax,'on'); xlabel(ax,'k'); title(ax,'고유값 간격 (K를 직접 추정)');
ax = nexttile(tl);
semilogy(ax, 1:numel(phiSeq), phiSeq, 'o-', 'LineWidth',1.1); hold(ax,'on');
yline(ax, phiC, '--', '\surd(2\lambda_2)'); yline(ax, phiGap, ':', 'φ-gap');
hold(ax,'off'); grid(ax,'on'); xlabel(ax,'자른 순서'); ylabel(ax,'\phi (log)');
title(ax,'φ 순서로 끝까지 잘랐을 때의 컷 비용');
ax = nexttile(tl); hold(ax,'on');
for c = find(startsWith(names,'K5') | startsWith(names,'K6'))'
    [lg, kc] = aliveCurve(dgs{c}.tree);
    stairs(ax, lg, kc, 'LineWidth',1.2, 'DisplayName', names{c});
end
hold(ax,'off'); grid(ax,'on'); legend(ax,'Location','northwest','Interpreter','none');
xlabel(ax,'log \lambda  (\lambda = 1/임베딩 거리)'); ylabel(ax,'살아 있는 군집 수');
title(ax,'HDBSCAN 계층: 평평한 구간 = 오래 유지되는 K');
exportgraphics(f, fullfile(outDir,'gap.png'), 'Resolution', 120); close(f);
fprintf('\n저장 완료: %s\n', outDir);

%% 로컬 함수
function o = mergeOpts(a, b)
o = a; f = fieldnames(b);
for i = 1:numel(f), o.(f{i}) = b.(f{i}); end
end

function m = embedComponents(Y, ky)
N = size(Y,1);
[idxNN, distNN] = knnsearch(Y, Y, 'K', ky+1);
S = sparse(repmat((1:N)', ky, 1), reshape(idxNN(:,2:end),[],1), ...
           reshape(distNN(:,2:end),[],1), N, N);
m = max(conncomp(graph(max(S, S.'), 'upper')));
end

function s = scoreAll(lab, gGrp, nGrp, gMajor, iH, iU, W)
% GLB 부품 라벨 기준 채점 + 군집 conductance
v  = lab > 0;
Ct = accumarray([lab(v) gGrp(v)], 1, [max(lab) nGrp]);
s.purity = sum(max(Ct,[],2)) / sum(Ct(:));
fr = zeros(numel(gMajor),1);
for q = 1:numel(gMajor), fr(q) = coverCount(Ct(:,gMajor(q))); end
s.frag = mean(fr);
s.fH = coverCount(Ct(:,iH)); s.fU = coverCount(Ct(:,iU));
[~, kw] = max(sum(Ct,2) - max(Ct,[],2));   % 다른 부품이 가장 많이 섞인 군집
s.worstN = sum(Ct(kw,:));
e = evalSegmentation(uint32(lab), gGrp); s.ari = e.ARI;
d = full(sum(W,2)); volTot = sum(d); K = max(lab); phi = zeros(K,1);
for k = 1:K
    S = lab == k; volS = sum(d(S));
    phi(k) = (volS - full(sum(sum(W(S,S))))) / max(min(volS, volTot - volS), eps);
end
s.phiMed = median(phi);
end

function n = coverCount(col)
% 이 부품 점의 90%를 덮는 데 필요한 군집 수
col = sort(col, 'descend'); n = find(cumsum(col) >= 0.9*sum(col), 1);
end

function [lg, kc] = aliveCurve(tr)
% 압축 계층에서 lambda 에 따른 살아 있는 군집 수 (루트 제외)
b = tr.birth(2:end); d = tr.death(2:end);
lv = unique([b; d]); lv = lv(lv > 0 & isfinite(lv));
kc = arrayfun(@(x) nnz(b <= x & d > x), lv);
lg = log(lv);
end

function drawSeg(ax, P, lab, pal, ttl)
hold(ax,'on');
m0 = lab == 0;
if any(m0), scatter3(ax, P(m0,1), P(m0,2), P(m0,3), 3, [0.6 0.55 0.55], 'filled', 'MarkerFaceAlpha',0.5); end
for k = 1:max(lab)
    m = lab == k;
    scatter3(ax, P(m,1), P(m,2), P(m,3), 5, pal(mod(k-1,size(pal,1))+1,:), 'filled');
end
hold(ax,'off'); axis(ax,'equal','tight'); grid(ax,'on'); view(ax,[0 -1 0]);
xlabel(ax,'X'); zlabel(ax,'Z'); title(ax, ttl, 'Interpreter','none', 'FontSize', 9);
end
