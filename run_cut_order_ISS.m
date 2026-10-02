% ISS: 컷 순서의 품질 vs 정지 기준의 품질 — 가중치는 C8a 로 고정
%
% (1) 순서: 간선 점수별로 끝까지(maxT 컷) 자르며 컷 수마다 채점한다.
%     → 어디서 멈추든 상관없이 "어떤 점수가 부품 경계를 먼저 찾는가"를 본다.
%        lensize : 간선 길이 (양쪽 >= minSize)
%        zahn    : 길이 / 국소 스케일
%        ward    : n_A n_B/(n_A+n_B) * ||mu_A - mu_B||^2   (거리 + 양쪽 점 개수)
%        wardvol : 위와 같되 부피(차수 합) 가중
%        phi_len : conductance 작은 순, 후보 = 긴 간선 상위 10%
%        phi_all : conductance 작은 순, 후보 = 전체
% (2) 정지: 각 점수 수열에서 log 점프가 1·2·3번째로 큰 곳에서 멈췄을 때의 점수.
%     phi 계열은 sqrt(2*lambda_2) 기준도 함께 표시.
% (3) 전제 확인: Fiedler 벡터가 이론처럼 두 값에 몰리고 0 에서 갈리는가.
%
% 평가는 GLB 부품 라벨로만 하며 방법 안에는 쓰지 않는다.
% 출력: out_cutting_ISS/cutorder_<yyyyMMdd_HHmmss>/  summary.csv, results.mat,
%       curves.png, scores.png, fiedler.png

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
% 국소화된 모드(점 몇 개짜리 조각에 몰린 고유벡터)는 임베딩에서 제외한다.
% 법선 가중은 일부 간선을 1e-10 수준까지 죽여 거의 끊어진 조각을 만들고, 그 조각이
% 작은 고유값을 차지한다 (걸러내지 않으면 15차원 중 7개, lambda_2 도 인공물).
minEff = 60;
[idx, Neff] = selectModes(U, nComp, 15, minEff);
Yuse = Y(:, idx);
phiC = sqrt(2*lam(idx(1)));                  % 첫 번째 정상 모드의 고유값 기준
phiAbs = 0.005;                              % 무차원 절대 기준
N    = size(W,1);
fprintf('임베딩 모드: %s (건너뛴 국소 모드 %d개), sqrt(2*lambda)=%.4f\n', ...
    mat2str(idx), idx(end) - nComp - numel(idx), phiC);

GT = load(fullfile('out_cutting_ISS','gt','gt_labels.mat'));
assert(isequal(GT.Pg, Pg), 'GT 점 집합이 다릅니다');
[gu, ~, gGrp] = unique(GT.groupId);
gMajor = find(accumarray(gGrp,1) >= 60);

%% 2) 순서별로 끝까지 자르고 컷 수마다 채점
maxT = 60;
base = struct('ky',40, 'minSize',60, 'maxCuts',maxT);
orders = {
    'lensize', struct('cutMethod','lensize'),                                 false
    'zahn',    struct('cutMethod','zahn'),                                    false
    'ward',    struct('cutMethod','ward'),                                    false
    'wardvol', struct('cutMethod','wardvol'),                                 false
    'phi_len', struct('cutMethod','phiorder','phiMax',1,'lenQuantile',0.90),  true
    'phi_all', struct('cutMethod','phiorder','phiMax',1,'candidates','all'),  true
};
nO = size(orders,1);
R = struct([]);
for o = 1:nO
    op = base; f = fieldnames(orders{o,2});
    for i = 1:numel(f), op.(f{i}) = orders{o,2}.(f{i}); end
    [~, dg] = segmentEmbedding(Yuse, W, op);
    ce = dg.cutEdges; Ekeep = dg.T_all.Edges.EndNodes; nT = size(ce,1);
    ari = zeros(nT+1,1); pur = ari; frg = ari; Kt = ari;
    labT = zeros(N, nT+1, 'uint16');
    for t = 0:nT
        Et  = [Ekeep; ce(t+1:end, 1:2)];
        lab = conncomp(graph(Et(:,1), Et(:,2), [], N))';
        [pur(t+1), frg(t+1), ari(t+1)] = scoreGT(lab, gGrp, numel(gu), gMajor);
        Kt(t+1) = max(lab); labT(:,t+1) = lab;
    end
    sc = ce(:,4);
    if orders{o,3}, jump = diff(log(sc)); else, jump = -diff(log(sc)); end   % t번째 컷 뒤의 점프
    [~, jo] = sort(jump, 'descend');
    R(o).name = orders{o,1}; R(o).isPhi = orders{o,3};
    R(o).score = sc; R(o).ari = ari; R(o).purity = pur; R(o).frag = frg; R(o).K = Kt;
    R(o).gapT = jo(1:3)';                         % 점프 1·2·3위: 이만큼 자르고 정지
    R(o).tCheeger = NaN; R(o).tAbs = NaN;
    if orders{o,3}                                % phi 가 처음 기준을 넘기 직전까지
        tc = find(sc > phiC, 1) - 1;   if isempty(tc), tc = nT; end
        ta = find(sc > phiAbs, 1) - 1; if isempty(ta), ta = nT; end
        R(o).tCheeger = tc; R(o).tAbs = ta;
    end
    R(o).labels = labT;
end

%% 3) 표
rows = struct([]);
fprintf('\n%-8s | %-16s | %-14s %-14s %-14s | %-14s %s\n', '순서', '최고 ARI (K)', 'gap 1위', 'gap 2위', 'gap 3위', 'sqrt(2*l)', 'phi<=0.005');
for o = 1:nO
    [bA, bi] = max(R(o).ari);
    g = R(o).gapT;
    r = struct('order', R(o).name, 'bestARI', bA, 'bestK', R(o).K(bi), ...
               'gap1_K', R(o).K(g(1)+1), 'gap1_ARI', R(o).ari(g(1)+1), ...
               'gap2_K', R(o).K(g(2)+1), 'gap2_ARI', R(o).ari(g(2)+1), ...
               'gap3_K', R(o).K(g(3)+1), 'gap3_ARI', R(o).ari(g(3)+1), ...
               'cheeger_K', NaN, 'cheeger_ARI', NaN, 'abs_K', NaN, 'abs_ARI', NaN, ...
               'ARI_K10', R(o).ari(min(10,end)), 'ARI_K14', R(o).ari(min(14,end)), ...
               'ARI_K20', R(o).ari(min(20,end)), 'ARI_K30', R(o).ari(min(30,end)));
    cs = '-';
    if R(o).isPhi
        tc = R(o).tCheeger; r.cheeger_K = R(o).K(tc+1); r.cheeger_ARI = R(o).ari(tc+1);
        ta = R(o).tAbs;     r.abs_K     = R(o).K(ta+1); r.abs_ARI     = R(o).ari(ta+1);
        cs = sprintf('%.3f (K=%2d)   %.3f (K=%2d)', r.cheeger_ARI, r.cheeger_K, r.abs_ARI, r.abs_K);
    end
    rows = [rows; r]; %#ok<AGROW>
    fprintf('%-8s | %.3f (K=%2d)     | %.3f (K=%2d)   %.3f (K=%2d)   %.3f (K=%2d)   | %s\n', ...
        r.order, r.bestARI, r.bestK, r.gap1_ARI, r.gap1_K, r.gap2_ARI, r.gap2_K, r.gap3_ARI, r.gap3_K, cs);
end
T = struct2table(rows);

%% 4) Fiedler 벡터가 두 값에 몰리는가
y2 = Y(:, idx(1));                           % 첫 번째 정상(비국소) 모드
d  = full(sum(W,2)); vTot = sum(d);
A  = y2 >= 0; vA = sum(d(A)); vB = vTot - vA;
idealA =  sqrt(vB / (vA * vTot));            % 이상적인 2분할이면 점들이 이 두 값에 몰린다
idealB = -sqrt(vA / (vB * vTot));
nearIdeal = mean(abs(y2 - idealA) < 0.25*abs(idealA) | abs(y2 - idealB) < 0.25*abs(idealB));
nearZero  = mean(abs(y2) < 0.25*min(abs(idealA), abs(idealB)));
fprintf('\nFiedler: 부호 분할 %d / %d 점, 이상값 %.4f / %.4f\n', nnz(A), nnz(~A), idealA, idealB);
fprintf('  이상값 ±25%% 안에 있는 점 %.1f%%, 0 근처(작은 이상값의 25%% 이내) 점 %.1f%%\n', ...
    100*nearIdeal, 100*nearZero);

%% 5) 저장 + 그림
outDir = fullfile('out_cutting_ISS', ['cutorder_' char(datetime('now','Format','yyyyMMdd_HHmmss'))]);
if ~exist(outDir,'dir'), mkdir(outDir); end
writetable(T, fullfile(outDir,'summary.csv'));
save(fullfile(outDir,'results.mat'), 'R','T','Pg','lam','phiC','y2','idealA','idealB','-v7.3');

f = figure('Position',[30 30 1800 560]);
tl = tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');
flds = {'ari','purity','frag'}; ylab = {'ARI','순도 (부품을 덜 합침 ↑)','조각 수 (부품을 덜 쪼갬 ↓)'};
for p = 1:3
    ax = nexttile(tl); hold(ax,'on');
    for o = 1:nO
        plot(ax, R(o).K, R(o).(flds{p}), '-', 'LineWidth',1.4, 'DisplayName', R(o).name);
    end
    hold(ax,'off'); grid(ax,'on'); xlabel(ax,'군집 수 K'); ylabel(ax, ylab{p});
    if p == 1, legend(ax,'Location','southeast','Interpreter','none'); title(ax,'순서별 점수 곡선 (어디서 멈추든)'); end
end
exportgraphics(f, fullfile(outDir,'curves.png'), 'Resolution', 120); close(f);

f = figure('Position',[30 30 1800 900]);
tl = tiledlayout(f,2,3,'TileSpacing','compact','Padding','compact');
for o = 1:nO
    ax = nexttile(tl);
    semilogy(ax, 1:numel(R(o).score), R(o).score, 'o-', 'LineWidth',1.1); hold(ax,'on');
    g = R(o).gapT;
    for q = 1:3, xline(ax, g(q)+0.5, '--', sprintf('%d위', q)); end
    if R(o).isPhi, yline(ax, phiC, ':', '\surd(2\lambda)'); yline(ax, phiAbs, '-.', '0.005'); end
    hold(ax,'off'); grid(ax,'on'); xlabel(ax,'자른 순서'); ylabel(ax,'간선 점수 (log)');
    title(ax, R(o).name, 'Interpreter','none');
end
exportgraphics(f, fullfile(outDir,'scores.png'), 'Resolution', 110); close(f);

f = figure('Position',[30 30 1800 560]);
tl = tiledlayout(f,1,3,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl);
histogram(ax, y2, 80, 'EdgeColor','none'); hold(ax,'on');
xline(ax, 0, 'k-', '0'); xline(ax, idealA, 'r--', '이상값'); xline(ax, idealB, 'r--', '이상값');
hold(ax,'off'); grid(ax,'on'); xlabel(ax,'Fiedler 좌표 y_2'); ylabel(ax,'점 개수');
title(ax,'Fiedler 값 분포: 두 값에 몰리는가');
ax = nexttile(tl);
plot(ax, sort(y2), '-', 'LineWidth',1.3); hold(ax,'on'); yline(ax, 0, 'k-');
yline(ax, idealA, 'r--'); yline(ax, idealB, 'r--'); hold(ax,'off'); grid(ax,'on');
xlabel(ax,'점 (y_2 오름차순)'); ylabel(ax,'y_2'); title(ax,'정렬한 Fiedler 값: 계단인가 비탈인가');
ax = nexttile(tl);
scatter3(ax, Pg(:,1), Pg(:,2), Pg(:,3), 6, y2, 'filled');
axis(ax,'equal','tight'); grid(ax,'on'); view(ax,[0 -1 0]); colorbar(ax);
xlabel(ax,'X'); zlabel(ax,'Z'); title(ax,'점군 위의 y_2');
exportgraphics(f, fullfile(outDir,'fiedler.png'), 'Resolution', 120); close(f);
fprintf('\n저장 완료: %s\n', outDir);

%% 로컬 함수
function [purity, fragAvg, ariGT] = scoreGT(lab, gGrp, nGrp, gMajor)
lab = double(lab);
Ct = accumarray([lab gGrp], 1, [max(lab) nGrp]);
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
