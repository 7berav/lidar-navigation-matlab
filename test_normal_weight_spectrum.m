% 법선 가중치가 스펙트럼 모드 순서에 주는 영향 테스트
%
% 씬: test_occlusion_spectral 과 같은 모듈 A - connector - 모듈 B + 패널 (가림 없음)
%     + 단독 아웃라이어 40개 (원통 옆면/패널 윗면에서 법선 방향으로 0.15 튐)
%     + 아웃라이어 뭉치 1개 (A 옆면 위 0.12, 반경 0.04 공 안 20점)
%
% 설정 (mutual kNN k=30, self-tuning gamma=1 은 공통)
%   S0 거리만 (현재 파이프라인)
%   S1 + 법선 각도 가중치 exp(-(1-|n_i.n_j|)/tauN), tauN = 0.3 (약하게)
%   S2 + 법선 각도 가중치 tauN = 0.1 (강하게)
%   S3 S2 + 접평면 거리 가중치 exp(-(d_perp/sigPerp)^2), sigPerp = 0.03
%   S4 S3 과 같은 가중치, alpha = 0 (나머지는 alpha = 1)
%   법선은 kNN(20) PCA, 부호 없음 (|n_i.n_j| 사용)
%
% 모드 분류 (앞쪽 고유벡터마다 무엇을 자르는 모드인지)
%   part   : 파트 지시벡터(A, connector, B, 패널) 공간 capture > 0.8
%   face   : 면 지시벡터(옆면, 뚜껑, 패널 면) 공간 capture > 0.8  → 꺾인 곳(crease) 분리
%   out    : 단독 아웃라이어에 에너지 30% 이상
%   clump  : 아웃라이어 뭉치에 에너지 30% 이상
%   smooth : 해당 없음 (길이 방향 파동 같은 LB 형 모드)
%
% 출력: out_cutting_normal/<yyyyMMdd_HHmmss>/ modes.png, results.mat
rng(1);

rho = 100; noise = 0.01; kEig = 16; kNorm = 20;
partNames = {'A','연결부','B','패널'};

%% 1) 씬 + 아웃라이어
[P0, part, face, nTrue] = makeScene(rho);
nS = size(P0,1);
cand = find(ismember(face, [1 7 10]));
src  = cand(randperm(numel(cand), 40));
Pout = P0(src,:) + 0.15*nTrue(src,:);
v    = randn(20,3); v = v ./ vecnorm(v,2,2) .* (0.04*rand(20,1).^(1/3));
Pcl  = [-4 0 1.12] + v;
P    = [P0 + noise*randn(nS,3); Pout; Pcl];
kind = [zeros(nS,1); ones(40,1); 2*ones(20,1)];        % 0 표면, 1 단독 아웃라이어, 2 뭉치
part = [part; part(src); ones(20,1)];
face = [face; face(src); ones(20,1)];
fprintf('점 %d개 (표면 %d, 단독 아웃라이어 40, 뭉치 20)\n', size(P,1), nS);

%% 2) 법선 추정 품질
nEst = pcaNormals(P, kNorm);
angS = acosd(min(1, abs(sum(nEst(1:nS,:) .* nTrue, 2))));
angO = acosd(min(1, abs(sum(nEst(kind==1,:) .* nTrue(src,:), 2))));
creaseBand = false(nS,1);   % 뚜껑 테두리/접합부 근처 (진짜 법선과 추정 법선이 크게 다른 곳)
creaseBand(angS > 20) = true;
fprintf('추정 법선 오차 [deg] 중앙값: 표면 %.1f, 단독 아웃라이어 %.1f (원래 표면 법선 기준)\n', ...
    median(angS), median(angO));
fprintf('표면 점 중 법선 오차 20도 초과 (꺾인 곳 근처로 번진 띠): %.1f%%\n', 100*mean(creaseBand));

% 간선별 이웃 법선 각도: 꺾인 곳을 가로지르는 간선 vs 매끈한 곡면 안의 간선
W0 = graphN(P, nEst, 30, Inf, Inf);
[ei, ej] = find(triu(W0, 1));
ss = kind(ei) == 0 & kind(ej) == 0;
ei = ei(ss); ej = ej(ss);
angE = acosd(min(1, abs(sum(nEst(ei,:) .* nEst(ej,:), 2))));
fi = face(ei); fj = face(ej);
thinPanel = fi >= 10 & fj >= 10;                   % 패널 윗면/아랫면/테두리끼리는 두께 0.04 라 제외
crossE = fi ~= fj & ~thinPanel;
fprintf('이웃 법선 각도 중앙값 [deg]: 꺾인 곳 가로지르는 간선 %.1f (90도 꺾임), A 옆면 안 %.1f, 연결부 옆면 안 %.1f\n', ...
    median(angE(crossE)), median(angE(fi==1 & fj==1)), median(angE(fi==4 & fj==4)));

%% 3) 설정별 스펙트럼 + 모드 분류
% split = true: alpha 밀도 보정 q 를 법선 항이 없는 거리 그래프에서 구하고, 법선 항은 그 뒤에 곱함
cfg = struct('name',  {'S0 거리만', 'S1 법선 약 tau=0.3', 'S2 법선 강 tau=0.1', ...
                       'S3 법선 강 + 접평면', 'S4 S3, alpha=0', 'S5 S3, 밀도는 거리그래프에서'}, ...
             'tauN',  {Inf, 0.3, 0.1, 0.1, 0.1, 0.1}, ...
             'sigP',  {Inf, Inf, Inf, 0.03, 0.03, 0.03}, ...
             'alpha', {1, 1, 1, 1, 0, 1}, ...
             'split', {false, false, false, false, false, true});
Wdist = graphN(P, nEst, 30, Inf, Inf);
surf = kind == 0;
Qp = orth(indic(part(surf)));
Qf = orth(indic(face(surf)));
res = struct([]);
for c = 1:numel(cfg)
    W = graphN(P, nEst, 30, cfg(c).tauN, cfg(c).sigP);
    nComp = max(conncomp(graph(W, 'upper')));
    if cfg(c).split
        [L, Dis] = normalizeSplit(W, Wdist, cfg(c).alpha);
    else
        [L, Dis] = normalizeGraph(W, cfg(c).alpha);
    end
    [U, lam] = embedSpectral(L, kEig);
    Y = Dis * U;
    fprintf('\n=== %s  (연결성분 %d) ===\n', cfg(c).name, nComp);
    fprintf(['   i    lambda      분류     파트capt 면capt  원주비  단독E   뭉치E   ' ...
             '에너지 최대 파트   파트별 평균 [A c B P]\n']);
    labs = cell(kEig,1);
    tubeSide = surf & ismember(face, [1 4 7]);
    for i = 2:kEig
        y  = Y(:,i);
        e  = U(:,i).^2 / sum(U(:,i).^2);     % D-가중 에너지 (u = D^{1/2} y). 저차수 점의 수치오차 증폭 배제
        eO = sum(e(kind==1)); eC = sum(e(kind==2));
        ys = y(surf) - mean(y(surf)); ys = ys / norm(ys);
        cP = norm(Qp.'*ys)^2; cF = norm(Qf.'*ys)^2;
        cr = circRatio(y, P(:,1), face, tubeSide);
        ep = accumarray(part(surf), e(surf), [4 1]).';
        [~, pBest] = max(ep);
        mp = accumarray(part(surf), y(surf), [4 1], @mean).'; mp = mp / max(abs(mp));
        if     eC > 0.3, lb = 'clump';
        elseif eO > 0.3, lb = 'out';
        elseif cP > 0.8, lb = 'part';
        elseif cF > 0.8, lb = 'face';
        elseif cr > 0.5, lb = 'width';
        else,            lb = 'smooth';
        end
        labs{i} = lb;
        fprintf('  %2d  %10.3e   %-7s  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f   %-4s %4.0f%%   [%+.1f %+.1f %+.1f %+.1f]\n', ...
            i, lam(i), lb, cP, cF, cr, eO, eC, partNames{pBest}, 100*ep(pBest), mp);
    end
    res(c).name = cfg(c).name; res(c).lambda = lam; res(c).labels = labs; %#ok<SAGROW>
    res(c).Y = Y; res(c).nComp = nComp;                                  %#ok<SAGROW>
end

%% 4) 그림: 설정별 앞쪽 6개 비상수 모드
outDir = fullfile('out_cutting_normal', char(datetime('now','Format','yyyyMMdd_HHmmss')));
if ~exist(outDir,'dir'), mkdir(outDir); end
show = [1 4 5 6];
f = figure('Visible','off', 'Position',[40 40 1900 300*numel(show)]);
tl = tiledlayout(numel(show), 6, 'TileSpacing','compact', 'Padding','compact');
title(tl, '법선 가중치별 앞쪽 모드 (큰 검은 테두리 점 = 아웃라이어)');
for rr = 1:numel(show)
    R = res(show(rr));
    for i = 2:7
        nexttile; y = R.Y(:,i); cl = max(abs(y(surf)));
        scatter3(P(surf,1), P(surf,2), P(surf,3), 1.5, y(surf), 'filled'); hold on;
        o = ~surf;
        scatter3(P(o,1), P(o,2), P(o,3), 14, y(o), 'filled', 'MarkerEdgeColor', 'k');
        hold off; axis equal off; view([0.35 -1 0.7]); caxis([-cl cl]); colormap(gca, bwr());
        title(sprintf('%s | u_{%d} %s  \\lambda=%.1e', R.name, i, R.labels{i}, R.lambda(i)), 'FontSize', 7);
    end
end
exportgraphics(f, fullfile(outDir,'modes.png'), 'Resolution', 110);
res = rmfield(res, 'Y');
save(fullfile(outDir,'results.mat'), 'res', 'cfg', 'angS', 'angO');
fprintf('\n저장: %s\n', outDir);

%% ===================== local functions =====================
function [P, part, face, nT] = makeScene(rho)
% face: 1 A옆, 2 A +x뚜껑(가까운), 3 A -x뚜껑(먼), 4 연결부옆, 7 B옆, 8 B +x뚜껑(먼), 9 B -x뚜껑(가까운)
%       10 패널윗면, 11 패널아랫면, 12 패널테두리
Rx = [0 0 1; 0 1 0; -1 0 0];
cyl = struct('r',{1.0,0.4,0.7}, 'h',{3.0,0.5,2.5}, 't',{[-3.5 0 0],[0 0 0],[3 0 0]}, ...
             'caps',{true,false,true});
P = zeros(0,3); part = zeros(0,1); face = zeros(0,1); nT = zeros(0,3);
for i = 1:3
    c = cyl(i);
    aSide = 2*pi*c.r*2*c.h; aCap = c.caps*2*pi*c.r^2;
    N = round(rho*(aSide + aCap));
    isCap = rand(N,1) < aCap/(aSide + aCap);
    th = 2*pi*rand(N,1);
    rr = c.r*ones(N,1); rr(isCap) = c.r*sqrt(rand(nnz(isCap),1));
    z  = c.h*(2*rand(N,1) - 1);
    sg = sign(rand(nnz(isCap),1) - 0.5); sg(sg==0) = 1; z(isCap) = c.h*sg;
    nl = [cos(th), sin(th), zeros(N,1)]; nl(isCap,:) = [zeros(nnz(isCap),2), sg];
    f  = 3*(i-1) + ones(N,1); f(isCap) = 3*(i-1) + 2 + (z(isCap) < 0);
    P  = [P; [rr.*cos(th), rr.*sin(th), z]*Rx.' + c.t]; %#ok<AGROW>
    nT = [nT; nl*Rx.']; part = [part; i*ones(N,1)]; face = [face; f]; %#ok<AGROW>
end
a = 1.2; b = 2.5; t = 0.02; ctr = [3.0, 0.7 + b, 0];
areas = [4*a*b, 4*a*b, 4*b*t, 4*b*t, 4*a*t, 4*a*t];
N  = round(rho*sum(areas));
fc = discretize(rand(N,1), [0 cumsum(areas)/sum(areas)]);
Q  = [(2*rand(N,1)-1)*a, (2*rand(N,1)-1)*b, (2*rand(N,1)-1)*t];
Q(fc==1,3) = t; Q(fc==2,3) = -t; Q(fc==3,1) = a; Q(fc==4,1) = -a; Q(fc==5,2) = b; Q(fc==6,2) = -b;
nl = zeros(N,3); nl(fc==1,3) = 1; nl(fc==2,3) = -1; nl(fc==3,1) = 1; nl(fc==4,1) = -1;
nl(fc==5,2) = 1; nl(fc==6,2) = -1;
f = 10*ones(N,1); f(fc==2) = 11; f(fc>=3) = 12;
P = [P; Q + ctr]; nT = [nT; nl]; part = [part; 4*ones(N,1)]; face = [face; f];
end

function n = pcaNormals(P, k)
idx = knnsearch(P, P, 'K', k);
N = size(P,1); n = zeros(N,3);
for i = 1:N
    X = P(idx(i,:),:); X = X - mean(X,1);
    [V, D] = eig(X.'*X); [~, m] = min(diag(D)); n(i,:) = V(:,m).';
end
end

function W = graphN(P, nEst, k, tauN, sigP)
% buildGraph 와 같은 mutual kNN + self-tuning (gamma=1, tauIdx=3) 에 법선 항을 곱함
N = size(P,1);
[idx, dist] = knnsearch(P, P, 'K', k+1);
idx = idx(:,2:end); dist = dist(:,2:end);
rows = repmat((1:N)', k, 1); cols = idx(:); di = dist(:);
A = sparse(rows, cols, 1, N, N); M = A & A.';
mask = M(sub2ind([N N], rows, cols));
rows = rows(mask); cols = cols(mask); di = di(mask);
tau = dist(:,3);
w = exp(-di.^2 ./ (tau(rows).*tau(cols)));
if isfinite(tauN)
    cs = abs(sum(nEst(rows,:).*nEst(cols,:), 2));
    w  = w .* exp(-(1 - cs)/tauN);
end
if isfinite(sigP)
    dv = P(cols,:) - P(rows,:);
    dp = max(abs(sum(nEst(rows,:).*dv, 2)), abs(sum(nEst(cols,:).*dv, 2)));
    w  = w .* exp(-(dp/sigP).^2);
end
W = sparse(rows, cols, w, N, N);
W = max(W, W.');
end

function [L, Dis] = normalizeSplit(W, Wq, alpha)
% normalizeGraph 와 같되, 밀도 q 는 법선 항 없는 거리 그래프 Wq 에서 구함.
% 법선 항이 줄인 간선을 "점이 드문 곳"으로 오인해 다시 키우는 것을 막기 위함.
N  = size(W,1);
q  = full(sum(Wq,2));
Dq = spdiags(q.^(-alpha), 0, N, N);
W  = Dq * W * Dq;
d  = full(sum(W,2));
Dis = spdiags(1./sqrt(d + eps), 0, N, N);
L  = speye(N) - Dis * W * Dis;
L  = (L + L.')/2;
end

function cr = circRatio(y, x, face, sel)
% 원통 옆면에서 축 방향 구간(40개) 안의 분산 / 전체 분산. 0 = 길이 방향 모드, 1 = 원주 방향 모드
num = 0; den = 0;
for f = [1 4 7]
    m = sel & face == f;
    if nnz(m) < 50, continue; end
    xv = x(m); yv = y(m);
    b = discretize(xv, linspace(min(xv), max(xv), 41));
    for j = 1:40
        s = yv(b == j);
        if numel(s) > 3, num = num + sum((s - mean(s)).^2); end
    end
    den = den + sum((yv - mean(y(sel))).^2);
end
cr = num / max(den, eps);
end

function H = indic(g)
u = unique(g); H = double(g(:) == u(:).');
end

function cm = bwr(n)
if nargin < 1, n = 256; end
x = linspace(-1, 1, n).';
cm = [min(1, 1+x), 1-abs(x), min(1, 1-x)];
end
