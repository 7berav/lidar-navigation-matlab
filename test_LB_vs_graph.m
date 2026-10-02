% 그래프 라플라시안(alpha=1) 고유벡터 vs 라플라스-벨트라미(LB) 고유함수 직접 비교
%
% 비교 순서 (고유벡터는 부호와 축퇴 고유공간 안의 회전이 임의라서 값끼리 바로 비교할 수 없음)
%   1) 모드 짝찾기 : 그래프 고유벡터 y_i 가 각 LB 고유공간(축퇴 쌍이면 2차원)에 담기는 비율 = capture
%   2) 고유값      : 짝지어진 모드끼리 lambda_graph = c * mu_LB 의 스케일 c 하나만 맞추고 상대오차
%   3) 모드 모양   : y_i 와, y_i 를 짝지어진 LB 고유공간에 투영한 함수를 나란히 그림
%
% 케이스
%   T1 열린 원통 (모듈 A 크기 r=1, L=6)          LB 해석해 cos/sin(m th) cos(n pi s/L)
%   T2 T1 을 축에 수직 방향에서 가린 반원통 띠      LB 해석해 cos(k pi u/W) cos(n pi s/L), 노이만
%   T3 모듈 A - connector - 모듈 B 체인 (뚜껑 포함, 패널 제외, test_occlusion_spectral 과 같은 치수)
%      회전면이므로 LB 고유함수 = g(s) cos/sin(m th), g 는 자오선 1D 유한요소로 계산
%
% 출력: out_cutting_LB/<yyyyMMdd_HHmmss>/
%   summary.png     capture 행렬 + 고유값 비교 (T1, T2, T3)
%   T1_modes.png    원통 펼친 좌표 (r*th, s) 에서 그래프 모드 | LB 모드
%   T2_modes.png    가린 띠 펼친 좌표에서 그래프 모드 | LB 모드
%   T3_modes.png    체인 3D 에서 그래프 모드 | LB 모드
%   T3_profile.png  체인 m=0 모드를 자오선 길이 s 에 대해 그린 1D 비교
%   results.mat
rng(1);

%% 0) 파라미터 (test_occlusion_spectral 과 동일)
rho    = 100;     % 표면 점 밀도 [pts / unit^2]
noise  = 0.01;
alpha  = 1;
kEig   = 20;
hprExp = 2.5;
gOpts  = struct('k',30, 'gamma',1.0, 'tauIdx',3, 'minCompSize',50);
idxShow = 2:9;    % 그림에 보일 그래프 모드 번호 (u_1 은 상수)

outDir = fullfile('out_cutting_LB', char(datetime('now','Format','yyyyMMdd_HHmmss')));
if ~exist(outDir,'dir'), mkdir(outDir); end

%% T1) 열린 원통 전체
r = 1; L = 6;
N1 = round(rho*2*pi*r*L);
th = 2*pi*rand(N1,1);
s  = L*rand(N1,1);
P1 = [r*cos(th), r*sin(th), s - L/2] + noise*randn(N1,3);   % 축 = z

G1 = graphModes(P1, gOpts, alpha, kEig);
k1 = G1.idx;
[Phi, mu, lab, grp] = tubeModes(th(k1), s(k1), r, L);
T1 = struct('G', G1, 'M', matchModes(G1.Y, G1.lambda, Phi, mu, grp, lab), ...
            'u', r*th(k1), 's', s(k1));
printMatch('T1 열린 원통 (가림 없음)', T1.M);

%% T2) 축에 수직 방향에서 가린 반원통 띠
C   = [0 25 0];
vis = find(hprVisible(P1, C, hprExp));
G2  = graphModes(P1(vis,:), gOpts, alpha, kEig);
k2  = vis(G2.idx);
thc = mod(th(k2) - pi/2 + pi, 2*pi) - pi;  % 시선 방향(+y, th=pi/2) 기준으로 펼침 (2pi 경계 끊김 방지)
lo  = prctile(thc, 0.5);
hi  = prctile(thc, 99.5);
W   = r*(hi - lo);                        % 띠 너비 (원주 방향 호 길이)
u2  = r*(thc - lo);
[Phi, mu, lab, grp] = stripModes(u2, s(k2), W, L);
T2 = struct('G', G2, 'M', matchModes(G2.Y, G2.lambda, Phi, mu, grp, lab), ...
            'u', u2, 's', s(k2), 'W', W);
fprintf('\n[T2] 가시 %d / %d 점, 띠 너비 W = %.3f (반원주 pi*r 의 %.0f%%)\n', ...
    numel(k2), numel(k1), W, 100*W/(pi*r));
printMatch('T2 가린 반원통 띠', T2.M);

axialTab = axialCompare(T1.M, T2.M, L, 'self-tuning 폭 (파이프라인 그대로, gamma=1)');

% 대조: 고정 폭 가우시안 (gamma=0 -> 모든 점이 같은 폭). 이론상 가려도 축 방향 고유값이 같아야 함
gFix = gOpts; gFix.gamma = 0;
G1f = graphModes(P1, gFix, alpha, kEig);
[Phi, mu, lab, grp] = tubeModes(th(G1f.idx), s(G1f.idx), r, L);
M1f = matchModes(G1f.Y, G1f.lambda, Phi, mu, grp, lab);
G2f = graphModes(P1(vis,:), gFix, alpha, kEig);
k2f = vis(G2f.idx);
thf = mod(th(k2f) - pi/2 + pi, 2*pi) - pi;
[Phi, mu, lab, grp] = stripModes(r*(thf - lo), s(k2f), W, L);
M2f = matchModes(G2f.Y, G2f.lambda, Phi, mu, grp, lab);
axialTabFix = axialCompare(M1f, M2f, L, '고정 폭 (gamma=0)');
fprintf('스케일 c: self-tuning 전체 %.4g / 가림 %.4g  |  고정 폭 전체 %.4g / 가림 %.4g\n', ...
    T1.M.c, T2.M.c, M1f.c, M2f.c);

%% T3) A - connector - B 체인 (회전면)
S = chainSegments();
[P3c, seg3, par3, th3] = sampleChain(S, rho);
P3 = P3c + noise*randn(size(P3c));
G3 = graphModes(P3, gOpts, alpha, kEig);
k3 = G3.idx;
F  = chainFEM(S, 0.01, 0:3, 12);
[Phi, mu, lab, grp, info] = chainModes(F, seg3(k3), par3(k3), th3(k3));
T3 = struct('G', G3, 'M', matchModes(G3.Y, G3.lambda, Phi, mu, grp, lab), ...
            'P', P3(k3,:), 'seg', seg3(k3), 'par', par3(k3), 'info', info, 'F', F);
printMatch('T3 체인 A-connector-B', T3.M);
fprintf('\n스케일 c 비교 (같은 밀도, 같은 kNN 이면 비슷해야 함): T1 %.4g  T2 %.4g  T3 %.4g\n', ...
    T1.M.c, T2.M.c, T3.M.c);

%% 그림
drawT1 = @(v, cl) drawUnrolled(T1.u, T1.s, v, cl);
f = modeFigure(T1.G, T1.M, idxShow, drawT1, 'T1 열린 원통: 그래프 모드 | LB 모드 (펼친 좌표)');
exportgraphics(f, fullfile(outDir,'T1_modes.png'), 'Resolution', 110);

drawT2 = @(v, cl) drawUnrolled(T2.u, T2.s, v, cl);
f = modeFigure(T2.G, T2.M, idxShow, drawT2, 'T2 가린 반원통 띠: 그래프 모드 | LB 모드 (펼친 좌표)');
exportgraphics(f, fullfile(outDir,'T2_modes.png'), 'Resolution', 110);

drawT3 = @(v, cl) draw3(T3.P, v, cl);
f = modeFigure(T3.G, T3.M, idxShow, drawT3, 'T3 체인: 그래프 모드 | LB 모드');
exportgraphics(f, fullfile(outDir,'T3_modes.png'), 'Resolution', 110);

f = profileFigure(T3, S, F);
exportgraphics(f, fullfile(outDir,'T3_profile.png'), 'Resolution', 110);

f = figure('Visible','off', 'Position',[50 50 1650 900]);
tl = tiledlayout(2, 3, 'TileSpacing','compact', 'Padding','compact');
cases = {T1, T2, T3}; names = {'T1 열린 원통', 'T2 가린 띠', 'T3 체인'};
for q = 1:3
    nexttile(tl, q);   plotCapture(cases{q}.M, names{q});
    nexttile(tl, q+3); plotEigen(cases{q}.M, names{q});
end
exportgraphics(f, fullfile(outDir,'summary.png'), 'Resolution', 110);

T1s = rmfield(T1, 'G'); T2s = rmfield(T2, 'G'); T3s = rmfield(T3, 'G');
save(fullfile(outDir,'results.mat'), 'T1s', 'T2s', 'T3s', 'axialTab', 'axialTabFix', 'S', ...
     'rho', 'noise', 'alpha', 'kEig', 'gOpts');
fprintf('\n저장: %s\n', outDir);

%% ===================== local functions =====================
function G = graphModes(P, gOpts, alpha, kEig)
[W, gi]   = buildGraph(P, gOpts);
[Lm, Dis] = normalizeGraph(W, alpha);
[U, lam]  = embedSpectral(Lm, kEig);
G = struct('idx', gi.idxKeep, 'Y', Dis*U, 'lambda', lam, ...
           'nComp', max(conncomp(graph(W,'upper'))));
end

function [Phi, mu, lab, grp] = tubeModes(th, s, r, L)
% 열린 원통 옆면, 양 끝 노이만: cos/sin(m th) cos(n pi s / L), mu = (m/r)^2 + (n pi/L)^2
list = zeros(0,3);
for m = 0:6
    for n = 0:14, list(end+1,:) = [m n (m/r)^2 + (n*pi/L)^2]; end %#ok<AGROW>
end
list = sortrows(list, 3); list = list(1:25,:);
Phi = []; mu = []; grp = []; lab = cell(1, size(list,1));
for g = 1:size(list,1)
    m = list(g,1); n = list(g,2); ax = cos(n*pi*s/L);
    if m == 0, Fg = ax; else, Fg = [cos(m*th).*ax, sin(m*th).*ax]; end
    Phi = [Phi, Fg]; mu = [mu; repmat(list(g,3), size(Fg,2), 1)]; %#ok<AGROW>
    grp = [grp; repmat(g, size(Fg,2), 1)];                           %#ok<AGROW>
    lab{g} = sprintf('(m=%d,n=%d)', m, n);
end
end

function [Phi, mu, lab, grp] = stripModes(u, s, W, L)
% 직사각형 띠 [0,W] x [0,L], 네 변 노이만: cos(k pi u/W) cos(n pi s/L)
list = zeros(0,3);
for k = 0:8
    for n = 0:14, list(end+1,:) = [k n (k*pi/W)^2 + (n*pi/L)^2]; end %#ok<AGROW>
end
list = sortrows(list, 3); list = list(1:30,:);
Phi = zeros(numel(u), size(list,1)); mu = list(:,3); grp = (1:size(list,1)).';
lab = cell(1, size(list,1));
for g = 1:size(list,1)
    Phi(:,g) = cos(list(g,1)*pi*u/W) .* cos(list(g,2)*pi*s/L);
    lab{g} = sprintf('(k=%d,n=%d)', list(g,1), list(g,2));
end
end

function M = matchModes(Y, lam, Phi, mu, grp, lab)
% capture(i,g) = || Q_g' y_i ||^2 / ||y_i||^2,  Q_g = LB 고유공간 g 의 정규직교 기저
nY = size(Y,2); nGr = max(grp);
Q = cell(nGr,1); muG = zeros(nGr,1);
for g = 1:nGr
    Q{g} = orth(Phi(:, grp==g));
    muG(g) = mu(find(grp==g, 1));
end
cap = zeros(nY, nGr);
for i = 1:nY
    y = Y(:,i) / norm(Y(:,i));
    for g = 1:nGr, cap(i,g) = norm(Q{g}.'*y)^2; end
end
[best, gi] = max(cap, [], 2);
muI = muG(gi);
% 근접 고유값 묶음: 정렬한 mu 에서 이웃과 상대 간격이 tol 미만이면 같은 묶음.
% 이산화 오차가 이 간격보다 크면 그래프는 묶음 안의 모드를 섞으므로, 묶음 단위 capture 도 본다.
tol = 0.12;
[ms, o] = sort(muG); clu = zeros(nGr,1); cid = 1; clu(o(1)) = cid;
for t = 2:nGr
    if ms(t) < 1e-9 || (ms(t) - ms(t-1))/ms(t) > tol, cid = cid + 1; end
    clu(o(t)) = cid;
end
capClu = zeros(nY,1);
for i = 1:nY
    Qc = orth(Phi(:, ismember(grp, find(clu == clu(gi(i))))));
    capClu(i) = norm(Qc.'*(Y(:,i)/norm(Y(:,i))))^2;
end
ok  = best > 0.8 & muI > 1e-9;
c   = (lam(ok).'*muI(ok)) / (muI(ok).'*muI(ok));      % lambda = c*mu 최소제곱
M = struct('cap',cap, 'best',best, 'capClu',capClu, 'clu',clu, 'gi',gi, 'muG',muG, ...
           'muI',muI, 'ok',ok, 'c',c, 'lam',lam, 'labG',{lab}, 'labI',{lab(gi)}, 'Q',{Q});
end

function tab = axialCompare(M1, M2, L, name)
% 같은 축 방향 모드 (m=0,n) 와 (k=0,n) 의 그래프 고유값을 원값 그대로 비교
fprintf('\n=== 축 방향 모드, 전체 원통 vs 가린 띠 : %s ===\n', name);
fprintf('  n    mu_LB     lambda_full    lambda_half    half/full\n');
tab = nan(4,4);
for n = 1:4
    i1 = find(reshape(strcmp(M1.labI, sprintf('(m=0,n=%d)', n)), [], 1) & M1.ok, 1);
    i2 = find(reshape(strcmp(M2.labI, sprintf('(k=0,n=%d)', n)), [], 1) & M2.ok, 1);
    if isempty(i1) || isempty(i2), fprintf('  %d  (짝 없음)\n', n); continue; end
    tab(n,:) = [(n*pi/L)^2, M1.lam(i1), M2.lam(i2), M2.lam(i2)/M1.lam(i1)];
    fprintf('  %d  %7.4f   %11.4e    %11.4e    %6.3f\n', n, tab(n,:));
end
end

function printMatch(name, M)
fprintf('\n=== %s : 스케일 c = %.4g ===\n', name, M.c);
fprintf('   i   lambda_graph   lambda/c   LB 모드          mu_LB    상대오차  capture  근접묶음\n');
for i = 1:numel(M.lam)
    if M.muI(i) > 1e-9
        e = sprintf('%+7.1f%%', 100*(M.lam(i)/M.c - M.muI(i))/M.muI(i));
    else
        e = '      -';
    end
    fprintf('  %2d   %11.4e   %8.4f   %-15s %8.4f   %s   %6.3f   %6.3f\n', ...
        i, M.lam(i), M.lam(i)/M.c, M.labI{i}, M.muI(i), e, M.best(i), M.capClu(i));
end
end

function S = chainSegments()
% 자오선 단면 (x, rho) 의 선분들. cap 은 파라미터 = rho, side 는 파라미터 = x
% 뚜껑은 연결부 안쪽까지 꽉 찬 원판 (test_occlusion_spectral 의 샘플링과 동일) -> J1, J2 에서 T자 접합
mk = @(a,b,p0,p1,ty) struct('a',a, 'b',b, 'p0',p0, 'p1',p1, 'type',ty);
S = [ mk('cAf','rAf', [-6.5 0],   [-6.5 1],   'cap')    % 1 A 먼쪽 뚜껑
      mk('rAf','rAn', [-6.5 1],   [-0.5 1],   'side')   % 2 A 옆면
      mk('rAn','J1',  [-0.5 1],   [-0.5 0.4], 'cap')    % 3 A 가까운 뚜껑 바깥 고리
      mk('J1','cAn',  [-0.5 0.4], [-0.5 0],   'cap')    % 4 A 가까운 뚜껑 안쪽 원판 (가지)
      mk('J1','J2',   [-0.5 0.4], [ 0.5 0.4], 'side')   % 5 connector 옆면
      mk('J2','cBn',  [ 0.5 0.4], [ 0.5 0],   'cap')    % 6 B 가까운 뚜껑 안쪽 원판 (가지)
      mk('J2','rBn',  [ 0.5 0.4], [ 0.5 0.7], 'cap')    % 7 B 가까운 뚜껑 바깥 고리
      mk('rBn','rBf', [ 0.5 0.7], [ 5.5 0.7], 'side')   % 8 B 옆면
      mk('rBf','cBf', [ 5.5 0.7], [ 5.5 0],   'cap') ]; % 9 B 먼쪽 뚜껑
end

function [P, segId, par, th] = sampleChain(S, rho)
P = zeros(0,3); segId = zeros(0,1); par = zeros(0,1); th = zeros(0,1);
for q = 1:numel(S)
    p0 = S(q).p0; p1 = S(q).p1;
    if strcmp(S(q).type, 'cap')
        r0 = min(p0(2),p1(2)); r1 = max(p0(2),p1(2));
        n  = round(rho*pi*(r1^2 - r0^2));
        rr = sqrt(r0^2 + rand(n,1)*(r1^2 - r0^2));
        xx = p0(1)*ones(n,1); pp = rr;
    else
        x0 = min(p0(1),p1(1)); x1 = max(p0(1),p1(1)); R = p0(2);
        n  = round(rho*2*pi*R*(x1 - x0));
        xx = x0 + (x1 - x0)*rand(n,1); rr = R*ones(n,1); pp = xx;
    end
    t = 2*pi*rand(n,1);
    P = [P; xx, rr.*cos(t), rr.*sin(t)];   %#ok<AGROW>
    segId = [segId; q*ones(n,1)];          %#ok<AGROW>
    par = [par; pp]; th = [th; t];         %#ok<AGROW>
end
end

function F = chainFEM(S, h, mList, nPer)
% 회전면 LB: phi = g(s) e^{i m th}
%   약형식  int ( rho g'^2 + m^2/rho g^2 ) ds = mu int rho g^2 ds
% 자오선 트리 위 1D 선형 유한요소. 접합점(J1, J2)은 노드 공유로 자동 처리.
% m >= 1 이면 축 위 노드(rho = 0)에서 g = 0.
nS = numel(S);
endNames = [{S.a}, {S.b}];
endPos   = [vertcat(S.p0); vertcat(S.p1)];
[uNames, iu] = unique(endNames, 'stable');
X = endPos(iu, :);
E = zeros(0,2); segNodes = cell(nS,1); segPar = cell(nS,1);
for q = 1:nS
    p0 = S(q).p0; p1 = S(q).p1;
    nE = max(4, ceil(norm(p1 - p0)/h));
    t  = (1:nE-1).' / nE;
    base = size(X,1);
    X = [X; p0 + t.*(p1 - p0)]; %#ok<AGROW>
    nodes = [find(strcmp(uNames, S(q).a)); base + (1:nE-1).'; find(strcmp(uNames, S(q).b))];
    segNodes{q} = nodes;
    if strcmp(S(q).type, 'cap'), segPar{q} = X(nodes,2); else, segPar{q} = X(nodes,1); end
    E = [E; nodes(1:end-1), nodes(2:end)]; %#ok<AGROW>
end
nN = size(X,1);
axisNodes = find(X(:,2) < 1e-12);

gp = [0.5 - sqrt(15)/10, 0.5, 0.5 + sqrt(15)/10];  gw = [5 8 5]/18;
nEl = size(E,1);
F = struct('m', mList, 'X', X, 'E', E, 'segNodes', {segNodes}, 'segPar', {segPar}, ...
           'mu', {cell(1,numel(mList))}, 'g', {cell(1,numel(mList))});
for mi = 1:numel(mList)
    m = mList(mi);
    I = zeros(4*nEl,1); J = I; Kv = I; Mv = I;
    for e = 1:nEl
        a = E(e,1); b = E(e,2); qa = X(a,:); qb = X(b,:); ell = norm(qb - qa);
        Ke = zeros(2); Me = zeros(2);
        for g = 1:3
            t = gp(g); rq = (1-t)*qa(2) + t*qb(2);
            Nv = [1-t, t]; dN = [-1, 1]/ell;
            Ke = Ke + gw(g)*ell*( rq*(dN.'*dN) + m^2/max(rq,1e-12)*(Nv.'*Nv) );
            Me = Me + gw(g)*ell*rq*(Nv.'*Nv);
        end
        ii = 4*(e-1) + (1:4);
        I(ii) = [a; b; a; b]; J(ii) = [a; a; b; b];
        Kv(ii) = Ke(:); Mv(ii) = Me(:);
    end
    K = sparse(I, J, Kv, nN, nN); Mm = sparse(I, J, Mv, nN, nN);
    free = true(nN,1); if m > 0, free(axisNodes) = false; end
    [V, D] = eig(full(K(free,free)), full(Mm(free,free)));
    [d, o] = sort(real(diag(D))); V = real(V(:,o));
    gfull = zeros(nN, nPer); gfull(free,:) = V(:,1:nPer);
    F.mu{mi} = d(1:nPer); F.g{mi} = gfull;
end
end

function gv = evalProfile(F, mi, j, seg, par)
gv = zeros(numel(seg),1);
for q = 1:numel(F.segNodes)
    pts = (seg == q);
    if ~any(pts), continue; end
    xp = F.segPar{q}; yp = F.g{mi}(F.segNodes{q}, j);
    if xp(end) < xp(1), xp = flipud(xp); yp = flipud(yp); end
    gv(pts) = interp1(xp, yp, par(pts), 'linear', 'extrap');
end
end

function [Phi, mu, lab, grp, info] = chainModes(F, seg, par, th)
list = zeros(0,4);   % [m, mi, j, mu]
for mi = 1:numel(F.m)
    for j = 1:numel(F.mu{mi}), list(end+1,:) = [F.m(mi), mi, j, F.mu{mi}(j)]; end %#ok<AGROW>
end
list = sortrows(list, 4); list = list(1:25,:);
Phi = []; mu = []; grp = []; lab = cell(1, size(list,1));
for g = 1:size(list,1)
    m = list(g,1); gv = evalProfile(F, list(g,2), list(g,3), seg, par);
    if m == 0, Fg = gv; else, Fg = [gv.*cos(m*th), gv.*sin(m*th)]; end
    Phi = [Phi, Fg]; mu = [mu; repmat(list(g,4), size(Fg,2), 1)]; %#ok<AGROW>
    grp = [grp; repmat(g, size(Fg,2), 1)];                           %#ok<AGROW>
    jm = nnz(list(1:g,1) == m);          % 같은 m 안에서 몇 번째인지 (상수 포함)
    lab{g} = sprintf('m=%d #%d', m, jm);
end
info = list;
end

function vis = hprVisible(P, C, gexp)
% Hidden Point Removal (Katz, Tal, Basri 2007)
p = P - C; n = vecnorm(p, 2, 2);
R = max(n) * 10^gexp;
f = p + 2*(R - n) .* p ./ n;
K = convhulln([f; 0 0 0]);
v = unique(K(:)); v(v > size(P,1)) = [];
vis = false(size(P,1), 1); vis(v) = true;
end

%% ---------- 그림 ----------
function cm = bwr(n)
if nargin < 1, n = 256; end
x = linspace(-1, 1, n).';
cm = [min(1, 1+x), 1-abs(x), min(1, 1-x)];
end

function drawUnrolled(u, s, v, cl)
scatter(u, s, 3, v, 'filled'); axis equal tight; caxis([-cl cl]); colormap(gca, bwr());
xlabel('원주 방향 r\theta'); ylabel('축 방향 s'); set(gca, 'FontSize', 7);
end

function draw3(P, v, cl)
scatter3(P(:,1), P(:,2), P(:,3), 2, v, 'filled'); axis equal; view([0.35 -1 0.6]);
caxis([-cl cl]); colormap(gca, bwr()); set(gca, 'FontSize', 7);
end

function f = modeFigure(G, M, idxShow, drawFn, figTitle)
n = numel(idxShow); nc = 4; nr = ceil(2*n/nc);
f = figure('Visible','off', 'Position',[50 50 1500 280*nr]);
tl = tiledlayout(nr, nc, 'TileSpacing','compact', 'Padding','compact');
title(tl, figTitle);
for t = 1:n
    i = idxShow(t); y = G.Y(:,i); Q = M.Q{M.gi(i)}; yL = Q*(Q.'*y);
    cl = max(abs(y));
    nexttile; drawFn(y, cl);
    title(sprintf('그래프 u_{%d}   \\lambda/c = %.3f', i, M.lam(i)/M.c), 'FontSize', 8);
    nexttile; drawFn(yL, cl);
    title(sprintf('LB %s   \\mu = %.3f   capture %.2f (근접묶음 %.2f)', ...
          M.labI{i}, M.muI(i), M.best(i), M.capClu(i)), 'FontSize', 8);
end
end

function plotCapture(M, name)
nG = min(numel(M.labG), 16);
imagesc(M.cap(:, 1:nG)); caxis([0 1]); colormap(gca, flipud(gray)); colorbar;
set(gca, 'XTick', 1:nG, 'XTickLabel', M.labG(1:nG), 'FontSize', 7); xtickangle(60);
ylabel('그래프 고유벡터 번호 i'); title([name ' : capture (검정 = 1)']);
end

function plotEigen(M, name)
i2 = (1:numel(M.lam)).' >= 2;
ok = M.ok & i2; bad = ~M.ok & i2;
mx = 1.1*max([M.muI(ok); M.lam(ok)/M.c]);
plot([0 mx], [0 mx], 'k-'); hold on;
plot(M.muI(ok), M.lam(ok)/M.c, 'o', 'MarkerFaceColor', [0.2 0.4 0.8], 'MarkerEdgeColor', 'none');
if any(bad), plot(M.muI(bad), M.lam(bad)/M.c, 'x', 'Color', [0.85 0.3 0.2], 'LineWidth', 1.2); end
hold off; grid on; axis([0 mx 0 mx]); axis square;
xlabel('\mu_{LB}'); ylabel('\lambda_{graph} / c');
title(sprintf('%s : 고유값 (c = %.3g)', name, M.c));
end

function f = profileFigure(T3, S, F)
% m=0 모드: 그래프 고유벡터 값을 자오선 길이 s 에 대해 찍고, LB g(s) 를 선으로 겹침
% 주 경로 = 세그먼트 1 2 3 5 7 8 9 (안쪽 원판 4, 6 은 가지라서 제외)
path = [1 2 3 5 7 8 9];
cum = 0; s0 = nan(numel(S),1);
for q = path, s0(q) = cum; cum = cum + norm(S(q).p1 - S(q).p0); end
par0 = @(q) S(q).p0(1 + strcmp(S(q).type, 'cap'));   % cap -> rho, side -> x
sPt = nan(numel(T3.seg),1);
for q = path
    m = T3.seg == q; sPt(m) = s0(q) + abs(T3.par(m) - par0(q));
end
M = T3.M;
sel = find(M.ok & T3.info(M.gi,1) == 0 & (1:numel(M.lam)).' >= 2);
sel = sel(1:min(6, end));
f = figure('Visible','off', 'Position',[50 50 1500 260*ceil(numel(sel)/2)]);
tl = tiledlayout(ceil(numel(sel)/2), 2, 'TileSpacing','compact', 'Padding','compact');
title(tl, 'T3 체인 m=0 모드: 점 = 그래프 고유벡터, 선 = LB 고유함수 (자오선 길이 s)');
bnd = s0([2 3 5 7 8 9]);
for t = 1:numel(sel)
    i = sel(t); y = T3.G.Y(:,i);
    mi = T3.info(M.gi(i),2); j = T3.info(M.gi(i),3);
    phi = evalProfile(F, mi, j, T3.seg, T3.par);
    a = (phi.'*y) / (phi.'*phi);
    nexttile; hold on;
    keep = ~isnan(sPt);
    scatter(sPt(keep), y(keep), 3, [0.3 0.45 0.8], 'filled', 'MarkerFaceAlpha', 0.35);
    for q = path
        nd = F.segNodes{q}; sN = s0(q) + abs(F.segPar{q} - par0(q));
        plot(sN, a*F.g{mi}(nd, j), 'r-', 'LineWidth', 1.4);
    end
    for b = bnd.', xline(b, ':', 'Color', [0.5 0.5 0.5]); end
    hold off; grid on; xlim([0 cum]);
    title(sprintf('그래프 u_{%d} vs LB %s   capture %.2f', i, M.labI{i}, M.best(i)), 'FontSize', 8);
    xlabel('s  (A뚜껑 | A옆면 | A고리 | connector | B고리 | B옆면 | B뚜껑)', 'FontSize', 7);
end
end
