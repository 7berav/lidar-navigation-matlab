function [labels, info] = cutTreeGreedy(Yemb, W, K, opts)
% 임베딩 MST 의 간선을 "후보"로만 쓰고, 어느 간선을 자를지는 원 그래프 NCut 으로 정한다
% (segmentEmbedding.m 과 독립 구현)
%
%   1) 임베딩 공간 kNN 간선 + 원 그래프 W 의 간선을 합쳐 길이 = 임베딩 거리인 그래프를 만들고
%      최소신장 포레스트를 구한다. W 간선을 합치는 이유: 임베딩이 퇴화(같은 좌표에 몰림)해도
%      트리의 연결성이 W 의 연결성보다 나빠지지 않게 하기 위함.
%   2) K 개 조각이 될 때까지 트리 간선을 하나씩 제거한다.
%        score 'ncut'   : 제거했을 때 K-way NCut = sum_i cut(S_i)/vol(S_i) 의 증가가 가장 작은 간선
%                         (양쪽 부피가 모두 minVolFrac 이상인 것만)
%        score 'length' : 가장 긴 간선부터 (기존 MST 방식, 비교 기준)
%
% 입력
%   Yemb : N x m 임베딩 좌표
%   W    : N x N 점수 계산에 쓸 원 그래프 (대칭, 대각 0)
%   K    : 목표 조각 수
%   opts : .score      'ncut' | 'length'             (기본 'ncut')
%          .ky         임베딩 kNN 이웃 수             (기본 40)
%          .minVolFrac ncut: 조각 최소 부피 비율      (기본 0.01)
%          .minSize    length: 이보다 작은 조각은 이웃 조각에 붙임 (기본 60)
%
% 출력
%   labels : N x 1 uint32
%   info   : .nTree (포레스트 성분 수), .cuts [u v 길이 점수], .E (트리 간선), .alive (안 잘린 간선)

if nargin < 4, opts = struct(); end
score   = getfield_def(opts, 'score', 'ncut');
ky      = getfield_def(opts, 'ky', 40);
minVolF = getfield_def(opts, 'minVolFrac', 0.01);
minSize = getfield_def(opts, 'minSize', 60);

N  = size(W,1);
ky = min(ky, N-2);

% --- 1) 후보 트리 ---
[idxNN, distNN] = knnsearch(Yemb, Yemb, 'K', ky+1);
[wi, wj] = find(triu(W, 1));
I = [repmat((1:N)', ky, 1); wi];
J = [reshape(idxNN(:,2:end), [], 1); wj];
V = [reshape(distNN(:,2:end), [], 1); vecnorm(Yemb(wi,:) - Yemb(wj,:), 2, 2)] + 1e-15;
keep = I ~= J;
Sg = sparse(I(keep), J(keep), V(keep), N, N);
Sg = max(Sg, Sg.');
Tf = minspantree(graph(Sg, 'upper'), 'Type', 'forest');
E  = Tf.Edges.EndNodes; EL = Tf.Edges.Weight;
nTree = max(conncomp(Tf));
nCuts = max(K - nTree, 0);
cuts  = zeros(0, 4);
alive = true(size(E,1), 1);

% --- 2) 절단 ---
switch lower(score)
    case 'length'
        [~, o] = sort(EL, 'descend');
        e = o(1:min(nCuts, numel(o)));
        alive(e) = false;
        cuts = [E(e,:), EL(e), EL(e)];
    case 'ncut'
        d = full(sum(W,2)); volTot = sum(d); minVol = minVolF * volTot;
        [ei, ej, ew] = find(triu(W, 1));
        Eid = sparse([E(:,1); E(:,2)], [E(:,2); E(:,1)], [1:size(E,1), 1:size(E,1)].', N, N);
        Lg  = ceil(log2(N + 2));
        for c = 1:nCuts
            T = graph(E(alive,1), E(alive,2), [], N);
            [par, ord, rootOf, depth] = rootForest(T, N);
            p = par; p(p == 0) = N + 1;
            up = zeros(N+1, Lg+1); up(:,1) = [p; N+1];
            for k = 2:Lg+1, up(:,k) = up(up(:,k-1), k-1); end

            same = rootOf(ei) == rootOf(ej);
            A    = accumarray(lcaVec(ei(same), ej(same), up, depth, Lg), ew(same), [N+1 1]);
            dout = accumarray([ei(~same); ej(~same)], [ew(~same); ew(~same)], [N 1]);
            ac = d - 2*A(1:N); av = d; ao = dout;          % 서브트리 합: 컷 / 부피 / 다른 조각으로 나가는 가중치
            for k = N:-1:1
                v = ord(k); pv = par(v);
                if pv > 0
                    ac(pv) = ac(pv) + ac(v); av(pv) = av(pv) + av(v); ao(pv) = ao(pv) + ao(v);
                end
            end
            volS = av; cutS = ac;
            volC = av(rootOf); cutC = ao(rootOf);
            volB = volC - volS;
            cutB = cutC - ao + (cutS - ao);                % cut(C\S) = cut(C) - w(S,밖) + w(S,C\S)
            delta = cutS ./ max(volS, eps) + cutB ./ max(volB, eps) - cutC ./ max(volC, eps);
            delta(par == 0 | volS < minVol | volB < minVol) = Inf;
            [dmin, v] = min(delta);
            if ~isfinite(dmin), break; end
            e = full(Eid(v, par(v)));
            alive(e) = false;
            cuts(end+1,:) = [E(e,:), EL(e), dmin]; %#ok<AGROW>
        end
    otherwise
        error('cutTreeGreedy: 알 수 없는 score "%s"', score);
end

lab = conncomp(graph(E(alive,1), E(alive,2), [], N)).';
if strcmpi(score, 'length')                                % 잔가지 조각은 이웃에 붙임
    sz = accumarray(lab, 1);
    lab(sz(lab) < minSize) = 0;
    lab = attachZero(lab, W);
end
[~, ~, lab] = unique(lab);
labels = uint32(lab);
info = struct('nTree', nTree, 'cuts', cuts, 'E', E, 'alive', alive);   % E(alive,:) = 남은 트리 간선
end

function [par, ord, rootOf, depth] = rootForest(T, N)
% 성분별로 루트를 잡아 BFS 순서·부모·루트·깊이를 만든다 (루트의 부모 = 0, 깊이 = 1)
comp = conncomp(T); par = zeros(N,1); ord = zeros(N,1); rootOf = zeros(N,1); pos = 0;
[~, first] = unique(comp, 'stable');
for r = first(:)'
    ed = bfsearch(T, r, 'edgetonew');
    if ~isempty(ed), par(ed(:,2)) = ed(:,1); nodes = [r; ed(:,2)]; else, nodes = r; end
    ord(pos+1:pos+numel(nodes)) = nodes; pos = pos + numel(nodes);
    rootOf(nodes) = r;
end
depth = zeros(N+1, 1);
for k = 1:N
    v = ord(k);
    if par(v) > 0, depth(v) = depth(par(v)) + 1; else, depth(v) = 1; end
end
end

function l = lcaVec(a, b, up, depth, Lg)
% binary lifting 으로 노드쌍의 최소공통조상 (같은 트리 안의 쌍만 넣을 것)
a = a(:); b = b(:);
sw = depth(a) < depth(b); t = a(sw); a(sw) = b(sw); b(sw) = t;
df = depth(a) - depth(b);
for k = 0:Lg
    m = bitand(df, 2^k) > 0; a(m) = up(a(m), k+1);
end
l = a; ne = a ~= b;
for k = Lg:-1:0
    m = ne & (up(a, k+1) ~= up(b, k+1));
    a(m) = up(a(m), k+1); b(m) = up(b(m), k+1);
end
l(ne) = up(a(ne), 1);
end

function lab = attachZero(lab, W)
lab = double(lab); [u, ~, ic] = unique(lab);
if u(1) == 0, lab = ic - 1; else, lab = ic; end          % 0 은 0 으로 유지, 나머지는 1..K
K = max(lab);
if K == 0, lab(:) = 1; return; end
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
end

function v = getfield_def(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
