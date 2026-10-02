function [labels, diag_out] = segmentEmbedding(Yembed, W, opts)
% 고유벡터 임베딩 공간에서 MST 기반 분할
% (ISS_normcut_4_modified.m 294~426행의 함수화 + conductance 판정 추가)
%
% 입력
%   Yembed : N x m 임베딩 좌표 (예: DinvSqrt*U 의 열 일부)
%   W      : N x N 원 그래프 유사도 행렬 (conductance 판정용)
%   opts   : .ky        임베딩 kNN 이웃 수               (기본 40)
%            .cutMethod 'length' | 'conductance' | 'phiorder' | 'hdbscan' |
%                       'lensize' | 'zahn' | 'ward' | 'wardvol'
%                                                        (기본 'conductance')
%            .scoreMin  lensize/zahn/ward: 최고 점수가 이보다 작으면 정지 (기본 0)
%            .coreK     zahn: 국소 스케일로 쓸 이웃 순번          (기본 10)
%            .qThr      length: 컷 분위수 / conductance: 후보 분위수 (기본 0.96)
%            .phiMax    conductance/phiorder 채택 임계값 (기본 0.05)
%            .minSize   최소 군집 크기, 미달 라벨=0      (기본 60)
%            .maxCuts   최대 컷 수 (eigengap K 연동: K-1) (기본 Inf)
%            .lenQuantile  phiorder 전용: 후보 간선 길이 분위수 (기본 0.90,
%                          0이면 포레스트 간선 전체를 후보로)
%            .candidates   phiorder 전용 후보 규칙               (기본 'length')
%                          'length'   : lenQuantile 이상 긴 간선
%                          'junction' : 트리의 유의미한 분기점(가지 3개 이상이
%                                       모두 minSize 이상)에 붙은 간선
%                          'all'      : 크기 조건만
%            .lenPower     phiorder 전용: 길이와 conductance 를 함께 봄 (기본 0)
%                          phi <= phiMax 인 간선 중 phi/(길이/중앙길이)^lenPower
%                          가 최소인 것을 자른다. 0 이면 phi 만 본다.
%            .mreachK      임베딩 kNN 거리를 상호도달거리로 바꿈, 0=끔 (기본 0)
%                          d'(a,b) = max(core_a, core_b, d(a,b)), core = mreachK번째
%                          이웃 거리. 성긴 다리를 통한 사슬(chaining)을 억제.
%            .assignNoise  잡음점(라벨 0)을 W로 가장 강하게 이어진 군집에 붙임
%                          (기본: hdbscan 은 true, length/conductance 는 false)
%
% 출력
%   labels   : N x 1 uint32 (0 = 버려진 점)
%   diag_out : .cutEdges  실제 잘린 간선 [u v len phi]
%              .candEdges 후보 간선 [u v len]
%              .T_all     최종 포레스트 graph 객체
%
% cutMethod 이론 근거 (가이드 §4.5):
%   'length'      : MST 긴 간선 절단 (현행 baseline, 사슬에 취약)
%   'conductance' : 후보 간선 제거 시 원 그래프 W에서의 conductance
%                   phi = cut(S,S~)/min(vol S, vol S~) 를 계산해
%                   phi <= phiMax 인 컷만 채택 (Cheeger 근거화)
%   'phiorder'    : 길이 대신 phi 가 작은 컷부터 자름. 매 단계 포레스트의 모든
%                   간선에 대해 phi 를 다시 계산하고, 양쪽 조각이 minSize 이상인
%                   것 중 phi 최소를 선택. phi > phiMax 이면 중단.
%                   lenQuantile 로 후보를 긴 간선으로 제한하면(기본 0.90)
%                   가늘고 긴 부품을 가로로 자르는 경향이 줄어든다.
%   간선 점수 계열 — 양쪽 조각이 모두 minSize 이상인 MST 간선 중 점수가 가장 큰 것부터
%   자른다 (maxCuts 또는 scoreMin 으로 정지). 간선 e 를 빼면 그 성분이 A, B 로 갈린다.
%   'lensize'     : 점수 = 간선 길이. 'length' 와 달리 잔가지에 컷을 낭비하지 않는다.
%   'zahn'        : 점수 = 길이 / 양 끝점의 국소 스케일(coreK 번째 이웃 거리 평균).
%                   주변보다 유난히 긴 간선 (Zahn 1971 의 inconsistent edge).
%   'ward'        : 점수 = n_A n_B/(n_A+n_B) * ||mu_A - mu_B||^2.
%                   거리(두 조각 중심 사이)와 양쪽 점 개수를 함께 본다. 이 간선으로
%                   나눴을 때 줄어드는 군집 내 분산 = MST 로 제한한 2-means 분할.
%   'wardvol'     : 'ward' 에서 점 개수 대신 부피(차수 합)를, 평균 대신 차수 가중
%                   평균을 쓴다 (normalized cut 의 이완 문제와 같은 가중).
%   'hdbscan'     : 탐욕적으로 하나씩 자르지 않고, MST를 단일연결 계층으로 보고
%                   최소 크기 minSize 로 압축한 뒤 안정도(excess of mass)가 최대인
%                   군집 조합을 한 번에 고른다 (Campello et al. 2013). 군집 수는
%                   계층에서 오래 살아남는 군집으로 저절로 정해진다.

if nargin < 3, opts = struct(); end
ky      = getfield_def(opts, 'ky', 40);
method  = getfield_def(opts, 'cutMethod', 'conductance');
qThr    = getfield_def(opts, 'qThr', 0.96);
phiMax  = getfield_def(opts, 'phiMax', 0.05);
minSize = getfield_def(opts, 'minSize', 60);
maxCuts = getfield_def(opts, 'maxCuts', Inf);
lenQ    = getfield_def(opts, 'lenQuantile', 0.90);
cand    = getfield_def(opts, 'candidates', 'length');
mreachK = getfield_def(opts, 'mreachK', 0);
aNoise  = getfield_def(opts, 'assignNoise', []);   % 기본: hdbscan 만 켬 (기존 결과 보존)
scoreMin = getfield_def(opts, 'scoreMin', 0);
lenPow   = getfield_def(opts, 'lenPower', 0);
coreK    = getfield_def(opts, 'coreK', 10);

N = size(Yembed,1);

% --- 1) 임베딩 공간 kNN 그래프 → 성분별 MST (포레스트) ---
Mdl = createns(Yembed, 'NSMethod','kdtree', 'Distance','euclidean');
[idxNN, distNN] = knnsearch(Mdl, Yembed, 'K', ky+1);
idxNN  = idxNN(:,2:end);
distNN = distNN(:,2:end);

I = repmat((1:N)', ky, 1);
J = idxNN(:);
V = distNN(:);
if mreachK > 0                          % 상호도달거리 (HDBSCAN의 robust single linkage)
    core = distNN(:, min(mreachK, ky));
    V = max(V, max(core(I), core(J)));
end
S = sparse(I, J, V, N, N);
S = max(S, S.');
G_knn = graph(S, 'upper');

comp = conncomp(G_knn);
m = max(comp);

% 포레스트 간선을 전역 인덱스로 수집
E_global = zeros(0,2);
W_global = zeros(0,1);
for i = 1:m
    Vi = find(comp == i);
    Ti = minspantree(subgraph(G_knn, Vi), 'Method','sparse');
    E  = Ti.Edges.EndNodes;
    E_global = [E_global; [Vi(E(:,1))', Vi(E(:,2))']]; %#ok<AGROW>
    W_global = [W_global; Ti.Edges.Weight];            %#ok<AGROW>
end

if strcmpi(method, 'phiorder')
    [labels, diag_out] = cutPhiOrdered(E_global, W_global, W, N, ...
                                       phiMax, minSize, maxCuts, lenQ, cand, lenPow);
    return
end
if any(strcmpi(method, {'lensize','zahn','ward','wardvol'}))
    [labels, diag_out] = cutTreeScore(E_global, W_global, Yembed, W, N, lower(method), ...
                                      minSize, maxCuts, scoreMin, coreK);
    return
end
if strcmpi(method, 'hdbscan')
    if isempty(aNoise), aNoise = true; end
    [labels, diag_out] = cutHDBSCAN(E_global, W_global, W, N, minSize, aNoise);
    return
end

% --- 2) 컷 후보: 분위수 이상 긴 간선, 긴 것부터 ---
tauLen = quantile(W_global, qThr);
candMask = W_global >= tauLen;
candIdx  = find(candMask);
[~, ord] = sort(W_global(candIdx), 'descend');
candIdx  = candIdx(ord);
candEdges = [E_global(candIdx,:), W_global(candIdx)];

% 포레스트 인접 리스트 (컷 반영하며 갱신)
adjF = sparse(E_global(:,1), E_global(:,2), 1:size(E_global,1), N, N);
adjF = adjF + adjF.';
removed = false(size(E_global,1),1);

dW = full(sum(W,2));            % 원 그래프 차수 (conductance 볼륨용)
volTotal = sum(dW);

% conductance용 간선 리스트 전계산:
%   cut(S,S~) = vol(S) - 2*assoc(S,S) - sum(diag W in S)
% 항등식으로 계산하면 W(Smask,~Smask) 부분행렬 추출 없이 O(nnz) 논리연산만 필요
[euW, evW, ewW] = find(triu(W, 1));
dgW = full(diag(W));

cutEdges = zeros(0,4);
for c = 1:numel(candIdx)
    if size(cutEdges,1) >= maxCuts, break; end
    e = candIdx(c);
    u = E_global(e,1); v = E_global(e,2);

    switch lower(method)
        case 'length'
            phi = NaN;
            accept = true;
        case 'conductance'
            % 간선 e를 제거했을 때 u쪽 서브트리 S를 BFS로 수집
            Smask = subtreeMask(adjF, removed, u, e, N);
            volS  = sum(dW(Smask));
            volSc = volTotal - volS;
            assocS = sum(ewW(Smask(euW) & Smask(evW)));
            cutVal = volS - 2*assocS - sum(dgW(Smask));
            phi = cutVal / max(min(volS, volSc), eps);
            accept = (phi <= phiMax);
        otherwise
            error('segmentEmbedding: 알 수 없는 cutMethod "%s"', method);
    end

    if accept
        removed(e) = true;
        cutEdges = [cutEdges; u, v, W_global(e), phi]; %#ok<AGROW>
    end
end

% --- 3) 라벨: 컷 반영된 포레스트의 연결 성분 ---
Ekeep = E_global(~removed, :);
T_all = graph(Ekeep(:,1), Ekeep(:,2), W_global(~removed), N);
labels_all = uint32(conncomp(T_all))';

% --- 4) 최소 크기 필터 → 연속 라벨 재맵핑, 미달=0 ---
K_now = double(max(labels_all));
sz = accumarray(double(labels_all), 1, [K_now, 1]);
keep = find(sz >= minSize);
labels = zeros(N,1,'uint32');
if ~isempty(keep)
    remap = zeros(K_now,1,'uint32'); remap(keep) = 1:numel(keep);
    labels = remap(labels_all);
end
covRaw = mean(labels > 0);
if ~isempty(aNoise) && aNoise
    labels = uint32(assignNoiseByW(labels, W));
end

diag_out = struct('cutEdges', cutEdges, 'candEdges', candEdges, ...
                  'T_all', T_all, 'coverageRaw', covRaw);
end

function Smask = subtreeMask(adjF, removed, src, eSkip, N)
% 포레스트에서 간선 eSkip을 막은 채 src에서 도달 가능한 노드 집합
Smask = false(N,1);
stack = src;
Smask(src) = true;
while ~isempty(stack)
    u = stack(end); stack(end) = [];
    [nbrs, ~, eids] = find(adjF(:,u));   % 대칭이므로 열 인덱싱(CSC에서 빠름)
    for t = 1:numel(nbrs)
        v = nbrs(t); eid = eids(t);
        if eid == eSkip || removed(eid) || Smask(v), continue; end
        Smask(v) = true;
        stack(end+1) = v; %#ok<AGROW>
    end
end
end

function [labels, diag_out] = cutPhiOrdered(E, EL, W, N, phiMax, minSize, maxCuts, lenQ, cand, lenPow)
% phi 가 작은 컷부터 자르는 규칙 ('phiorder').
%
% 매 단계 포레스트의 모든 간선 (v, parent(v)) 에 대해
%   cut(S_v) = sum_{u in sub(v)} ( d_u - 2*A_u ),   A_u = LCA가 u인 W 간선 가중치 합
% 을 서브트리 합으로 한 번에 구한다. W 간선의 양 끝이 모두 sub(v) 에 있을 필요충분조건이
% "그 간선의 LCA가 sub(v) 안에 있다" 이므로 성립한다 (W 대칭, 대각 0 가정).
d = full(sum(W,2)); volTot = sum(d);
[wi, wj, ww] = find(triu(W,1));
if nargin < 9,  cand = 'length'; end
if nargin < 10, lenPow = 0; end
lenRef = median(EL);
lenMin = 0;
if strcmpi(cand, 'length') && lenQ > 0, lenMin = quantile(EL, lenQ); end
Lsp   = sparse([E(:,1); E(:,2)], [E(:,2); E(:,1)], [EL; EL], N, N);
alive = true(size(E,1),1);
Lg    = ceil(log2(N+2));
cutEdges = zeros(0,4);
while size(cutEdges,1) < maxCuts
    T = graph(E(alive,1), E(alive,2), [], N);
    [par, ord] = rootForest(T, N);
    p = par; p(p==0) = N+1;
    U = zeros(N+1, Lg+1); U(:,1) = [p; N+1];
    for k = 2:Lg+1, U(:,k) = U(U(:,k-1), k-1); end
    depth = zeros(N+1,1); rootOf = zeros(N,1);
    for k = 1:N
        v = ord(k);
        if par(v) > 0, depth(v) = depth(par(v)) + 1; rootOf(v) = rootOf(par(v));
        else,          depth(v) = 1;                 rootOf(v) = v;
        end
    end
    A  = accumarray(lcaVec(wi, wj, U, depth, Lg), ww, [N+1 1]);
    ac = d - 2*A(1:N); av = d; an = ones(N,1);
    for k = N:-1:1
        v = ord(k); pv = par(v);
        if pv > 0
            ac(pv) = ac(pv) + ac(v); av(pv) = av(pv) + av(v); an(pv) = an(pv) + an(v);
        end
    end
    valid = par > 0 & an >= minSize & (an(rootOf) - an) >= minSize;   % 유의미한 간선
    switch lower(cand)
        case 'length'
            if lenMin > 0
                iv = find(valid);
                valid(iv) = full(Lsp(sub2ind([N N], iv, par(iv)))) >= lenMin;
            end
        case 'junction'
            % 유의미한 간선만 남긴 트리에서 차수 3 이상인 노드 = 자연스러운 분기점
            sigDeg = accumarray(par(valid), 1, [N 1]) + double(valid);
            isJ = sigDeg >= 3;
            iv = find(valid);
            valid(iv) = isJ(iv) | isJ(par(iv));
        case 'all'
        otherwise
            error('segmentEmbedding: 알 수 없는 candidates "%s"', cand);
    end
    phi = ac ./ max(min(av, volTot - av), eps);
    phi(~valid) = inf;
    if lenPow > 0
        % 길이와 conductance 를 함께: phi <= phiMax 인 간선 중 phi/(길이/기준길이)^p 최소.
        % 같은 phi 라면 임베딩에서 더 멀리 떨어진 쪽을 먼저 자른다.
        ok = find(phi <= phiMax);
        if isempty(ok), break; end
        le = full(Lsp(sub2ind([N N], ok, par(ok))));
        [~, j] = min(phi(ok) ./ (le / lenRef).^lenPow);
        v = ok(j); pmin = phi(v);
    else
        [pmin, v] = min(phi);
        if ~isfinite(pmin) || pmin > phiMax, break; end
    end
    e = find(alive & ((E(:,1)==v & E(:,2)==par(v)) | (E(:,2)==v & E(:,1)==par(v))), 1);
    alive(e) = false;
    cutEdges = [cutEdges; E(e,1), E(e,2), EL(e), pmin]; %#ok<AGROW>
end
T_all  = graph(E(alive,1), E(alive,2), EL(alive), N);
labels = uint32(conncomp(T_all))';
K_now  = double(max(labels));
sz     = accumarray(double(labels), 1, [K_now 1]);
keep   = find(sz >= minSize);
lab    = zeros(N,1,'uint32');
if ~isempty(keep)
    remap = zeros(K_now,1,'uint32'); remap(keep) = 1:numel(keep);
    lab = remap(labels);
end
labels = lab;
diag_out = struct('cutEdges', cutEdges, 'candEdges', [E(EL >= lenMin, :), EL(EL >= lenMin)], ...
                  'T_all', T_all);
end

function [par, ord] = rootForest(T, N)
% 포레스트를 성분별로 루트를 잡아 BFS 순서와 부모 배열을 만든다 (루트의 부모 = 0)
comp = conncomp(T); par = zeros(N,1); ord = zeros(N,1); pos = 0;
[~, first] = unique(comp, 'stable');
for r = first(:)'
    ed = bfsearch(T, r, 'edgetonew');
    if ~isempty(ed), par(ed(:,2)) = ed(:,1); nodes = [r; ed(:,2)]; else, nodes = r; end
    ord(pos+1:pos+numel(nodes)) = nodes; pos = pos + numel(nodes);
end
end

function l = lcaVec(a, b, U, depth, Lg)
% binary lifting 으로 노드쌍의 LCA를 한 번에 계산 (다른 트리면 가상 루트 N+1)
sw = depth(a) < depth(b); t = a(sw); a(sw) = b(sw); b(sw) = t;
df = depth(a) - depth(b);
for k = 0:Lg
    m = bitand(df, 2^k) > 0; a(m) = U(a(m), k+1);
end
l = a; ne = a ~= b;
for k = Lg:-1:0
    m = ne & (U(a,k+1) ~= U(b,k+1));
    a(m) = U(a(m), k+1); b(m) = U(b(m), k+1);
end
l(ne) = U(a(ne), 1);
end

function [labels, diag_out] = cutTreeScore(E, EL, Y, W, N, mode, minSize, maxCuts, scoreMin, coreK)
% 간선 점수가 가장 큰 것부터 자르는 규칙 (lensize / zahn / ward / wardvol).
% 후보는 양쪽 조각이 모두 minSize 이상인 간선. cutEdges = [u v len score].
Lsp   = sparse([E(:,1); E(:,2)], [E(:,2); E(:,1)], [EL; EL], N, N);
alive = true(size(E,1),1);
isWard = any(strcmp(mode, {'ward','wardvol'}));
if isWard
    if strcmp(mode, 'wardvol'), wt = full(sum(W,2)); else, wt = ones(N,1); end
    YW = Y .* wt;
elseif strcmp(mode, 'zahn')
    [~, dk] = knnsearch(Y, Y, 'K', coreK+1);
    core = dk(:, end);
end
cutEdges = zeros(0,4);
while size(cutEdges,1) < maxCuts
    T = graph(E(alive,1), E(alive,2), [], N);
    [par, ord] = rootForest(T, N);
    rootOf = zeros(N,1);
    for k = 1:N
        v = ord(k);
        if par(v) > 0, rootOf(v) = rootOf(par(v)); else, rootOf(v) = v; end
    end
    an = ones(N,1);
    if isWard, sw = wt; sy = YW; end
    for k = N:-1:1
        v = ord(k); pv = par(v);
        if pv > 0
            an(pv) = an(pv) + an(v);
            if isWard, sw(pv) = sw(pv) + sw(v); sy(pv,:) = sy(pv,:) + sy(v,:); end
        end
    end
    iv = find(par > 0 & an >= minSize & (an(rootOf) - an) >= minSize);
    if isempty(iv), break; end
    len = full(Lsp(sub2ind([N N], iv, par(iv))));
    switch mode
        case 'lensize'
            sc = len;
        case 'zahn'
            sc = len ./ max(0.5*(core(iv) + core(par(iv))), eps);
        otherwise                                   % ward / wardvol
            rt  = rootOf(iv);
            nA  = sw(iv); nR = sw(rt); nB = nR - nA;
            muA = sy(iv,:) ./ nA;
            muB = (sy(rt,:) - sy(iv,:)) ./ nB;
            sc  = (nA .* nB ./ nR) .* sum((muA - muB).^2, 2);
    end
    [smax, j] = max(sc);
    if smax < scoreMin, break; end
    v = iv(j);
    e = find(alive & ((E(:,1)==v & E(:,2)==par(v)) | (E(:,2)==v & E(:,1)==par(v))), 1);
    alive(e) = false;
    cutEdges = [cutEdges; E(e,1), E(e,2), EL(e), smax]; %#ok<AGROW>
end
T_all  = graph(E(alive,1), E(alive,2), EL(alive), N);
labels = uint32(conncomp(T_all))';
K_now  = double(max(labels));
sz     = accumarray(double(labels), 1, [K_now 1]);
keep   = find(sz >= minSize);
lab    = zeros(N,1,'uint32');
if ~isempty(keep)
    remap = zeros(K_now,1,'uint32'); remap(keep) = 1:numel(keep);
    lab = remap(labels);
end
labels = lab;
diag_out = struct('cutEdges', cutEdges, 'candEdges', zeros(0,3), 'T_all', T_all);
end

function [labels, diag_out] = cutHDBSCAN(E, EL, W, N, m, assignNoise)
% MST 를 단일연결 계층으로 보고 최소 크기 m 으로 압축한 뒤,
% 안정도 합이 최대인 군집 조합을 고른다 (Campello, Moulavi, Sander 2013).
%   lambda = 1/거리. 군집 C 의 안정도 = sum_{p in C} (lambda_p - lambda_birth(C))
%   (lambda_p = 점 p 가 C 를 떠나는 수준). 자식 안정도 합보다 크면 C 를 고른다.

% 1) 덴드로그램: 짧은 간선부터 합침 (union-find)
[w, o] = sort(EL, 'ascend'); E = E(o,:);
M   = 2*N;
lft = zeros(M,1); rgt = zeros(M,1); dst = zeros(M,1); sz = zeros(M,1); sz(1:N) = 1;
uf  = (1:N)'; node = (1:N)'; nNode = N;
for t = 1:size(E,1)
    a = E(t,1); while uf(a) ~= a, uf(a) = uf(uf(a)); a = uf(a); end
    b = E(t,2); while uf(b) ~= b, uf(b) = uf(uf(b)); b = uf(b); end
    if a == b, continue; end
    nNode = nNode + 1;
    lft(nNode) = node(a); rgt(nNode) = node(b); dst(nNode) = w(t);
    sz(nNode)  = sz(node(a)) + sz(node(b));
    uf(b) = a; node(a) = nNode;
end
rt = zeros(N,1);
for x = 1:N, a = x; while uf(a) ~= a, a = uf(a); end, rt(x) = a; end
rs  = unique(rt);
top = node(rs(1));
for q = 2:numel(rs)                                   % 여러 성분이면 거리 Inf 로 묶음
    nNode = nNode + 1;
    lft(nNode) = top; rgt(nNode) = node(rs(q)); dst(nNode) = Inf;
    sz(nNode)  = sz(top) + sz(node(rs(q))); top = nNode;
end

% 2) 잎 순서: 각 노드의 잎이 leafOrder(lo:hi) 로 연속하게
leafOrder = zeros(N,1); lo = zeros(nNode,1); hi = zeros(nNode,1);
stack = top; pos = 0;
while ~isempty(stack)
    v = stack(end); stack(end) = [];
    if v <= N
        pos = pos + 1; leafOrder(pos) = v; lo(v) = pos; hi(v) = pos;
    else
        stack(end+1) = rgt(v); stack(end+1) = lft(v); %#ok<AGROW>
    end
end
for v = N+1:nNode, lo(v) = lo(lft(v)); hi(v) = hi(rgt(v)); end

% 3) 압축 계층: 두 자식이 모두 m 이상이면 군집 분기, 아니면 작은 쪽이 떨어져 나감
cPar = 0; cBirth = 0; cDeath = NaN; cStab = 0; nC = 1;     % 군집 1 = 루트
ptC = zeros(N,1); ptLam = zeros(N,1);
stack = [top 1];
while ~isempty(stack)
    v = stack(end,1); c = stack(end,2); stack(end,:) = [];
    if v <= N
        ptC(v) = c; ptLam(v) = cBirth(c); continue
    end
    l = lft(v); r = rgt(v); lam = 1 / max(dst(v), eps);
    if sz(l) >= m && sz(r) >= m
        cDeath(c) = lam;
        cStab(c)  = cStab(c) + sz(v) * (lam - cBirth(c));
        for ch = [l r]
            nC = nC + 1;
            cPar(nC) = c; cBirth(nC) = lam; cDeath(nC) = NaN; cStab(nC) = 0;
            stack(end+1,:) = [ch nC]; %#ok<AGROW>
        end
    else
        for ch = [l r]
            if sz(ch) >= m
                stack(end+1,:) = [ch c]; %#ok<AGROW>
            else
                idx = leafOrder(lo(ch):hi(ch));
                ptC(idx) = c; ptLam(idx) = lam;
                cStab(c) = cStab(c) + sz(ch) * (lam - cBirth(c));
            end
        end
        if sz(l) < m && sz(r) < m, cDeath(c) = lam; end
    end
end

% 4) 안정도 최대 조합 (excess of mass): 아래에서 위로, 루트는 제외
isLeafC = true(nC,1); isLeafC(cPar(cPar > 0)) = false;
sel = false(nC,1); best = zeros(nC,1); childSum = zeros(nC,1);
for c = nC:-1:2
    if isLeafC(c) || cStab(c) >= childSum(c)
        sel(c) = true; best(c) = cStab(c);
    else
        best(c) = childSum(c);
    end
    childSum(cPar(c)) = childSum(cPar(c)) + best(c);
end
hasSelAnc = false(nC,1);                                   % 선택된 군집의 자손은 해제
for c = 2:nC
    p = cPar(c);
    hasSelAnc(c) = p > 1 && (hasSelAnc(p) || sel(p));
    if hasSelAnc(c), sel(c) = false; end
end

% 5) 점 라벨 = 떨어져 나온 군집의 가장 가까운 선택 조상
selAnc = zeros(nC,1);
for c = 2:nC
    if sel(c), selAnc(c) = c; else, selAnc(c) = selAnc(cPar(c)); end
end
lab = selAnc(ptC);
if ~any(lab > 0), lab(:) = 1; end                         % 분기가 없으면 전체 하나
[~, ~, ic] = unique(lab(lab > 0));
labs = zeros(N,1); labs(lab > 0) = ic;
covRaw = mean(labs > 0);
if assignNoise, labs = assignNoiseByW(labs, W); end
labels = uint32(labs);

diag_out = struct('cutEdges', zeros(0,4), 'candEdges', zeros(0,3), 'T_all', [], ...
                  'coverageRaw', covRaw, 'nClusters', double(max(labs)), ...
                  'tree', struct('parent', cPar(:), 'birth', cBirth(:), 'death', cDeath(:), ...
                                 'stability', cStab(:), 'selected', sel));
end

function lab = assignNoiseByW(lab, W)
% 라벨 0 인 점을 W 가중치 합이 가장 큰 이웃 군집에 붙임 (바깥쪽부터 반복)
lab = double(lab); K = max(lab);
for it = 1:100
    un = find(lab == 0);
    if isempty(un), break; end
    [ii, jj, ww] = find(W(un, :));
    ok = lab(jj) > 0;
    if ~any(ok), break; end
    A = accumarray([ii(ok), lab(jj(ok))], ww(ok), [numel(un) K]);
    [mx, bestK] = max(A, [], 2);
    upd = mx > 0;
    if ~any(upd), break; end
    lab(un(upd)) = bestK(upd);
end
end

function v = getfield_def(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
