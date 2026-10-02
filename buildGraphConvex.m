function [W, info] = buildGraphConvex(P, S, opts)
% mutual kNN + self-tuning 그래프에 "오목 접합만 약화" 하는 법선 가중을 곱한다
% (buildGraph.m 과 같은 기본 그래프를 쓰되 독립 구현 — buildGraph 가 바뀌어도 영향 없음)
%
% 배경: 법선이 다르면 무조건 약화하는 방식(buildGraph 의 normalSigma, mode='unsigned')은
%   상자 모서리·원기둥 뚜껑 같은 볼록 모서리도 끊어서 한 부품을 면 단위로 쪼갠다.
%   부품은 대체로 볼록하고 부품 사이 접합은 오목하므로 (LCCP, Stein et al. 2014),
%   바깥 방향으로 정렬한 법선으로 간선을 볼록/오목으로 나눠 볼록 간선은 그대로 둔다.
%
% 입력
%   P    : N x 3 점군
%   S    : N x 3 각 점을 관측한 센서 위치 (법선을 센서 쪽 = 바깥쪽으로 정렬). [] 이면 'convex' 불가
%   opts : .k           최근접 이웃 수                         (기본 30)
%          .tauIdx      self-tuning 스케일로 쓸 이웃 순번      (기본 3)
%          .minCompSize 이보다 작은 연결성분 제거, 0=끔        (기본 0)
%          .mode        'none' | 'unsigned' | 'convex'         (기본 'convex')
%                       unsigned : 모든 간선에 exp(-(1-|n_i.n_j|)^2/sigma^2)   (= C8a)
%                       convex   : 위 감쇠를 주되, 볼록 간선은 면제 (양 끝이 모두 한쪽 면만
%                                  관측된 점일 때만). 얇은 판처럼 양면이 섞인 점은 방향이
%                                  의미 없으므로 각도만으로 감쇠한다.
%          .sigma       법선 감쇠 계수                          (기본 0.2)
%          .floor       감쇠 하한 (convex 기본 0.02, unsigned 기본 0)
%
% 출력
%   W    : sparse 대칭 유사도 (성분 필터 후)
%   info : .idxKeep, .normals (정렬된 법선), .twoSided (양면 점), .penFrac (감쇠된 간선 비율)

if nargin < 3, opts = struct(); end
k       = getfield_def(opts, 'k', 30);
tauIdx  = getfield_def(opts, 'tauIdx', 3);
minComp = getfield_def(opts, 'minCompSize', 0);
mode    = getfield_def(opts, 'mode', 'convex');
sigma   = getfield_def(opts, 'sigma', 0.2);
wFloor  = getfield_def(opts, 'floor', 0.02 * strcmpi(mode, 'convex'));

N = size(P,1);
[idx, dist] = knnsearch(P, P, 'K', k+1);
idx = idx(:,2:end); dist = dist(:,2:end);

rows = repmat((1:N)', k, 1); cols = idx(:); di = dist(:);
A = sparse(rows, cols, 1, N, N);
M = A & A.';
mask = full(M(sub2ind([N N], rows, cols)));
rows = rows(mask); cols = cols(mask); di = di(mask);

tau = dist(:, tauIdx);
w = exp(-(di.^2) ./ (tau(rows) .* tau(cols) + eps));

nrm = []; two = false(N,1); penFrac = 0;
if ~strcmpi(mode, 'none')
    nrm = pcaNormals(P, idx);
    if ~isempty(S)                                  % 센서 쪽(바깥)으로 정렬
        flip = sum(nrm .* (S - P), 2) < 0;
        nrm(flip,:) = -nrm(flip,:);
    end
    dotN = sum(nrm(rows,:) .* nrm(cols,:), 2);
    pen  = exp(-((1 - abs(dotN)).^2) ./ sigma^2);
    if strcmpi(mode, 'convex')
        assert(~isempty(S), 'buildGraphConvex: convex 모드는 센서 위치 S 가 필요합니다');
        % 양면 점: 이웃 중 반대 방향 법선이 20% 이상 (얇은 판의 앞뒤가 섞여 찍힌 경우)
        kk = min(15, k);
        opp = false(N, kk);
        for j = 1:kk, opp(:,j) = sum(nrm(idx(:,j),:) .* nrm, 2) < -0.7; end
        two = mean(opp, 2) >= 0.2;
        % 볼록: 간선을 따라 법선이 벌어진다  (n_j - n_i).(p_j - p_i) > 0
        kappa  = sum((nrm(cols,:) - nrm(rows,:)) .* (P(cols,:) - P(rows,:)), 2);
        exempt = kappa > 0 & ~two(rows) & ~two(cols);
        pen(exempt) = 1;
    end
    pen = max(pen, wFloor);
    penFrac = mean(pen < 0.5);
    w = w .* pen;
end

W = sparse(rows, cols, w, N, N);
W = max(W, W.');

idxKeep = (1:N)';
if minComp > 0
    comp = conncomp(graph(W, 'upper')).';
    sz   = accumarray(comp, 1);
    idxKeep = find(sz(comp) >= minComp);
    W = W(idxKeep, idxKeep);
end
info = struct('idxKeep', idxKeep, 'twoSided', two(idxKeep), 'penFrac', penFrac);
if ~isempty(nrm), info.normals = nrm(idxKeep,:); end
end

function nrm = pcaNormals(P, idx)
% 이웃 공분산의 최소 고유벡터 = 국소 법선 (부호 미정)
N = size(P,1); nrm = zeros(N,3);
for i = 1:N
    Q = P(idx(i,:), :); Q = Q - mean(Q,1);
    [V, D] = eig(Q.'*Q);
    [~, o] = min(diag(D));
    nrm(i,:) = V(:,o).';
end
end

function v = getfield_def(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
