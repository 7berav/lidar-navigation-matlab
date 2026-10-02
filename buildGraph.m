function [W, info] = buildGraph(P, opts)
% mutual kNN 그래프 구성 + 자기조율(self-tuning) 가중치
% (ISS_normcut_4_modified.m 62~134행 로직의 함수화)
%
% 입력
%   P    : N x 3 점군
%   opts : struct (생략 가능)
%          .k           최근접 이웃 수                  (기본 30)
%          .gamma       tau 지수 (1=순수 self-tuning)   (기본 1.0)
%          .tauIdx      sigma_i로 쓸 이웃 순번          (기본 3)
%          .minCompSize 이보다 작은 연결성분 제거, 0=끔 (기본 0)
%          .normalSigma 법선 일치도 가중 계수, 0=끔          (기본 0)
%          .normalK     법선 추정에 쓸 이웃 수               (기본 k)
%          .normalFloor 법선 감쇠의 하한 (0~1)               (기본 0)
%                       감쇠 계수를 max(계수, normalFloor) 로 제한한다. 하한이 없으면
%                       법선이 크게 다른 간선이 1e-10 수준까지 떨어져 사실상 끊어진
%                       작은 조각이 생기고, 그 조각에 국소화된 고유벡터가 작은
%                       고유값들을 차지해 임베딩 차원을 낭비한다.
%
% 법선 가중 (normalSigma > 0):
%   국소 PCA 법선 n_i (공분산의 최소 고유벡터) 를 구해
%       W_ij <- W_ij * exp( -(1 - |n_i' n_j|)^2 / normalSigma^2 )
%   곡률이 급하게 바뀌는 경계를 지나는 간선을 약화시킨다.
%   PCA 법선은 부호가 정해지지 않으므로 내적에 절댓값을 쓴다.
%
% 출력
%   W    : N2 x N2 sparse 대칭 유사도 행렬 (성분 필터 후)
%   info : .idxKeep  원본 P에서 살아남은 행 인덱스 (N2 x 1)
%          .rows,.cols,.w  mutual kNN 간선 목록 (새 번호 기준)
%          .G        graph 객체
%          .normals  N x 3 국소 PCA 법선 (normalSigma > 0 일 때만)

if nargin < 2, opts = struct(); end
k       = getfield_def(opts, 'k', 30);
gamma   = getfield_def(opts, 'gamma', 1.0);
tauIdx  = getfield_def(opts, 'tauIdx', 3);
minComp = getfield_def(opts, 'minCompSize', 0);
nSigma  = getfield_def(opts, 'normalSigma', 0);
nK      = getfield_def(opts, 'normalK', k);
nFloor  = getfield_def(opts, 'normalFloor', 0);

N = size(P,1);

% 1) kNN
Mdl = createns(P, 'NSMethod','kdtree', 'Distance','euclidean');
[idx, dist] = knnsearch(Mdl, P, 'K', k+1);
idx  = idx(:,2:end);                    % 자기 자신 제거
dist = dist(:,2:end);

rows_all = repmat((1:N)', k, 1);
cols_all = idx(:);
di_all   = dist(:);

% 2) mutual kNN 마스크
A = sparse(rows_all, cols_all, 1, N, N);
M = A & A.';                            % 상호 이웃만 유지
mask = M(sub2ind([N,N], rows_all, cols_all));
rows = rows_all(mask);
cols = cols_all(mask);
di   = di_all(mask);

% 3) self-tuning 가중치 (Zelnik-Manor–Perona 변형, gamma로 등방성 혼합)
sigma_i = median(dist(:));
tau = dist(:, tauIdx);
w_sim = exp( -(di.^2) ./ ( (tau(rows).^gamma) .* (tau(cols).^gamma) ...
                           .* (sigma_i.^(2*(1-gamma))) + eps ) );

% 3-2) 법선 일치도 가중 (옵션)
nrm = [];
if nSigma > 0
    nrm = pcaNormals(P, idx(:, 1:min(nK, size(idx,2))));
    dotN  = abs(sum(nrm(rows,:) .* nrm(cols,:), 2));
    w_sim = w_sim .* max(exp( -((1 - dotN).^2) ./ (nSigma^2) ), nFloor);
end

G = graph(rows, cols, w_sim, N);
G = simplify(G, 'min');
W = adjacency(G, 'weighted');
W = max(W, W.');                        % 수치적 비대칭 보정

% 4) 작은 연결성분 제거 (옵션)
idxKeep = (1:N)';
if minComp > 0
    comp = conncomp(G);
    compSize = accumarray(comp(:), 1);
    keepNodes = compSize(comp(:)) >= minComp;
    idxKeep = find(keepNodes);
    W = W(idxKeep, idxKeep);
    old2new = zeros(N,1); old2new(idxKeep) = 1:numel(idxKeep);
    edgeMask = keepNodes(rows) & keepNodes(cols);
    rows = old2new(rows(edgeMask));
    cols = old2new(cols(edgeMask));
    w_sim = w_sim(edgeMask);
    G = graph(W, 'upper');
end

info = struct('idxKeep', idxKeep, 'rows', rows, 'cols', cols, ...
              'w', w_sim, 'G', G);
if ~isempty(nrm), info.normals = nrm; end
end

function nrm = pcaNormals(P, idx)
% 이웃 공분산의 최소 고유벡터 = 국소 법선 (부호는 정해지지 않음)
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
