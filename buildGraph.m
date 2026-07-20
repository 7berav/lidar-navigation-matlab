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
%
% 출력
%   W    : N2 x N2 sparse 대칭 유사도 행렬 (성분 필터 후)
%   info : .idxKeep  원본 P에서 살아남은 행 인덱스 (N2 x 1)
%          .rows,.cols,.w  mutual kNN 간선 목록 (새 번호 기준)
%          .G        graph 객체

if nargin < 2, opts = struct(); end
k       = getfield_def(opts, 'k', 30);
gamma   = getfield_def(opts, 'gamma', 1.0);
tauIdx  = getfield_def(opts, 'tauIdx', 3);
minComp = getfield_def(opts, 'minCompSize', 0);

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
end

function v = getfield_def(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
