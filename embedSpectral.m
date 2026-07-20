function [U, lambda] = embedSpectral(L, kEig, opts)
% 라플라시안의 최소 고유값 kEig개에 대한 고유벡터 계산
%
% 입력
%   L    : sparse 대칭 라플라시안 (normalizeGraph 출력)
%   kEig : 고유쌍 개수 (기본 45)
%   opts : .tol (기본 1e-4), .maxIter (기본 500)
%
% 출력
%   U      : N x kEig 고유벡터 (고유값 오름차순 정렬)
%   lambda : kEig x 1 고유값 (오름차순)

if nargin < 2 || isempty(kEig), kEig = 45; end
if nargin < 3, opts = struct(); end
tol     = getfield_def(opts, 'tol', 1e-4);
maxIter = getfield_def(opts, 'maxIter', 500);

[U, lamb] = eigs(L, kEig, 'smallestabs', ...
                 'Tolerance', tol, 'MaxIterations', maxIter);
lambda = diag(lamb);
[lambda, ord] = sort(lambda, 'ascend');
U = U(:, ord);
end

function v = getfield_def(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
