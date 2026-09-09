function [K, gaps, diag_out] = selectEigenGap(lambda, opts)
% eigengap 휴리스틱으로 군집 수 K 자동 선택
%
%   K* = argmax_{kMin <= k <= kMax} (lambda_{k+1} - lambda_k)
%
% 근거: 이상적으로 K개의 성분/군집이면 lambda_1..lambda_K ≈ 0 이고
% lambda_{K+1}부터 크게 뛴다 (perturbation 논거, von Luxburg 튜토리얼 §8.3).
%
% 입력
%   lambda : kEig x 1 고유값 (오름차순, embedSpectral 출력)
%   opts   : .kMin     탐색 하한                          (기본 2)
%            .kMax     탐색 상한                          (기본 numel(lambda)-1)
%            .relative true면 상대 갭 (λ_{k+1}-λ_k)/λ_{k+1} 사용 (기본 false)
%
% 출력
%   K        : 추정 군집 수 (임베딩 차원 = K, MST 목표 컷 수 = K-1로 쓰면 됨)
%   gaps     : (numel(lambda)-1) x 1 갭 벡터 (진단/플롯용)
%   diag_out : .kRange 탐색 구간, .gapAtK 선택된 갭 크기

if nargin < 2, opts = struct(); end
kMin = getfield_def(opts, 'kMin', 2);
kMax = getfield_def(opts, 'kMax', numel(lambda)-1);
rel  = getfield_def(opts, 'relative', false);

lambda = lambda(:);
assert(numel(lambda) >= 3, 'selectEigenGap: 고유값이 최소 3개 필요');
kMax = min(kMax, numel(lambda)-1);
kMin = max(2, min(kMin, kMax));

gaps = diff(lambda);
if rel
    gaps = gaps ./ max(lambda(2:end), eps);
end

[gapAtK, iRel] = max(gaps(kMin:kMax));
K = kMin + iRel - 1;

diag_out = struct('kRange', [kMin kMax], 'gapAtK', gapAtK);
end

function v = getfield_def(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
