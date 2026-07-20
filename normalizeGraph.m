function [L, DinvSqrt, d] = normalizeGraph(W, alpha)
% Coifman–Lafon alpha-정규화 + 대칭 정규화 라플라시안
%
%   W_tilde = W ./ (q_i^alpha * q_j^alpha),  q_i = sum_j W_ij  (KDE)
%   L = I - D~^{-1/2} W~ D~^{-1/2}
%
%   alpha = 0   : 밀도 반영 그대로 (기존 ISS_normcut_4_modified와 동일)
%   alpha = 0.5 : Fokker–Planck 극한
%   alpha = 1   : Laplace–Beltrami 극한 (밀도 아티팩트 소거) ← 권장 기본
%
% 출력
%   L        : N x N sparse 대칭 라플라시안
%   DinvSqrt : D~^{-1/2} (고유벡터를 random-walk 좌표로 되돌릴 때 사용:
%              y = DinvSqrt * u)
%   d        : 정규화 후 차수벡터

if nargin < 2, alpha = 1; end
N = size(W,1);

if alpha > 0
    q = full(sum(W,2));
    Dq = spdiags(q.^(-alpha), 0, N, N);
    W = Dq * W * Dq;                     % alpha-정규화
end

d        = full(sum(W,2));
DinvSqrt = spdiags(1./sqrt(d + eps), 0, N, N);
L        = speye(N) - DinvSqrt * W * DinvSqrt;
L        = (L + L.')/2;                  % 수치 대칭화
end
