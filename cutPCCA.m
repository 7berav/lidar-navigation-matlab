function [labels, chi, vertIdx] = cutPCCA(Psi)
% PCCA+ (inner simplex 초기해): 고유벡터 임베딩에서 simplex 꼭짓점 K 개를 찾고
% 무게중심 좌표가 가장 큰 꼭짓점에 각 점을 할당한다 (Weber & Galliat 2002, Deuflhard & Weber 2005)
%
% 군집이 K 개로 잘 나뉘는 그래프에서 L_rw 의 앞쪽 고유벡터 K 개로 만든 점들은
% (K-1)-simplex 위에 놓이고 꼭짓점이 각 군집의 중심부에 대응한다. 파라미터가 없다.
%
% 입력
%   Psi : N x K, L_rw 의 앞쪽 K 개 고유벡터 (첫 열은 상수 벡터. psi = D^{-1/2} u)
%
% 출력
%   labels  : N x 1 uint32
%   chi     : N x K 소속도 (행 합 1, 꼭짓점 행은 단위벡터)
%   vertIdx : K x 1 꼭짓점으로 뽑힌 점 번호

[N, K] = size(Psi);
if K <= 1
    labels = ones(N, 1, 'uint32'); chi = ones(N, 1); vertIdx = 1; return
end

% inner simplex: 가장 먼 점부터 차례로 뽑고, 뽑힌 방향을 사영으로 제거
X = Psi;
vertIdx = zeros(K, 1);
[~, vertIdx(1)] = max(sum(X.^2, 2));
X = X - X(vertIdx(1), :);
for j = 2:K
    [~, vertIdx(j)] = max(sum(X.^2, 2));
    v = X(vertIdx(j), :) / norm(X(vertIdx(j), :));
    X = X - (X * v.') * v;
end

V = Psi(vertIdx, :);
if rcond(V) > 1e-12
    chi = Psi / V;
else
    chi = Psi * pinv(V);
end
[~, lab] = max(chi, [], 2);
[~, ~, lab] = unique(lab);                 % 빈 군집 제거 후 연속 번호
labels = uint32(lab);
end
