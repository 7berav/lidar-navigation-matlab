function [beta, error] = regressionFourthOrder(inputMatrix, symbterm, lambda)
    % Ridge regression (N-normalized):
    %   argmin  (1/N)||Phi*beta - 1||^2 + lambda*||beta||^2
    %   =>  (Phi'Phi + lambda*N*I) beta = Phi'*1
    %
    % lambda = 0 (default) => ordinary least squares (OLS)
    % lambda > 0           => ridge regression (L2 penalty on beta)
    %
    % N 스케일링: Phi'Phi ~ O(N) 이므로 lambda*N 을 곱해야
    % lambda 의 의미가 점 개수 N 에 무관하게 일정.
    %
    % 페널티는 beta 추정에만 적용하고 RANSAC 잔차/스코어 계산은 변경 없음.

    if nargin < 3 || isempty(lambda)
        lambda = 0;
    end

    Temp = calculateFourthOrder(inputMatrix, symbterm);
    N    = size(Temp, 1);
    nT   = size(Temp, 2);

    % Ridge normal equation: (Phi'Phi + lambda*N*I) beta = Phi'*1
    A    = Temp.' * Temp + lambda * N * eye(nT);
    beta = A \ (Temp.' * ones(N, 1));

    error = ones(N, 1) - Temp * beta;
end
