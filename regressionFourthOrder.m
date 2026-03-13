function [beta, error] = regressionFourthOrder(inputMatrix, symbterm, lambda)
    % Ridge regression: beta = (Phi'*Phi + lambda*I) \ (Phi'*1)
    % lambda = 0 (default) => ordinary least squares (OLS)
    % lambda > 0           => ridge regression (L2 penalty on beta)
    %
    % Score Sj = mean(exp(-r^2/(2*sigma^2))) 은 호출부에서 그대로 사용.
    % 페널티는 beta 추정에만 적용하고 RANSAC 잔차/스코어 계산은 변경 없음.

    if nargin < 3 || isempty(lambda)
        lambda = 0;
    end

    Temp = calculateFourthOrder(inputMatrix, symbterm);
    N    = size(Temp, 1);
    nT   = size(Temp, 2);

    % Ridge normal equation: (Phi'Phi + lambda*I) beta = Phi'*1
    A    = Temp.' * Temp + lambda * eye(nT);
    beta = A \ (Temp.' * ones(N, 1));

    error = ones(N, 1) - Temp * beta;
end
