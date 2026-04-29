function [beta, error] = regressionFourthOrder(inputMatrix, symbterm, lambda, order, s)
    % Ridge regression (N-normalized):
    %   argmin  (1/N)||Phi*beta - 1||^2 + beta'*D*beta
    %   =>  (Phi'Phi + N*D) beta = Phi'*1
    %
    % 인수
    %   lambda : 스칼라 (전역 페널티 강도, 기본값 0 = OLS)
    %   order  : 다항식 최고 차수 (선택, 주면 Sobolev 가중 ridge 적용)
    %   s      : Sobolev 차수 (선택, 기본값 1 = H^1)
    %
    % 동작 모드
    %   order 없음 → D = lambda * I          (등방성 ridge)
    %   order 있음 → D = lambda * diag(w)    (Sobolev 가중 ridge)
    %     w_i = [l_i*(l_i+1)]^s  (level l 항의 Laplace-Beltrami 고유값)
    %     l=0(상수항): w=0, l=2: w=6^s, l=4: w=20^s, l=6: w=42^s
    %
    % N 스케일링: Phi'Phi ~ O(N) 이므로 N*D 를 더해야
    % lambda 의 의미가 점 개수 N 에 무관하게 일정.

    if nargin < 3 || isempty(lambda), lambda = 0; end
    if nargin < 4, order = []; end
    if nargin < 5 || isempty(s), s = 1; end

    Temp = calculateFourthOrder(inputMatrix, symbterm);
    N    = size(Temp, 1);
    nT   = size(Temp, 2);

    % 페널티 대각 벡터 w 구성
    if lambda == 0
        P = zeros(nT);                         % OLS

    elseif isempty(order)
        P = lambda * N * eye(nT);              % 등방성 ridge

    else
        % Sobolev 가중 ridge: level l → [l*(l+1)]^s
        levels = order:-2:0;
        w = zeros(nT, 1);
        idx = 1;
        for l = levels
            nk = 2*l + 1;
            w(idx : idx+nk-1) = (l*(l+1))/order/(order+1);  % l=0 → 0 (상수항 면제)
            idx = idx + nk;
        end
        P = lambda * N * diag(w);
    end

    % Ridge normal equation: (Phi'Phi + P) beta = Phi'*1
    A    = Temp.' * Temp + P;
    beta = A \ (Temp.' * ones(N, 1));

    error = ones(N, 1) - Temp * beta;
end
