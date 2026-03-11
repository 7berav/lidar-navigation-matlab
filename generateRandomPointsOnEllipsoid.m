function P = generateRandomPointsOnEllipsoid(N, a, b, c)


    arguments
        N (1,1) {mustBeInteger, mustBePositive}
        a (1,1) {mustBePositive}
        b (1,1) {mustBePositive}
        c (1,1) {mustBePositive}
    end

    P = zeros(N,3);
    count = 0;

    % --- 1) 수락확률 상계 p_max 추정 (프리패스) ---
    M = 20000;                
    u0 = randn(M,3); u0 = u0 ./ vecnorm(u0,2,2);
    ux = u0(:,1); uy = u0(:,2); uz = u0(:,3);

    D = (ux.^2)/(a^2) + (uy.^2)/(b^2) + (uz.^2)/(c^2);
    S = (ux.^2)/(a^4) + (uy.^2)/(b^4) + (uz.^2)/(c^4);
    w = (D.^2) ./ sqrt(S);      % 비정규화 수락 가중치 ∝ dA 보정

    pmax = max(w) * 1.10;       % 10% 안전계수

    % --- 2) 본 샘플링 (거절법) ---
    while count < N
        batch = max(4*(N-count), 256);   % 벡터화 배치
        u = randn(batch,3); u = u ./ vecnorm(u,2,2);
        ux = u(:,1); uy = u(:,2); uz = u(:,3);

        D = (ux.^2)/(a^2) + (uy.^2)/(b^2) + (uz.^2)/(c^2);
        S = (ux.^2)/(a^4) + (uy.^2)/(b^4) + (uz.^2)/(c^4);
        w = (D.^2) ./ sqrt(S);
        accProb = w / pmax;

        r = rand(batch,1);
        keep = r < accProb;

        if ~any(keep), continue; end

        u_keep = u(keep,:);
        Dk = D(keep);

        % 광선-타원체 교점 (표면)
        X = u_keep ./ sqrt(Dk);

        k = size(X,1);
        take = min(k, N-count);
        P(count+1:count+take, :) = X(1:take, :);
        count = count + take;
    end
end