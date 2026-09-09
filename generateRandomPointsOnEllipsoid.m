function P = generateRandomPointsOnEllipsoid(N, a, b, c)


    arguments
        N (1,1) {mustBeInteger, mustBePositive}
        a (1,1) {mustBePositive}
        b (1,1) {mustBePositive}
        c (1,1) {mustBePositive}
    end

    P = zeros(N,3);
    count = 0;

    % 면적균일 ray-casting 의 야코비안:
    %   X = u/sqrt(D),   D = Σ u_i²/a_i²,   S = Σ u_i²/a_i⁴
    %   r = |X| = 1/sqrt(D),  법선 n ∝ (u_i/a_i²),  n̂·û = D/sqrt(S)
    %   dA/dΩ = r² / (n̂·û) = (1/D) / (D/sqrt(S)) = **sqrt(S)/D²**
    % (2026-09-01 수정) 이전 구현은 이 역수인 D²/sqrt(S) 를 써서 면적 보정이 뒤집혀
    % 있었다. 그 결과 축비와 무관하게 좌표 분산이 거의 같아졌다 — 5:1:1 타원체도
    % 공분산 고유값비 1.00, 심지어 var(z) > var(x) (장축의 분산이 더 작음).
    % 즉 기하학적으로는 타원면 위의 점이 맞지만 2차 모멘트가 등방이라 자세(주축)
    % 평가에 쓸 수 없는 점군이었다. 수정 후 5:1:1 의 고유값비 17.3 으로 정상화.
    % 기존 산출물에는 ellipsoid 가 쓰인 적이 없어 소급 영향 없음.

    % --- 1) 수락확률 상계 p_max 추정 (프리패스) ---
    M = 20000;
    u0 = randn(M,3); u0 = u0 ./ vecnorm(u0,2,2);
    ux = u0(:,1); uy = u0(:,2); uz = u0(:,3);

    D = (ux.^2)/(a^2) + (uy.^2)/(b^2) + (uz.^2)/(c^2);
    S = (ux.^2)/(a^4) + (uy.^2)/(b^4) + (uz.^2)/(c^4);
    w = sqrt(S) ./ (D.^2);      % 비정규화 수락 가중치 ∝ dA/dΩ

    pmax = max(w) * 1.10;       % 10% 안전계수

    % --- 2) 본 샘플링 (거절법) ---
    while count < N
        batch = max(4*(N-count), 256);   % 벡터화 배치
        u = randn(batch,3); u = u ./ vecnorm(u,2,2);
        ux = u(:,1); uy = u(:,2); uz = u(:,3);

        D = (ux.^2)/(a^2) + (uy.^2)/(b^2) + (uz.^2)/(c^2);
        S = (ux.^2)/(a^4) + (uy.^2)/(b^4) + (uz.^2)/(c^4);
        w = sqrt(S) ./ (D.^2);
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