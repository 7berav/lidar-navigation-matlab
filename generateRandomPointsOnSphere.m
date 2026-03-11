function P = generateRandomPointsOnSphere(N)


    R = 1;
    u = randn(N,3);
    u = u ./ vecnorm(u,2,2);   % 각 행을 L2 정규화 (단위벡터)
    P = R * u;                 % 반지름 스케일
end