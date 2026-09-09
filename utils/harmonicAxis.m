function [ax, gapRatio, l2frac, Q] = harmonicAxis(beta, order, Funcs)
%HARMONICAXIS  Major axis from the l=2 block of a Fischer-basis fit.
%
%   [ax, gapRatio, l2frac, Q] = harmonicAxis(beta, order, Funcs)
%
%   Fischer 기저는 구면조화 분해다: order N 이면 levels = N:-2:0 이고 레벨 l 의 항이
%   2l+1 개다. 그중 l=2 블록(5개)이 사중극(quadrupole) 성분이고, 대칭 무자취 3x3
%   텐서 Q 하나와 일대일이다. Q 의 고유벡터가 주축이다.
%
%   등위면을 샘플링해서 PCA 를 거는 방식(principalAxis)과 달리 **계수만** 쓴다.
%   외삽면이 가려진 쪽으로 길게 뻗어도 그 로브가 2차 모멘트를 지배하는 일이 없다.
%
%   ax        주축 단위벡터. 부호는 임의이므로 acosd(min(1,abs(dot(a,b)))) 로 비교한다.
%   gapRatio  |e1|/|e2| (내림차순 절대 고유값). principalAxis 의 gapRatio 와 같은 역할.
%   l2frac    ||beta_l2|| / ||beta||. 사중극 정보가 아예 없으면(cube, sphere) 0 에 가깝다.
%             자세 유효성 게이트는 gapRatio 가 아니라 이 값으로 거는 것이 맞다.
%   Q         3x3 대칭 무자취 텐서.
%
%   주축 선택: 등위면은 f=1 이므로 어떤 방향으로 f 가 작을수록 표면이 멀다.
%   따라서 장축은 사중극 기여가 가장 음인 방향 = 최소 고유값의 고유벡터다.

    persistent MAPS
    if isempty(MAPS), MAPS = containers.Map('KeyType','double','ValueType','any'); end

    lv  = order:-2:0;
    idx = find(lv == 2, 1);
    if isempty(idx)
        ax = [NaN NaN NaN].'; gapRatio = NaN; l2frac = NaN; Q = nan(3); return;
    end
    i0 = sum(2*lv(1:idx-1) + 1) + 1;      % l=2 블록 시작 인덱스
    b2 = beta(i0:i0+4);
    l2frac = norm(b2) / max(norm(beta), eps);

    % 기저 -> 이차형식 맵은 order 마다 한 번만 만든다
    if ~isKey(MAPS, order)
        MAPS(order) = buildMap(order, Funcs, i0);
    end
    M = MAPS(order);                       % 6 x 5, 열마다 [Qxx Qyy Qzz Qxy Qxz Qyz]
    q = M * b2(:);
    Q = [q(1) q(4) q(5); q(4) q(2) q(6); q(5) q(6) q(3)];
    Q = Q - trace(Q)/3 * eye(3);           % 무자취로 강제(수치 잔여 제거)

    [V, D] = eig(Q, 'vector');
    [~, ord] = sort(abs(D), 'descend');
    dA = D(ord); V = V(:, ord);
    gapRatio = abs(dA(1)) / max(abs(dA(2)), eps);

    [~, imin] = min(D);                    % 가장 음인 방향 = 표면이 가장 먼 방향
    ax = V(:, find(ord == imin, 1));
    if isempty(ax), ax = V(:,1); end
    ax = ax / norm(ax);
end

function M = buildMap(order, Funcs, i0)
% l=2 기저함수 5개를 단위구 위에서 이차형식에 맞춰 푼다.
% Phi_m(x) = r^(order-2) * H_2m(x) 이고 r=1 위에서 H_2m 은 순수 이차형식이다.
    rng(0);
    P = randn(400,3);
    P = P ./ vecnorm(P,2,2);               % 단위구
    Phi = Funcs(P);                        % 400 x nT
    A = Phi(:, i0:i0+4);                   % 400 x 5
    x = P(:,1); y = P(:,2); z = P(:,3);
    D = [x.^2, y.^2, z.^2, 2*x.*y, 2*x.*z, 2*y.*z];   % 400 x 6
    M = D \ A;                             % 6 x 5
end
