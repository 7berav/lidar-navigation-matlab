function P = makeShapePoints(shapeName, N, aspect)
%MAKESHAPEPOINTS  Sample N points on the surface of a named synthetic shape.
%   P = makeShapePoints(shapeName, N) returns an N-by-3 point set.
%   P = makeShapePoints(shapeName, N, aspect) stretches the body z axis by
%   `aspect` before re-normalising, giving an elongated variant.
%
%   Supported names, ordered by symmetry order (= difficulty for attitude
%   recovery from the fitted polynomial):
%     'sphere'    - full rotational symmetry; 2nd order suffices (sanity floor)
%     'cylinder'  - continuous symmetry about one axis; only the long axis is
%                   physically meaningful
%     'ellipsoid' - three distinct 2nd moments; axes recoverable at 2nd order
%     'hexagon'   - D6h hexagonal prism; 2nd moments give the long axis only,
%                   the in-plane 6-fold angle needs 6th order
%     'cube'      - Oh; 2nd moments are exactly isotropic, so PCA and the
%                   quadratic coefficients carry NO orientation. 4th order is
%                   required. This is why the basis is 4th order.
%
%   Shared by Sec 3.1 and Sec 3.2 -- do not fork.

%   Elongated variants ('*_long') exist because attitude recoverability is
%   bounded by how elongated the object actually is, not only by occlusion.
%   The stock cylinder has axis_gap ~1.3 and the hexagonal prism ~2.1, which is
%   barely above the gapMin=1.15 validity gate -- so a large theta_axis there
%   may say more about the object than about the estimator. The '_long'
%   variants raise the gap well clear of the gate and separate those two causes.
%   They also match real RPOD targets (rocket bodies, stacked buses) better.

    if nargin < 3 || isempty(aspect), aspect = 1; end

    switch lower(shapeName)
        case {'hexagon','hexagonal_prism'}
            P = generateRandomPointsOnHexagonPrism(N);
        case 'hexagon_long'
            P = generateRandomPointsOnHexagonPrism(N); aspect = 3;
        case 'cube'
            P = generateRandomPointsOnCube(N);
        case 'cylinder'
            P = generateRandomPointsOnCylinder(N);
        case 'cylinder_long'
            P = generateRandomPointsOnCylinder(N); aspect = 3;
        case 'sphere'
            P = generateRandomPointsOnSphere(N);
        case 'ellipsoid'
            P = generateRandomPointsOnEllipsoid(N, 1.2, 0.9, 0.7);
        case 'ellipsoid_long'
            P = generateRandomPointsOnEllipsoid(N, 1.2, 0.9, 0.7); aspect = 3;
        otherwise
            error('makeShapePoints:unknownShape', 'Unknown shapeName: %s', shapeName);
    end

    if aspect ~= 1
        % z 를 늘린 뒤 **최대 반경을 1.5 로 재정규화**한다. 그냥 늘리기만 하면 회전 후
        % 등위면 격자(fit 프레임 ±1.9)를 벗어나 표면이 잘린다.
        P(:,3) = P(:,3) * aspect;
        P = P / max(abs(P(:))) * 1.5;
    end
end
