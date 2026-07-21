function points = generateRandomPointsOnCylinder(N)
%GENERATERANDOMPOINTSONCYLINDER Uniform points on a closed unit cylinder.
%   Radius = 1, z in [-1,1]. Samples are uniform with respect to surface
%   area: side area 4*pi (probability 2/3), two caps total area 2*pi
%   (probability 1/3).

    arguments
        N (1,1) {mustBeInteger,mustBePositive}
    end

    points = zeros(N,3);
    isSide = rand(N,1) < 2/3;

    % Curved side: uniform azimuth and height.
    nSide = nnz(isSide);
    theta = 2*pi*rand(nSide,1);
    points(isSide,:) = [cos(theta), sin(theta), 2*rand(nSide,1)-1];

    % End caps: uniform disk area requires radius sqrt(U).
    capIdx = find(~isSide);
    nCap = numel(capIdx);
    theta = 2*pi*rand(nCap,1);
    radius = sqrt(rand(nCap,1));
    zSign = 2*(rand(nCap,1) >= 0.5)-1;
    points(capIdx,:) = [radius.*cos(theta), radius.*sin(theta), zSign];
end
