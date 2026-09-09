function P = meshSurfacePoints(beta, PhiGrid, GX, GY, GZ, center, nSamples)
%MESHSURFACEPOINTS  Sample points uniformly on the implicit level set f=1.
%   P = meshSurfacePoints(beta, PhiGrid, GX, GY, GZ, center, nSamples)
%
%   Evaluates f = PhiGrid*beta on the grid (GX,GY,GZ), extracts the f=1
%   isosurface, then draws nSamples points with probability proportional to
%   triangle area (area-weighted uniform surface sampling). Vertices are
%   shifted by `center` to undo the fit-frame translation.
%
%   PhiGrid must be the design matrix evaluated at [GX(:) GY(:) GZ(:)] with
%   the same basis that produced beta (see calculateFourthOrder).
%
%   Returns zeros(0,3) when the level set is empty or degenerate; callers
%   must treat that as a failed fit (see symmetricMetrics' emptyPenalty).
%
%   Shared by Sec 3.1 and Sec 3.2 -- do not fork.

    F  = reshape(PhiGrid*beta, size(GX));
    fv = isosurface(GX, GY, GZ, F, 1);
    if isempty(fv.vertices) || isempty(fv.faces)
        P = zeros(0,3); return;
    end

    V = fv.vertices + center;
    A = V(fv.faces(:,1),:); B = V(fv.faces(:,2),:); C = V(fv.faces(:,3),:);
    area  = 0.5*vecnorm(cross(B-A, C-A, 2), 2, 2);
    valid = isfinite(area) & area > eps;
    if ~any(valid)
        P = zeros(0,3); return;
    end
    triPool = find(valid);
    pick = randsample(triPool, nSamples, true, area(valid));
    u = sqrt(rand(nSamples,1)); v = rand(nSamples,1);
    P = (1-u).*A(pick,:) + (u.*(1-v)).*B(pick,:) + (u.*v).*C(pick,:);
end
