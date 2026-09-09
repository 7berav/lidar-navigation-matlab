function [ax, gapRatio] = principalAxis(P)
%PRINCIPALAXIS  Dominant axis of a point set via covariance eigendecomposition.
%   [ax, gapRatio] = principalAxis(P) returns the unit eigenvector of cov(P)
%   with the largest eigenvalue, and gapRatio = lambda1/lambda2 (>= 1).
%
%   gapRatio is a VALIDITY FLAG, not a quality score. The dominant axis is only
%   meaningful when the two leading 2nd moments are well separated:
%     gapRatio >> 1  -- elongated (cylinder, hexagonal prism): axis meaningful
%     gapRatio -> 1  -- 2nd moments degenerate: axis is arbitrary noise
%
%   A cube has EXACTLY isotropic 2nd moments (gapRatio == 1 up to sampling
%   noise), so no covariance-based method -- PCA on the raw cloud or on the
%   fitted level set -- can recover its orientation. That is a property of the
%   object, not of the estimator, and is precisely why a 4th-order basis is
%   needed. Callers must gate on gapRatio and report NaN rather than a number.
%
%   Sign is arbitrary (ax and -ax are the same axis); compare axes with
%   acosd(min(1,abs(dot(a,b)))) so the result lands in [0,90] degrees:
%   0 deg = aligned (antiparallel counts as aligned), 90 deg = perpendicular.
%
%   REFERENCE LEVEL: two uniformly random 3-D axes give |cos| ~ U[0,1], hence
%   a median angle of 60 deg (mean 57.3, IQR 41.5-75.5). Always report angles
%   against this line -- above ~60 deg the estimate carries no orientation
%   information at all, and well above it the axis is systematically wrong
%   rather than merely noisy. Without the reference an angle of "50 deg" reads
%   as a moderate error when it is in fact near chance.

    if size(P,1) < 3
        ax = [NaN NaN NaN].'; gapRatio = NaN; return;
    end

    C = cov(P);
    [Vec, D] = eig(C, 'vector');
    [d, ord] = sort(D, 'descend');
    Vec = Vec(:, ord);

    ax = Vec(:,1);
    ax = ax / norm(ax);
    gapRatio = d(1) / max(d(2), eps);
end
