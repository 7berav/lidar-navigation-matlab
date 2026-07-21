function [d, dPQ, dQP] = hausdorffDistance(P, Q, pctl)
%HAUSDORFFDISTANCE  Symmetric (percentile) Hausdorff distance between point sets.
%   d = hausdorffDistance(P, Q) returns the symmetric Hausdorff distance
%       d = max( h(P,Q), h(Q,P) ),   h(A,B) = max_i min_j ||a_i - b_j||
%   i.e. the worst-case nearest-neighbor error over both directions.
%
%   d = hausdorffDistance(P, Q, pctl) replaces the max over nearest-neighbor
%   distances with the pctl-th percentile (0 < pctl <= 100). This is the
%   robust Hausdorff variant: pctl = 100 is the classical max; pctl = 95
%   ignores the single worst point and is recommended for noisy LiDAR clouds.
%
%   [d, dPQ, dQP] = ... also returns the two directed distances:
%       dPQ = h(P,Q)  (P -> Q, completeness),   dQP = h(Q,P)  (Q -> P, accuracy).
%
%   Mirrors chamferDistance (mean -> max/percentile of the same per-point
%   nearest-neighbor distances). NN distances via nnDist.

    if nargin < 3 || isempty(pctl), pctl = 100; end

    nnPQ = nnDist(P, Q);
    nnQP = nnDist(Q, P);

    if pctl >= 100
        dPQ = max(nnPQ);
        dQP = max(nnQP);
    else
        dPQ = prctile(nnPQ, pctl);
        dQP = prctile(nnQP, pctl);
    end
    d = max(dPQ, dQP);
end
