function [cd, hd] = symmetricMetrics(Pgt, Pest, hdPctl, emptyPenalty)
%SYMMETRICMETRICS  Symmetric Chamfer + percentile Hausdorff between two point sets.
%   [cd, hd] = symmetricMetrics(Pgt, Pest, hdPctl, emptyPenalty)
%
%   cd = 0.5*(mean_i min_j||p_i-q_j|| + mean_j min_i||q_j-p_i||)
%   hd = max( prctile(d_{gt->est}, hdPctl), prctile(d_{est->gt}, hdPctl) )
%
%   emptyPenalty is returned for both when either set is empty (a degenerate
%   fit that produces no isosurface must be penalised, not skipped). Callers
%   use the GT bounding-box diagonal.
%
%   Shared by Sec 3.1 (lcurve_demo) and Sec 3.2 (ransac_*) so that cd_occ /
%   cd_all are computed by literally the same code in both sections -- this
%   is what makes the 3.1 -> 3.2 loss-decomposition ladder a valid comparison.
%   Do not fork this function.
%
%   NN distances via nnDist (knnsearch or brute-force fallback).

    if isempty(Pgt) || isempty(Pest)
        cd = emptyPenalty; hd = emptyPenalty; return;
    end
    d1 = nnDist(Pgt, Pest);
    d2 = nnDist(Pest, Pgt);
    cd = 0.5*(mean(d1) + mean(d2));
    hd = max(prctile(d1, hdPctl), prctile(d2, hdPctl));
end
