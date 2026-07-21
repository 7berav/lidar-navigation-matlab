function d = chamferDistance(P, Q)
%CHAMFERDISTANCE  Symmetric Chamfer distance between two point sets.
%   d = chamferDistance(P, Q) returns the (symmetric) Chamfer distance
%   between point clouds P (Np-by-3) and Q (Nq-by-3):
%       d = 0.5 * ( mean_i min_j ||p_i - q_j|| + mean_j min_i ||q_j - p_i|| )
%
%   Nearest-neighbor distances via nnDist (knnsearch or brute-force fallback).

    d = 0.5 * (mean(nnDist(P, Q)) + mean(nnDist(Q, P)));
end
