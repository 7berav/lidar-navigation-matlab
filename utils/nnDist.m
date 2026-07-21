function dmin = nnDist(A, B)
%NNDIST  Nearest-neighbor Euclidean distance from each point in A to set B.
%   dmin = nnDist(A, B) returns, for each row of A (Na-by-3), the distance to
%   its closest point in B (Nb-by-3). Uses knnsearch when available, else a
%   block-wise brute-force fallback.
%
%   Shared helper for chamferDistance and hausdorffDistance.

    try
        [~, dmin] = knnsearch(B, A);           % Statistics Toolbox
    catch
        nA   = size(A, 1);
        dmin = zeros(nA, 1);
        sqB  = sum(B.^2, 2).';
        blk  = 2048;
        for i0 = 1:blk:nA
            i1 = min(i0 + blk - 1, nA);
            D2 = sum(A(i0:i1,:).^2, 2) + sqB - 2 * (A(i0:i1,:) * B.');
            dmin(i0:i1) = sqrt(max(min(D2, [], 2), 0));
        end
    end
end
