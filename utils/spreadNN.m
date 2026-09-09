function s = spreadNN(P)
%SPREADNN  Density-normalised mean nearest-neighbour spacing of a point set.
%   s = spreadNN(P) = mean(distance to nearest OTHER point) * sqrt(size(P,1)).
%
%   The sqrt(n) factor removes the trivial dependence on point count: for a
%   fixed 2-D surface, mean NN spacing scales as sqrt(Area/n), so s is a
%   density-free proxy for "how much surface area does this set cover".
%
%   Purpose (Sec 3.2, fn_clump): comparing spreadNN of the rejected true
%   inliers against spreadNN of all true inliers tells whether a given recall
%   loss is spread evenly over the surface or concentrated in one region:
%     ratio ~ 1  -- rejected points scattered uniformly (mild, CD barely moves)
%     ratio < 1  -- rejected points clumped, a whole patch was dropped (severe)
%   Without the sqrt(n) normalisation the ratio would just track the count.
%
%   Returns NaN for fewer than 2 points.

    n = size(P,1);
    if n < 2, s = NaN; return; end

    try
        [~, d] = knnsearch(P, P, 'K', 2);   % Statistics Toolbox
        dnn = d(:,2);                       % column 1 is the point itself
    catch
        dnn  = zeros(n,1);
        sqP  = sum(P.^2, 2).';
        blk  = 2048;
        for i0 = 1:blk:n
            i1 = min(i0+blk-1, n);
            D2 = sum(P(i0:i1,:).^2, 2) + sqP - 2*(P(i0:i1,:)*P.');
            D2(sub2ind(size(D2), 1:(i1-i0+1), i0:i1)) = inf;   % mask self
            dnn(i0:i1) = sqrt(max(min(D2, [], 2), 0));
        end
    end

    s = mean(dnn) * sqrt(n);
end
