function D = make_digraph_from_ugraph(G)
    [i,j,w] = find(adjacency(G,'weighted'));
    s = [i; j]; t = [j; i]; c = [w; w];
    D = digraph(s, t, c, numnodes(G));
end

