function A = adj_undirCycle(n)
    d = ones(n,1);
    diags = [d, d, d, d];
    v = [-(n-1) -1 1 n-1];
    A = spdiags(diags, v, n, n);
end