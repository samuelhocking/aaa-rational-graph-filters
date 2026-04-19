function A = adj_dirCycle(n)
    d = ones(n,1);
    diags = [d, d];
    v = [-(n-1) 1];
    A = spdiags(diags, v, n, n);
end