function A = adj_Kn(n)
    A = ones(n,n) - diag(ones(n,1));
end