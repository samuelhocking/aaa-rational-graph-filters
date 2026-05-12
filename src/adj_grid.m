function A = adj_grid(nrow,ncol)
    arguments
        nrow
        ncol
    end
    N = nrow*ncol;
    linear_idx = @(r,c) (r-1)*ncol+c;
    A = zeros(N);
    for i = 1:ncol
        for j = 1:nrow
            disp([i,j])
            linidx = linear_idx(i,j);
            if j == 1
                A(linidx, linear_idx(i,j+1)) = 1;
            elseif j == ncol
                A(linidx, linear_idx(i,j-1)) = 1;
            else
                A(linidx, linear_idx(i,j+1)) = 1;
                A(linidx, linear_idx(i,j-1)) = 1;
            end
            if i == 1
                A(linidx, linear_idx(i+1,j)) = 1;
            elseif i == nrow
                A(linidx, linear_idx(i-1,j)) = 1;
            else
                A(linidx, linear_idx(i+1,j)) = 1;
                A(linidx, linear_idx(i-1,j)) = 1;
            end
        end
    end
end