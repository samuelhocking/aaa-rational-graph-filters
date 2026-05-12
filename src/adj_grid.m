function A = adj_grid(nrow,ncol)
    arguments
        nrow
        ncol
    end
    N = nrow*ncol;
    linear_idx = @(r,c) (r-1)*ncol+c;
    rr = [];
    cc = [];
    vv = [];
    % A = zeros(N);
    for i = 1:nrow
        for j = 1:ncol
            disp([i,j])
            linidx = linear_idx(i,j);
            if j == 1
                disp("j=1")
                rr(end+1) = linidx;
                cc(end+1) = linear_idx(i,j+1);
                vv(end+1) = 1;
                % A(linidx, linear_idx(i,j+1)) = 1;
            elseif j == ncol
                disp("j=ncol")
                rr(end+1) = linidx;
                cc(end+1) = linear_idx(i,j-1);
                vv(end+1) = 1;
                % A(linidx, linear_idx(i,j-1)) = 1;
            else
                disp("j=else")
                rr(end+1) = linidx;
                cc(end+1) = linear_idx(i,j+1);
                vv(end+1) = 1;
                rr(end+1) = linidx;
                cc(end+1) = linear_idx(i,j-1);
                vv(end+1) = 1;
                % A(linidx, linear_idx(i,j+1)) = 1;
                % A(linidx, linear_idx(i,j-1)) = 1;
            end
            if i == 1
                disp("i=1")
                rr(end+1) = linidx;
                cc(end+1) = linear_idx(i+1,j);
                vv(end+1) = 1;
                % A(linidx, linear_idx(i+1,j)) = 1;
            elseif i == nrow
                disp("i=nrow")
                rr(end+1) = linidx;
                cc(end+1) = linear_idx(i-1,j);
                vv(end+1) = 1;
                % A(linidx, linear_idx(i-1,j)) = 1;
            else
                disp("i=else")
                rr(end+1) = linidx;
                cc(end+1) = linear_idx(i+1,j);
                vv(end+1) = 1;
                rr(end+1) = linidx;
                cc(end+1) = linear_idx(i-1,j);
                vv(end+1) = 1;
                % A(linidx, linear_idx(i+1,j)) = 1;
                % A(linidx, linear_idx(i-1,j)) = 1;
            end
        end
    end
    A = sparse(rr,cc,vv,N,N);
end