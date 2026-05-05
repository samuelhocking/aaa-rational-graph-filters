function y = applyrLx(x,L,pol,res,opts)
    arguments
        x
        L
        pol
        res
        opts.a0 = 0
    end
    n = length(x);
    npol = length(pol);
    y = opts.a0*x;
    for j = 1:npol
        y = y + res(j) * ((L-pol(j)*eye(n)) \ x);
        % Lj = L-pol(j)*eye(n);
        % [Qj,Rj] = qr(Lj);
        % yj(:,j) = Rj \ (Qj' * x);
        % yj(:,j) = Lj \ x;
    end
    % y = opts.a0*x + yj * res;
end