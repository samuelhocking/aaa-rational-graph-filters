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
    yj = zeros(n,npol);
    for j = 1:npol
        Lj = L-pol(j)*eye(n);
        [Qj,Rj] = qr(Lj);
        yj(:,j) = Rj \ (Qj' * x);
    end
    y = opts.a0 + yj * res;
end