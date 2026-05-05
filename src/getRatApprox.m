function [r, pol, res, zer] = getRatApprox(h, lam_max, opts)
    arguments
        h
        lam_max
        opts.npts = 10000
        opts.eps = 1e-12
        opts.tol = 1e-13
        opts.cutoff_lo = false
        opts.cutoff_hi = false
    end
    cutoff_lo = opts.cutoff_lo;
    cutoff_hi = opts.cutoff_hi;
    eps = opts.eps;
    tol = opts.tol;

    zz = linspace(-lam_max,lam_max,opts.npts).';
    if isnumeric(cutoff_lo)
        zz = [zz; cutoff_lo; cutoff_lo-eps; cutoff_lo+eps];
    end
    if isnumeric(cutoff_hi)
        zz = [zz; cutoff_hi; cutoff_hi-eps; cutoff_hi+eps];
    end
    zz = unique(zz);

    hh = h(zz);         % compute reference data for filtered sample points

    [~, pol, res, zer, ~, ~, ~, ~] = aaa(hh, zz, 'tol', tol);   % AAA algorithm

    I = find(imag(pol) == 0);   % find indices of real poles
    pol(I) = []; res(I) = [];   % remove real poles

    res = (1./(zz - pol.'))\hh;                         % compute LS residues

    r = @(zz) real(pfeval(zz, pol, res));     % get function handle for partial fraction rational filter
end