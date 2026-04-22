function r = pfeval(zz, pol, res, opts)
    arguments
        zz
        pol
        res
        opts.a0 = 0
    end
    zz = zz(:);
    pol = pol(:);
    res = res(:);
    r = opts.a0 + (1./(zz - pol.'))*res;
end