function I = expClust(m, opts);
    % produces logspaced points on [a,b] clustered near a
    % reverse = true clusters them near b 
    arguments
        m
        opts.a = 0
        opts.b = 1
        opts.reverse = false
        opts.exponent = 1
    end
    if ~opts.reverse
        I = (opts.b-opts.a)*logspace(-14,0,m).^opts.exponent+opts.a;
    else
        I = -(opts.b-opts.a)*logspace(-14,0,m).^opts.exponent+opts.b;
    end
    I = I(:);
end