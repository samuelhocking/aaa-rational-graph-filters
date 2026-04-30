function h = idealHPF(cutoff)
    h = @(lambda) 1 .* (lambda >= cutoff);
end