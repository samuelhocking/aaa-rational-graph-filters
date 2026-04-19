function h = idealHPF(cutoff)
    h = @(lambda) lambda .* (lambda >= cutoff);
end