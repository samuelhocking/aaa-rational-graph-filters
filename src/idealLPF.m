function h = idealLPF(cutoff)
    h = @(lambda) lambda .* (lambda <= cutoff);
end