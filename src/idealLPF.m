function h = idealLPF(cutoff)
    h = @(lambda) 1 .* (lambda <= cutoff);
end