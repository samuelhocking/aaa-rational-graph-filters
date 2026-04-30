function h = idealBPF(low,high)
    h = @(lambda) 1 .* ((lambda >= low) & (lambda <= high));
end