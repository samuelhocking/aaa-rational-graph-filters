function h = idealBPF(low,high)
    h = @(lambda) lambda .* ((lambda >= low) & (lambda <= high));
end