function A = build_pixel_graph(img, sigma_r)
    %BUILD_PIXEL_GRAPH  Undirected 4-neighbor graph with Gaussian intensity weights.
    %
    %   A = build_pixel_graph(img, sigma_r)
    %
    %   img       H-by-W grayscale (double). Edge weight between neighbors i,j is
    %             exp(-(I_i - I_j)^2 / (2*sigma_r^2)).
    %   sigma_r   intensity scale for similarity (larger => weaker similarity penalty).

    arguments
        img (:,:) double
        sigma_r (1,1) double {mustBePositive}
    end

    [H, W] = size(img);
    n = H * W;
    nnz_est = 4 * n;
    ri = zeros(nnz_est, 1);
    ci = zeros(nnz_est, 1);
    vv = zeros(nnz_est, 1);
    k = 0;
    inv_den = 1 / (2 * sigma_r^2);

    for r = 1:H
        for c = 1:W
            i = (r - 1) * W + c;
            val_i = img(r, c);
            if c < W
                j = i + 1;
                val_j = img(r, c + 1);
                w = exp(-(val_i - val_j)^2 * inv_den);
                k = k + 1;
                ri(k) = i;
                ci(k) = j;
                vv(k) = w;
            end
            if r < H
                j = i + W;
                val_j = img(r + 1, c);
                w = exp(-(val_i - val_j)^2 * inv_den);
                k = k + 1;
                ri(k) = i;
                ci(k) = j;
                vv(k) = w;
            end
        end
    end

    ri = ri(1:k);
    ci = ci(1:k);
    vv = vv(1:k);
    A = sparse([ri; ci], [ci; ri], [vv; vv], n, n);
end
