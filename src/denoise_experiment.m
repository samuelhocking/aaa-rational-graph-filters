%DENOISE_EXPERIMENT  Graph low-pass denoising on a similarity-weighted pixel grid.
%
%   Builds a checkerboard test image, adds Gaussian noise, forms weights from the
%   noisy intensities, applies an ideal low-pass filter h(L)x via eigendecomposition,
%   and reports RMSE. Optional rational/polynomial variants can be plugged in later.

clearvars -except colors LW lw MS ms fs cmap

H = 48;
W = 48;
[RR, CC] = meshgrid(1:W, 1:H);
block = 12;
clean = double(mod(floor((RR - 1) / block) + floor((CC - 1) / block), 2) == 0);

sigma_n = 0.25;
rng(42);
noise = sigma_n * randn(H, W);
noisy = clean + noise;

sigma_r = 0.2;
A = build_pixel_graph(noisy, sigma_r);
L = laplacian(A);
n = H * W;

[V, lam] = eig(full(L), "vector");
lam = real(lam);

dd = full(sum(A, 2));
lam_max = max(lam);

cutoff_frac = 0.15;
cutoff = cutoff_frac * lam_max;
h = idealLPF(cutoff);

x = noisy(:);
x_clean = clean(:);
yh = V * (h(lam) .* (V.' * x));

rmse_noisy = sqrt(mean((x - x_clean).^2));
rmse_filt = sqrt(mean((yh - x_clean).^2));

fprintf("RMSE noisy: %.4f\n", rmse_noisy);
fprintf("RMSE filtered (ideal LPF): %.4f\n", rmse_filt);

figure("Name", "Graph denoising");
tiledlayout(1, 3, "Padding", "compact");
nexttile;
imagesc(clean);
axis image off;
title("Clean");
colorbar;

nexttile;
imagesc(noisy);
axis image off;
title("Noisy");
colorbar;

nexttile;
imagesc(reshape(yh, H, W));
axis image off;
title("Filtered");
colorbar;
