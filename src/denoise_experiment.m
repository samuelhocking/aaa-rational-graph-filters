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

cutoff_lo = 0+1e-2;
cutoff_hi = cutoff;
alpha_r = (2/3)/(cutoff_hi-cutoff_lo);
beta_r = -1/2 * (cutoff_lo+cutoff_hi) * alpha_r;

% load approxParams
% pol_r = (pol - beta_r)./alpha_r;
% res_r = res/alpha_r;
% r = @(zz) real(pfeval(zz, pol_r, res_r));

[r, pol_r, res_r, zer_r] = getRatApprox(h, lam_max, npts=10000, cutoff_lo=cutoff_lo, cutoff_hi=cutoff_hi);

rL = @(x) real(applyrLx(x,L,pol_r,res_r));


pcf = chebfun(h, [0,lam_max]);
cheb_poly_coeffs = chebcoeffs(pcf);  
poly_deg = 2*length(pol_r)-1;
cheb_poly_coeffs_deg = cheb_poly_coeffs(1:poly_deg+1);
alpha_p = 2/(lam_max-0);
beta_p = -1/2 * (lam_max+0) * alpha_p;
L_p = alpha_p * L + beta_p*eye(n);
pL = @(x) applypLx(x,L_p,cheb_poly_coeffs_deg);
p = @(zz) chebpolyval(flip(cheb_poly_coeffs_deg), alpha_p*zz + beta_p);

yr = rL(x);
yp = pL(x);

rmse_noisy = sqrt(mean((x - x_clean).^2));
rmse_filt = sqrt(mean((yh - x_clean).^2));
rmse_poly = sqrt(mean((yp - x_clean).^2));
rmse_rat = sqrt(mean((yr - x_clean).^2));

fprintf("RMSE noisy: %.4f\n", rmse_noisy);
fprintf("RMSE filtered (ideal LPF): %.4f\n", rmse_filt);
fprintf("RMSE polynomial: %.4f\n", rmse_poly);
fprintf("RMSE rational: %.4f\n", rmse_rat);


f1 = figure("Name", "Graph denoising");
tl = tiledlayout(1, 4, "TileSpacing", "tight", "Padding", "compact");
c_min = min([noisy(:); yh(:); yr(:); yp(:)]);
c_max = max([noisy(:); yh(:); yr(:); yp(:)]);

nexttile;
imagesc(noisy);
axis image off;
title(sprintf("Noisy \n(RMSE = $%.4f$)", rmse_noisy))
% colorbar;

nexttile;
imagesc(reshape(yh, H, W));
axis image off;
title(sprintf("Exact filter \n(RMSE = $%.4f$)", rmse_filt))
% colorbar;

nexttile;
imagesc(reshape(yp, H, W));
axis image off;
title(sprintf("Polynomial \n(RMSE = $%.4f$)", rmse_poly))
% colorbar;

nexttile;
imagesc(reshape(yr, H, W));
axis image off;
title(sprintf("Rational \n(RMSE = $%.4f$)", rmse_rat))
% title("Filtered \\ RMSE=$$");
% colorbar;

clim([c_min c_max]);
cb = colorbar;
cb.Layout.Tile = 'east';

fontsize(fs, "points")

exportgraphics(f1, "./figures/denoising_panel.png", 'Resolution', 300);


f2 = figure();
f2.Position(3:4) = [400 400];
spy(L);
title("Sparsity pattern of $L$")
fontsize(fs, "points")

exportgraphics(f2, "./figures/denoising_sparsity.png", 'Resolution', 300);


