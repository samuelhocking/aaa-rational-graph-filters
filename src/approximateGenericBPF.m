% CS 2252 Final Project
% Script to approximate a generic band pass filter
% The generic approximation can be shifted and scaled to form an LPF, BPF, or HPF

%% Initial Setup

clear all
close all

set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter', 'latex')
% set(groot, 'defaultFigurePosition', [100 200 800 200]); % [x y width height]

colors = get(groot,'defaultAxesColorOrder');

LW = "LineWidth"; lw = 1.2;
MS = "MarkerSize"; ms = 15;
fs = 12;
cmap = 'parula';

warning('off', 'MATLAB:legend:IgnoringExtraEntries');

%%

plottingflag = 0;

%% Plan

% Given Laplacian L and input x filtered signal y is
% y = h(L)x

% 0. Construct G and L
% 1. Declare ideal filter h(lambda)
% 2. Compute L=VDV^{-1} and exact application y = Vh(D)V^{-1}x
% 3a. Compute rational approximation r(lambda) ~= h(lambda)
% 3b. Compute approximate application y ~= r(L)x
% 4a. Compute polynomial approximation p(lambda) ~= h(lambda)
% 4b. Compute approximate application y ~= p(L)x
% 5. Plotting

%% 1. Declare ideal filter h(lambda)

% generic BPF     
h_generic = idealBPF(-1/2,1/2);

%% 3a. Compute rational approximation r(lambda) ~= h(lambda)

% generic upper bound for lambda is 2*max(degree)
m_equi = 50000;
zz = linspace(-100,100,m_equi).';                     % basic equispaced grid
eps = 1e-12;                                % define perturbation
zz = [zz; -1; 1; -1-eps; -1+eps; 1-eps; 1+eps];    % manually add points very close to the jumps
zz = unique(zz);                                    % remove duplicates
m = length(zz);

m = length(zz);

hh = h_generic(zz);         % compute reference data for filtered sample points

tol = 1e-13;
[r, pol, res, zer, zj, fj, wj, errvec] = aaa(hh, zz, 'tol', tol);   % AAA algorithm

% I = find((imag(pol) == 0) .* (real(pol) >= 0) .* (real(pol) <= 2*max(dd)));  % remove real poles in approximation domain
% I = find((imag(pol) == 0) .* (real(pol) >= min(zz)) .* (real(pol) <= max(zz)));  % remove real poles in approximation domain
I = find(imag(pol) == 0);   % find indices of real poles
pol(I) = []; res(I) = [];   % remove real poles

AA = [1./(zz - pol.')];             % set up least-squares matrix for preset poles
bb = hh;                            % RHS is filtered sample points

aa = AA\bb;                         % compute LS residues
res = aa;

% rpf = @(zz) pfeval(zz, pol, aa(2:end), a0=aa(1));
rpf = @(zz) real(pfeval(zz, pol, res));     % get function handle for partial fraction rational filter
% rpf = @(zz) pfeval(zz, pol, res);     % get function handle for partial fraction rational filter

%% 4a. Compute polynomial approximation p(lambda) ~= h(lambda)

p = chebfun(h_generic, [-1,1]);              % compute polynomial approximation to filter
% p0 = chebfun(h, [0 2*max(dd)], 'splitting', 'on');
% p = polyfit(p0,40);
cheb_poly_coeffs = chebcoeffs(p);                         % get chebyshev coefficients of polynomial approx.

% % for fair comparison, equate number of free params
% % a degree n rational has 2n params, deg n poly has n+1
% % cheb polynomials are orthogonal so can truncate coeffs freely
deg = 2*length(pol)-1;

cc = flip(cheb_poly_coeffs);
cc_deg = flip(cheb_poly_coeffs(1:deg+1));

p_full = @(x) chebpolyval(cc, x);
p = @(x) chebpolyval(cc_deg, x);

%%

tt = linspace(-1,1,1000).';
hh = h_generic(tt);
hr = rpf(tt);
hp = p(tt);
hpfull = p_full(tt);
figure();
subplot(211); hold on;
plot(tt,hh);
plot(tt,hr);
plot(tt,hp);
plot(tt,hpfull);
legend(["exact" "rational" "polynomial" "full polynomial"])

subplot(212); hold on;
semilogy(tt,abs(hh-hr));
semilogy(tt,abs(hh-hp));
semilogy(tt,abs(hh-hpfull));
legend(["rational" "polynomial" "full polynomial"])

%%

save('approxParams.mat', 'pol', 'res', 'cheb_poly_coeffs')