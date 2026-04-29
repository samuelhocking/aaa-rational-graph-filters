% CS 2252 Final Project

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

plottingflag = 1;

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

%% 0. Construct G and L

n = 300;

% A = adj_Kn(n);
A = adj_undirCycle(n);
% A = adj_dirCycle(n);

% degree vector
dd = full(sum(A,2));

L = laplacian(A);
[V,lam] = eig(full(L),"vector");

% set signal either directly (x)
% x = 1/n*ones(n,1);                    % uniform signal

% or set signal indirectly by structuring the eigenvector expansion coeffs
% w = 1-1/n*(1:n).'                     % linearly decaying mode weights
% w = 1./(1:n).';                       % decay like 1/n
% w = exp(-(1:n).');                    % exponential decay
% w = 1 - (1/n*(1:n).').^2;             % quadratic decay
w = 1 - (1/n*(0:(n-1)).').^(1/2);       % square root decay
% w = 1/n*ones(n,1);                    % uniform spectrum
x = V * w;                              % assemble expansion
y = L * x;                              % compute unfiltered output

%% 1. Declare ideal filter h(lambda)

cutoff = 1;         
h = idealLPF(cutoff);
cutoff_idx = find(lam > cutoff,1);

% h = idealHPF(cutoff);
% cutoff_idx = find(lam < cutoff,1);

% h = idealBPF(cutoff,cutoff+1);
% cutoff_idx = find(lam < cutoff,1);

%% 2. Compute L=VDV^{-1} and exact application y = Vh(D)V^{-1}x

Lh = V * diag(h(lam)) * V.';            % compute exactly filtered Laplacian h(L)
yh = Lh * x;                            % compute exactly filtered signal

VD = V * diag(lam);                     % compute weighted eigenvectors
VhD = V * diag(h(lam));                 % compute filter weighted eigenvectors

%% 3a. Compute rational approximation r(lambda) ~= h(lambda)

% generic upper bound for lambda is 2*max(degree)
m = 1000;
% zz = linspace(0,2*max(dd),m-1).';         % linearly spaced points

zz = [                                      % two intervals of exponentially clustered points near cutoff
    expClust(m,a=0,b=cutoff,reverse=true,exponent=1/2);
    expClust(m,a=cutoff,b=2*max(dd),exponent=1/2)
];

m = length(zz);

hh = h(zz);         % compute reference data for filtered sample points

tol = 1e-13;
[r, pol, res, zer, zj, fj, wj, errvec] = aaa(hh, zz, 'tol', tol);   % AAA algorithm
% I = find(wj == 0);
% zj(I) = []; wj(I) = []; fj(I) = [];
[pol, res, zer] = prz(zj, fj, wj);      % compute poles, residues, and zeros of r

% I = find(imag(pol) == 0);                                                 % remove all real poles
I = find((imag(pol) == 0) .* (real(pol) >= 0) .* (real(pol) <= 2*max(dd)))  % remove real poles in approximation domain
pol(I) = []; res(I) = [];

% AA = [ones(m,1) 1./(zz - pol.')];
AA = [1./(zz - pol.')];             % set up least-squares matrix for preset poles
bb = hh;                            % RHS is filtered sample points
aa = AA\bb;                         % compute LS residues

% rpf = @(zz) pfeval(zz, pol, aa(2:end), a0=aa(1));
rpf = @(zz) pfeval(zz, pol, aa(1:end));     % get function handle for partial fraction rational filter

%% 3b. Compute approximate application y ~= r(L)x

% yr = applyrLx(x,L,pol,aa(2:end),a0=aa(1));
yr = applyrLx(x,L,pol,aa);                  % compute r(L)

%% 4a. Compute polynomial approximation p(lambda) ~= h(lambda)

p = chebfun(h, [0 2*max(dd)]);              % compute polynomial approximation to filter
% p0 = chebfun(h, [0 2*max(dd)], 'splitting', 'on');
% p = polyfit(p0,40);
cc = chebcoeffs(p);                         % get chebyshev coefficients of polynomial approx.

%% 4b. Compute approximate application y ~= p(L)x

% deg = length(cc)-1;
% deg = length(pol)-1;
deg = 2*length(pol)-1;                      % for fair comparison, equate number of free params
                                            % a degree n rational has 2n params, deg n poly has n+1
                                            % cheb polynomials are orthogonal so can truncate coeffs freely

Ls = (2*L - (2*max(dd))*eye(n))/(2*max(dd));    % scaled Laplacian for numerical stability

vprev = x;                                      % evaluate p(L) by chebyshev recurrence
vcurr = Ls*x;
yp = cc(1)*vprev + cc(2)*vcurr;
for j = 3:(deg+1)
    vnext = 2*Ls*vcurr - vprev;
    yp = yp + cc(j)*vnext;
    vprev = vcurr;
    vcurr = vnext;
end

%% 5. Plotting
% fig1: 

%%
    f1 = figure();
    subplot(311); hold on;
    scatter(1:n, lam, ms, LW, lw)
    scatter(1:n, h(lam), 0.5*ms, 'filled', LW, lw)
    xline(cutoff_idx-1/2,":r",LW,lw)
    xlabel("$j$")
    legend(["$\lambda_j$" "$h(\lambda_j)$" "cutoff"],"Location","best")

    subplot(312); hold on;
    semilogy(1:n, abs(V.' * x), LW, lw)
    % scatter(1:n, h(lam), 0.5*ms, 'filled', LW, lw)
    % xline(cutoff_idx-1/2,":r",LW,lw)
    xlabel("$j$")
    ylabel("$V^Tx$")
    % legend(["$\lambda_j$" "$h(\lambda_j)$" "cutoff"],"Location","best")
    
    subplot(313); hold on;
    yyaxis left
    plot(x,LW,lw);
    ylabel("$x$")
    yyaxis right; hold on;
    plot(y,"-",LW,lw);
    plot(yh,".-",LW,lw);
    ylabel("$y$")
    xlabel("$j$")
    legend(["$x$" "$y$" "$y_h$"],"Location","best")
    %%
    figure();
    subplot(311); hold on;
    plot(yh,"-",LW,lw);
    plot(yr,"-",LW,lw);
    plot(yp,"-",LW,lw);
    xlim([0 n])
    ylabel("Output")
    legend(["Exact" "Rational" "Polynomial"])
    
    subplot(312); hold on;
    semilogy(abs(yh-yr),"-",LW,lw);
    semilogy(abs(yh-yp),"-",LW,lw);
    xlim([0 n])
    ylabel("Error")
    legend(["Rational" "Polynomial"])

    subplot(313); hold on;
    plot(real(yh-yr),"-",LW,lw);
    plot(real(yh-yp),"-",LW,lw);
    xlim([0 n])
    ylabel("Signed error")
    legend(["Rational" "Polynomial"])
    
    %%
    figure()
    subplot(211); hold on;
    plot(lam,h(lam),LW,lw);
    plot(lam,rpf(lam),LW,lw);
    plot(lam,p(lam),LW,lw);
    xlim([0 2*max(dd)])
    xlabel("$\lambda$")
    ylabel("Response")
    legend(["Exact" "Rational" "Polynomial"])
    
    subplot(212); hold on;
    semilogy(lam,abs(h(lam)-rpf(lam)),LW,lw);
    semilogy(lam,abs(h(lam)-p(lam)),LW,lw);
    xlim([0 2*max(dd)])
    xlabel("$\lambda$")
    ylabel("Error")
    legend(["Rational" "Polynomial"])

%%
if plottingflag
    f3 = figure();
    tcl = tiledlayout(1,2, "TileSpacing", "tight", "Padding", "loose");
    ax1 = nexttile; imagesc(V * diag(lam)); axis square
    ax2 = nexttile; imagesc(V * diag(h(lam))); axis square

    % Add shared colorbar
    clim([min([VD VhD],[],"all") max([VD VhD],[],"all")])
    cb = colorbar(ax2);
    cb.Layout.Tile = 'east';
    linkaxes([ax1, ax2], 'xy');
    fontsize(fs, "points")

    f4 = figure();
    tcl = tiledlayout(1,2, "TileSpacing", "tight", "Padding", "loose");
    ax3 = nexttile; imagesc(L); axis square
    ax4 = nexttile; imagesc(Lh); axis square

    % Add shared colorbar
    clim([min([L Lh],[],"all") max([L Lh],[],"all")])
    cb = colorbar(ax4);
    cb.Layout.Tile = 'east';
    linkaxes([ax3, ax4], 'xy');
    fontsize(fs, "points")

    figure();
    hold on;
    scatter(zz, 0*zz, 5, 'filled');
    scatter(real(pol), imag(pol), ms);

end