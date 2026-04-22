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

% x = 1/n*ones(n,1);
% x = sin(4*pi*(0:(n-1))/n).^4.';
% x = tanh(sin(2*pi*40.1*(1:n)/n)).';
% w = 1-1/n*(1:n).'
% w = 1./(1:n).';
% w = exp(-(1:n).');
% w = 1 - (1/n*(1:n).').^2;
w = 1 - (1/n*(0:(n-1)).').^(1/2);
% w = 1/n*ones(n,1);
x = V * w;
y = L * x;

%% 1. Declare ideal filter h(lambda)

cutoff = 1;
% h = idealLPF(cutoff);
% cutoff_idx = find(lam > cutoff,1);
% h = idealHPF(cutoff);
% cutoff_idx = find(lam < cutoff,1);
h = idealBPF(cutoff,cutoff+1);
cutoff_idx = find(lam < cutoff,1);

%% 2. Compute L=VDV^{-1} and exact application y = Vh(D)V^{-1}x

Lh = V * diag(h(lam)) * V.';
yh = Lh * x;

VD = V * diag(lam);
VhD = V * diag(h(lam));

%% 3a. Compute rational approximation r(lambda) ~= h(lambda)

% generic upper bound for lambda is 2*max(degree)
m = 2000;
zz = linspace(0,2*max(dd),m-1).';
% zz = [zz; cutoff];
m = length(zz);

% λ(t)=λ∗+(λmax​−λ∗)⋅sign(t)∣t∣1/2,t∈[−1,1]

% tt = linspace(-1,1,m-1);
% zz = (cutoff + (2*max(dd)-cutoff) .* sign(tt) .* sqrt(abs(tt))).';

hh = h(zz);

tol = 1e-13;
% [r, pol, res, zer, zj, fj, wj, errvec] = aaa(hh, zz, 'tol', tol, 'lawson', 0, 'cleanup', 0);
[r, pol, res, zer, zj, fj, wj, errvec] = aaa(hh, zz, 'tol', tol);
% I = find(wj == 0);
% zj(I) = []; wj(I) = []; fj(I) = [];
[pol, res, zer] = prz(zj, fj, wj);

I = find(imag(pol) == 0);
pol(I) = []; res(I) = [];

% AA = [ones(m,1) 1./(zz - pol.')];
AA = [1./(zz - pol.')];
bb = hh;
aa = AA\bb;

% rpf = @(zz) pfeval(zz, pol, aa(2:end), a0=aa(1));
rpf = @(zz) pfeval(zz, pol, aa(1:end));

%% 3b. Compute approximate application y ~= r(L)x

% yr = applyrLx(x,L,pol,aa(2:end),a0=aa(1));
yr = applyrLx(x,L,pol,aa);

%% 4a. Compute polynomial approximation p(lambda) ~= h(lambda)

p = chebfun(h, [0 2*max(dd)]);
% p0 = chebfun(h, [0 2*max(dd)], 'splitting', 'on');
% p = polyfit(p0,40);
cc = chebcoeffs(p);

%% 4b. Compute approximate application y ~= p(L)x

% deg = length(pol);
deg = 2*length(pol);

Ls = (2*L - (2*max(dd))*eye(n))/(2*max(dd));

vprev = x;
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
    subplot(211); hold on;
    plot(yh,"-",LW,lw);
    plot(yr,"-",LW,lw);
    plot(yp,"-",LW,lw);
    xlim([0 n])
    ylabel("Output")
    legend(["Exact" "Rational" "Polynomial"])
    
    subplot(212); hold on;
    semilogy(abs(yh-yr),"-",LW,lw);
    semilogy(abs(yh-yp),"-",LW,lw);
    xlim([0 n])
    ylabel("Error")
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