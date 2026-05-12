% CS 2252 Final Project

%% Initial Setup

clear all
close all

set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter', 'latex')

colors = get(groot,'defaultAxesColorOrder');

LW = "LineWidth"; lw = 1.2;
MS = "MarkerSize"; ms = 15;
fs = 12;
cmap = 'parula';

warning('off', 'MATLAB:legend:IgnoringExtraEntries');

%%

% load pol, res from generic BPF
load approxParams

%% Form G and L
n = 500;

A = adj_undirCycle(n);

% degree vector
dd = full(sum(A,2));
dmax = max(dd);

L = laplacian(A);

%% Declare h

% LPF example
cutoff_lo = 3;
cutoff_hi = 4;  % UB on spectrum is 2*dmax

h = idealBPF(cutoff_lo,cutoff_hi);

% approximate p
pcf = chebfun(h, [0,2*dmax]);
cheb_poly_coeffs = chebcoeffs(pcf);  

alpha_r = (2/3)/(cutoff_hi-cutoff_lo);
beta_r = -1/2 * (cutoff_lo+cutoff_hi) * alpha_r;

alpha_p = 2/(2*dmax-0);
beta_p = -1/2 * (2*dmax) * alpha_p;

pol_r = (pol - beta_r)./alpha_r;
res_r = res/alpha_r;

r = @(zz) real(pfeval(zz, pol_r, res_r));
rL = @(x) real(applyrLx(x,L,pol_r,res_r));

poly_deg = 2*length(pol)-1;
cheb_poly_coeffs_deg = cheb_poly_coeffs(1:poly_deg+1);
L_p = alpha_p * L + beta_p*eye(n);

p = @(zz) chebpolyval(flip(cheb_poly_coeffs_deg), alpha_p*zz + beta_p);
pL = @(x) applypLx(x,L_p,cheb_poly_coeffs_deg);

%%

tt = linspace(0,2*dmax,5000).';
hh = h(tt);
hr = r(tt);
hp = p(tt);

figure();
subplot(211); hold on;
plot(tt,hh);
plot(tt,hr);
plot(tt,hp);
legend(["exact" "rational" "polynomial"])

subplot(212); hold on;
plot(tt,abs(hh-hr));
plot(tt,abs(hh-hp));
yscale log
legend(["rational" "polynomial"])

%%
[V,lam] = eig(full(L),"vector");

% set signal either directly (x)
% x = 1/n*ones(n,1);                    % uniform signal

% or set signal indirectly by structuring the eigenvector expansion coeffs
% w = 1-1/n*(0:n-1).';                     % linearly decaying mode weights
% w = 1./(1:n).';                       % decay like 1/n
% w = (1./(1:n).').^2;                       % decay like 1/n^2
% w = exp(-(1:n).');                    % exponential decay
% w = 1 - (1/n*(1:n).').^2;             % quadratic decay
% w = 1 - (1/n*(0:(n-1)).').^(1/2);       % square root decay
% w = 1/n*ones(n,1);                    % uniform spectrum
w = ones(n,1);                    % uniform spectrum
x = V * w;                              % assemble expansion
y = L * x;                              % compute unfiltered output

Lh = V * diag(h(lam)) * V.';            % compute exactly filtered Laplacian h(L)
hDVT = diag(h(lam)) * V.';
yh = Lh * x;                            % compute exactly filtered signal

yr = rL(x);
yp = pL(x);

%%

f2 = figure();
f2.Position(4) = 600;

tl = tiledlayout(f2,3,1,"TileSpacing","compact","Padding","loose");

ax1 = nexttile(tl); hold on;

% plot(x, LW, lw);
scatter(lam, V.' * x, ms, LW, lw);
scatter(lam, hDVT * x, 0.5*ms, 'filled', LW, lw);
% ax1.XTick = [];
% ax1.XTickLabel = [];
xlabel("Spectrum, $\lambda$")
% ylabel("$V^Tx$")
legend(["$V^Tx$" "$h(D)V^Tx$"],"Location","best","NumColumns",2)
grid on;

ax2 = nexttile(tl); hold on;
plot(yh,LW,lw)
plot(yr,LW,lw)
plot(yp,LW,lw)
legend(["Exact" "Rational" "Polynomial"],"Location","best","NumColumns",3)
ylabel("Response")
% xlabel("Coordinate index")
% ax2.XTick = [];
ax2.XTickLabel = [];
grid on;

ax3 = nexttile(tl); hold on;
plot(abs(yh-yr),LW,lw)
plot(abs(yh-yp),LW,lw)
yscale log
legend(["Rational" "Polynomial"],"Location","best","NumColumns",2)
ylabel("Error")
xlabel("Coordinate index")
grid on
ax3.YMinorGrid = 'off';

fontsize(fs,"points")
linkaxes([ax2 ax3],"x")

% exportgraphics(f2, "./figures/exp_LPF_lindecay.png", 'Resolution', 300);
% exportgraphics(f2, "./figures/exp_BPF_unif.png", 'Resolution', 300);
% exportgraphics(f2, "./figures/exp_HPF_unif.png", 'Resolution', 300);
