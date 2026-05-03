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

% load pol, res, cheb_poly_coeffs from BPF with window [-1 1]
load approxParams

%% Form G and L
n = 500;

% A = adj_Kn(n);
A = adj_undirCycle(n);
% A = adj_dirCycle(n);

% degree vector
dd = full(sum(A,2));
dmax = max(dd);

L = laplacian(A);

%% Declare h

% LPF example
cutoff_lo = 1;
cutoff_hi = 2;  % UB on spectrum is 2*dmax

h = idealBPF(cutoff_lo,cutoff_hi);

% approximate p
% p = chebfun(h, [0,2*dmax], 2*length(pol)+1);
p = chebfun(h, [0,2*dmax]);
cheb_poly_coeffs = chebcoeffs(p);  

alpha_r = (2/3)/(cutoff_hi-cutoff_lo);
beta_r = -1/2 * (cutoff_lo+cutoff_hi) * alpha_r;

alpha_p = 2/(2*dmax-0);
beta_p = -1/2 * (2*dmax) * alpha_p;

pol_r = (pol - beta_r)./alpha_r;
res_r = res/alpha_r;
L_r = alpha_r * L + beta_r*eye(n);
L_p = alpha_p * L + beta_p*eye(n);

r = @(zz) real(pfeval(zz, pol_r, res_r));
rL = @(x) applyrLx(x,L_r,pol_r,res_r);

poly_deg = 2*length(pol)-1;
cheb_poly_coeffs_deg = cheb_poly_coeffs(1:poly_deg+1);

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


