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

%% 1. Declare ideal filter h(lambda)

% generic BPF     
cutoff_lo = -1/3;
cutoff_hi = 1/3;
h_generic = idealBPF(cutoff_lo,cutoff_hi);

%% 2. Compute rational approximation r(lambda) ~= h(lambda)

% generic upper bound for lambda is 2*max(degree)
m_equi = 50000;
zz = linspace(-1,1,m_equi).';                     % basic equispaced grid
eps = 1e-12;                                % define perturbation
zz = [zz; cutoff_lo; cutoff_hi; cutoff_lo-eps; cutoff_lo+eps; cutoff_hi-eps; cutoff_hi+eps];    % manually add points very close to the jumps
zz = unique(zz);                                    % remove duplicates
m = length(zz);

hh = h_generic(zz);         % compute reference data for filtered sample points

tol = 1e-13;
[r, pol, res, zer, zj, fj, wj, errvec] = aaa(hh, zz, 'tol', tol);   % AAA algorithm

I = find(imag(pol) == 0);   % find indices of real poles
pol(I) = []; res(I) = [];   % remove real poles

AA = [1./(zz - pol.')];             % set up least-squares matrix for preset poles
bb = hh;                            % RHS is filtered sample points

res = AA\bb;                         % compute LS residues

rpf = @(zz) real(pfeval(zz, pol, res));     % get function handle for partial fraction rational filter

%%

tt = zz;
hh = h_generic(tt);
hr = rpf(tt);
figure();
subplot(211); hold on;
plot(tt,hh);
plot(tt,hr);
legend(["exact" "rational"])

subplot(212); hold on;
semilogy(tt,abs(hh-hr));
legend(["rational"])

%%

save('approxParams.mat', 'pol', 'res')