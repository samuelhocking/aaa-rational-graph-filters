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

%% 0.

n = 100;

% A = adj_Kn(n);
A = adj_undirCycle(n);
% A = adj_dirCycle(n);

L = laplacian(A);
[V,lam] = eig(full(L),"vector");

f1 = figure(1); hold on;
plot(lam, '.', MS, ms);

x = randn(n,1);
y = L * x;

%% 1.

h = idealLPF(0.5);
yh = V * diag(h(lam)) * V.' * x;

f2 = figure(2); hold on;
plot(x,LW,lw);
plot(y,LW,lw);
plot(yh,LW,lw);
legend(["x" "y" "y_h"])