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

cutoff_lo = -1/3;
cutoff_hi = 1/3;
h_generic = idealBPF(cutoff_lo,cutoff_hi);

f1 = figure();
f1.Position(4) = 200;

subplot(131)
tt = linspace(-1,1,1000).';
plot(tt, h_generic(tt), LW,lw)
ylim([-0.1 1.1])
grid on
xlabel("$\tilde{\lambda}$")
ylabel("$h(\tilde{\lambda})$")
title("Parent indicator function")


subplot(132)
h_LPF = idealBPF(0,1);
tt = linspace(0,4,1000).';
plot(tt, h_LPF(tt), LW,lw)
grid on
% axis equal
ylim([-0.1 1.1])
xlabel("$\lambda$")
ylabel("$h(\lambda)$")
title("Ideal low-pass filter on $[0,1]$")

subplot(133)
h_HPF = idealBPF(2.5,4);
tt = linspace(0,4,1000).';
plot(tt, h_HPF(tt), LW,lw)
grid on
% axis equal
ylim([-0.1 1.1])
xlabel("$\lambda$")
ylabel("$h(\lambda)$")
title("Ideal high-pass filter on $[2.5,4]$")
fontsize(fs,"points")

exportgraphics(f1, "../figures/indicator_gallery.png", 'Resolution', 300);
% exportgraphics(f1, "../figures/indicator_gallery.png", 'Resolution', 300, 'Padding', 'figure');