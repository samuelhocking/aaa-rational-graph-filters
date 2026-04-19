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

% 1. Declare ideal filter h(lambda)
% 2. Compute L=VDV^{-1} and exact application y = Vh(D)V^{-1}x
% 3a. Compute rational approximation r(lambda) ~= h(lambda)
% 3b. Compute approximate application y ~= r(L)x
% 4a. Compute polynomial approximation p(lambda) ~= h(lambda)
% 4b. Compute approximate application y ~= p(L)x

%%