%% - FIGURE 1B D-PRIME - %%
% computes d-prime for switch vs. no-switch trials.
%
%%% --- REQUIRED INPUTS --- %%%
% - dprime.mat 

%% - HOUSE KEEPING - %%

clear
close all;
clc;

%% - LOAD DATA - %%

% load data related to figure
load dprime.mat

%% - PLOT D-PRIME - %%

figure;

cmap = viridis(100);

color = [cmap(80,:); cmap(10,:)];

raincloud_pairplot(dprime, color);

set(gca, 'ytick', []);
xlabel('dprime');

h = gca; h.YAxis.Visible = 'off';

set(gca, 'fontsize', 13);

box off
