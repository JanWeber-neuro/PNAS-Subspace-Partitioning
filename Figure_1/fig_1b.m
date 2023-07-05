%% - FIGURE 1B RT & ACCURACY - %%
% plots group level reaction time and accuracy for switch vs. no-switch
% trials.
%
%%% --- REQUIRED INPUTS --- %%%
% - behav_rt_acc.mat 

%% - HOUSE KEEPING - %%

clear
close all;
clc;

%% - LOAD DATA - %%

% load data related to figure
load behav_rt_acc.mat;

%% - PLOT REACTION TIME - %%

cmap = viridis(100);

figure;

color = [cmap(80,:); cmap(10,:)];

raincloud_pairplot(RT, color);

set(gca, 'ytick', []);
xlabel('RT [s]');

h = gca; h.YAxis.Visible = 'off';

set(gca, 'fontsize', 13);

box off

%% - PLOT ACCURACY - %%

figure;

color = [cmap(80,:); cmap(10,:)];

raincloud_pairplot(Acc, color);

set(gca, 'ytick', []);
xlabel('Accuracy');

h = gca; h.YAxis.Visible = 'off';

set(gca, 'fontsize', 13);

box off

xlim([0.3 1.2]);

