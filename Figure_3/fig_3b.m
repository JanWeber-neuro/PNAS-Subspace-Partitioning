%% - FIGURE 3B - %%
% This code plots the cumulative explained variance for history and context
% subspaces as in Fig. 3B.
%
%%% --- REQUIRED INPUTS --- %%%
% - subspace_dimensionality.mat

%%

clear
ft_defaults;

%% - LOAD DATA - %%

% load cumulative explained variance for context and history subspace
load('subspace_dimensionality.mat');

%% - PLOT CUMULATIVE SUM FOR CONTEXT AND HISTORY DIMENSIONS - %%

nPCs = size(ev_context,2);

figure;

% colors
cmap       = viridis(20);
color{1}   = cmap(3,:);
color{2}   = cmap(15,:);

% - [PLOT CONTEXT] - %
a = plot(1:nPCs, ev_context', 'Color', color{1}, 'LineWidth', 2);
for jj = 1:size(ev_context,1)
    a(jj).Color(4) = 0.2;
end
hold on
scatter(1:nPCs, mean(ev_context,1), 80, 'filled', 'MarkerFaceColor', color{1}, 'MarkerEdgeColor', 'none');

hold on

% - [PLOT HISTORY] - %
a = plot(1:nPCs, ev_history', 'Color', color{2}, 'LineWidth', 2);
for jj = 1:size(ev_history,1)
    a(jj).Color(4) = 0.2;
end
hold on
scatter(1:nPCs, mean(ev_history,1), 80, 'filled', 'MarkerFaceColor', color{2}, 'MarkerEdgeColor', 'none');

% figure settings
set(gca, 'fontsize', 13, 'linewidth', 1.5);

xticks(1:nPCs);
xticklabels({'PC1', 'PC2', 'PC3', 'PC4', 'PC5'});
xtickangle(45);
xlim([0 nPCs+1]);
ylabel('Cumulative %EV');

box off

set(gcf, 'Position',[520 382 294 248]);

