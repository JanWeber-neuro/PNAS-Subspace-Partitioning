%% - FIGURE 4C - %%
% This code plots the correlation between subspace alignment and behavioral
% switch cost as in Fig. 4C.
%
%%% --- REQUIRED INPUTS --- %%%
% - subspace_alignment.mat
% - behav_switch_cost_rt.mat

%%

clear
ft_defaults;

%% - LOAD DATA - %%

% load cumulative explained variance for context and history subspace
load('subspace_alignment.mat');

% load switch cost
load('behav_switch_cost_rt.mat');

%% - PLOT CORRELATION SUBSPACE ALIGNMENT AND BEHAVIOR - %%

figure;

% get current behavioral vector
tmpbehav = DiffRT;
        
% colormap to use
cmap = cbrewer('seq', 'YlOrRd', length(all_rho), 'PCHIP');

% sort behavior for color coding
[behav_sort,sortidx] = sort(tmpbehav, 'ascend');
         
% scatterplot %EV context vs. %EV history
scatter(behav_sort, all_rho(sortidx), 1500, cmap, ...
    'Marker', '.', 'MarkerFaceAlpha', 1, 'MarkerEdgeAlpha', 1);

hold on

% fit uncertainty estimates
mdl = fitlm(behav_sort, all_rho(sortidx));
xpred = linspace(min(behav_sort), max(all_rho(sortidx)), 1000)';
[ypred, yci] = predict(mdl, xpred);

opts={'EdgeColor', 'none',...
      'FaceColor', 'k',...
      'FaceAlpha',0.2};

fill_between(xpred, yci(:,1), yci(:,2),1,opts{:});
hold on
plot(xpred, ypred, 'k');
        
ylabel('subspace alignment');
xlabel('switch cost');

ylim([-0.7 1]);
xlim([-0.005 0.042]);

set(gca, 'fontsize', 13, 'linewidth', 1.5);
box off

set(gcf, 'Position', [480 425 232 183]);        

% - COMPUTE THE CORRELATION COEFFICIENT  - %

% fit uncertainty estimates
[r,p] = corr(behav_sort, all_rho(sortidx), 'rows', 'complete', 'type', 'Spearman');

n = size(behav_sort,1) - sum(isnan(all_rho(sortidx)));

% step 1: perform fisher z
zr = log((1+r)/(1-r)) / 2;

% step 2: find upper and lower bounds
LB = zr - (1.96/sqrt(n-3));

UB = zr + (1.96/sqrt(n-3));

% step 3: get confidence estimates
CI = [(exp(2*LB)-1) / (exp(2*LB)+1); (exp(2*UB)-1) / (exp(2*UB)+1)];

fprintf('\nrho = %.3f, p = %.3f, CI = [%.3f %.3f]\n', r, p, CI(1), CI(2));

