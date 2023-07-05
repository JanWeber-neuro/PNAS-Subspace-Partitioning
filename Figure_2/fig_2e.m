%% - FIGURE 2E - %%
% This code plots the correlation traces between neural information and
% behavioral switch cost as in Fig. 2E.
%
%%% --- REQUIRED INPUTS --- %%%
% - w2_roi.mat
% - bootstrapped_corr.mat

%%

clear
close all;
clc;

%% - MAIN SETTINGS - %%

% load w2 data within ROIs
load("w2_roi.mat")
nseg     = numel(fieldnames(w2));
segNames = fieldnames(w2);

% supporting vars
factorNames = {'Cue', 'Past Trial Type'};
nFactor     = numel(factorNames);

% load bootstrapped correlation coefficients
load('bootstrapped_corr.mat');

%% - CONTRAST CORRELATION COEFFICIENTS BETWEEN CONTEXT AND HISTORY FOR PFC - %%

% extract the number of participants with neural information in region x
% about variable y
nsub_context = [];
nsub_history = [];

nsub_context.pfc = size(squeeze(w2.br.Prefrontal(1,:,:)),2) - ...
    sum(all(isnan(squeeze(w2.br.Prefrontal(1,:,:)))),2);
nsub_history.pfc = size(squeeze(w2.br.Prefrontal(2,:,:)),2) - ...
    sum(all(isnan(squeeze(w2.br.Prefrontal(2,:,:)))),2);

% - [ACCURACY] - %

% test for significant differences in correlation coefficients over time
pvals_pfc_acc = NaN(1, length(w2.br.time));

for Itime = 1:length(pvals_pfc_acc)

    pvals_pfc_acc(Itime) = compare_correlation_coefficients(boot_rho.acc.br.Prefrontal.Context(Itime), ...
        boot_rho.acc.br.Prefrontal.History(Itime), nsub_context.pfc, nsub_history.pfc);

end%Itime

[~, ~, ~, pvals_pfc_acc] = fdr_bh(pvals_pfc_acc);

%% - CONTRAST CORRELATION COEFFICIENTS BETWEEN CONTEXT AND HISTORY FOR MOTOR CORTEX - %%

% extract the number of participants with neural information in region x
% about variable y
nsub_context.m1 = size(squeeze(w2.br.Motor(1,:,:)),2) - ...
    sum(all(isnan(squeeze(w2.br.Motor(1,:,:)))),2);
nsub_history.m1 = size(squeeze(w2.br.Motor(2,:,:)),2) - ...
    sum(all(isnan(squeeze(w2.br.Motor(2,:,:)))),2);

% - [ACCURACY] - %

% test for significant differences in correlation coefficients over time
pvals_m1_acc = NaN(1, length(w2.br.time));

for Itime = 1:length(pvals_m1_acc)

    pvals_m1_acc(Itime) = compare_correlation_coefficients(boot_rho.acc.br.Motor.Context(Itime), ...
        boot_rho.acc.br.Motor.History(Itime), nsub_context.m1, nsub_history.m1);

end%Itime

[~, ~, ~, pvals_m1_acc] = fdr_bh(pvals_m1_acc);

%% - PLOT CORRELATION TRACES - %%

% region of interests
roi2use = {'Prefrontal', 'Motor'};

% location for scatters
locscatter = [-1.2 -1.1];

% factor names
fcnames = fieldnames(boot_rho.acc.br.Motor);

% colors
cmap       = viridis(20);
color{1}   = cmap(3,:);
color{2}   = cmap(15,:);

for Iroi = 1:numel(roi2use)
    
    % current ROI
    tmpROI = roi2use{Iroi};
  
    if strcmp(tmpROI, 'Prefrontal')
        sigCorrDiff.acc.pfc = find(pvals_pfc_acc < 0.05);
        tmpCorrDiff = sigCorrDiff.acc.pfc;
    else
        sigCorrDiff.acc.m1  = find(pvals_m1_acc < 0.05);
        tmpCorrDiff = sigCorrDiff.acc.m1;
    end
    
    figure;

    for Ifactor = 1:numel(fcnames)

        % get current factor name and index
        tmpfactor = fcnames{Ifactor};
  
        % get current mean correlation and confidence interval
        tmp_avg = boot_rho.acc.br.(tmpROI).(tmpfactor);
        tmp_ci  = [boot_ci.acc.br.(tmpROI).(tmpfactor)(1,:) - tmp_avg; ...
                   tmp_avg - boot_ci.acc.br.(tmpROI).(tmpfactor)(2,:)];
        
        a = shadedErrorBar(w2.br.time, tmp_avg, tmp_ci);
        a.mainLine.Color = color{Ifactor}; a.mainLine.LineWidth = 3;
        a.patch.FaceColor = color{Ifactor}; a.patch.EdgeColor = color{Ifactor}; a.patch.FaceAlpha = 0.5;
        a.edge(1).Color = color{Ifactor}; a.edge(1).Color(4) = 0.1; a.edge(1).LineWidth = eps;
        a.edge(2).Color = color{Ifactor}; a.edge(2).Color(4) = 0.1; a.edge(2).LineWidth = eps;
        
        hold on

        % - [MARK SIGNIFICANT CLUSTERS] - %
    
        % check if there is a negative cluster
        if isfield(true_clustermass.acc.br, tmpROI) && ...
                isfield(true_clustermass.acc.br.(tmpROI), tmpfactor) && ...
                    isfield(true_clustermass.acc.br.(tmpROI).(tmpfactor), 'neg')
            
            % compute p-val
            pval_neg = NaN(length(true_clustermass.acc.br.(tmpROI).(tmpfactor).neg),1);
            for jj = 1:size(pval_neg,1)
                pval_neg(jj) = sum(clustermass_surro.acc.br.(tmpROI).(tmpfactor).neg > ...
                    abs(true_clustermass.acc.br.(tmpROI).(tmpfactor).neg(jj))) / nshuffle;
            end
    
        end

        % check if there is a positive cluster
        if isfield(true_clustermass.acc.br, tmpROI) && ...
                    isfield(true_clustermass.acc.br.(tmpROI), tmpfactor) && ...
                    isfield(true_clustermass.acc.br.(tmpROI).(tmpfactor), 'pos')
            
            % compute p-val
            pval_pos = NaN(length(true_clustermass.acc.br.(tmpROI).(tmpfactor).pos),1);
            for jj = 1:size(pval_pos,1)
                pval_pos(jj) = sum(clustermass_surro.acc.br.(tmpROI).(tmpfactor).pos > ...
                    abs(true_clustermass.acc.br.(tmpROI).(tmpfactor).pos(jj))) / nshuffle;
            end
            
        end
        
        % mark significant timepoints from negative clusters
        if true_island_neg.acc.br.(tmpROI).(tmpfactor).NumObjects > 0
            for Icluster = 1:true_island_neg.acc.br.(tmpROI).(tmpfactor).NumObjects
                if pval_neg(Icluster) < 0.05
                    sigpnts = true_island_neg.acc.br.(tmpROI).(tmpfactor).PixelIdxList{Icluster};
                    plot(w2.br.time(sigpnts), ones(1, length(sigpnts)) * locscatter(Ifactor), ...
                       'Color', color{Ifactor}, 'LineWidth', 2);
                end
            end
        end
    
        % mark significant timepoints from positive clusters
        if true_island_pos.acc.br.(tmpROI).(tmpfactor).NumObjects > 0
            for Icluster = 1:true_island_pos.acc.br.(tmpROI).(tmpfactor).NumObjects
                if pval_pos(Icluster) < 0.05
                    sigpnts = true_island_pos.acc.br.(tmpROI).(tmpfactor).PixelIdxList{Icluster};
                    plot(w2.br.time(sigpnts), ones(1, length(sigpnts)) * locscatter(Ifactor), ...
                       'Color', color{Ifactor}, 'LineWidth', 2);
                end
            end
        end   
        
        % mark significant difference between correlations (context vs. history)
        plot(w2.br.time(tmpCorrDiff), ones(1, length(tmpCorrDiff)) * -0.9, ...
           'Color', 'k', 'LineWidth', 2);

    end%Ifactor

    % figure settings
    xlabel('Time [s]');
    ylabel('rho [acc]');

    xlim([w2.br.time(1) w2.br.time(end)]);
    ylim([-1.3 1]);

    set(gca, 'fontsize', 13, 'linewidth', 2);
    box off
    set(gcf, 'Position', [523 453 230 191]);
   
end%Iroi


