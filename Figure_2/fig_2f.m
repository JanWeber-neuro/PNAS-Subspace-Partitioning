%% - FIGURE 2F - %%
% This code plots the correlation traces between context- and
% history-dependent neural information for prefrontal and motor cortex as
% in Fig. 2F.
%
%%% --- REQUIRED INPUTS --- %%%
% - w2_roi.mat

%%

clear
close all;
clc

%% - LOAD DATA - %%

% load omega squared data
load("w2_roi.mat");

%% - GET TRACES FOR CONTEXT AND PAST TRIAL TYPE - %%

% variables of interest
vars      = {'Cue', 'Past Trial Type'};

% get ROIs of interest
roi2use   = {'Prefrontal', 'Motor'};

% critical alpha for thresholding
critalpha = 0.05;

% get index of variable
varsidx = find(ismember(w2.br.factorNames, vars));

% preallocate memory for correlation coefficient and pvalues
true_rho  = struct();
true_pval = struct();

% initialize true clustermass
true_clustermass = struct();

% initialize surrogate clustermass
clustermass_surro = struct();

for Iroi = 1:numel(roi2use)

    % get current ROI
    tmpROI = roi2use{Iroi};

    % get the traces for the variables
    trace_1 = squeeze(w2.br.(tmpROI)(varsidx(1),:,:))';
    trace_2 = squeeze(w2.br.(tmpROI)(varsidx(2),:,:))';
      
    for Itime = 1:length(w2.br.time)
    
        [r,p] = corr(trace_1(:,Itime), trace_2(:,Itime), 'rows', 'complete', 'type', 'Spearman');
        true_rho.(tmpROI)(Itime)  = r;
        true_pval.(tmpROI)(Itime) = p;
    
    end

    % - [COMPUTE POSITIVE CLUSTER] - %

    % pass variable
    tmp_rho = true_rho.(tmpROI);

    % threshold at critical alpha
    tmp_rho(true_pval.(tmpROI) > critalpha) = 0;

    % set positive r values to zero for two-sided testing
    tmp_rho(tmp_rho < 0) = 0;

    true_island_pos.(tmpROI)      = bwconncomp(tmp_rho);

    true_labelisland_pos.(tmpROI) = labelmatrix(true_island_pos.(tmpROI));
    
    % get the clustermass
    if true_island_pos.(tmpROI).NumObjects > 0

        true_clustermass.(tmpROI).pos = NaN(true_island_pos.(tmpROI).NumObjects,1);

        for Icluster = 1:true_island_pos.(tmpROI).NumObjects
            
            % extract clustermass
            true_clustermass.(tmpROI).pos(Icluster) = ...
                sum(abs(tmp_rho(true_island_pos.(tmpROI).PixelIdxList{Icluster})));

        end%Icluster

    end

    % - [COMPUTE NEGATIVE CLUSTER] - %

    % pass variable
    tmp_rho = true_rho.(tmpROI);

    % threshold at critical alpha
    tmp_rho(true_pval.(tmpROI) > critalpha) = 0;

    % set positive r values to zero for two-sided testing
    tmp_rho(tmp_rho > 0) = 0;

    true_island_neg.(tmpROI)      = bwconncomp(tmp_rho);

    true_labelisland_neg.(tmpROI) = labelmatrix(true_island_neg.(tmpROI));

    % get the clustermass
    if true_island_neg.(tmpROI).NumObjects > 0

        true_clustermass.(tmpROI).neg = NaN(true_island_neg.(tmpROI).NumObjects,1);

        for Icluster = 1:true_island_neg.(tmpROI).NumObjects
            
            % extract clustermass
            true_clustermass.(tmpROI).neg(Icluster) = ...
                sum(abs(tmp_rho(true_island_neg.(tmpROI).PixelIdxList{Icluster})));

        end%Icluster

    end
    
    %% - COMPUTE A SURROGATE DISTRIBITION USING BLOCK SWAPPING PROCEDURE - %%
    
    % number of random shuffles
    nshuffle = 1000;

    clustermass_surro.(tmpROI).pos = NaN(nshuffle,1);
    clustermass_surro.(tmpROI).neg = NaN(nshuffle,1);

    for Ishuffle = 1:nshuffle

        if mod(Ishuffle,100) == 0
            fprintf('-- running permutation %d/%d\n', Ishuffle, nshuffle);
        end
        
        % randomly cut the data at some point (block swapping)
        rndcut = randsample(length(w2.br.time),1);

        while rndcut == 1
            rndcut = randsample(length(w2.br.time),1);
        end
            
        % compute the correlation over time
        surro_rho  = NaN(1, length(w2.br.time));
        surro_pval = NaN(1, length(w2.br.time));

        % perform block swapping for trace 1
        tmp_trace_1 = [trace_1(:,rndcut+1:end) trace_1(:,1:rndcut)];
        
        % corr over time
        for Itime = 1:length(w2.br.time)
        
            [r,p] = corr(tmp_trace_1(:,Itime), trace_2(:,Itime), 'rows', 'complete', 'type', 'Spearman');
            surro_rho(Itime)  = r;
            surro_pval(Itime) = p;
        
        end

        % - [COMPUTE POSITIVE CLUSTER] - %
    
        % pass variable
        tmp_rho = surro_rho;
    
        % threshold at critical alpha
        tmp_rho(surro_pval > critalpha) = 0;
    
        % set positive r values to zero for two-sided testing
        tmp_rho(tmp_rho < 0) = 0;
    
        surro_island_pos      = bwconncomp(tmp_rho);
    
        surro_labelisland_pos = labelmatrix(surro_island_pos);
            
        % get the clustermass
        if surro_island_pos.NumObjects > 0
    
            tmpsurro_cluster = NaN(surro_island_pos.NumObjects,1);
    
            for Icluster = 1:surro_island_pos.NumObjects
                
                % extract clustermass
                tmpsurro_cluster(Icluster) = ...
                    sum(abs(tmp_rho(surro_island_pos.PixelIdxList{Icluster})));
    
            end%Icluster

            % get maximum cluster
            clustermass_surro.(tmpROI).pos(Ishuffle) = max(tmpsurro_cluster);
    
        else

            clustermass_surro.(tmpROI).pos(Ishuffle) = 0;

        end

        % - [COMPUTE NEGATIVE CLUSTER] - %
    
        % pass variable
        tmp_rho = surro_rho;
    
        % threshold at critical alpha
        tmp_rho(surro_pval > critalpha) = 0;
    
        % set positive r values to zero for two-sided testing
        tmp_rho(tmp_rho > 0) = 0;
    
        surro_island_neg      = bwconncomp(tmp_rho);
    
        surro_labelisland_neg = labelmatrix(surro_island_neg);
            
        % get the clustermass
        if surro_island_neg.NumObjects > 0
    
            tmpsurro_cluster = NaN(surro_island_neg.NumObjects,1);
    
            for Icluster = 1:surro_island_neg.NumObjects
                
                % extract clustermass
                tmpsurro_cluster(Icluster) = ...
                    sum(abs(tmp_rho(surro_island_neg.PixelIdxList{Icluster})));
    
            end%Icluster

            % get maximum cluster
            clustermass_surro.(tmpROI).neg(Ishuffle) = max(tmpsurro_cluster);
    
        else

            clustermass_surro.(tmpROI).neg(Ishuffle) = 0;

        end

    end%Ishuffle

end%Iroi
        
%% - PLOT TRACES - %%

% preallocate memory for bootstrapped correlations
boot_rho = struct();
boot_ci  = struct(); % confidence intervals

% number of iterations
niter = 500;

% sampling freq
fsample = 512;

% smoothing window
smoothwin = round(fsample*0.01);

fprintf('\n\nBootstrapping Correlations...\n');

for Iroi = 1:numel(roi2use)
    
    % get current ROI
    tmpROI = roi2use{Iroi};

    % get the traces for the variables
    trace_1 = squeeze(w2.br.(tmpROI)(varsidx(1),:,:))';
    trace_2 = squeeze(w2.br.(tmpROI)(varsidx(2),:,:))';
          
    % IDs to use (remove subjects with all NaNs)
    idx2use = setdiff(1:size(trace_1,1), find(all(isnan(trace_1),2)));
          
    % - [BOOTSTRAPPING WITH REPLACEMENT] - %

    % preallocate memory for bootstrapped correlations 
    tmp_rho = NaN(niter, size(trace_1,2));

    for Ishuffle = 1:niter

        % randomly sample with replacement
        bootsample = randsample(idx2use, floor(length(idx2use)*0.8), 'true');

        for Itime = 1:length(w2.br.time)
        
            [r,~] = corr(trace_1(idx2use,Itime), trace_2(idx2use,Itime), 'rows', 'complete', 'type', 'Spearman');
            tmp_rho(Ishuffle,Itime)  = r;

        end%Itime

    end%Ishuffle

    % smooth data 
    tmp_rho = smoothdata(tmp_rho, 2, 'movmean', [smoothwin smoothwin]);

    % - [FISHER Z TRANSFORMATION TO COMPUTE CONFIDENCE INTERVALS] - %
     
    x    = nanmean(tmp_rho); % current corr
    z    = 0.5*log((1+x)./(1-x)); % fisher z transformation
    ci   = [z + 1.96 * sqrt(1/(niter-3)); z - 1.96 * sqrt(1/(niter-3))]; % confidence intervals in z-space
    ci_r = [(exp(2*ci(1,:))-1) ./ (exp(2*ci(1,:))+1); (exp(2*ci(2,:))-1) ./ (exp(2*ci(2,:))+1)]; % back transformation into r-space  

    boot_ci.(tmpROI)  = ci_r; clear ci_r
    boot_rho.(tmpROI) = x; clear x
        
end%Iroi

%% - CONTRAST CORRELATION COEFFICIENTS BETWEEN PFC & MOTOR CORTEX OVER TIME - %%

% extract the number of participants with neural information per region
nsub_m1  = size(squeeze(w2.br.Motor(1,:,:)),2) - ...
    sum(all(isnan(squeeze(w2.br.Motor(1,:,:)))),2);
nsub_pfc = size(squeeze(w2.br.Prefrontal(1,:,:)),2) - ...
    sum(all(isnan(squeeze(w2.br.Prefrontal(1,:,:)))),2);

% test for significant differences in correlation coefficients over time
pvals_uncorr = NaN(1, length(w2.br.time));

for Itime = 1:length(pvals_uncorr)

    pvals_uncorr(Itime) = compare_correlation_coefficients(boot_rho.Prefrontal(Itime), ...
        boot_rho.Motor(Itime), nsub_pfc, nsub_m1);

end%Itime

[~, ~, ~, pvals_corr] = fdr_bh(pvals_uncorr, 0.05);

% get timepoints where the correlation coefficients are significantly
% different between motor cortex and pfc
SigCorrDiff = find(pvals_corr < 0.05);

fprintf('\nminimum p-value for correlation difference between PFC and M1 = %.3f\n', min(pvals_corr));

%% - PLOT TRACES - %%
    
% colors
color{1}   = 'r';
color{2}   = 'b';

% location for scatters
locscatter = -0.9;

figure;

for Iroi = 1:numel(roi2use)
    
    % current ROI
    tmpROI = roi2use{Iroi};
  
    % get current mean correlation and confidence interval
    tmp_avg = boot_rho.(tmpROI);
    tmp_ci  = [boot_ci.(tmpROI)(1,:) - tmp_avg; ...
               tmp_avg - boot_ci.(tmpROI)(2,:)];
    
    a = shadedErrorBar(w2.br.time, tmp_avg, tmp_ci);
    a.mainLine.Color = color{Iroi}; a.mainLine.LineWidth = 3;
    a.patch.FaceColor = color{Iroi}; a.patch.EdgeColor = color{Iroi}; a.patch.FaceAlpha = 0.5;
    a.edge(1).Color = color{Iroi}; a.edge(1).Color(4) = 0.1; a.edge(1).LineWidth = eps;
    a.edge(2).Color = color{Iroi}; a.edge(2).Color(4) = 0.1; a.edge(2).LineWidth = eps;
    
    hold on

    % - [MARK SIGNIFICANT CLUSTERS] - %

    % check if there is a negative cluster
    if isfield(true_clustermass.(tmpROI), 'neg')
        
        % compute p-val
        pval_neg = NaN(length(true_clustermass.(tmpROI).neg),1);
        for jj = 1:size(pval_neg,1)
            pval_neg(jj) = sum(clustermass_surro.(tmpROI).neg > ...
                abs(true_clustermass.(tmpROI).neg(jj))) / nshuffle;
        end

    end

    % check if there is a positive cluster
    if isfield(true_clustermass.(tmpROI), 'pos')

        % compute p-val
        pval_pos = NaN(length(true_clustermass.(tmpROI).pos),1);
        for jj = 1:size(pval_pos,1)
            pval_pos(jj) = sum(clustermass_surro.(tmpROI).pos > ...
                abs(true_clustermass.(tmpROI).pos(jj))) / nshuffle;
        end
        
    end
    
    % mark significant timepoints from negative clusters
    if true_island_neg.(tmpROI).NumObjects > 0
        for Icluster = 1:true_island_neg.(tmpROI).NumObjects
            if pval_neg(Icluster) < 0.05
                sigpnts = true_island_neg.(tmpROI).PixelIdxList{Icluster};
                plot(w2.br.time(sigpnts), ones(1, length(sigpnts)) * locscatter, ...
                    'Color', color{Iroi}, 'LineWidth', 2);
            end
        end
    end

    % mark significant timepoints from positive clusters
    if true_island_pos.(tmpROI).NumObjects > 0
        for Icluster = 1:true_island_pos.(tmpROI).NumObjects
            if pval_pos(Icluster) < 0.05
                sigpnts = true_island_pos.(tmpROI).PixelIdxList{Icluster};
                plot(w2.br.time(sigpnts), ones(1, length(sigpnts)) * locscatter, ...
                    'Color', color{Iroi}, 'LineWidth', 2);
            end
        end
    end   

    hold on

end

% mark significant difference between correlations (context vs. history)
plot(w2.br.time(SigCorrDiff), ones(1, length(SigCorrDiff)) * -0.8, ...
   'Color', 'k', 'LineWidth', 2);

% figure settings
xlabel('Time [s]');
ylabel('rho');

xlim([w2.br.time(1) w2.br.time(end)]);
ylim([-1 1]);

yline(0, 'k', '--', 'LineWidth', 1.5)

set(gca, 'fontsize', 13, 'linewidth', 2);
box off
set(gcf, 'Position', [480 414 355 243]);
