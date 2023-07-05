%% - FIGURE 2A/B - %%
% This code plots the traces of context- and history-dependent neural
% information for prefrontal and motor cortex as in Fig. 2A/B
%
%%% --- REQUIRED INPUTS --- %%%
% - w2_roi.mat
% - ROI_sub_clean.mat

%%

clear
close all;
clc;

%% - LOAD DATA - %%

% load omega-square + ROI 
load('w2_roi.mat')
load('ROIsub_clean.mat')

%% - MAIN SETTINGS - %%

% supporting vars
[~, ~, nSub] = size(ROIsub.br.nElec); 

% find factor indices
factoridx = [find(strcmp(w2.br.factorNames, 'Cue')) ...
             find(strcmp(w2.br.factorNames, 'Past Trial Type'))];

% # factors
nFactor  = length(factoridx);

% find ROI indices
ROIidx   = [find(strcmp(w2.br.ROInames, 'Prefrontal')) ...
    find(strcmp(w2.br.ROInames, 'Motor'))];

% # ROIs
nROI     = length(ROIidx);

% segment names
segNames = {'br'};
nseg     = numel(segNames);
nTime.br = numel(w2.br.time);

%% - GET TEMPORAL CLUSTERS - %%

% preallocate memory for observed clusters
obsCluster = cell(nseg, nROI, nFactor);

for iSeg = 1:nseg
    
    tmpseg = segNames{iSeg};

    % initialize tvalues (factor x ROI x time)
    tvalues.(tmpseg) = NaN(nFactor,nROI,nTime.(tmpseg));

    for iROI = 1:nROI
        
        % get current ROI index
        tmp_ROIidx = ROIidx(iROI);

        for iFactor = 1:nFactor
            
            tmp_factoridx = factoridx(iFactor);

            % -- [OBSERVED DISTRIBUTION] -- %

            % sample selection
            encoding  = squeeze(w2.(tmpseg).(ROIsub.(tmpseg).labels(tmp_ROIidx))(tmp_factoridx,:,:));
            nencoding = squeeze(w2.(tmpseg).baseline.(ROIsub.(tmpseg).labels(tmp_ROIidx))(tmp_factoridx,:,:));
            
            idxout = unique([find(all(isnan(encoding),1)) find(all(isnan(nencoding),1))]);
            encoding(:,idxout)  = [];
            nencoding(:,idxout) = [];

            [~,~,~, statstemp] = ttest(double(encoding)', double(nencoding)'); 
            tvalues.(tmpseg)(iFactor,iROI, :) = statstemp.tstat;  
                        
            % get threshold for clustering
            df = size(encoding,2)-1; % sample size - 1
            criticalT = abs(tinv(.05,df));
            
            % remove t-values below threshold
            tvalues_temp = squeeze(tvalues.(tmpseg)(iFactor,iROI, :));
            tvalues_temp(abs(tvalues_temp) < criticalT) = 0; 
            
            % find clusters
            island = bwconncomp(tvalues_temp); 
            
            % sum tvalues in the clusters
            obsCluster{iSeg,iROI,iFactor}.Idxlist = island.PixelIdxList;
            for ii = 1:island.NumObjects
                cluster_tvalues = tvalues_temp(island.PixelIdxList{ii});
                obsCluster{iSeg,iROI,iFactor}.sumtval(ii) = sum(cluster_tvalues);
            end
            
        end%iFactor
        
    end%iROI

end%iSeg

%% - PERMUTATION TEST - %%

echo off;
warning off;

% -- [PERMUTATION] -- % 

% number of iterations for permutation testing
nShuffle = 1000;

for iSeg = 1:nseg

    tmpseg = segNames{iSeg};

    % preallocate memory shuffled cluster and t-values
    randCluster_maxtval.(tmpseg) = NaN(nROI,nFactor,nShuffle);
    tvalues_rand.(tmpseg)        = NaN(nROI, nFactor, nTime.(tmpseg));

    for iROI = 1:nROI

        % get current ROI index
        tmp_ROIidx = ROIidx(iROI);
            
        for iFactor = 1:nFactor

            tmp_factoridx = factoridx(iFactor);

            fprintf('Segment %d/%d, ROI %d/%d, Factor %d/%d\n', iSeg, nseg, iROI, nROI, iFactor, nFactor)
    
            temp = NaN(nShuffle,1);
            
            for iSample = 1:nShuffle
                
                % sample selection
                encoding  = squeeze(w2.(tmpseg).(ROIsub.(tmpseg).labels(tmp_ROIidx))(tmp_factoridx,:,:));
                nencoding = squeeze(w2.(tmpseg).baseline.(ROIsub.(tmpseg).labels(tmp_ROIidx))(tmp_factoridx,:,:));
                
                idxout = unique([find(all(isnan(encoding),1)) find(all(isnan(nencoding),1))]);
                encoding(:,idxout)  = [];
                nencoding(:,idxout) = [];
        
                % random permutation
                [~, n]     = size(encoding);
                shufdata   = cat(2,encoding,nencoding); % concatenate data to shuffle then
                
                rndshuffle = randperm(size(shufdata,2));

                shuf_c1    = shufdata(:,rndshuffle(1:end/2));
                shuf_c2    = shufdata(:,rndshuffle(end/2+1:end));

                % contrast
                [~,~,~, statstemp] = ttest(double(shuf_c1)', double(shuf_c2)');
                tvalues_rand.(tmpseg)(iFactor,iROI, :) = statstemp.tstat; 
                
                % get threshold for clustering
                df = size(encoding,2)-1; % sample size - 1
                criticalT = abs(tinv(.05, df));
                
                % remove t-values below threshold
                tvalues_temp = squeeze(tvalues_rand.(tmpseg)(iFactor,iROI, :));
                tvalues_temp(abs(tvalues_temp) < criticalT) = 0; 
                
                % find clusters
                island = bwconncomp(tvalues_temp); 
                                        
                % take sum of tvalues within every cluster
                if island.NumObjects > 0
                    
                    cluster_tvalues_rand = NaN(island.NumObjects,1);
                    for ii = 1:island.NumObjects
                        cluster_tvalues_rand(ii) = sum(abs(tvalues_temp(island.PixelIdxList{ii})));
                    end
                   
                    % select only the highest value
                    temp(iSample) = max(cluster_tvalues_rand,[],'all','omitnan');
                                                               
                end  
                
            end%iSample
            
            % save max tvalue of every iteration
            randCluster_maxtval.(tmpseg)(iROI,iFactor,:) = temp;
    
        end%iFactor

    end%iROI
    
end%iSeg

fprintf('\n')
   
%% - COMPUTE P-VALUE FOR EACH CLUSTER - %%

for iSeg = 1:nseg

    tmpseg = segNames{iSeg};

    for iROI = 1:nROI

        for iFactor = 1:nFactor
            
            if isfield(obsCluster{iSeg, iROI, iFactor}, 'sumtval')

                % observed values for cluster
                obstval  = obsCluster{iSeg, iROI, iFactor}.sumtval;
    
                % shuffled max sum(t) per permutation
                randtval = squeeze(randCluster_maxtval.(tmpseg)(iROI,iFactor, :));
                randtval(isnan(randtval)) = [];
    
                % compute p-value based on null distribution
                pval = sum(randtval > obstval) / nShuffle;
                obsCluster{iSeg, iROI, iFactor}.pval = pval;  

            else

                obsCluster{iSeg,iROI,iFactor}.pval = NaN;

            end
                
        end%iFactor

    end%iROI

end%iSeg

%% - GET SIGNIFICANT CLUSTERS - %%

% critical alpha for clusters
critalpha = 0.05;

for iSeg = 1:nseg

    tmpseg = segNames{iSeg};

    % ROI x Factor x Time
    sigClusters = zeros(nROI,nFactor,nTime.(tmpseg));
    
    for iROI = 1:nROI

        for iFactor = 1:nFactor
            
            my_data = obsCluster{iSeg,iROI, iFactor}.Idxlist;
            idx     = [];
            
            for jj = 1:numel(my_data) % select cluster
                if obsCluster{iSeg,iROI,iFactor}.pval(jj) < 0.05
                    idx = cat(1, idx, my_data{1,jj}); % concatenate clusters idx
                end
            end
            
            sigClusters(iROI, iFactor, idx) = 1;
            
        end%iFactor

    end%iROI

    % pass over to structure
    w2.(tmpseg).sigClusters = sigClusters;
    
end%iSeg

%% - PLOT %EV FOR PFC (HISTORY AND CONTEXT) - %%

close all;

plot_colors  = {'r', 'b'};

% - [PFC CONTEXT] -%

figure;

tmpData  = squeeze(w2.br.Prefrontal(1,:,:))'*100;
tmpStats = logical(squeeze(w2.br.sigClusters(1,1,:)));

a = shadedErrorBar(w2.br.time, mean(tmpData, 1, 'omitnan'), ...
  std(tmpData, 1, 'omitnan') / sqrt(size(tmpData,1) - ...
  sum(all(isnan(tmpData),2))));
a.mainLine.Color = plot_colors{1}; a.mainLine.LineWidth = 3;
a.patch.FaceColor = plot_colors{1}; a.patch.EdgeColor = plot_colors{1}; a.patch.FaceAlpha = 0.5;
a.edge(1).Color = plot_colors{1}; a.edge(1).Color(4) = 0.1; a.edge(1).LineWidth = eps;
a.edge(2).Color = plot_colors{1}; a.edge(2).Color(4) = 0.1; a.edge(2).LineWidth = eps;

hold on 

tmpBaseline = squeeze(w2.br.baseline.Prefrontal(1,:,:))'*100;
a = shadedErrorBar(w2.br.time, mean(tmpBaseline, 1, 'omitnan'), ...
  std(tmpBaseline, 1, 'omitnan') / sqrt(size(tmpBaseline,1) - ...
  sum(all(isnan(tmpBaseline),2))));
a.mainLine.Color = [0.3 0.3 0.3]; a.mainLine.LineWidth = 3;
a.patch.FaceColor = [0.3 0.3 0.3]; a.patch.EdgeColor = [0.3 0.3 0.3]; a.patch.FaceAlpha = 0.5;
a.edge(1).Color = [0.3 0.3 0.3]; a.edge(1).Color(4) = 0.1; a.edge(1).LineWidth = eps;
a.edge(2).Color = [0.3 0.3 0.3]; a.edge(2).Color(4) = 0.1; a.edge(2).LineWidth = eps;


xlabel('Time [s]');
ylabel('Neural Information [%EV]');

set(gca, 'fontsize', 13, 'linewidth', 1.5, 'xlim', ...
    [w2.br.time(1) w2.br.time(end)]);

hold on

% get minimum and maximum value
minval = min(mean(tmpBaseline, 1, 'omitnan'));

% scatter stats
x1 = w2.br.time(tmpStats);
y1 = ones(1,length(x1))*minval-abs(minval*4);
plot(x1, y1, 'Color', 'k', 'LineWidth', 2);

ylim([-0.5 3.5]);

box off

set(gcf, 'Position', [560 427 268 212]);

% - [PFC HISTORY] -%

figure;

tmpData  = squeeze(w2.br.Prefrontal(2,:,:))'*100;
tmpStats = logical(squeeze(w2.br.sigClusters(1,2,:)));

a = shadedErrorBar(w2.br.time, mean(tmpData, 1, 'omitnan'), ...
  std(tmpData, 1, 'omitnan') / sqrt(size(tmpData,1) - ...
  sum(all(isnan(tmpData),2))));
a.mainLine.Color = plot_colors{1}; a.mainLine.LineWidth = 3;
a.patch.FaceColor = plot_colors{1}; a.patch.EdgeColor = plot_colors{1}; a.patch.FaceAlpha = 0.5;
a.edge(1).Color = plot_colors{1}; a.edge(1).Color(4) = 0.1; a.edge(1).LineWidth = eps;
a.edge(2).Color = plot_colors{1}; a.edge(2).Color(4) = 0.1; a.edge(2).LineWidth = eps;

hold on 

tmpBaseline = squeeze(w2.br.baseline.Prefrontal(2,:,:))'*100;
a = shadedErrorBar(w2.br.time, mean(tmpBaseline, 1, 'omitnan'), ...
  std(tmpBaseline, 1, 'omitnan') / sqrt(size(tmpBaseline,1) - ...
  sum(all(isnan(tmpBaseline),2))));
a.mainLine.Color = [0.3 0.3 0.3]; a.mainLine.LineWidth = 3;
a.patch.FaceColor = [0.3 0.3 0.3]; a.patch.EdgeColor = [0.3 0.3 0.3]; a.patch.FaceAlpha = 0.5;
a.edge(1).Color = [0.3 0.3 0.3]; a.edge(1).Color(4) = 0.1; a.edge(1).LineWidth = eps;
a.edge(2).Color = [0.3 0.3 0.3]; a.edge(2).Color(4) = 0.1; a.edge(2).LineWidth = eps;

xlabel('Time [s]');
ylabel('Neural Information [%EV]');

set(gca, 'fontsize', 13, 'linewidth', 1.5, 'xlim', ...
    [w2.br.time(1) w2.br.time(end)]);

hold on

% get minimum and maximum value
minval = min(mean(tmpBaseline, 1, 'omitnan'));

% scatter stats
x1 = w2.br.time(tmpStats);
y1 = ones(1,length(x1))*minval-abs(minval*1);
plot(x1, y1, 'Color', 'k', 'LineWidth', 2);

ylim([-0.5 3.5]);

box off

set(gcf, 'Position', [560 427 268 212]);

%% - PLOT %EV FOR MOTOR CORTEX (HISTORY AND CONTEXT) - %%

close all;

% - [MOTOR CONTEXT] -%

figure;

tmpData  = squeeze(w2.br.Motor(1,:,:))'*100;
tmpStats = logical(squeeze(w2.br.sigClusters(2,1,:)));

a = shadedErrorBar(w2.br.time, mean(tmpData, 1, 'omitnan'), ...
  std(tmpData, 1, 'omitnan') / sqrt(size(tmpData,1) - ...
  sum(all(isnan(tmpData),2))));
a.mainLine.Color = plot_colors{2}; a.mainLine.LineWidth = 3;
a.patch.FaceColor = plot_colors{2}; a.patch.EdgeColor = plot_colors{2}; a.patch.FaceAlpha = 0.5;
a.edge(1).Color = plot_colors{2}; a.edge(1).Color(4) = 0.1; a.edge(1).LineWidth = eps;
a.edge(2).Color = plot_colors{2}; a.edge(2).Color(4) = 0.1; a.edge(2).LineWidth = eps;

hold on 

tmpBaseline = squeeze(w2.br.baseline.Motor(1,:,:))'*100;
a = shadedErrorBar(w2.br.time, mean(tmpBaseline, 1, 'omitnan'), ...
  std(tmpBaseline, 1, 'omitnan') / sqrt(size(tmpBaseline,1) - ...
  sum(all(isnan(tmpBaseline),2))));
a.mainLine.Color = [0.3 0.3 0.3]; a.mainLine.LineWidth = 3;
a.patch.FaceColor = [0.3 0.3 0.3]; a.patch.EdgeColor = [0.3 0.3 0.3]; a.patch.FaceAlpha = 0.5;
a.edge(1).Color = [0.3 0.3 0.3]; a.edge(1).Color(4) = 0.1; a.edge(1).LineWidth = eps;
a.edge(2).Color = [0.3 0.3 0.3]; a.edge(2).Color(4) = 0.1; a.edge(2).LineWidth = eps;

xlabel('Time [s]');
ylabel('Neural Information [%EV]');

set(gca, 'fontsize', 13, 'linewidth', 1.5, 'xlim', ...
    [w2.br.time(1) w2.br.time(end)]);

hold on

% get minimum and maximum value
minval = min(mean(tmpBaseline, 1, 'omitnan'));

% scatter stats
x1 = w2.br.time(tmpStats);
y1 = ones(1,length(x1))*minval-abs(minval*2);
plot(x1, y1, 'Color', 'k', 'LineWidth', 2);

box off

set(gcf, 'Position', [560 427 268 212]);

ylim([-0.5 3.5]);

% - [MOTOR HISTORY] -%

figure;

tmpData  = squeeze(w2.br.Motor(2,:,:))'*100;
tmpStats = logical(squeeze(w2.br.sigClusters(2,2,:)));

a = shadedErrorBar(w2.br.time, mean(tmpData, 1, 'omitnan'), ...
  std(tmpData, 1, 'omitnan') / sqrt(size(tmpData,1) - ...
  sum(all(isnan(tmpData),2))));
a.mainLine.Color = plot_colors{2}; a.mainLine.LineWidth = 3;
a.patch.FaceColor = plot_colors{2}; a.patch.EdgeColor = plot_colors{2}; a.patch.FaceAlpha = 0.5;
a.edge(1).Color = plot_colors{2}; a.edge(1).Color(4) = 0.1; a.edge(1).LineWidth = eps;
a.edge(2).Color = plot_colors{2}; a.edge(2).Color(4) = 0.1; a.edge(2).LineWidth = eps;

hold on 

tmpBaseline = squeeze(w2.br.baseline.Motor(2,:,:))'*100;
a = shadedErrorBar(w2.br.time, mean(tmpBaseline, 1, 'omitnan'), ...
  std(tmpBaseline, 1, 'omitnan') / sqrt(size(tmpBaseline,1) - ...
  sum(all(isnan(tmpBaseline),2))));
a.mainLine.Color = [0.3 0.3 0.3]; a.mainLine.LineWidth = 3;
a.patch.FaceColor = [0.3 0.3 0.3]; a.patch.EdgeColor = [0.3 0.3 0.3]; a.patch.FaceAlpha = 0.5;
a.edge(1).Color = [0.3 0.3 0.3]; a.edge(1).Color(4) = 0.1; a.edge(1).LineWidth = eps;
a.edge(2).Color = [0.3 0.3 0.3]; a.edge(2).Color(4) = 0.1; a.edge(2).LineWidth = eps;

xlabel('Time [s]');
ylabel('Neural Information [%EV]');

set(gca, 'fontsize', 13, 'linewidth', 1.5, 'xlim', ...
    [w2.br.time(1) w2.br.time(end)]);

hold on

% get minimum and maximum value
minval = min(mean(tmpBaseline, 1, 'omitnan'));

% scatter stats
x1 = w2.br.time(tmpStats);
y1 = ones(1,length(x1))*minval-abs(minval*1);
plot(x1, y1, 'Color', 'k', 'LineWidth', 2);

ylim([-0.5 3.5]);

box off

set(gcf, 'Position', [560 427 268 212]);
