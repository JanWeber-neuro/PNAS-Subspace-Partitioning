%% raincloud_pairplot plots paired data as 1-d scatter (raindrops), linking lines, and density plots (clouds)
% Use like: h = raincloud_pairplot(data, colours, plot_top_to_bottom)
% Where 'data' is a n x 2 matrix of paired data
% And 'colours' is a 2 x 3 matrix containing two RGB values
% By default data is plotted left to right, plot_top_to_bottom rotates it
% Example:
% data(:,1)    = 0.2 + rand(30,1);
% data(:,2)    = data(:,1) + rand(30,1) - 0.5;
% cl(1,:)      = [0.1 0.3 0.7]; 
% cl(2,:)      = [0.4 0.8 1];
% figure('pos',[0 0 300 600]); h = raincloud_pairplot(data, cl);

function h = raincloud_pairplot(data, colours, plot_top_to_bottom)

%% optional input args

% if user doesn't specify colours, use these
if nargin < 2
    colours(1,:) = [0.1 0.4 0.7];
    colours(2,:) = [0.4 0.7 1.0];
end

if nargin < 3
    plot_top_to_bottom = 0;
end

%% data;
hold on

sz = size(data);

xpos = zeros(sz);
xpos(:,1) = 1;
xpos(:,2) = 2;

%% lines

h.line = line(data', xpos', 'Color', [0.2 0.2 0.2], 'LineWidth', 2);

%% dots

h.dots1 = scatter(data(:,1), xpos(:,1));
h.dots1.MarkerEdgeColor = 'none';
h.dots1.MarkerFaceColor = colours(1,:);
h.dots1.SizeData = 120;

h.dots2 = scatter(data(:,2), xpos(:,2));
h.dots2.MarkerEdgeColor = 'none';
h.dots2.MarkerFaceColor = colours(2,:);
h.dots2.SizeData = 120;

%% densities

n_bins = 200;

[ks{1}, x{1}]   = ksdensity(data(:,1), 'NumPoints', n_bins, 'bandwidth', []);
[ks{2}, x{2}]   = ksdensity(data(:,2), 'NumPoints', n_bins, 'bandwidth', []);

q               = (1: n_bins-1)';
faces           = [q, q+1, q+n_bins+1, q+n_bins];
ks_offsets      = [0.8 2.2];
max_height      = 0.4;

% scale the density plots so they fill the empty space
ks{1} = ks{1} .* (max_height / max(ks{1}));
ks{2} = ks{2} .* (max_height / max(ks{2}));

% remember to plot this one upside down!
vertices{1} = [x{1}', -ks{1}' + ks_offsets(1); x{1}', ones(n_bins, 1) * ks_offsets(1)];
vertices{2} = [x{2}', ks{2}'  + ks_offsets(2); x{2}', ones(n_bins, 1) * ks_offsets(2)]; 

h.patch1 = patch('Faces', faces, 'Vertices', vertices{1}, 'FaceVertexCData', colours(1, :), 'FaceColor', 'flat', 'EdgeColor', 'none', 'FaceAlpha', 1);
h.patch2 = patch('Faces', faces, 'Vertices', vertices{2}, 'FaceVertexCData', colours(2, :), 'FaceColor', 'flat', 'EdgeColor', 'none', 'FaceAlpha', 1);

%% 

set(gca,'YLim',[0.3 2.7]);

%% rotation options
% note, because it's easiest to think about things this way, the whole plot
% is actually arranged BOTTOM-TO-TOP, then either flipped or rotated here
% depending on the value of 'plot_top_to_bottom'.

if plot_top_to_bottom
    % flip everything so '1' is on the top of the y-axis
    axis ij
else 
    % default is to plot left-to-right, which requires the plot to be rotated
    view([90 -90]);
end


