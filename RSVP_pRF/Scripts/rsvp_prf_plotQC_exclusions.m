function [fh] = rsvp_prf_plotQC_exclusions(ii_trial,ii_cfg)
%RSVP_PRF_PLOTQC_EXCLUSIONS Plot which TRs are flagged for saccades
%
%   fh = rsvp_prf_plotQC_exclusions(ii_trial,ii_cfg)
%
%   Plots ii_trial{run}.tr_excl, made by rsvp_prf_proc: a TR is flagged
%   (red) if a saccade started in it. Trial-level exclusion codes
%   (ii_trial{run}.excl_trial) are not plotted. Saves, next to the EDFs:
%   - <run>_TR_matrix.png: TR x sweep matrix for each run
%   - <folder>_BarHoriz_matrix.png: horizontal sweeps from all runs, aligned
%     by bar location (left to right), with the number of included and
%     excluded TRs at each location
%   - <folder>_BarVert_matrix.png: the same for vertical sweeps (top to
%     bottom)
%
% Sweep directions (used for labels and for grouping sweeps by direction)
% are read per run from ii_trial{ii}.sweep_dir (added by
% rsvp_prf_collectTrials from ii_cfg.task_timing), so RDK runs, whose
% directions come from p.seq, are labeled correctly. If that field is
% missing, ii_cfg.task_timing.sweep_dir is used for every run, and failing
% that the RSVP direction cycle (L2R, T2B, R2L, B2T).



% Get params
n_bars = ii_cfg.rsvp_prf_params.n_bars;
tr_per_bar = ii_cfg.rsvp_prf_params.tr_per_bar;
tr_per_sweep = n_bars * tr_per_bar;


for ii = 1:size(ii_trial,1)

    % sweeps in this run and their directions
    n_trials = size(ii_trial{ii}(1).tr_excl, 2);
    sweep_dir = get_sweep_dirs(ii_trial{ii}(1), ii_cfg, n_trials);

    % Set up matrix for plot
    excl_matrix = ii_trial{ii}(1).tr_excl';
    excl_matrix = flip(excl_matrix);
    excl_matrix(end+1,:) = zeros(tr_per_sweep,1)';
    excl_matrix(:,end+1) = zeros(n_trials + 1,1);

    % Set color of trs
    color_idx = [2 3 4 5];
    for i = 1:size(excl_matrix,2)
        these_boxes = find(excl_matrix(:,i) == 0);
        tr_col = mod(i, tr_per_bar);
        if tr_col == 0
            tr_col = tr_per_bar;
        end
        tr_col = tr_col + 1;
        excl_matrix(these_boxes,i) = tr_col;
    end


    % Plot matrix
    figure('visible','off');
    x = [repmat(1:(tr_per_sweep + 1), (n_trials + 1), 1)];
    y = [repmat(1:(n_trials + 1), (tr_per_sweep + 1), 1)]';
    map = [1 0 0; 0.8 0.8 0.8; 0.7 0.7 0.7; 0.55 0.55 0.55; 0.4 0.4 0.4;];
    c = excl_matrix;
    fh = pcolor(x,y,c);
    colormap(map);
    set(gca, 'CLim', [1 tr_per_bar + 1]); % fixed range, so 1 (excluded) is red even if no TRs are excluded
    fh.EdgeColor = [0.95 0.95 0.95];

    % Set figure params
    fig_height = n_trials * 50/1440; % to normalized units
    fig_width  = tr_per_sweep * 50/2560;
    fig_center = [0.5 0.5];
    fig_pos = [fig_center - [fig_width fig_height]*0.5 fig_width fig_height];
    set(gcf,'Units','Normalized');
    set(gcf,'Position',fig_pos);
    set(fh, 'facealpha', 0.5);

    % Draw boxes around timepoints
    set(fh, 'LineWidth', 5);
    yline(1, 'LineWidth', 5);
    yline(n_trials + 1, 'LineWidth', 5);
    xline(1, 'LineWidth', 5);
    for i = 1:(tr_per_sweep + 1)
        if mod(i-1,tr_per_bar) == 0
            xline(i, 'LineWidth', 5);
        end
    end


    % Y axis: one label per sweep, centered on its row (rows run from the
    % last sweep at the bottom to sweep 1 at the top)
    ylabels = cell(1, n_trials);
    for i = 1:n_trials
        ylabels{i} = ['Trial ' num2str(i) ' - ' sweep_dir{i} '  '];
    end
    set(gca, 'YTick', (1:n_trials) + 0.5, 'YTickLabel', flip(ylabels), 'fontsize', 16);
    % X axis
    new_xticks = [1:tr_per_sweep];
    set(gca,'XTick', new_xticks);
    % Titles
    title('TRs Excluded for Saccades','FontSize',16);
    xlabel('TR', 'FontSize',16);

    %file_loc = find(ii_cfg.edf_file == '/');
    [run_path, run_name] = fileparts(ii_trial{ii}(1).edf_file);
    saveas(fh, fullfile(run_path, [run_name '_TR_matrix.png']));
    clear excl_matrix;

end

%% Align via space
horiz_trial_labels = {};
vert_trial_labels = {};
horiz_matrix = [];
vert_matrix = [];
% Group sweeps by axis, ordered in space (left->right, top->bottom). R2L and
% B2T sweeps are reversed bar by bar, keeping the TR order within each bar.
for i = 1:size(ii_trial,1)
    excl_matrix = ii_trial{i}(1).tr_excl;
    sweep_dir = get_sweep_dirs(ii_trial{i}(1), ii_cfg, size(excl_matrix, 2));
    for j = 1:size(excl_matrix, 2)
        switch sweep_dir{j}
            case 'L2R'
                horiz_matrix(end+1,:) = excl_matrix(:,j);
                horiz_trial_labels{end+1} = ['\newline' 'Run ' num2str(i) ', Trial ' num2str(j) ' - L2R'];
            case 'R2L'
                horiz_matrix(end+1,:) = flip_bars(excl_matrix(:,j), tr_per_bar);
                horiz_trial_labels{end+1} = ['\newline' 'Run ' num2str(i) ', Trial ' num2str(j) ' - R2L'];
            case 'T2B'
                vert_matrix(end+1,:) = excl_matrix(:,j);
                vert_trial_labels{end+1} = ['\newline' ' Run ' num2str(i) ', Trial ' num2str(j) ' - T2B'];
            case 'B2T'
                vert_matrix(end+1,:) = flip_bars(excl_matrix(:,j), tr_per_bar);
                vert_trial_labels{end+1} = ['\newline' ' Run ' num2str(i) ', Trial ' num2str(j) ' - B2T'];
            otherwise
                error('rsvp_prf_plotQC_exclusions:unknownDir','Unknown sweep direction ''%s'' (run %i, sweep %i)', sweep_dir{j}, i, j);
        end
    end
end
clear excl_matrix x y map c fh;
% REMEMBER to flip Y trials

%% Horizontal Alignment

% Set up matrix for plot
horiz_matrix = flip(horiz_matrix);
% Get n included & excluded trials for each bar location
incl_sum = sum(horiz_matrix == 0);
excl_sum = sum(horiz_matrix);
horiz_incl = [];
horiz_excl = [];
sum_idc = {};
for i = 1:tr_per_bar
    sum_idc{i} = i:tr_per_bar:tr_per_sweep;
end
horiz_incl = incl_sum(sum_idc{1});
horiz_excl = excl_sum(sum_idc{1});
for i = 2:tr_per_bar
    horiz_incl = horiz_incl + incl_sum(sum_idc{i});
    horiz_excl = horiz_excl + excl_sum(sum_idc{i});
end
horiz_sub_xlabels = [horiz_incl; horiz_excl];
horiz_sub_xlabels = horiz_sub_xlabels(:)';
% Add filler row & column for plotting
horiz_matrix(end+1,:) = zeros(tr_per_sweep,1)';
horiz_matrix(:,end+1) = zeros(size(horiz_matrix,1),1);


% Set color of trs
for i = 1:size(horiz_matrix,2)
    these_boxes = find(horiz_matrix(:,i) == 0);
    tr_col = mod(i, tr_per_bar);
    if tr_col == 0
        tr_col = tr_per_bar;
    end
    tr_col = tr_col + 1;
    horiz_matrix(these_boxes,i) = tr_col;
end


% Plot matrix (in a new figure: resizing the last run's figure doesn't
% always take effect before it's saved)
figure('visible','off');
subplot(size(horiz_matrix,1) + 4, 1, [1:size(horiz_matrix,1)]);
x = [repmat(1:(tr_per_sweep + 1), size(horiz_matrix,1), 1)];
y = [repmat(1:size(horiz_matrix,1), tr_per_sweep + 1, 1)]';
map = [1 0 0; 0.8 0.8 0.8; 0.7 0.7 0.7; 0.6 0.6 0.6; 0.5 0.5 0.5;];
c = horiz_matrix;
fh = pcolor(x,y,c);
colormap(map);
set(gca, 'CLim', [1 tr_per_bar + 1]); % fixed range, so 1 (excluded) is red even if no TRs are excluded
fh.EdgeColor = [0.95 0.95 0.95];


% Set figure params
fig_height = size(horiz_matrix,1) * 50/1440; % to normalized units
fig_width  = tr_per_sweep * 50/2560;
fig_center = [0.5 0.5];
fig_pos = [fig_center - [fig_width fig_height]*0.5 fig_width fig_height];
set(gcf,'Units','Normalized');
set(gcf,'Position',fig_pos);
set(fh, 'facealpha', 0.5);


% Draw boxes around timepoints
set(fh, 'LineWidth', 5);
yline(1, 'LineWidth', 5);
yline(size(horiz_matrix,1), 'LineWidth', 5);
xline(1, 'LineWidth', 5);
for i = 1:(tr_per_sweep + 1)
    if mod(i-1,tr_per_bar) == 0
        xline(i, 'LineWidth', 5);
    end
end


% Y axis
yticks = [1:size(horiz_matrix,1)];
set(gca, 'YTick', yticks);
horiz_trial_labels{end+1} = '';
horiz_trial_labels = flip(horiz_trial_labels);
set(gca, 'YTickLabel', horiz_trial_labels, 'fontsize', 16);

% X axis
new_xticks = [2:tr_per_bar:tr_per_sweep];
set(gca,'XTick', new_xticks);
xlabels = 1:n_bars;
set(gca, 'XTickLabel', xlabels);

% Titles
title('TRs Excluded, Horizontal Alignment','FontSize',16);
xlabel('Bar Location', 'FontSize',16);
ax_main = gca;

% Create subplot
subplot(size(horiz_matrix,1) + 4, 1, [size(horiz_matrix,1) + 3, size(horiz_matrix,1) + 4]);
%summary_matrix = [repmat([1 0],1,n_bars); zeros(tr_per_sweep,1)'];
summary_matrix = [repmat([1 0],1,n_bars); zeros(2 * n_bars, 1)'];
summary_matrix(:,end+1) = zeros(size(summary_matrix,1),1);

% Set up matrix
sub_x = [repmat(1:(2 * n_bars + 1),2,1)];
sub_y = [repmat(1,1, 2 * n_bars + 1); repmat(2,1, 2 * n_bars + 1);];
sub_c = summary_matrix;
fh_sub = pcolor(sub_x,sub_y,sub_c);
colormap(map);
fh_sub.EdgeColor = [0.95 0.95 0.95];
set(fh_sub, 'facealpha', 0.5);

% Draw boxes around timepoints
set(fh_sub, 'LineWidth', 5);
yline(1, 'LineWidth', 5);
yline(2, 'LineWidth', 5);
xline(1, 'LineWidth', 5);
for i = 1:(tr_per_sweep + 1)
    if mod(i-1,2) == 0
        xline(i, 'LineWidth', 5);
    end
end
for i = 1:(2 * n_bars)
    text(i + 0.5 ,1.5,num2str(horiz_sub_xlabels(i)), 'FontSize', 12, 'HorizontalAlignment', 'center');
end
set(gca,'xtick',[]);
set(gca,'ytick',[]);
xlabel('Included TRs, Excluded TRs', 'FontSize',16);
ax_sub = gca;

% Make room left of the main panel for the run/trial labels: move the left
% edge of both panels right by the space the labels need, so the summary
% panel's columns stay aligned with the main panel's. (Done in normalized
% units: this figure can be taller than the screen, and setting its size
% in pixels gets clamped to the screen.)
drawnow;
set([ax_main ax_sub], 'Units', 'normalized');
extra = ax_main.TightInset(1) + 0.01 - ax_main.Position(1);
if extra > 0
    ax_main.Position = ax_main.Position + [extra 0 -extra 0];
    ax_sub.Position = ax_sub.Position + [extra 0 -extra 0];
end


%file_loc = find(ii_cfg.edf_file == '/');
saveas(fh, session_fig_name(ii_trial, 'BarHoriz_matrix.png'));

clear incl_sum excl_sum sum_idc these_boxes;



%% Vertical Alignment


% Set up matrix for plot
vert_matrix = flip(vert_matrix);
% Get n included & excluded trials for each bar location
incl_sum = sum(vert_matrix == 0);
excl_sum = sum(vert_matrix);
vert_incl = [];
vert_excl = [];
sum_idc = {};
for i = 1:tr_per_bar
    sum_idc{i} = i:tr_per_bar:tr_per_sweep;
end
vert_incl = incl_sum(sum_idc{1});
vert_excl = excl_sum(sum_idc{1});
for i = 2:tr_per_bar
    vert_incl = vert_incl + incl_sum(sum_idc{i});
    vert_excl = vert_excl + excl_sum(sum_idc{i});
end
vert_sub_xlabels = [vert_incl; vert_excl];
vert_sub_xlabels = vert_sub_xlabels(:)';
% Flip for vertical alignment
vert_matrix = vert_matrix';
% Add filler row & column for plotting
vert_matrix = flip(vert_matrix,2);
vert_matrix(end+1,:) = zeros(size(vert_matrix,2),1);
vert_matrix(:,end+1) = zeros(tr_per_sweep + 1,1)';
vert_matrix = flip(vert_matrix,2);

% Set color of trs
vert_matrix = vert_matrix';
for i = 1:size(vert_matrix,2)
    these_boxes = find(vert_matrix(:,i) == 0);
    tr_col = mod(i, tr_per_bar);
    if tr_col == 0
        tr_col = tr_per_bar;
    end
    tr_col = tr_col + 1;
    vert_matrix(these_boxes,i) = tr_col;
end
vert_matrix = flip(vert_matrix',2);


% Plot matrix
figure('visible','off');
subplot(1, size(vert_matrix,2) + 4, [1:size(vert_matrix,2)]);
y = [repmat(1:(tr_per_sweep + 1), size(vert_matrix,2), 1)]';
x = [repmat(1:size(vert_matrix,2), tr_per_sweep + 1, 1)];
map = [1 0 0; 0.8 0.8 0.8; 0.7 0.7 0.7; 0.6 0.6 0.6; 0.5 0.5 0.5;];
c = vert_matrix;
fh_vert = pcolor(x,y,c);
colormap(map);
set(gca, 'CLim', [1 tr_per_bar + 1]); % fixed range, so 1 (excluded) is red even if no TRs are excluded
fh_vert.EdgeColor = [0.95 0.95 0.95];


% Set figure params
fig_height = size(vert_matrix,2) * 50/1440; % to normalized units
fig_width  = tr_per_sweep * 100/2560;
fig_center = [0.5 0.5];
fig_pos = [fig_center - [fig_width fig_height]*0.5 fig_width fig_height];
set(gcf,'Units','Normalized');
set(gcf,'Position',fig_pos);
set(fh_vert, 'facealpha', 0.5);


% Draw boxes around timepoints
set(fh_vert, 'LineWidth', 5);
xline(1, 'LineWidth', 5);
xline(size(vert_matrix,2), 'LineWidth', 5);
yline(1, 'LineWidth', 5);
for i = 1:(tr_per_sweep + 1)
    if mod(i-1,tr_per_bar) == 0
        yline(i, 'LineWidth', 5);
    end
end


% X axis: one label per sweep, centered on its column
vert_trial_labels = strtrim(strrep(vert_trial_labels, '\newline', ''));
set(gca, 'XTick', (1:numel(vert_trial_labels)) + 0.5, 'XTickLabel', vert_trial_labels, 'fontsize', 16);
xtickangle(90);

% Y axis: rows run from the top of the screen (bar location 1) at y = 1, so
% reverse the axis to draw location 1 at the top
set(gca, 'YDir', 'reverse');
new_yticks = [2:tr_per_bar:tr_per_sweep];
set(gca,'YTick', new_yticks);
ylabels = 1:n_bars;
set(gca, 'YTickLabel', ylabels);
% Titles
ylabel('Bar Location', 'FontSize',16);
title('TRs Excluded, Vertical Alignment','FontSize',16);
ax_main = gca;


% Create subplot
subplot(1, size(vert_matrix,2) + 4, [size(vert_matrix,2) + 3, size(vert_matrix,2) + 4]);
summary_matrix = summary_matrix';

% Set up matrix
sub_x = [repmat(1,1,2* n_bars + 1); repmat(2,1, 2* n_bars + 1);]';
sub_y = [repmat(1:(2 * n_bars + 1),2,1)]';
sub_c = summary_matrix;
fh_sub = pcolor(sub_x,sub_y,sub_c);
set(gca, 'YDir', 'reverse'); % bar location 1 at the top, as in the main panel
colormap(map);
fh_sub.EdgeColor = [0.95 0.95 0.95];
set(fh_sub, 'facealpha', 0.5);


% Draw boxes around timepoints
set(fh_sub, 'LineWidth', 5);
xline(1, 'LineWidth', 5);
xline(2, 'LineWidth', 5);
yline(1, 'LineWidth', 5);
for i = 1:(2 * n_bars + 1)
    if mod(i-1,2) == 0
        yline(i, 'LineWidth', 5);
    end
end
for i = 1:(2 * n_bars)
    text(1.4, i + 0.5, num2str(vert_sub_xlabels(i)), 'FontSize', 14);
end
set(gca,'xtick',[]);
set(gca,'ytick',[]);
% Titles
ylabel('Included TRs, Excluded TRs', 'FontSize',16);
ax_sub = gca;

% Make room below the main panel for the rotated run labels: grow the
% figure by the space they need and move both panels up by the same
% amount, so the summary panel's rows stay aligned with the main panel's
drawnow;
set([gcf ax_main ax_sub], 'Units', 'pixels');
pos_main = ax_main.Position;
pos_sub = ax_sub.Position;
extra_px = ax_main.TightInset(2) + 10 - pos_main(2);
if extra_px > 0
    set(gcf, 'Position', get(gcf, 'Position') + [0 0 0 extra_px]);
    ax_main.Position = pos_main + [0 extra_px 0 0];
    ax_sub.Position = pos_sub + [0 extra_px 0 0];
end

saveas(fh_sub, session_fig_name(ii_trial, 'BarVert_matrix.png'));




return

end


function sweep_dir = get_sweep_dirs(trial, ii_cfg, n_sweeps)
% Direction label of each sweep in one run (cellstr, n_sweeps x 1)
if isfield(trial, 'sweep_dir')
    sweep_dir = trial.sweep_dir;
elseif isfield(ii_cfg, 'task_timing')
    warning('rsvp_prf_plotQC_exclusions:noRunDirs', ...
        ['ii_trial has no sweep_dir (collected before it was added); using ii_cfg.task_timing.sweep_dir ' ...
         'for every run, which assumes all runs used the same sweep sequence']);
    sweep_dir = ii_cfg.task_timing.sweep_dir;
else
    rsvp_dirs = {'L2R'; 'T2B'; 'R2L'; 'B2T'};
    sweep_dir = rsvp_dirs(mod((1:n_sweeps).' - 1, 4) + 1);
end
if numel(sweep_dir) < n_sweeps
    error('rsvp_prf_plotQC_exclusions:dirCount', '%i sweep directions for %i sweeps', numel(sweep_dir), n_sweeps);
end
end


function fn = session_fig_name(ii_trial, suffix)
% File name for a figure combining all runs: saved in the first run's
% folder and named after that folder (e.g. Subjects/AY_RF1_behav/
% AY_RF1_behav_BarHoriz_matrix.png)
run_path = fileparts(ii_trial{1}(1).edf_file);
[~, folder_name] = fileparts(run_path);
fn = fullfile(run_path, [folder_name '_' suffix]);
end


function v = flip_bars(v, tr_per_bar)
% Reverse the order of bar positions in a sweep's TR vector, keeping the
% order of TRs within each bar
v = reshape(v, tr_per_bar, []);
v = reshape(fliplr(v), [], 1);
end

