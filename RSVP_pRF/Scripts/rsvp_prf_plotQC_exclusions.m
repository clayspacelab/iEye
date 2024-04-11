function [fh] = rsvp_prf_plotQC_exclusions(ii_trial,ii_cfg,which_excl,fig_visible)
% ii_plotQC_exclusions Plots 



% Get params
n_bars = ii_cfg.rsvp_prf_params.n_bars;
tr_per_bar = ii_cfg.rsvp_prf_params.tr_per_bar;
tr_per_sweep = n_bars * tr_per_bar;
n_trials = ii_cfg.rsvp_prf_params.n_trials;


for ii = 1:size(ii_trial,1)

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


    % Y axis
    yticks = [1:n_trials];
    ylabels = {};
    sweep_dir = ['L2R';'T2B';'R2L';'B2T';];
    for i = 1:n_trials
        ylab = ['\newline' '\newline' 'Trial ' num2str(i) ' - ' sweep_dir(moditer(i,4),:) '  ' ];
        ylabels{end+1} = ylab;
    end
    ylabels{end+1} = '';
    ylabels = flip(ylabels);
    yticks(yticks);
    yticklabels(ylabels);
    y = get(gca,'YTickLabel');
    set(gca, 'YTickLabel', y, 'fontsize', 16);
    % X axis
    new_xticks = [1:tr_per_sweep];
    set(gca,'XTick', new_xticks);
    % Titles
    title('TRs Excluded for Saccades','FontSize',16);
    xlabel('TR', 'FontSize',16);

    %file_loc = find(ii_cfg.edf_file == '/');
    filename_pre = ii_trial{ii}(1).edf_file(1:end-4);
    saveas(fh, [filename_pre '_TR_matrix.png']);
    clear excl_matrix;

end

%% Align via space
horiz_trial_labels = {};
vert_trial_labels = {};
horiz_matrix = [];
vert_matrix = [];
for i = 1:size(ii_trial,1)
    excl_matrix = ii_trial{i}(1).tr_excl;
    for j = 1:size(excl_matrix, 2)
        rem = mod(j,4);
        if mod(rem,2) ~= 0
            if mod(rem,3) == 1
                horiz_matrix(end+1,:) = excl_matrix(:,j);
                horiz_trial_labels{end+1} = ['\newline' 'Run ' num2str(i) ', Trial ' num2str(j) ' - L2R'];
            else
                horiz_matrix(end+1,:) = flip(excl_matrix(:,j));
                horiz_trial_labels{end+1} = ['\newline' 'Run ' num2str(i) ', Trial ' num2str(j) ' - R2L'];
            end
        else
            if mod(rem,4) == 2
                vert_matrix(end+1,:) = excl_matrix(:,j);
                vert_trial_labels{end+1} = ['\newline' ' Run ' num2str(i) ', Trial ' num2str(j) ' - T2B'];
            else
                vert_matrix(end+1,:) = flip(excl_matrix(:,j));
                vert_trial_labels{end+1} = ['\newline' ' Run ' num2str(i) ', Trial ' num2str(j) ' - B2T'];
            end
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


% Plot matrix
subplot(size(horiz_matrix,1) + 4, 1, [1:size(horiz_matrix,1)]);
%figure;
x = [repmat(1:(tr_per_sweep + 1), size(horiz_matrix,1), 1)];
y = [repmat(1:size(horiz_matrix,1), tr_per_sweep + 1, 1)]';
map = [1 0 0; 0.8 0.8 0.8; 0.7 0.7 0.7; 0.6 0.6 0.6; 0.5 0.5 0.5;];
c = horiz_matrix;
fh = pcolor(x,y,c);
colormap(map);
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
    text(i + 0.3 ,1.5,num2str(horiz_sub_xlabels(i)), 'FontSize', 14);
end
set(gca,'xtick',[]);
set(gca,'ytick',[]);
xlabel('Included TRs, Excluded TRs', 'FontSize',16);


%file_loc = find(ii_cfg.edf_file == '/');
filename_pre = ii_trial{ii}(1).edf_file(1:end-6);
saveas(fh, [filename_pre 'BarHoriz_matrix.png']);

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


% X axis
xticks = [1:size(vert_matrix,2)];
set(gca, 'xTick', xticks);
vert_trial_labels{end+1} = '';
%vert_trial_labels = flip(vert_trial_labels);
set(gca, 'XTickLabel', vert_trial_labels, 'fontsize', 16);
xtickangle(90);

% Y axis
new_yticks = [2:tr_per_bar:tr_per_sweep];
set(gca,'YTick', new_yticks);
ylabels = 1:n_bars;
ylabels = flip(ylabels);
set(gca, 'YTickLabel', ylabels);
% Titles
ylabel('Bar Location', 'FontSize',16);
title('TRs Excluded, Vertical Alignment','FontSize',16);


% Create subplot
subplot(1, size(vert_matrix,2) + 4, [size(vert_matrix,2) + 3, size(vert_matrix,2) + 4]);
summary_matrix = summary_matrix';

% Set up matrix
sub_x = [repmat(1,1,2* n_bars + 1); repmat(2,1, 2* n_bars + 1);]';
sub_y = [repmat(1:(2 * n_bars + 1),2,1)]';
sub_c = summary_matrix;
fh_sub = pcolor(sub_x,sub_y,sub_c);
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


filename_pre = ii_trial{ii}(1).edf_file(1:end-6);
saveas(fh_sub, [filename_pre 'BarVert_matrix.png']);




return

end



