function [ f_han ] = rsvp_prf_plottimeseries( ii_data, ii_cfg, which_chans, varargin )
%ii_plottimeseries Plot full timeseries in ii_data
%   By default, will try to show trial boundaries (vertical lines), XDAT
%   epoch codes (colored bar), and current selections. 
%
% ii_plottimeseries(ii_data,ii_cfg) will try to plot the most verbose
% possible depiction of ii_data. If trials have been defined, they'll be
% delineated by light gray vertical lines. If saccades have been
% identified, they'll be indicated by black triangles. If fixations have
% been identified, they'll be drawn as black lines over viewed channels.
% Default viewed channels are {'X','Y'}
%
% ii_plottimeseries(ii_data,ii_cfg,which_chans) plots only the channels
% within which_chans (cell array of strings, or a string)
%
% f_han = ii_plottimeseries(...) returns a figure handle to plot generated
%
% several plotting arguments available (use as many or few, in any order):
% ii_plottimeseries(..., 'noselections') disables drawing of selections
% ii_plottimeseries(..., 'nosaccades') disables drawing of saccades
% ii_plottimeseries(..., 'notrials') disables drawing of trial boundaries
% ii_plottimeseries(..., 'noxdat') disables plotting of epoch colorbar
% ii_plottimeseries(..., 'nofixations') disables plotting of fixations
% ii_plottimeseries(..., 'YLim', [min max]) manually sets y limits (before
% colorbar, if that's used)
% ii_plottimeseries(..., 'nofigure') disables visible plotting of the
% figure; useful when saving out figure from outside script
%
% All plot options above are automatically set based on fields in ii_cfg,
% ii_data. 

% Tommy Sprague, 8/17/17
% TODO: all plotting options!!!!

% a few constants for plotting:
TRIAL_BORDER_COLOR = [0.4 0.4 0.4];
SELECTION_COLOR = [0.6 0 0]; % DARK RED
SELECTION_ALPHA = 0.25;

SACCADE_MARKER = 'v';
SACCADE_MARKERSIZE = 4;
SACCADE_COLOR = [1.0 1.0 1.0] * 0.3;

FIXATION_COLOR = [0 0 0];

XDAT_COLORMAP = 'lines';

% 80% wide, 0.25 tall as wide
scr_size = get(0,'Screensize');
FIG_POSITION = [[0.05 0.5-0.22*0.9]*scr_size(3) [0.9 0.9*0.22]*scr_size(3)];





if nargin < 3 || isempty(which_chans)
    which_chans = {'X','Y','TarX','TarY'}; %Tarx TarY added on jun17 
end

if ~iscell(which_chans), which_chans = {which_chans}; end

% TODO: make sure which_chans exists

%which_chans = {'X','Y','TarX','TarY', 'PupilZ'};
chan_colors = lines(length(which_chans));
if size(chan_colors,1) >= 4 % swap colors of channels 3 and 4 (stim X, Y)
    tmp = chan_colors(3,:);
    chan_colors(3,:) = chan_colors(4,:);
    chan_colors(4,:) = tmp;
end
myt = (1:length(ii_data.(which_chans{1}))).' * 1/ii_cfg.hz;


% Timing and position of every bar step (RSVP or RDK task); see
% rsvp_prf_extract_task_timing. Saved by rsvp_prf_proc; recomputed for files
% preprocessed before it was.
if isfield(ii_cfg, 'task_timing')
    task_timing = ii_cfg.task_timing;
else
    task_timing = rsvp_prf_extract_task_timing(ii_data, ii_cfg, ii_cfg.rsvp_prf_params);
end

% Set screen edges for graph (stimulus area, deg)
stim_max_deg = max(task_timing.stim_extent_deg);
myy = [-stim_max_deg/2, stim_max_deg/2];

% Stimulus trace: bar center on the axis the bar moves along during each bar
% step; NaN (not drawn) when there's no sweep and on the axis that doesn't move
TarX = nan(size(myt));
TarY = nan(size(myt));
for rr = 1:numel(task_timing.sweep)
    idx = max(task_timing.start_idx(rr), 1):min(task_timing.end_idx(rr), length(myt));
    if ismember(task_timing.dir{rr}, {'L2R','R2L'})
        TarX(idx) = task_timing.tar_x(rr);
    else
        TarY(idx) = task_timing.tar_y(rr);
    end
end

% break the trace between sweeps (RDK sweeps run back to back, with no ITI)
sweep_on = task_timing.start_idx(task_timing.step == 1);
sweep_on = sweep_on(sweep_on >= 1 & sweep_on <= length(myt));
TarX(sweep_on) = NaN;
TarY(sweep_on) = NaN;

ii_data.TarX = TarX;
ii_data.TarY = TarY;






% draw the channels

if ismember('nofigure',varargin)
    f_han = figure('visible','off');
else
    f_han = figure;
end
hold on; p_han = [];
for cc = 1:length(which_chans)
    if length(ii_data.(which_chans{cc})) > length(myt)
        ii_data.(which_chans{cc}) = ii_data.(which_chans{cc})(1:length(myt));
    end
    p_han(end+1) = plot(myt,ii_data.(which_chans{cc}),'-','LineWidth',1.5,'Color',chan_colors(cc,:));
    
    % if there's a corresponding fixation channel and we don't say no, plot
    % the fixation on top
    if ismember(sprintf('%s_fix',which_chans{cc}),fieldnames(ii_data)) && ~ismember('nofixations',varargin)
        plot(myt,ii_data.(sprintf('%s_fix',which_chans{cc})),'-','Color',FIXATION_COLOR);
    end
    if cc == 4
     legend(p_han, {'X','Y','Stim X','Stim Y'}, 'AutoUpdate', 'off')
    end
end

 %added jun 17 
% hmm...this isn't working, keeps adding other elements. I just want
% which_chans lines...
% if ~ismember('nolegend',varargin)
%     legend(p_han,which_chans,'location','eastoutside');
% end

% optimize y-lims
%myy = get(gca,'YLim');





% if trials defined and user doesn't override, draw trial boundaries
% TODO: put these in the back?
if ismember('tcursel',fieldnames(ii_cfg)) && ~ismember('notrials',varargin)
    
    
    for tt = 1:size(ii_cfg.tcursel,1)
        
        plot(myt(ii_cfg.tcursel(tt,1)) * [1 1], myy, '-', 'Color', TRIAL_BORDER_COLOR );
        
        if tt == size(ii_cfg.tcursel,1)
            plot(myt(ii_cfg.tcursel(tt,2)) * [1 1], myy, '-', 'Color', TRIAL_BORDER_COLOR );
        end
    end    
end



% draw saccade markers (if possible, not forbidden)
if ismember('saccades',fieldnames(ii_cfg)) && ~ismember('nosaccades',varargin)
    
    for ss = 1:size(ii_cfg.saccades,1)
        
        plot(mean(myt(ii_cfg.saccades(ss,:))),myy(2)-0.25,SACCADE_MARKER,'MarkerSize',SACCADE_MARKERSIZE,'Color',SACCADE_COLOR,'MarkerFaceColor',SACCADE_COLOR);
        
    end
end




% Select trs to be excluded, mark as red transluscent patch
for i = 1:size(ii_cfg.bad_tr_idc,1)
    p1 = ii_cfg.bad_tr_idc(i,1);
    p2 = ii_cfg.bad_tr_idc(i,2);
    if ii_cfg.bad_tr_idc(i,2) > length(myt)
        p2 = length(myt);
    end
    v1 = (1/ii_cfg.hz)*[-0.5 0.5 0.5 -0.5]+myt([p1 p2 p2 p1]).'; 
    patch(v1,[myy(1) myy(1) myy(2) myy(2)],SELECTION_COLOR,'linewidth',1,'FaceAlpha',SELECTION_ALPHA,'EdgeAlpha',0);
    clear p1 p2 v1;
end



% Adaptation for rsvp prf data
% Draw patches for each TR within each bar step, one gray level per TR
tr_names = fieldnames(ii_cfg.tr_bricks);
tr_names = tr_names(~cellfun(@isempty, regexp(tr_names, '^tr_\d+$')));
tr_colors = num2cell(repmat(linspace(0.8, 0.4, max(numel(tr_names), 2)).', 1, 3), 2);
for i = 1:numel(tr_names)
    tr_name = tr_names{i};
    for j = 1:size(ii_cfg.tr_bricks.(tr_name), 1)
        p1 = ii_cfg.tr_bricks.(tr_name)(j,1);
        p2 = ii_cfg.tr_bricks.(tr_name)(j,2);
        if ii_cfg.tr_bricks.(tr_name)(j,2) > length(myt)
            p2 = length(myt);
        end
        if ii_cfg.tr_bricks.(tr_name)(j,1) > length(myt)
            p1 = length(myt);
        end
        v1 = (1/ii_cfg.hz)*[-0.5 0.5 0.5 -0.5]+myt([p1 p2 p2 p1]).'; 
        patch(v1,[myy(1) myy(1) myy(2) myy(2)],tr_colors{i},'linewidth',1,'FaceAlpha',SELECTION_ALPHA,'EdgeAlpha',0);
        clear p1 p2 v1;
    end
end



%{
% draw XDAT ribbon
if ismember('XDAT',fieldnames(ii_data)) && ~ismember('noxdat',varargin)
    
    image(myt,myy(1)+0.5,ii_data.XDAT.');
    colormap(XDAT_COLORMAP);
    
end
%}



% make figure big (80% screen width, 0.25*width tall)
set(gcf,'Position',FIG_POSITION)

% add a title
[~, fn_name, fn_ext] = fileparts(ii_cfg.edf_file);
fn_fortitle = [fn_name fn_ext];
set(gcf,'NumberTitle','off','Name',fn_fortitle);


% fill up figure better, optimize display
set(gca,'Position',[0.05 0.2 0.9 0.75],'TickDir','out','LineWidth',1.5,'FontSize',14,'TickLength',[0.0075 0.0075]);

% label axes, etc
xlabel('Time (s)');
ylabel('Distance from Center (Degrees Visual Angle)');
xlim([myt(1) myt(end)]);

end

