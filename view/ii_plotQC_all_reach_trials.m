function fhr = ii_plotQC_all_reach_trials(ii_trial_reach,ii_cfg,which_excl,fig_visible)
% ii_plotQC_alltrials Plots all trials, indicating response epoch, whether
% trial was excluded, why, etc. All timecourses (X and Y channels
% separately). 
%
% Usage:
%    ii_plotQC_alltrials(ii_trial_reach), where ii_trial_reach is a single run from
%    ii_scoreMGS.m, plots the timecourse of all trials (except for ITI,
%    last epoch), marking the response epoch, and highlights i_reach and
%    f_reach in thicker lines. Inidicates whether trial is excluded for any
%    particular reason. when ii_trial_reach is a run-concatenated variable
%    (output from ii_combineruns), loops over all runs and fills a figure
%    with each run as best as possible.
%
%    ii_plotQC_alltrials(ii_trial_reach,which_excl) uses the exclusion criteria
%    defined in which_excl rather than the default set (any) 
%
%    ii_plotQC_alltrials(ii_trial_reach,[],fig_visible) sets whether the plotted
%    figure(s) will be set to visible or whether they'll be hidden
%
%    fh = ii_plotQC_alltrials(...) returns handles to all figures
%
% Tommy Sprague, 6/12/2018

%% setup
% all possible exclusions and their string labels
all_excl = [11 12 13 20 21 22];
excl_labels = {'drift','calibration','delay fixation','no i_reach','bad i_reach','i_reach err'};
excl_labels_abbrev = {'d','c','fb','xi','bi','ei'};

SAMPLING_PER = 1/ii_cfg.hz;


if nargin < 3 || isempty(which_excl)
    which_excl = all_excl;
end

if nargin < 4 || isempty(fig_visible)
    fig_visible = 1;
end

if fig_visible == 0
    fig_arg = {'Visible','off'};
else
    fig_arg = {'Visible','on'};
end

% very rare we have more than 48 trials per run in any MGS experiment
MAX_NROWS = 8;
MAX_NCOLS = 6;

% how big to make the figure? look up the monitor size and start with that
tmp = get(groot,'MonitorPositions');
fig_size = tmp(1,:); clear tmp;

if isfield(ii_trial_reach,'r_num')
    ru = unique(ii_trial_reach.r_num);
else
    ru = 1;
end

% load plotting params
plot_params = ii_loadplotparams; % where we save things like EXCL_COLOR, etc...
% (note: can override them here if interested...)

fhr = [];

%% loop over runs
for rr = 1:length(ru)
    
    fhr(end+1) = figure(fig_arg{:},'Position',fig_size);
    
    % which trials are we considering right now?
    if isfield(ii_trial_reach,'r_num')
        thisidx = find(ii_trial_reach.r_num==ru(rr));
    else
        thisidx = 1:size(ii_trial_reach.i_reach,1);
    end
    
    % for this run, actual number of rows, cols determined based on using all rows
    nrows = min(MAX_NROWS,length(thisidx));
    ncols = min(MAX_NCOLS,ceil(length(thisidx)/nrows));
    
    for tt = 1:length(thisidx)
        subplot(nrows,ncols,tt);
        hold on;
        
        myt = (1:length(ii_trial_reach.rX{thisidx(tt)})) * SAMPLING_PER;
        
        
        % when did response epoch start? (XDAT=4)
        this_resp_start = myt(find(ismember(ii_trial_reach.XDAT{thisidx(tt)},ii_trial_reach.params.resp_epoch),1,'first'));
        this_resp_end   = myt(find(ismember(ii_trial_reach.XDAT{thisidx(tt)},ii_trial_reach.params.resp_epoch),1,'last'));
        
        %control for catch trials where there will be no resp epoch. geh
        %7/2/19
        if isempty(this_resp_start)
            this_resp_start = myt(find(ismember(ii_trial_reach.XDAT{thisidx(tt)},10),1,'first'));
            this_resp_end = myt(find(ismember(ii_trial_reach.XDAT{thisidx(tt)},10),1,'last'));
        else
        end

        % TARG positions
        plot([myt(1) myt(end)],[1 1]*ii_trial_reach.targ(thisidx(tt),1),'k--');
        plot([myt(1) myt(end)],[1 1]*ii_trial_reach.targ(thisidx(tt),2),'k--');
        plot([myt(1) myt(end)],[0 0],'k-');
        
        % plot start of response epoch
        plot(this_resp_start*[1 1],plot_params.MAXECC*[-1 1],'-','Color',[0.5 0.5 0.5]);
        plot(this_resp_end*[1 1],  plot_params.MAXECC*[-1 1],'-','Color',[0.5 0.5 0.5]);

        % X, Y eye positions
        plot(myt,ii_trial_reach.(ii_trial_reach.params.score_chans{1}){thisidx(tt)},'-','LineWidth',1,'Color',plot_params.RAW_COLORS(1,:));
        plot(myt,ii_trial_reach.(ii_trial_reach.params.score_chans{2}){thisidx(tt)},'-','LineWidth',1,'Color',plot_params.RAW_COLORS(2,:));

        
        % overlay the primary reach
        i_reach_t = (1:size(ii_trial_reach.i_reach_trace{thisidx(tt)},1)) * SAMPLING_PER + this_resp_start + ii_trial_reach.i_reach_rt(thisidx(tt));
        if ~isempty(ii_trial_reach.i_reach_trace{thisidx(tt)})
            plot(i_reach_t,ii_trial_reach.i_reach_trace{thisidx(tt)}(:,1),'-','LineWidth',2,'Color',plot_params.REACH_COLORS(1,:));
            plot(i_reach_t,ii_trial_reach.i_reach_trace{thisidx(tt)}(:,2),'-','LineWidth',2,'Color',plot_params.REACH_COLORS(1,:));
        end
        
        % and the final reach, if it exists
        if ii_trial_reach.n_reach(thisidx(tt)) >= 2
            f_reach_t = (1:size(ii_trial_reach.f_reach_trace{thisidx(tt)},1)) * SAMPLING_PER + this_resp_start + ii_trial_reach.f_reach_rt(thisidx(tt));
            if ~isempty(ii_trial_reach.f_reach_trace{thisidx(tt)})
                plot(f_reach_t,ii_trial_reach.f_reach_trace{thisidx(tt)}(:,1),'-','LineWidth',1.5,'Color',plot_params.REACH_COLORS(2,:));
                plot(f_reach_t,ii_trial_reach.f_reach_trace{thisidx(tt)}(:,2),'-','LineWidth',1.5,'Color',plot_params.REACH_COLORS(2,:));
            end
            
        end

        % mark if trial excluded

        % if actually exlcude (if ismember this trial, which_excl), red
        % if satisfies an exclusion criterion (even if we don't exclude),
        % italics
        
        if isfield(ii_trial_reach,'r_num')
            this_txt = sprintf('r%02.f, %02.f',ru(rr),ii_trial_reach.t_num(thisidx(tt)));
        else
            this_txt = sprintf('Trial %02.f',tt);
        end
        
        % add any/all exclusion criteria, w/ an asterisk if we use it to
        % exclude
        for ee = 1:length(ii_trial_reach.excl_trial{thisidx(tt)})
            if ee == 1
                this_txt = strcat(this_txt,'; ');
            end
            
            % which all_excl index?
            %tmpidx = find(ii_trial_reach.excl_trial{thisidx(tt)}(ee)==all_excl);
            
            this_txt = strcat(this_txt,excl_labels_abbrev{ii_trial_reach.excl_trial{thisidx(tt)}(ee)==all_excl});
            
            % if this is a useful exclusion criterion, add an asterisk
            if ismember(ii_trial_reach.excl_trial{thisidx(tt)}(ee),which_excl)
                this_txt = strcat(this_txt,'*');
            end
            
            % if not the last one, add a comma
            if ee~=length(ii_trial_reach.excl_trial{thisidx(tt)})
                this_txt = strcat(this_txt,', ');
            end
        end
        
        this_angle = 'normal';
        this_color = [0 0 0];
        if ~isempty(ii_trial_reach.excl_trial{thisidx(tt)})
            this_angle = 'italic';
            if any(ismember(ii_trial_reach.excl_trial{thisidx(tt)},which_excl))
                this_color = plot_params.EXCL_COLOR;
            end
        end
        
        text(0,plot_params.MAXECC,this_txt,'FontSize',10,'FontAngle',this_angle,'Color',this_color);
        
        
        this_trial_end = myt(find(ismember(ii_trial_reach.XDAT{thisidx(tt)},ii_trial_reach.params.resp_epoch(end)+1),1,'last'));
        
        %control for catch trials, which do not have same end epoch
        if isempty(this_trial_end)
            this_trial_end = 16;
        else
        end
        xlim([0 this_trial_end]);
        ylim(plot_params.MAXECC*[-1 1]);
        
        axis off;
        
        
        
        hold off;
    end
end

return