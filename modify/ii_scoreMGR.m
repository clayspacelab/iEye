function [ii_trial_reach,ii_cfg] = ii_scoreMGR(ii_data,ii_cfg,ii_reach,targ_coords,resp_epoch,fix_epoch,excl_criteria,save_chans,score_mode,align_to)
%ii_scoreMGS Extract typical parameters from each trial of an MGS dataset:
% - primary reach: first reach, exceeding some amplitude threshold,
%   after a response cue in the direction of one or more target positions
% - final reach: endpoint of final eye position before a feedback or
%   return-to-fixation stimulus appears
% - RT: time from beginning of go cue until primary reach is detected
%
% Then, knowing these paramters & features of the trial like target
% position(s), we also save out an aligned version of all reach endpoints
% and traces. 
%
% ii_trial_reach contains one entry per trial of each of several fields:
% - i_reach: n_trials x 2, ALIGNED X Y of primary reach endpoint
% - f_reach: n_trials x 2, ALIGNED X Y of final reach endpoint (note:
%   final reachade & primary reach endpoint can be same value when only a
%   single reach is detected)
% - n_reach: n_trials x 1, number of reach detected after and including  
%    primary reach(if no primary reach, this will be 0)
% - excl_trial: structure listing several exclusion criterion applying to
%   drift correction, epochs that are considered 'fixation', calibration
%   criteria, etc. 
% - i/f_reach_trace: cell array (n_trials x 1) with trace of initial/final
%   reachade through time, beginning at i_reach_rt
% - i/f_reach_err: euclidean distance between first target and endpoint of
%   initial/final reach
% - targ: n_trials x 2, extracted (or input) target coordinates for X,Y on 
%   each trial
%
% Because ii_data, ii_reach are not modified, we won't return them here.
%
% We extract the same fields as were extracted in ii_extractreachades
%
% Usage:
%  [ii_trial_reach, ii_cfg] = ii_scoreMGS(ii_data,ii_cfg,ii_reach) uses default
%  reach scoring parameters to extract primary and final eye positions
%  for a standard MGS-style task. Be default, looks for target coordinate
%  in TarX, TarY channels, extracts X,Y,Pupil and XDAT, and identifies
%  resp_epoch as first one following a the longest non-final epoch in a
%  trial (fixation epochs defined as all those previous). [TODO]
%
%  [ii_trial_reach, ii_cfg] = ii_scoreMGS(ii_data,ii_cfg,ii_reach,targ_coords)
%  uses coordinates specified by targ_coords instead of default (TarX,TarY
%  channels of ii_data). targ_coords can be a cell array of 2 strings, in
%  which case those channels are looked up in ii_data, a set of n_trials x
%  2n coordiantes, in which case pairs of columns are considered X,Y pairs
%  per trial (and we determine which of the n coords was chosen), or an array
%  of 2n numbers, in which case ii_cfg.trialinfo(:,targ_coords) will act as
%  X,Y pairs)
%
%  [ii_trial_reach, ii_cfg] = ii_scoreMGS(ii_data,ii_cfg,ii_reach,[],resp_epoch) only
%  looks at reach beginning and ending with XDAT within resp_epoch 
%  (integer or array of integers)
%
%  [ii_trial_reach, ii_cfg] = ii_scoreMGS(ii_data,ii_cfg,ii_reach,[],[],fix_epoch)
%  only considers fixation breaks when XDAT within fix_epoch (integer or
%  array of integers)
%
%  [ii_trial_reach, ii_cfg] = ii_scoreMGS(ii_data,ii_cfg,ii_reach,[],[],[],excl_criteria)
%  uses values stored within excl_criteria to determine primary reach
%  amplitude/duration threshold, error threshold, drift correction
%  threshold, and delay-period fixation break threshold (see below for
%  exact field names [TODO])
%
%  [ii_trial_reach, ii_cfg] = ii_scoreMGS(ii_data,ii_cfg,ii_reach,[],[],[],[],save_chans)
%  saves trial timecourse for each trial from each channel listed in
%  save_chans rather than the default set ('X','Y',Pupil and XDAT; note
%  that X/Y are derived from those used when extracting reachades)
%
% [ii_trial_reach, ii_cfg] =
% ii_scoreMGS(ii_data,ii_cfg,ii_reach,[],[],[],[],[],score_mode) determines
% whether a scored reach must begin AND end during indicated epoch(s)
% ('strict'), or just begin ('lenient'). Default is 'strict'.
%
% [ii_trial_reach, ii_cfg] =
% ii_scoreMGS(ii_data,ii_cfg,ii_reach,[],[],[],[],[],[],align_to) aligns
% i_reach_raw and f_reach_raw to coordinate in align_to (1x2). if absent,
% aligns only based on polar angle, DOES NOT adjust ecc!!! If input for X 
% is nan, also does not adjust ecc. (default: [NaN 0])
%
% Examples:
% TODO
%
% NOTES:
% - at present, we require at least one target coordinate to be defined; in
%   principle, that's not necessary - we can add functionality for that later
% - also, only a single target coordinate for now: ii_nearestcoord could be
%   called from here, but let's assume it's already been called
% - asusming coordinate pairs: X,Y (no support for 1D or 3D, etc, coords) -
%   this is consistent with ii_extractreachades
% - for aligning all targets, do we want to scale in polar coords, or
%   translate in cartesian? [discuss!]
%
% TODO:
% - add a i_reach, f_reach idx to ii_reach - it's zero if not a scored
%   reachade; trial number if it's either of those (useful for plotting main
%   sequences later on)
%
% TRIAL EXCLUSION CODES:
% first-digit:
% - 1 - trial-level exclusion (bad drift correction [11], calibration [12], or delay-
%       fixation break [13]
% - 2 - primary reach exclusion (no primary reach [20]; too small/short [21], large error [22])
%
%
% Tommy Sprague, 6/11/2018 - adapted from individual study
% _extractsaccadeData scripts (e.g., MGSMap_extractsaccadeData)

%% sanitize inputs & set default values

% just not enough inputs
if nargin < 3
    error('iEye:ii_scoreMGR:insufficientInputs','At minimum, to score MGR data you must submit ii_data,ii_cfg, and ii_reach (extracted reach data)');
end

% if targ_coords empty, use TarX, TarY
if nargin < 4 || isempty(targ_coords)
    targ_coords = {'TarX','TarY'}; %continuing to use XDAT TarX TarY coords 
end

% if no response epoch defined:
% deduce which epoch to use for response...first, we look for delay,
% which we define as the longest epoch >= 1 and <= max(epoch(within-tiral))
% (this assumes that the ITI is the last epoch...)
if nargin < 5 || isempty(resp_epoch)
    tmp_iti = max(ii_data.XDAT(ii_cfg.trialvec~=0));
    tmp_delay = mode(ii_data.XDAT(ii_cfg.trialvec~=0 & ii_data.XDAT~=tmp_iti));
    resp_epoch = tmp_delay + 1;
    clear tmp_delay tmp_iti;
end

% by default, use all epochs up until resp_epoch
if nargin < 6 || isempty(fix_epoch)
    fix_epoch = 1:(resp_epoch-1);
end

% these are default exclusion criteria, excl_criteria input only needs a
% subset of these, if a field is missing it'll be overwritten w/ a default
% value

excl_default.i_dur_thresh = 150; % must be shorter than 150 ms
excl_default.i_amp_thresh = 5;   % must be longer than 5 dva [if FIRST reachade in denoted epoch is not at least this long and at most this duration, drop the trial]
excl_default.i_err_thresh = 5;   % i_reach must be within this many DVA of target pos to consider the trial

excl_default.drift_thresh = 2.5;     % if drift correction norm is this big or more, drop
excl_default.delay_fix_thresh = 2.5; % if any fixation is this far from 0,0 during delay (epoch 3)


% if excl_criteria is empty, make it a struct and insert default into
% missing values
if nargin < 7 || isempty(excl_criteria)
    excl_criteria = struct(); % empty struct
end
    
excl_fields = fields(excl_default);
for ee = 1:length(excl_fields)
    if ~ismember(excl_fields{ee},fields(excl_criteria))
        excl_criteria.(excl_fields{ee}) = excl_default.(excl_fields{ee});
    end
end
% now, we have excl_criteria w/ all the wanted fields, and default values
% filling in those not provided originally
clear excl_default; % make sure we don't accidentally use this guy again


if nargin < 9 || isempty(score_mode)
    score_mode = 'strict';
elseif ~ismember(score_mode,{'strict','lenient'})
    error('iEye:ii_scoreMGR:invalidMode','Only strict and lenient allowed as score_mode inputs');
end

if nargin < 10 || isempty(align_to)
    align_to = [NaN 0];
end


%% Figure out which channels we're looking at
% Channels to extract is defined for ii_extractreachades, so we need to
% figure out what those fields were - can pull it out of field names in
% ii_reach (use $$_trace) - look for whatever comes before _trace (sorry, a
% bit ugly below, but cellfun's weren't working the way I wanted)
tmp_fields = fields(ii_reach);
tmp_idx = find(cellfun(@any,strfind(tmp_fields,'_trace')));
which_chans = cell(length(tmp_idx),1);
for cc = 1:length(tmp_idx)
    which_chans{cc} = tmp_fields{tmp_idx(cc)}(1:(strfind(tmp_fields{tmp_idx(cc)},'_trace')-1));
end
clear tmp_fields tmp_idx;


if nargin < 8 || isempty(save_chans)
    save_chans = {which_chans{:},'Pupil','XDAT','Velocity'};
end


%% NOW: start extracting data!

%% Build ii_trial_reach; fill with nothing
for chan_idx = 1:length(save_chans)
    ii_trial_reach.(save_chans{chan_idx}) = cell(ii_cfg.numtrials,1);
end

% also want to store the raw coordinates
ii_trial_reach.i_reach_raw = nan(ii_cfg.numtrials,2);
ii_trial_reach.f_reach_raw = nan(ii_cfg.numtrials,2);

% aligned coordinates
ii_trial_reach.i_reach = nan(ii_cfg.numtrials,2);
ii_trial_reach.f_reach = nan(ii_cfg.numtrials,2);


ii_trial_reach.i_reach_err = nan(ii_cfg.numtrials,1);
ii_trial_reach.f_reach_err = nan(ii_cfg.numtrials,1);

ii_trial_reach.n_reach = nan(ii_cfg.numtrials,1); % how many reachades are there total? (in each epoch maybe?)
ii_trial_reach.n_reach_epoch = nan(ii_cfg.numtrials,5); % TODO: fill this w/ n_trials x n_epochs (exclude ITI)

ii_trial_reach.i_reach_rt = nan(ii_cfg.numtrials,1); % latency from go cue to each of these
ii_trial_reach.f_reach_rt = nan(ii_cfg.numtrials,1);

ii_trial_reach.i_reach_trace = cell(ii_cfg.numtrials,1);
ii_trial_reach.f_reach_trace = cell(ii_cfg.numtrials,1);

% save some calibration, drift correction info for convenience
if isfield(ii_cfg,'calibrate')
    ii_trial_reach.calib_amt = ii_cfg.calibrate.amt;
    ii_trial_reach.calib_adj = ii_cfg.calibrate.adj;
    ii_trial_reach.calib_err = ii_cfg.calibrate.err;
end

if isfield(ii_cfg,'drift')
    ii_trial_reach.drift_amt = ii_cfg.drift.amt;
end

ii_trial_reach.excl_trial = cell(ii_cfg.numtrials,1);  % why is this trial excluded? each cell includes several markers


% let's copy over ii_cfg.trialinfo, if it exists
if isfield(ii_cfg,'trialinfo')
    ii_trial_reach.trialinfo = ii_cfg.trialinfo;
end

% add parameters used for extracting reachades: ii_trial_reach.params (as they
% were input/sanitized - but, for e.g., targ_coord, not updated [?])
ii_trial_reach.params.excl_criteria = excl_criteria;
ii_trial_reach.params.resp_epoch  = resp_epoch;
ii_trial_reach.params.fix_epoch   = fix_epoch;
ii_trial_reach.params.targ_coords = targ_coords;
ii_trial_reach.params.save_chans  = save_chans;
ii_trial_reach.params.score_mode  = score_mode;
ii_trial_reach.params.score_chans = which_chans; % the two channels that were used for scoring


%% extract targ_coords for each trial

targ_coords_extracted = nan(ii_cfg.numtrials,2);

% if targ_coords is n_trials x 2, then just extract from there
if size(targ_coords,1)==ii_cfg.numtrials && size(targ_coords,2) == 2
    targ_coords_extracted = targ_coords;

% if a cell, then we need to look in ii_data at these fields for each trial
% we assume that the target coordinate is given by the epoch immediately
% after the end of the response epoch (resp_epoch(end)+1)
% always takes the mode - so uses whatever is present during resp_epoch of
% this channel most commonly (to account for mis-sampling, etc)
elseif iscell(targ_coords)
    for tt = 1:ii_cfg.numtrials
        thisidx = ii_cfg.trialvec==tt & ismember(ii_data.XDAT,resp_epoch(end)+1);
        for cc = 1:length(targ_coords)
            targ_coords_extracted(tt,cc) = mode(ii_data.(targ_coords{cc})(thisidx));
        end
    end
    
% if 2 #'s, pull info out of those columns ii_cfg.trialinfo
elseif numel(targ_coords)==2
    if max(targ_coords) > size(ii_cfg.trialinfo,2)
        error('iEye:ii_scoreMGS:invalidInput','Column index submitted for targ_coords exceeds size of ii_cfg.trialinfo (%i columns)',size(ii_cfg.trialinfo,2));
    end
    targ_coords_extracted = ii_cfg.trialinfo(:,targ_coords);
else
    error('iEye:ii_scoreMGS:invalidInput','Invalid input for targ_coords: input should be a pair of channels (cell array of strings), list of coordinates (numtrials x 2), or pair of columns into ii_cfg.trialinfo');
end

% save these in ii_trial_reach as ii_trial_reach.targ
ii_trial_reach.targ = targ_coords_extracted;



% get idx reachades to look at (those that start [end] in resp_epoch), and
% corresponding trial index
if strcmpi(score_mode,'strict')
    which_reach = find(ismember(ii_reach.epoch_start,resp_epoch) & ismember(ii_reach.epoch_end,resp_epoch));
elseif strcmpi(score_mode,'lenient')
    % if lenient, just needs to start in this epoch
    which_reach = find(ismember(ii_reach.epoch_start,resp_epoch));
end

which_trials = ii_reach.trial_start(which_reach); % trial number for each reachade in which_reach

%% loop over trials
for tt = 1:ii_cfg.numtrials
    
    % save the data from each channel from each trial
    for chan_idx = 1:length(save_chans)
        ii_trial_reach.(save_chans{chan_idx}){tt} = ii_data.(save_chans{chan_idx})(ii_cfg.trialvec==tt);
    end
    
    % time that relevant epoch of trial started
    t_start = find(ii_cfg.trialvec==tt & ismember(ii_data.XDAT,resp_epoch) ,1,'first')/ii_cfg.hz; 
    
    %% score reachades: find primary/initial; final reachade
    
    % all reachades that started (and ended?) in resp_epoch on this trial
    this_i_reach = which_reach(which_trials==tt);
    
    % extract amplitude of each of these
    this_i_amp = ii_reach.amplitude(this_i_reach);
    this_i_dur = ii_reach.duration(this_i_reach);
    
    
    % if amplitude too small or duration too short (for ALL
    % reachades in response epoch
    
    % first, if no reachade detected (if i_reach is empty),
    % record as "20" - no identified reachades (whatosever!)
    if isempty(this_i_reach)
        ii_trial_reach.excl_trial{tt}(end+1) = 20; % no reachades identified in response epoch
    elseif ~any(this_i_dur<=excl_criteria.i_dur_thresh & this_i_amp>=excl_criteria.i_amp_thresh)  % if none of the reach pass the amplitude & duration criteria
        ii_trial_reach.excl_trial{tt}(end+1) = 21; % none of the identified reaches (>0) passed exclusion criteria for primary reachade
    else
        % ok, now we can use the first of the reachades that do
        % pass exclusion criteria as the 'primary' reachade
        
        % find the first of this_i_reach for which amp & dur
        % pass threshold
        this_i_reach = this_i_reach(this_i_amp>=excl_criteria.i_amp_thresh & this_i_dur <= excl_criteria.i_dur_thresh);
        
        % just index into the first element of this_i_reach... [TODO: use
        % which_chan to compute this....]
        ii_trial_reach.i_reach_raw(tt,:) = [ii_reach.rX_end(this_i_reach(1)) ii_reach.rY_end(this_i_reach(1))];
        ii_trial_reach.i_reach_rt(tt) = ii_reach.t(this_i_reach(1),1)-t_start;
        ii_trial_reach.i_reach_trace{tt} = [ii_reach.rX_trace{this_i_reach(1)} ii_reach.rY_trace{this_i_reach(1)}];
        
    end
    clear this_i_amp this_i_dur;
    
    % same for final reach
    this_f_reach = which_reach(find(which_trials==tt,1,'last'));
    if ~isempty(this_f_reach)
        ii_trial_reach.f_reach_raw(tt,:) = [ii_reach.rX_end(this_f_reach) ii_reach.rY_end(this_f_reach)];
        ii_trial_reach.f_reach_rt(tt) = ii_reach.t(this_f_reach,1)-t_start;
        ii_trial_reach.f_reach_trace{tt} = [ii_reach.rX_trace{this_f_reach} ii_reach.rY_trace{this_f_reach}];
    end
    
    % count reachades in this trial during relevant epoch AFTER
    % primary reachade, up to and including final reachade
    
    % first reachade: primary
    
    % last reachade: final
    
    % so, I think we just need the length of i_reach:f_reach
    if ~isempty(this_i_reach) && ~isempty(this_f_reach)
        % in most cases, we'll identify both i_reach & f_reach -
        % use this algo
        ii_trial_reach.n_reach(tt) = length(this_i_reach:this_f_reach);
    else
        % if either of them is undefined (or both), this will
        % give us the # of reachades: 1 if only i_reach or only
        % f_reach, 0 if neither
        ii_trial_reach.n_reach(tt) = sum([~isempty(this_i_reach) ~isempty(this_f_reach)]);
    end

    %% trial exclusions: find reasons we may want to exclude each trial
    
    % ~~~~~ FIRST: exclude based on trial-level features (see above)
    
    % note: only reject trials based on these criteria if those steps were
    % run!
    
    % DRIFT CORRECTION TOO BIG
    if isfield(ii_cfg,'drift')
        if sqrt(sum(ii_cfg.drift.amt(tt,:).^2)) > excl_criteria.drift_thresh
            ii_trial_reach.excl_trial{tt}(end+1) = 11;
        end
    end
    
    % CALIBRATION OUTSIDE OF RANGE
    if isfield(ii_cfg,'calibrate')
        if ii_cfg.calibrate.adj(tt)~=1
            ii_trial_reach.excl_trial{tt}(end+1) = 12;
        end
    end
    
    % DURING DELAY, FIXATION OUTSIDE OF RANGE
    
    % find fixations in this trial; epoch [TODO: make sure I'm using the
    % right channels here!!]
    this_fix_idx = ii_cfg.trialvec==tt & ismember(ii_data.XDAT,fix_epoch);
    if max(sqrt(ii_data.rX_fix(this_fix_idx).^2+ii_data.rY_fix(this_fix_idx).^2)) > excl_criteria.delay_fix_thresh
        ii_trial_reach.excl_trial{tt}(end+1) = 13;
    end
    
    % I_reach ERROR TOO HIGH (use targ_coords_extracted)
    if ~isempty(this_i_reach)
        if sqrt(sum((ii_trial_reach.i_reach_raw(tt,:)-targ_coords_extracted(tt,:)).^2,2)) > ...
                excl_criteria.i_err_thresh
            ii_trial_reach.excl_trial{tt}(end+1) = 22;
        end
    end
    
    clear this_i_reach this_f_reach;
    
    
end

%% compute error for initial, final reachade based on extracted targ coords
% error for initial (note: mean absolute error..)
ii_trial_reach.i_reach_err = sqrt(sum((ii_trial_reach.i_reach_raw-targ_coords_extracted).^2,2));

% error for final
ii_trial_reach.f_reach_err = sqrt(sum((ii_trial_reach.f_reach_raw-targ_coords_extracted).^2,2));


%% rotate all trials such that i_reach, f_reach aligned to targ

align_fields = {'i_reach','f_reach'}; % use these_raw, align them

for aa = 1:length(align_fields)
    
    % rotate all so that i_reach, f_reach are at known position... (align_to)
    [tmpth,tmpr] = cart2pol(ii_trial_reach.(sprintf('%s_raw',align_fields{aa}))(:,1),  ii_trial_reach.(sprintf('%s_raw',align_fields{aa}))(:,2));
    [adjth,adjr] = cart2pol(targ_coords_extracted(:,1),targ_coords_extracted(:,2)); % for now, we don't use radius here for anything...
    tmpth = tmpth - adjth;  % TODO: use align_to?
    [ii_trial_reach.(align_fields{aa})(:,1),ii_trial_reach.(align_fields{aa})(:,2)] = pol2cart(tmpth,tmpr);
    
    % adjust x if align_to is not nan [NOTE: this is a translation, not a
    % scaling]
    if ~isnan(align_to(1))
        ii_trial_reach.(align_fields{aa})(:,1) = ii_trial_reach.(align_fields{aa})(:,1) - (adjr-align_to(1));
    end
    clear tmpth tmpr adjth adjr;
    
end


return