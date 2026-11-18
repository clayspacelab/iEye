function [ii_trial,ii_cfg] = rsvp_prf_collectTrials(ii_data,ii_cfg,ii_sacc,fix_epoch,excl_criteria,save_chans)
%RSVP_PRF_COLLECTTRIALS Collect each trial's data and trial-level exclusion flags
%
%   [ii_trial,ii_cfg] = rsvp_prf_collectTrials(ii_data,ii_cfg,ii_sacc,fix_epoch)
%   [ii_trial,ii_cfg] = rsvp_prf_collectTrials(...,fix_epoch,excl_criteria)
%   [ii_trial,ii_cfg] = rsvp_prf_collectTrials(...,fix_epoch,excl_criteria,save_chans)
%
%   For the RSVP and RDK pRF tasks, one trial is one bar sweep (as defined
%   by ii_definetrial in rsvp_prf_preproc). Adapted from iEye's ii_scoreMGS,
%   without its target, response and saccade-scoring features.
%
% Inputs:
%   fix_epoch:     XDAT value(s) during which the subject should be
%                  fixating (required); used for exclusion code 13
%   excl_criteria: struct with any of these fields (default in brackets):
%                  - drift_thresh [2.5]: deg; flag the trial (code 11) if
%                    its drift correction was larger than this
%                  - fix_thresh [2.5]: deg; flag the trial (code 13) if any
%                    fixation during fix_epoch was farther than this from
%                    (0,0)
%   save_chans:    channels to save for each trial (default: the channels
%                  saccades were extracted from (from ii_sacc), plus Pupil,
%                  XDAT and Velocity)
%
% Output ii_trial:
%   <chan>:     n_trials x 1 cell, each trial's samples of each save_chans
%               channel
%   drift_amt:  n_trials x 2 drift correction applied to each trial (if
%               drift correction was run)
%   excl_trial: n_trials x 1 cell, exclusion codes flagged for each trial:
%               11 = drift correction larger than drift_thresh
%               13 = fixation during fix_epoch farther than fix_thresh
%                    from (0,0)
%   trialinfo:  copy of ii_cfg.trialinfo, if present
%   params:     fix_epoch, excl_criteria and save_chans used
%   edf_file, tr_excl, sweep_dir: run-level info used by
%               rsvp_prf_plotQC_exclusions
%
% ii_cfg is returned unchanged.
%
% See also rsvp_prf_preproc, rsvp_prf_plotQC_exclusions, ii_combineruns.

% Tommy Sprague, 6/11/2018 (ii_scoreMGS); adapted for RSVP/RDK pRF tasks

%% check inputs & set defaults
if nargin < 4 || isempty(fix_epoch)
    error('rsvp_prf_collectTrials:noFixEpoch', ...
        'fix_epoch (XDAT value(s) during which the subject should fixate) must be given');
end

excl_default.drift_thresh = 2.5; % deg
excl_default.fix_thresh   = 2.5; % deg
if nargin < 5 || isempty(excl_criteria)
    excl_criteria = struct();
end
unknown = setdiff(fieldnames(excl_criteria), fieldnames(excl_default));
if ~isempty(unknown)
    error('rsvp_prf_collectTrials:unknownCriteria', 'Unknown excl_criteria field(s): %s (allowed: %s)', ...
        strjoin(unknown, ', '), strjoin(fieldnames(excl_default), ', '));
end
default_fields = fieldnames(excl_default);
for ee = 1:numel(default_fields)
    if ~isfield(excl_criteria, default_fields{ee})
        excl_criteria.(default_fields{ee}) = excl_default.(default_fields{ee});
    end
end

% channels saccades were extracted from (ii_sacc has a <chan>_trace field
% for each), e.g. {'X','Y'}
sacc_fields = fieldnames(ii_sacc);
is_trace = endsWith(sacc_fields, '_trace');
which_chans = extractBefore(sacc_fields(is_trace), '_trace');

if nargin < 6 || isempty(save_chans)
    save_chans = [which_chans(:).' {'Pupil','XDAT','Velocity'}];
end

% fixation channels for the fixation check (made by ii_findfixations)
fix_chans = strcat(which_chans, '_fix');
if ~all(isfield(ii_data, fix_chans))
    error('rsvp_prf_collectTrials:noFixChans', 'ii_data is missing fixation channel(s) %s; run ii_findfixations first', ...
        strjoin(fix_chans(~isfield(ii_data, fix_chans)), ', '));
end


%% collect trials
n_trials = ii_cfg.numtrials;
for cc = 1:numel(save_chans)
    ii_trial.(save_chans{cc}) = cell(n_trials, 1);
end

has_drift = isfield(ii_cfg, 'drift');
if has_drift
    ii_trial.drift_amt = ii_cfg.drift.amt;
end

ii_trial.excl_trial = cell(n_trials, 1);

if isfield(ii_cfg, 'trialinfo')
    ii_trial.trialinfo = ii_cfg.trialinfo;
end

ii_trial.params.fix_epoch     = fix_epoch;
ii_trial.params.excl_criteria = excl_criteria;
ii_trial.params.save_chans    = save_chans;

% distance of the current fixation from (0,0) at every sample
fix_dist = zeros(size(ii_data.(fix_chans{1})));
for cc = 1:numel(fix_chans)
    fix_dist = fix_dist + ii_data.(fix_chans{cc}).^2;
end
fix_dist = sqrt(fix_dist);
in_fix_epoch = ismember(ii_data.XDAT, fix_epoch);

for tt = 1:n_trials
    in_trial = ii_cfg.trialvec == tt;

    % each trial's data from each channel
    for cc = 1:numel(save_chans)
        ii_trial.(save_chans{cc}){tt} = ii_data.(save_chans{cc})(in_trial);
    end

    % 11: drift correction too big
    if has_drift && sqrt(sum(ii_cfg.drift.amt(tt,:).^2)) > excl_criteria.drift_thresh
        ii_trial.excl_trial{tt}(end+1) = 11;
    end

    % 13: fixation away from center while subject should be fixating
    if any(fix_dist(in_trial & in_fix_epoch) > excl_criteria.fix_thresh)
        ii_trial.excl_trial{tt}(end+1) = 13;
    end
end


%% run-level info (for rsvp_prf_plotQC_exclusions)
ii_trial.edf_file = ii_cfg.edf_file;
ii_trial.tr_excl  = ii_cfg.tr_excl;
if isfield(ii_cfg, 'task_timing')
    ii_trial.sweep_dir = ii_cfg.task_timing.sweep_dir; % n_trials x 1 direction labels
end

end
