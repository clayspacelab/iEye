% rsvp_prf_analysis.m
%
% example analysis script for example RSVP fMRI dataset (Subjects/example_sub_*)
%
% Here's an example data analysis script that walks from raw EDFs for a
% RSVP task prf mapping session.


%make sure edf2asc is on the path
EDF_ASC_PATH = '/usr/local/bin';

PATH = getenv('PATH');
%kinda dumb path check
if ~contains(PATH,EDF_ASC_PATH)
    setenv('PATH',[PATH ':' EDF_ASC_PATH]);
    PATH = getenv('PATH');
end

% Configure paths & file locations
tmp = mfilename('fullpath'); tmp2 = strfind(tmp,filesep);
addpath(genpath(tmp(1:(tmp2(end-1)))));
root = tmp(1:(tmp2(end)-1));
root = strcat(root,'/Subjects/');


% Select subject to run
%edf_prefix = 'example_sub_1/';
edf_prefix = 'example_sub_2/';

edf_files = dir(fullfile(root,sprintf('%s*.edf',edf_prefix)));
%edf_files = dir(fullfile(root,sprintf('%s*_R5.edf',edf_prefix)));

%load prf task parameters from the rsvp_params.txt settings file.
%ALTERNATIVELY use rsvp_read_params.m to load task parameters from the
%*_summary.csv from the rsvp task output.
rsvp_params = rsvp_prf_loadparams();

% very first thing we want to do is define all parameters for processing
% (see ii_loadparams.m for default values)
% ifg_fn = 'home/jeff/docs/prf_ieye/iEye/prf/p_500hz.ifg';
config_root = tmp(1:(tmp2(end-1)-1));
ifg_fn = [config_root '/ifg/p_500hz.ifg'];
%ifg_fn = [config_root '/ifg/p_1000hz.ifg'];
ii_params = ii_loadparams; % load default set of analysis parameters, only change what we have to
ii_params.resolution = rsvp_params.screen_res; % gaze is recorded in screen pixels; ii_rescale centers on resolution/2
ii_params.valid_epochs =[1 2];
ii_params.trial_end_value = 2;   % XDAT value for trial end
ii_params.fix_epoch = [2];% XDAT values for fixation epoch
ii_params.drift_epoch = [2]; % XDAT values for drift correction
ii_params.drift_fixation_mode  = 'mode'; % how to choose which fixation?
% ii_params.calibrate_epoch = 1;   % XDAT value for when we calibrate (feedback stim)
% ii_params.calibrate_select_mode = 'last'; % how do we select fixation with which to calibrate?
% ii_params.calibrate_mode = 'scale'; % scale: trial-by-trial, rescale each trial; 'run' - run-wise polynomial fit
ii_params.blink_window = [200 200]; % how long before/after blink (ms) to drop?
ii_params.plot_epoch = [1 2];  % what epochs do we plot for preprocessing?
% ii_params.calibrate_limits = [2.5]; % when amount of adj exceeds this, don't actually calibrate (trial-wise); ignore trial for polynomial fitting (run)
ii_params.ppd = rsvp_params.ppd; %31.8578; % for scanner, 1280 x 1024 - convert pix to DVA
%ii_params.sacc_duration_thresh = 0.01;
ii_params.sacc_amplitude_thresh = 1.25;


% create empty cell array of all our trial data for combination later on
ii_trial = cell(length(edf_files),1);

for ff = 1:length(edf_files)
    
    % what is the output filename?
    preproc_fn = sprintf('%s/%s_preproc.mat',edf_files(ff).folder,edf_files(ff).name(1:(end-4)));
    
    
    % run preprocessing!
     [ii_data, ii_cfg, ii_sacc] = rsvp_prf_preproc(fullfile(edf_files(ff).folder,edf_files(ff).name),ifg_fn,preproc_fn,ii_params,rsvp_params,[],{});
    
    if ff == 1
        % plot some features of the data
        % (check out the docs for each of these; lots of options...)
        rsvp_prf_plottimeseries(ii_data,ii_cfg, {'X','Y','TarX','TarY'})%,'PupilZ'}); % pltos the full timeseries
        
        ii_plotalltrials(ii_data,ii_cfg, {'X','Y'},{},{},{},{},ii_params.show_plots); % plots each trial individually
        
        ii_plotalltrials2d(ii_data,ii_cfg, {'X','Y'}); % plots all trials, in 2d, overlaid on one another w/ fixations
    end
    
    % score trials
    % default parameters should work fine - but see docs for other
    % arguments you can/should give when possible
    [ii_trial{ff},ii_cfg] = rsvp_prf_collectTrials(ii_data,ii_cfg,ii_sacc,ii_params.fix_epoch); 
    
end

ii_sess = ii_combineruns(ii_trial);

%% look at the processed data, make sure it looks ok

% plot which TRs are flagged for saccades (rsvp_prf_proc flags a TR if a
% saccade started in it); figures are saved next to the EDFs:
% - one TR x sweep matrix per run (*_TR_matrix.png)
% - horizontal and vertical sweeps from all runs, aligned by bar location,
%   with the number of included/excluded TRs at each location
%   (*_BarHoriz_matrix.png, *_BarVert_matrix.png)
% trial-level exclusion codes (ii_trial{run}.excl_trial: 11 = drift
% correction too large, 13 = fixation away from center) are not plotted
fh_excl = rsvp_prf_plotQC_exclusions(ii_trial,ii_cfg);

