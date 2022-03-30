% example_analysis.m
%
% example analysis script for example fMRI dataset (exfmri*.edf)
%
% Here's an example data analysis script that walks from raw EDFs for a
% single-item WM task conducted in the scanner (12 s delay interval)
% through preprocessing, saccade extraction, MGS scoring, QC plotting, and
% some basic analyses/plots
%
% Tommy Sprague, 6/12/2018 - subj CC from MGSMap (sess 1, runs 1-3)


% Configure paths & file locations
tmp = mfilename('fullpath'); tmp2 = strfind(tmp,filesep);
addpath(genpath(tmp(1:(tmp2(end-1)))));
root = tmp(1:(tmp2(end)-1));
root = strcat(root,'/Subjects/');


% Select subject to run
%edf_prefix = 'example_sub_1/';
edf_prefix = 'example_sub_2/';

edf_files = dir(fullfile(root,sprintf('%s*.edf',edf_prefix)));


% very first thing we want to do is define all parameters for processing
% (see ii_loadparams.m for default values)
% ifg_fn = 'home/jeff/docs/prf_ieye/iEye/prf/p_500hz.ifg';
config_root = tmp(1:(tmp2(end-1)-1));
ifg_fn = [config_root '/ifg/p_500hz.ifg'];
%ifg_fn = [config_root '/ifg/p_1000hz.ifg'];
ii_params = ii_loadparams; % load default set of analysis parameters, only change what we have to
ii_params.resolution = [1080 1080];
ii_params.valid_epochs =[1 2];
ii_params.trial_end_value = 2;   % XDAT value for trial end
ii_params.drift_epoch = [1]; % XDAT values for drift correction
ii_params.calibrate_epoch = 1;   % XDAT value for when we calibrate (feedback stim)
ii_params.calibrate_select_mode = 'last'; % how do we select fixation with which to calibrate?
ii_params.calibrate_mode = 'scale'; % scale: trial-by-trial, rescale each trial; 'run' - run-wise polynomial fit
ii_params.blink_window = [200 200]; % how long before/after blink (ms) to drop?
ii_params.plot_epoch = [1 2];  % what epochs do we plot for preprocessing?
ii_params.calibrate_limits = [2.5]; % when amount of adj exceeds this, don't actually calibrate (trial-wise); ignore trial for polynomial fitting (run)
ii_params.ppd = 31.8578; % for scanner, 1280 x 1024 - convert pix to DVA
%ii_params.sacc_duration_thresh = 0.01;
ii_params.sacc_amplitude_thresh = 1.25;


% create empty cell array of all our trial data for combination later on
ii_trial = cell(length(edf_files),1);

for ff = 1:length(edf_files)
    
    % what is the output filename?
    preproc_fn = sprintf('%s/%s_preproc.mat',edf_files(ff).folder,edf_files(ff).name(1:(end-4)));
    
    
    % run preprocessing!
     [ii_data, ii_cfg, ii_sacc] = rsvp_prf_preproc(fullfile(edf_files(ff).folder,edf_files(ff).name),ifg_fn,preproc_fn,ii_params,[],{});%{'calibration'});
    
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
    [ii_trial{ff},ii_cfg] = rsvp_prf_collectTrials(ii_data,ii_cfg,ii_sacc); 
    
end

ii_sess = rsvp_prf_combineruns(ii_trial);

%% look at the processed data, make sure it looks ok

% based on what criteria shoudl we exclude trials?
which_excl = [12 31 32];

% first, plot an exclusion report over all runs
% - this will open a few figures. first one will be a 'dot plot', which
%   shows all exclusion criteria exceeded for a given trial by plotting a dot
%   above that criterion. if any dot for a trial is red, we're going to throw that
%   trial away. if black, we'll ignore that criterion for now.
% - second, plots a summary graph, showing % of trials excluded due to each
%   criterion, and how many overall will be excluded. If plotting
%   run-concatenated data, also show how this varies across runs
fh_excl = rsvp_prf_plotQC_exclusions(ii_trial,ii_cfg,which_excl);

