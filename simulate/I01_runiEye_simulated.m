function [ii_sess_pro, ii_sess_anti] = I01_runiEye_simulated(direct, taskMap, end_block)

if nargin < 3
    end_block = 10;
end

ii_params = ii_loadparams; % load default set of analysis parameters, only change what we have to
ii_params.valid_epochs = [1 2 3 4 5 6 7 8];
ii_params.trial_start_value = 1; %XDAT value for trial start
ii_params.trial_end_value = 8;   % XDAT value for trial end
ii_params.drift_epoch = [1 2 3 4]; % XDAT values for drift correction
ii_params.drift_fixation_mode  = 'mode';
ii_params.calibrate_epoch = 7; % XDAT value for when we calibrate (feedback stim)
ii_params.calibrate_select_mode = 'last'; % how do we select fixation with which to calibrate?
ii_params.calibrate_mode = 'run'; % scale: trial-by-trial, rescale each trial; 'run' - run-wise polynomial fit
ii_params.blink_thresh = 0.1;
ii_params.blink_window = [200 200]; % how long before/after blink (ms) to drop?
ii_params.plot_epoch = [4 5 7];  % what epochs do we plot for preprocessing?
ii_params.calibrate_limits = [2.5]; % when amount of adj exceeds this, don't actually calibrate (trial-wise); ignore trial for polynomial fitting (run)
ii_params.resolution = [960*2 540*2];
ii_params.ppd = 55 * tand(1) / 0.0264;
excl_criteria.i_dur_thresh = 850; % must be shorter than 150 ms
excl_criteria.i_amp_thresh = 2;   % must be longer than 5 dva [if FIRST saccade in denoted epoch is not at least this long and at most this duration, drop the trial]
excl_criteria.i_err_thresh = 10;   % i_sacc must be within this many DVA of target pos to consider the trial

excl_criteria.drift_thresh = 2.5;     % if drift correction norm is this big or more, drop
excl_criteria.delay_fix_thresh = 2.5; % if any fixation is this far from 0,0 during delay (epoch 3)

block_pro = 1;
block_anti = 1;
for block = 1:end_block
    disp(['Running block ' num2str(block, "%02d")])
    direct.block = [direct.day '/block' num2str(block,"%02d")];
    
    % Loading task, display and timeReport for the block
%     matFile_extract = dir(fullfile(direct.block, '*.mat'));
%     matFile = [direct.block filesep matFile_extract(end).name];
%     load(matFile);
%     parameters = matFile.parameters;
%     screen = matFile.screen;
%     timeReport = matFile.timeReport;
%     trialInfo{block} = timeReport;
    eyecond = taskMap(block).condition;
    
    edfFileName = parameters.edfFile;
    edfFile_original = [direct.block filesep edfFileName '.edf'];
    edf_block_fold = [direct.save_eyedata '/block' num2str(block, "%02d")];
    if ~exist(edf_block_fold, 'dir')
        mkdir(edf_block_fold);
    end

    edfFile = [edf_block_fold filesep edfFileName '.edf'];
    copyfile(edfFile_original, edfFile);
    % what is the output filename?
    preproc_fn = edfFile(1:end-4);
    % simulate data
    [ii_data, ii_cfg, real_error] = simulate_eyetracking(subjID, day, block);
    % run preprocessing!
    [ii_data, ii_cfg, ii_sacc, real_error] = I02_iipreproc_simulated(preproc_fn, ii_params);
        
    % score trials
    % default parameters should work fine - but see docs for other
    % arguments you can/should give when possible
    if ii_sacc.epoch_start == 5
        ii_sacc.epoch_start = 6;
    end
%     if block == 5
%         taskMap(block).stimVF = taskMap(block).stimVF(2:end);
%     end
    if strcmp(eyecond, 'pro')
        [ii_trial_pro{block_pro},~] = ii_scoreMGS(ii_data,ii_cfg,ii_sacc,[],6,[],excl_criteria,[],'lenient');
        ii_trial_pro{block_pro}.stimVF = taskMap(block).stimVF;
        block_pro = block_pro+1;
    elseif strcmp(eyecond, 'anti')
        [ii_trial_anti{block_anti},~] = ii_scoreMGS(ii_data,ii_cfg,ii_sacc, [], 6,[],excl_criteria,[],'lenient');
        ii_trial_anti{block_anti}.stimVF = taskMap(block).stimVF;
        block_anti = block_anti+1;
    end
end
%% Combining runs
disp('Combining runs')
disp(['Total trials = ', num2str(end_block * 40)])
% For pro trials
if ~exist("ii_trial_pro", "var")
    ii_sess_pro = [];
else
    ii_sess_pro = ii_combineruns(ii_trial_pro);
    disp(['Pro trials = ', num2str(size(ii_sess_pro.i_sacc_err, 1))])
    disp(['nan trials ii_sess_pro.i_sacc_err = ', num2str(sum(isnan(ii_sess_pro.i_sacc_err)))])
    disp(['nan trials ii_sess_pro.f_sacc_err = ', num2str(sum(isnan(ii_sess_pro.f_sacc_err)))])
    ii_sess_pro.break_fix = zeros(length(ii_sess_pro.excl_trial), 1);
    ii_sess_pro.prim_sacc = ones(length(ii_sess_pro.excl_trial), 1);
    ii_sess_pro.small_sacc = zeros(length(ii_sess_pro.excl_trial), 1);
    ii_sess_pro.large_error = zeros(length(ii_sess_pro.excl_trial), 1);

    for ii = 1:length(ii_sess_pro.excl_trial)
        % Flag trials with fixation breaks
        if sum(ii_sess_pro.excl_trial{ii} == 13) > 0
            ii_sess_pro.break_fix(ii) = 1;
        end
        % Flag trials with no primary saccades
        if sum(ii_sess_pro.excl_trial{ii} == 20) > 0
            ii_sess_pro.prim_sacc(ii) = 0;
        end
        % Flag trials with small or short saccades
        if sum(ii_sess_pro.excl_trial{ii} == 21) > 0
            ii_sess_pro.small_sacc(ii) = 1;
        end
        % Flag trials with large MGS errors
        if sum(ii_sess_pro.excl_trial{ii} == 22) > 0
            ii_sess_pro.large_error(ii) =  1;
        end
    end
end
% For anti trials
if ~exist("ii_trial_anti", "var")
    ii_sess_anti = [];
else
    ii_sess_anti = ii_combineruns(ii_trial_anti);
    disp(['Anti trials = ', num2str(size(ii_sess_anti.i_sacc_err, 1))])
    disp(['nan trials ii_sess_anti.i_sacc_err = ', num2str(sum(isnan(ii_sess_anti.i_sacc_err)))])
    disp(['nan trials ii_sess_anti.f_sacc_err = ', num2str(sum(isnan(ii_sess_anti.f_sacc_err)))])
    ii_sess_anti.break_fix = zeros(length(ii_sess_anti.excl_trial), 1);
    ii_sess_anti.prim_sacc = ones(length(ii_sess_anti.excl_trial), 1);
    ii_sess_anti.small_sacc = zeros(length(ii_sess_anti.excl_trial), 1);
    ii_sess_anti.large_error = zeros(length(ii_sess_anti.excl_trial), 1);

    for ii = 1:length(ii_sess_anti.excl_trial)
        % Flag trials with fixation breaks
        if sum(ii_sess_anti.excl_trial{ii} == 13) > 0
            ii_sess_anti.break_fix(ii) = 1;
        end
        % Flag trials with no primary saccades
        if sum(ii_sess_anti.excl_trial{ii} == 20) > 0
            ii_sess_anti.prim_sacc(ii) = 0;
        end
        % Flag trials with small or short saccades
        if sum(ii_sess_anti.excl_trial{ii} == 21) > 0
            ii_sess_anti.small_sacc(ii) = 1;
        end
        % Flag trials with large MGS errors
        if sum(ii_sess_anti.excl_trial{ii} == 22) > 0
            ii_sess_anti.large_error(ii) =  1;
        end
    end
end
end
