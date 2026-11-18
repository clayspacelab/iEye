function params = rdk_read_params(filename,tr)
%RSVP_READ_PARAMS Read the parameters from a RDK task .mat file
%
%   params = rdk_read_params(filename,tr)
%
%   TR must be specified because it's not given in the task output.
%
%   Because RDK task data is already loaded in directly as a struct from a
%   .mat file, no parsing is needed as with the RSVP data. Rather, this
%   function serves to align field names with the RSVP data for common
%   preprocessing and also strips off the response output so only task params
%   are saved with the eye data preproc


%param key for translation between rsvp and rdk
rdk_params = {'screen_width_deg','screen_height_deg','n_steps','viewing_distance','resolution'};
rsvp_params = {'screen_w_deg','screen_h_deg','n_bars','view_dist','screen_res'};

%load in params
F = load(filename);
params = F.p;
clear F

%add tr and derive tr_per_bar
params.tr = tr;
params.tr_per_bar = params.step_dur / params.tr;

%add explicit sweep variable as in RSVP task
params.n_trials = size(params.seq,1); %this is really the number of bar sweeps, but RSVP task calls it n_trials
n_behave_trials = params.n_trials * params.n_steps;

%derive stim bounds
params.stim_bounds = repmat(params.bar_extent * params.resolution(2),1,2);

%rename key variables to match rsvp for code interoperability
for pp=1:length(rdk_params)
    params = renameStructField(params,rdk_params{pp},rsvp_params{pp});
end

%remove behavioral variables/per behavioral trial variables
param_names = fieldnames(params);
param_lens = structfun(@(x) size(x,1),params);
params = rmfield(params,param_names(param_lens==n_behave_trials));

%label w/ task for preproc handling
params.task = 'rdk';

end