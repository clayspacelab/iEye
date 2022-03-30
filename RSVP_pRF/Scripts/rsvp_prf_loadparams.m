function [rsvp_params] = rsvp_prf_loadparams()


% Load parameters from RSVP pRF task
params_fid = fopen('rsvp_params.txt', 'r');

% Older version of rsvp_params.txt below
%{ 
rsvp_params_list = {'disp_units'; 'targ_rate'; 'n_trials'; 'n_bars'; 'tr'; 'tr_per_bar'; 'testing'; 'screen_res'; 
    'stim_bounds'; 'screen_height'; 'view_dist'; 'background_color'; 'targ_cooldown'; 'response_period'; 
    'constrained_set'; 'calibrate_targ'; 'bore_mask'; 'save_log'; 'stair_upper'; 'stair_lower'; 'eyetracker_ip'; 
    'pulse_cue'; 'extended_start'; 'fix_size'; 'calib_type';};
%}
rsvp_params_list = {'eyetracker_ip';'calib_type';'response_key';'pulse_cue';'background_color';'text_color';'fix_color';'tr';
    'tr_per_bar';'n_bars';'n_trials';'targ_rate';'targ_cooldown';'response_period';'response_delay';'show_feedback';
    'feedback_frames';'calibrate_targ';'extended_start';'bore_mask';'constrained_set';'stair_upper';'stair_lower';
    'fix_size';'save_log';'screen_res';'stim_bounds';'screen_height';'view_dist';'disp_units';'testing';};



rsvp_params = {};
i = 1;
for j = 1:93
    thisline = fgetl(params_fid);
    if mod(j,3) == 2
        rsvp_params{i,1} = thisline;
        i = i+1;
    end
end
rsvp_params = cell2struct(rsvp_params, rsvp_params_list);
rsvp_params.n_trials = str2num(rsvp_params.n_trials);
rsvp_params.n_bars = str2num(rsvp_params.n_bars);
rsvp_params.tr = str2double(rsvp_params.tr);
rsvp_params.tr_per_bar = str2num(rsvp_params.tr_per_bar);


% Get screen params, convert to DVA
mon_height = str2double(rsvp_params.screen_height);
mon_dist = str2double(rsvp_params.view_dist);
mon_vert_res = split(rsvp_params.screen_res, ',');
mon_vert_res = str2num(mon_vert_res{2});
mon_vert_deg = (180 / pi) * 2 * atan(mon_height / (2 * mon_dist));
ppd = mon_vert_res / mon_vert_deg;
rsvp_params.ppd = ppd;


end