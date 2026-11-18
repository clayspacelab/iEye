function [ii_cfg, ii_sacc] = rsvp_prf_proc(ii_data, ii_cfg, ii_sacc, rsvp_params)
% Load parameters file from rsvp_prf program and process data specific to
% rsvp_prf
%
% Works for both the RSVP task and the RDK task (rsvp_params from
% rdk_read_params). The timing and position of every bar step comes from
% rsvp_prf_extract_task_timing, which is saved in ii_cfg.task_timing.
%
% Outputs added:
%   ii_cfg.tr_bricks.tr_<j>: [start end] sample of TR j of each bar step
%       (one row per bar step, in time order)
%   ii_cfg.tr_bricks.sweep_n, .bar_loc: sweep number and bar step of each row
%   ii_cfg.tr_excl: (n_bars * tr_per_bar) x n_trials, 1 = TR flagged
%       (currently: a saccade started in it)
%   ii_cfg.bad_tr_idc: [start end] samples of each flagged TR
%   ii_sacc.sacc_tr, .sacc_bar_loc: TR (within bar step) and bar step in
%       which each saccade started; NaN if it started outside every bar step
%       (e.g. in the RSVP loading period between sweeps)

task_timing = rsvp_prf_extract_task_timing(ii_data, ii_cfg, rsvp_params);
tr_per_bar = task_timing.tr_per_bar;

% TR start/end indices for each bar step, stored in tr_bricks
for jj = 1:tr_per_bar
    tr_bricks.(sprintf('tr_%i', jj)) = [task_timing.tr_start(:, jj) task_timing.tr_end(:, jj)];
end
tr_bricks.sweep_n = task_timing.sweep;
tr_bricks.bar_loc = task_timing.step;


% Flag TRs in which a saccade starts
n_sacc = size(ii_sacc.idx, 1);
tr_excl = zeros(rsvp_params.n_bars * tr_per_bar, rsvp_params.n_trials);
sacc_tr = nan(n_sacc, 1);
sacc_bar_loc = nan(n_sacc, 1);
bad_tr_idc = zeros(0, 2);
for ii = 1:n_sacc
    sacc_start = ii_sacc.idx(ii, 1);
    rr = find(task_timing.start_idx <= sacc_start & task_timing.end_idx >= sacc_start, 1);
    if isempty(rr)
        continue; % not during a bar step (e.g. RSVP loading period)
    end
    jj = find(task_timing.tr_start(rr, :) <= sacc_start & task_timing.tr_end(rr, :) >= sacc_start, 1);
    sacc_tr(ii) = jj;
    sacc_bar_loc(ii) = task_timing.step(rr);
    bad_tr_idc(end+1, :) = [task_timing.tr_start(rr, jj) task_timing.tr_end(rr, jj)]; %#ok<AGROW>
    tr_excl(tr_per_bar * (task_timing.step(rr) - 1) + jj, task_timing.sweep(rr)) = 1;
end
[~, uniq_idc] = unique(bad_tr_idc(:, 1));
bad_tr_idc = bad_tr_idc(uniq_idc, :);


% Load data into stream
ii_sacc.sacc_tr = sacc_tr;
ii_sacc.sacc_bar_loc = sacc_bar_loc;
ii_cfg.tr_excl = tr_excl;
ii_cfg.bad_tr_idc = bad_tr_idc;
ii_cfg.tr_bricks = tr_bricks;
ii_cfg.task_timing = task_timing;
ii_cfg.rsvp_prf_params = rsvp_params;

end
