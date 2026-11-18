function task_timing = rsvp_prf_extract_task_timing(ii_data, ii_cfg, params, epoch_chan)
%RSVP_PRF_EXTRACT_TASK_TIMING Sample indices and bar position of every bar step
%
%   task_timing = rsvp_prf_extract_task_timing(ii_data, ii_cfg, params)
%   task_timing = rsvp_prf_extract_task_timing(ii_data, ii_cfg, params, epoch_chan)
%
%   Works out when each bar position of each sweep was on screen (as sample
%   indices into ii_data) and where the bar was (in degrees), for both the
%   RSVP task (rsvp_sweep.py) and the RDK task (bar_multSize_sameSeq_PTB.m).
%   Used by rsvp_prf_proc (TR bricks, saccade flagging) and
%   rsvp_prf_plottimeseries (trial structure and stimulus trace), and (via
%   rsvp_prf_collectTrials) rsvp_prf_plotQC_exclusions (sweep directions).
%
%   Sweep onsets: each sweep starts at the first sample of a run of
%   epoch_chan == onset epoch, which depends on params.task: 1 for 'rdk'
%   (after rdk_relabel_xdat, XDAT 1 = first bar position of each sweep), 2
%   otherwise (RSVP: XDAT 1 = loading period before each sweep, 2 = sweep).
%   Within a sweep, bar step k starts (k-1) * tr_per_bar * tr seconds after
%   the onset, and TR j of that step starts another (j-1) * tr seconds later. Both tasks hold each bar
%   position for exactly tr_per_bar TRs, so tr_per_bar must be a whole number.
%
%   Warns if a sweep's length in epoch_chan (onset to the next onset or to
%   the end of the sweep epochs) is off from n_bars * tr_per_bar * tr by
%   more than 0.1 TR.
%
%   Sweep numbers: for RDK runs processed with rdk_relabel_xdat, each sweep's
%   number is read from the sweep labels it saved (so sweeps lost from the
%   recording don't shift the numbering); otherwise sweeps are numbered in
%   order of appearance.
%
%   Bar positions (degrees, + = right/up, bar center):
%   - RSVP: directions cycle L2R, T2B, R2L, B2T; the bar is one image wide
%     (stim_bounds(2)/6 px) and its center runs from margin to margin of
%     stim_bounds in n_bars equal steps (as in rsvp_sweep.py).
%   - RDK: direction and width come from params.seq (col 1 = width in deg,
%     col 2 = 1/2/3/4 = down/up/left/right); the center runs over
%     +/-(bar_extent * screen_h_deg - width)/2 in n_bars equal steps (as in
%     bar_multSize_sameSeq_PTB.m).
%
%   Output task_timing (one row per bar step, in time order):
%     sweep, step       sweep number, bar position within the sweep (1..n_bars)
%     start_idx, end_idx  first/last sample of the step
%     tr_start, tr_end  [n_rows x tr_per_bar] first/last sample of each TR
%     tar_x, tar_y      bar center (deg); the coordinate that doesn't move is 0
%     dir               direction label ('L2R','R2L','T2B','B2T')
%     width_deg         bar width (deg)
%   plus per sweep (n_trials x 1, from the task params, whether or not the
%   sweep was found in the data):
%     sweep_dir         direction label of each sweep
%     sweep_width_deg   bar width of each sweep (deg)
%   plus scalars: tr_per_bar, n_bars, n_trials, samples_per_tr,
%   stim_extent_deg ([width height] of the stimulus area), onset_epoch.
%   Sample indices can run past the end of the data if the recording stopped
%   early; callers should clip them.
%
%   See also rsvp_prf_proc, rsvp_prf_plottimeseries, rdk_relabel_xdat.

if nargin < 4 || isempty(epoch_chan)
    epoch_chan = 'XDAT';
end

is_rdk = isfield(params, 'task') && strcmp(params.task, 'rdk');
hz = ii_cfg.hz;
tr = params.tr;
n_bars = params.n_bars;
n_trials = params.n_trials;

% ---- check timing params ----
tr_per_bar = round(params.tr_per_bar);
if abs(params.tr_per_bar - tr_per_bar) > 1e-6 || tr_per_bar < 1
    error('rsvp_prf_extract_task_timing:trPerBar', ...
        ['tr_per_bar is %g; each bar position must last a whole number of TRs ' ...
         '(for RDK, step_dur / tr = %g / %g)'], params.tr_per_bar, params.tr_per_bar * tr, tr);
end

% ---- sweep onsets ----
if is_rdk
    onset_epoch = 1;
else
    onset_epoch = 2;
end

xdat = ii_data.(epoch_chan)(:);
in_onset_epoch = xdat == onset_epoch;
onset_idx = find(in_onset_epoch & [true; ~in_onset_epoch(1:end-1)]);

% sweep numbers
sweep_id = (1:numel(onset_idx)).';
if is_rdk && isfield(ii_cfg, 'rdk_relabel')
    lab_field = sprintf('sweep_%s', epoch_chan);
    if isfield(ii_cfg.rdk_relabel, lab_field) && numel(ii_cfg.rdk_relabel.(lab_field)) == numel(xdat)
        sweep_id = ii_cfg.rdk_relabel.(lab_field)(onset_idx);
        sweep_id = sweep_id(:);
    end
end

ok = sweep_id >= 1 & sweep_id <= n_trials;
if any(~ok)
    warning('rsvp_prf_extract_task_timing:extraSweeps', ...
        'Found %i sweep onsets in %s but the task has %i sweeps; ignoring the extra/unlabeled ones', ...
        numel(onset_idx), epoch_chan, n_trials);
    onset_idx = onset_idx(ok);
    sweep_id = sweep_id(ok);
end
if numel(onset_idx) < n_trials
    error('rsvp_prf_extract_task_timing:missingSweeps', ...
        'Found %i of %i sweeps in %s (value %i); every sweep must be in the recording', ...
        numel(onset_idx), n_trials, epoch_chan, onset_epoch);
end

% ---- check sweep lengths ----
% each sweep lasts until the next sweep onset or until epoch_chan leaves the
% sweep epochs (RDK: 1-2, back to back; RSVP: 2, then the loading period)
if is_rdk
    sweep_epochs = [1 2];
else
    sweep_epochs = 2;
end
in_sweep = ismember(xdat, sweep_epochs);
sweep_end = zeros(size(onset_idx));
for ss = 1:numel(onset_idx)
    last = onset_idx(ss) - 1 + find(~in_sweep(onset_idx(ss):end), 1) - 1;
    if isempty(last), last = numel(xdat); end
    next_on = onset_idx(find(onset_idx > onset_idx(ss), 1));
    if ~isempty(next_on), last = min(last, next_on - 1); end
    sweep_end(ss) = last;
end
sweep_dur = (sweep_end - onset_idx + 1) / hz;                % s
nominal_dur = n_bars * tr_per_bar * tr;                      % s
off = abs(sweep_dur - nominal_dur) > 0.1 * tr;
if any(off)
    warning('rsvp_prf_extract_task_timing:sweepLength', ...
        ['Sweep(s) %s last %s s in %s, not %.2f s. The TR grid is timed from sweep ' ...
         'onset, so the last TR of a short sweep includes samples after the sweep ' ...
         '(saccades there are flagged), and the end of a long sweep has no TR ' ...
         '(saccades there are not flagged).'], ...
        mat2str(sweep_id(off).'), mat2str(round(sweep_dur(off).', 2)), epoch_chan, nominal_dur);
end

% ---- bar positions for each sweep (deg) ----
[pos, axis_moves, dir_lbl, width_deg, stim_extent_deg] = bar_positions(params, is_rdk);

% ---- assemble one row per bar step ----
n_sw = numel(onset_idx);
n_rows = n_sw * n_bars;
task_timing.sweep = reshape(repmat(sweep_id.', n_bars, 1), [], 1);
task_timing.step = repmat((1:n_bars).', n_sw, 1);

% time of each TR boundary relative to its sweep's onset, in samples
% (rounded from seconds so rounding error doesn't accumulate)
tr_bounds = round((0:(n_bars * tr_per_bar)) * tr * hz);      % 1 x (n_TR+1)
tr_on = tr_bounds(1:end-1);                                   % TR start offsets
tr_off = tr_bounds(2:end) - 1;                                % TR end offsets
on_rows = reshape(repmat(onset_idx.', n_bars, 1), [], 1);     % onset of each row's sweep

task_timing.tr_start = zeros(n_rows, tr_per_bar);
task_timing.tr_end = zeros(n_rows, tr_per_bar);
for jj = 1:tr_per_bar
    tr_k = (task_timing.step - 1) * tr_per_bar + jj;                % TR number within sweep
    task_timing.tr_start(:, jj) = on_rows + tr_on(tr_k).';
    task_timing.tr_end(:, jj) = on_rows + tr_off(tr_k).';
end
task_timing.start_idx = task_timing.tr_start(:, 1);
task_timing.end_idx = task_timing.tr_end(:, end);

task_timing.tar_x = zeros(n_rows, 1);
task_timing.tar_y = zeros(n_rows, 1);
task_timing.dir = cell(n_rows, 1);
task_timing.width_deg = zeros(n_rows, 1);
for rr = 1:n_rows
    ss = task_timing.sweep(rr);
    if axis_moves(ss) == 'x'
        task_timing.tar_x(rr) = pos(task_timing.step(rr), ss);
    else
        task_timing.tar_y(rr) = pos(task_timing.step(rr), ss);
    end
    task_timing.dir{rr} = dir_lbl{ss};
    task_timing.width_deg(rr) = width_deg(ss);
end

task_timing.sweep_dir = dir_lbl(:);
task_timing.sweep_width_deg = width_deg(:);
task_timing.tr_per_bar = tr_per_bar;
task_timing.n_bars = n_bars;
task_timing.n_trials = n_trials;
task_timing.samples_per_tr = tr * hz;
task_timing.stim_extent_deg = stim_extent_deg;
task_timing.onset_epoch = onset_epoch;

end


function [pos, axis_moves, dir_lbl, width_deg, stim_extent_deg] = bar_positions(params, is_rdk)
% pos: n_bars x n_trials bar-center positions (deg) along the axis that moves
n_bars = params.n_bars;
n_trials = params.n_trials;
pos = zeros(n_bars, n_trials);
axis_moves = repmat('x', 1, n_trials);
dir_lbl = cell(1, n_trials);
width_deg = zeros(1, n_trials);

if is_rdk
    % bar_multSize_sameSeq_PTB.m: whole sweep spans bar_extent * screen height
    extent = params.bar_extent * params.screen_h_deg;
    stim_extent_deg = [extent extent];
    rdk_dirs = {'T2B', 'B2T', 'R2L', 'L2R'}; % seq(:,2) = 1/2/3/4 = down/up/left/right
    for ss = 1:n_trials
        w = params.seq(ss, 1);
        d = params.seq(ss, 2);
        p = linspace(-extent/2 + w/2, extent/2 - w/2, n_bars);
        if d == 1 || d == 3
            p = fliplr(p);
        end
        pos(:, ss) = p.';
        if d == 1 || d == 2
            axis_moves(ss) = 'y';
        end
        dir_lbl{ss} = rdk_dirs{d};
        width_deg(ss) = w;
    end
else
    % rsvp_sweep.py: 6 images (image_h px each) per bar; bar center runs
    % between margins half an image inside stim_bounds
    ppd = params.ppd;
    sb = params.stim_bounds;
    image_h = sb(2) / 6;
    x_marg = sb(1)/2 - image_h/2;
    y_marg = sb(2)/2 - image_h/2;
    stim_extent_deg = sb / ppd;
    rsvp_dirs = {'L2R', 'T2B', 'R2L', 'B2T'};
    for ss = 1:n_trials
        d = rsvp_dirs{mod(ss - 1, 4) + 1};
        switch d
            case 'L2R', p = linspace(-x_marg, x_marg, n_bars);
            case 'R2L', p = linspace(x_marg, -x_marg, n_bars);
            case 'T2B', p = linspace(y_marg, -y_marg, n_bars); axis_moves(ss) = 'y';
            case 'B2T', p = linspace(-y_marg, y_marg, n_bars); axis_moves(ss) = 'y';
        end
        pos(:, ss) = p.' / ppd;
        dir_lbl{ss} = d;
        width_deg(ss) = image_h / ppd;
    end
end
end
