function [ii_data,ii_cfg] = rdk_relabel_xdat(ii_data,ii_cfg,rdk_params,epoch_chan)
%RDK_RELABEL_XDAT Relabel RDK sweep epochs to match the RSVP XDAT convention
%   [ii_data,ii_cfg] = rdk_relabel_xdat(ii_data,ii_cfg,rdk_params,epoch_chan)
%
%   The RDK task labels every sample of a bar sweep with the sweep number
%   (XDAT = 1:n_trials). The RSVP task labels each ITI 1 and each sweep 2,
%   which is what the iEye/RSVP preprocessing expects. The RDK task has no
%   ITI, so instead the first bar position of each sweep (the first
%   rdk_params.step_dur s) is labeled 1 and the rest of the sweep 2.
%
%   A sweep starts wherever epoch_chan changes to a sweep number. Samples
%   whose value is not a sweep number (e.g. 101/111 start/end codes, if not
%   already removed by ii_trim) are left unchanged. Run after ii_trim.
%
%   Checking for lost messages: sometimes xDAT messages are lost from the
%   EDF; e.g. if every even-numbered message is lost, sweep 2's samples keep
%   sweep 1's label. So first, the sweeps found in epoch_chan are compared
%   with the completed sweeps in rdk_params (those with both a sweep_start
%   and a sweep_end). If every completed sweep's onset is present, in order,
%   the labels are used as they are. They come straight from the task's
%   messages, so they are right even when the sample count doesn't match the
%   task clock (e.g. EyeLink sometimes records extra samples 1 ms apart
%   during tracking loss, which iEye counts as 2 ms each).
%
%   Otherwise a warning is given and the sweep labels are reconstructed from
%   the task's timing: onsets that survived in epoch_chan are kept, and
%   missing ones are filled in from sweep_start. Only completed sweeps are
%   reconstructed, and only samples already labeled with some sweep number
%   are changed. Sweep-labeled samples that don't fall in a completed sweep
%   (e.g. an aborted sweep) are set to 0, i.e. not a sweep. A completed
%   sweep whose samples were all trimmed (e.g. its own message and all
%   earlier ones were lost) can't be recovered.
%
%   For reconstruction, the task clock is aligned to samples using the
%   median offset of the surviving onsets. If any surviving onset is more
%   than TOL_MS from the task timing (e.g. because of a gap in the
%   recording, or many extra samples), the labels can't be reconstructed
%   reliably and an error is thrown.
%
% Inputs:
%   rdk_params:  output of rdk_read_params (uses n_trials, step_dur, and
%                sweep_start/sweep_end if present)
%   epoch_chan:  channel name containing epoch label (default 'XDAT')
%
%   Outputs saved in ii_cfg.rdk_relabel:
%   - orig_<epoch_chan>: the original labels
%   - sweep_<epoch_chan>: the sweep number labels that were used, after any
%     reconstruction
%   - reconstructed: true if the labels were reconstructed from timing
%
%   See also rdk_read_params, ii_trim, rsvp_prf_preproc.

TOL_MS = 20; % max allowed difference between an XDAT sweep onset and the task's sweep_start

if nargin < 4
    epoch_chan = 'XDAT';
end

orig_xdat = ii_data.(epoch_chan);
xdat = orig_xdat;

% check sweep labels against task timing, reconstruct if they disagree
reconstructed = false;
if all(isfield(rdk_params, {'sweep_start','sweep_end'}))
    [xdat, reconstructed] = check_sweep_timing(xdat, rdk_params, ii_cfg.hz, TOL_MS, epoch_chan);
else
    warning('rdk_relabel_xdat:noSweepTiming','rdk_params has no sweep_start/sweep_end; %s not checked against task timing', epoch_chan);
end

% samples in a sweep, and the first sample of each sweep
is_sweep = ismember(xdat, 1:rdk_params.n_trials);
is_onset = is_sweep & [true; diff(xdat) ~= 0];
onset_idx = find(is_onset);
if isempty(onset_idx)
    warning('rdk_relabel_xdat:noSweepsFound','No sweep epochs (1-%i) found in %s; leaving it unchanged', ...
        rdk_params.n_trials, epoch_chan);
    return;
end

% position of each sample within its sweep (1 = sweep onset)
sweep_num = max(cumsum(is_onset), 1); % samples before first sweep aren't used
pos_in_sweep = (1:numel(xdat)).' - onset_idx(sweep_num) + 1;

n_first = round(rdk_params.step_dur * ii_cfg.hz); % samples in first bar position

new_xdat = xdat;
new_xdat(is_sweep & pos_in_sweep <= n_first) = 1;
new_xdat(is_sweep & pos_in_sweep >  n_first) = 2;
ii_data.(epoch_chan) = new_xdat;

% add history to ii_cfg
if reconstructed
    recon_str = ' (sweep labels reconstructed from task timing)';
else
    recon_str = '';
end
ii_cfg.history{end+1} = sprintf('rdk_relabel_xdat - relabeled %i sweeps in channel %s%s: first %i samples (%.2f s) = 1, rest = 2 on %s', ...
    numel(onset_idx), epoch_chan, recon_str, n_first, rdk_params.step_dur, datestr(now,30));
ii_cfg.rdk_relabel.(sprintf('orig_%s',epoch_chan)) = orig_xdat;
ii_cfg.rdk_relabel.(sprintf('sweep_%s',epoch_chan)) = xdat;
ii_cfg.rdk_relabel.reconstructed = reconstructed;

end


function [xdat, reconstructed] = check_sweep_timing(xdat, rdk_params, hz, tol_ms, epoch_chan)
% Compare sweep onsets in xdat with the task's sweep_start/sweep_end (s, task
% clock); if they don't match, reconstruct the sweep labels from timing.

reconstructed = false;
n_sweeps  = rdk_params.n_trials;
start_s   = rdk_params.sweep_start(:);
end_s     = rdk_params.sweep_end(:);
completed = find(~isnan(start_s) & ~isnan(end_s)); % only these can be recovered

% sweep onsets found in xdat
is_sweep = ismember(xdat, 1:n_sweeps);
on_idx = find(is_sweep & [true; diff(xdat) ~= 0]);
on_val = xdat(on_idx);

if isequal(on_val, completed)
    return; % no messages lost: every completed sweep is labeled, in order
end

% align task clock to samples using the onsets of completed sweeps
use = ismember(on_val, completed);
if ~any(use)
    warning('rdk_relabel_xdat:cannotAlign','No onsets of completed sweeps found in %s; cannot check it against task timing', epoch_chan);
    return;
end
offset = median(on_idx(use) - start_s(on_val(use)) * hz);
pred_start = round(offset + start_s * hz); % predicted onset sample of each sweep
pred_end   = round(offset + end_s * hz);

resid_ms = (on_idx(use) - pred_start(on_val(use))) / hz * 1000;
if any(abs(resid_ms) > tol_ms)
    error('rdk_relabel_xdat:timingMismatch', ...
        ['Sweep onsets in %s disagree with the task''s sweep_start by up to %.0f ms (tolerance %i ms); ' ...
         'the recording may have gaps, so sweep labels can''t be reconstructed reliably'], ...
        epoch_chan, max(abs(resid_ms)), tol_ms);
end

% onsets that survived in xdat are exact, so only predict the missing ones
pred_start(on_val(use)) = on_idx(use);

% reconstruct: label k from sweep k's onset until the next completed sweep's
% onset (as the eyetracker would have), or until sweep_end if the next sweep
% wasn't completed; the last sweep's label runs to the end of the data
% (normally trimmed at the xDAT 111 end code)
n_samp = numel(xdat);
new_xdat = xdat;
new_xdat(is_sweep) = 0; % sweep-labeled samples not in a completed sweep
for k = completed.'
    if ismember(k+1, completed)
        last = pred_start(k+1) - 1;
    elseif k == n_sweeps
        last = n_samp;
    else
        last = pred_end(k);
    end
    next_onset = on_idx(find(on_idx > pred_start(k), 1)); % never run past a recorded onset
    if ~isempty(next_onset)
        last = min(last, next_onset - 1);
    end
    idx = max(pred_start(k), 1):min(last, n_samp);
    idx = idx(is_sweep(idx)); % only change samples that were labeled as a sweep
    new_xdat(idx) = k;
end

recovered = intersect(completed, new_xdat);
list_str = @(v) ['[' strtrim(sprintf('%i ', v)) ']'];
warning('rdk_relabel_xdat:sweepLabelsReconstructed', ...
    ['Sweep labels in %s don''t match the task''s sweep timing.\n' ...
     '  Completed sweeps: %s\n  Sweep onsets found in %s: %s\n' ...
     'Reconstructed sweep labels from sweep_start/sweep_end (surviving onsets within %.1f ms of prediction).\n' ...
     '  Sweeps recovered: %s\n  Completed sweeps with no samples to label (not recovered): %s\n' ...
     '  Sweep-labeled samples outside any completed sweep, set to 0: %i'], ...
    epoch_chan, list_str(completed), epoch_chan, list_str(on_val), max(abs(resid_ms)), ...
    list_str(recovered), list_str(setdiff(completed, recovered)), sum(is_sweep & new_xdat == 0));

xdat = new_xdat;
reconstructed = true;

end
