function [ii_cfg, ii_sacc] = rsvp_prf_proc(ii_data, ii_cfg, ii_sacc)
% Load parameters file from rsvp_prf program and process data specific to
% rsvp_prf


rsvp_params = rsvp_prf_loadparams();

% Get sample info from params
tr_per_sweep = rsvp_params.tr_per_bar * rsvp_params.n_bars;
samples_per_tr = rsvp_params.tr * ii_cfg.hz;
samples_per_bar = samples_per_tr * rsvp_params.tr_per_bar;


% Find indices for each epoch start and end
tr_sel = ii_cfg.sel * 0;
iti_swhere = find(diff(ii_data.XDAT) == -1);
task_ewhere = [iti_swhere' length(tr_sel)]';
iti_swhere = [1 iti_swhere']';
task_swhere = find(diff(ii_data.XDAT) == 1);
iti_ewhere = task_swhere - 1;
iti_idc = [iti_swhere iti_ewhere];
tr_idc = repmat(1:rsvp_params.tr_per_bar,1,rsvp_params.n_bars);
for i = 1:tr_per_sweep
    if i == 1
        task_idc = [task_swhere (task_swhere + samples_per_tr)];
    else
        task_idc = [task_idc (task_idc(:,size(task_idc,2)) + samples_per_tr)];
    end
end


% Assign tr number to selection vector
tr_names = {};
for i = 1:rsvp_params.tr_per_bar
    tr_name = strcat('tr_',num2str(i));
    tr_names{end+1} = tr_name;
    tr_bricks.(tr_names{end}) = [];
    these_trs = find(tr_idc == i);
    for j = 1:length(these_trs)
        this_tr = these_trs(j);
        tr_sel(task_idc(:,this_tr):task_idc(:,this_tr + 1)) = i;
        % Track indices of tr spans for plotting
        tr_bricks.(tr_names{end}) = [tr_bricks.(tr_names{end}); task_idc(:,this_tr) task_idc(:,this_tr + 1)];
    end
    tr_bricks.(tr_names{end}) = sort(tr_bricks.(tr_names{end}));
end




directions = repmat(['L2R'; 'T2B'; 'R2L'; 'B2T';], ceil(rsvp_params.n_trials / 4), 1);
sweep_n = [];
sweep_dir = [];
bar_loc = [];
for i = 1:rsvp_params.n_trials
    sweep_n = [sweep_n; repmat(i, rsvp_params.n_bars, 1);];
    sweep_dir = [sweep_dir; repmat(directions(i,:), rsvp_params.n_bars, 1);];
    bar_loc = [bar_loc 1:rsvp_params.n_bars];
end
bar_loc = bar_loc';



tr_excl = zeros(rsvp_params.n_bars * rsvp_params.tr_per_bar, rsvp_params.n_trials);


sacc_tr = [];
sacc_bar_loc = [];
sacc_sweep_dir = [];
bad_tr_idc = [];
break_now = false;
% For each saccade,
for i = 1:size(ii_cfg.saccades, 1)
    % Find indices of trs in same trial as saccade
    these_idc = find(sweep_n == ii_sacc.trial_start(i));
    % Find index where saccade starts
    sacc_start = ii_sacc.idx(i,1);
    % Mark sweep direction in which each saccade occurs
    sacc_sweep_dir = [sacc_sweep_dir; directions(ii_sacc.trial_start(i),:)];
    % Ignore saccade if made in ITI
    if ii_sacc.epoch_start(i) == 1
        sacc_tr = [sacc_tr; NaN;];
        sacc_bar_loc = [sacc_bar_loc; NaN;];
        continue;
    end
    % For each group of trs, 
    for j = 1:rsvp_params.tr_per_bar
        % Grab indices of trs from tr_bricks
        this_brick = ['tr_' num2str(j)];
        these_trs = tr_bricks.(this_brick)(these_idc,:);
        % For each tr in group,
        for k = 1:size(these_trs,1)
            if sacc_start >= these_trs(k,1) && sacc_start <= these_trs(k,2)
                sacc_tr = [sacc_tr; j;];
                sacc_bar_loc = [sacc_bar_loc; k;];
                bad_tr_idc = [bad_tr_idc; these_trs(k,:);];
                % Add exclusion info to tr exclusion vec
                excl_idx = rsvp_params.tr_per_bar * (k - 1) + j;
                tr_excl(excl_idx, ii_sacc.trial_start(i)) = 1;
                break_now = true;
                break;
            end
        end
        % Break loop once tr is found
        if break_now
            break_now = false;
            break;
        end
    end
end
[tmp1, uniq_idc, tmp2] = unique(bad_tr_idc(:,1));
bad_tr_idc = bad_tr_idc(uniq_idc,:);
clear tmp1 tmp2;



% Load data into stream
tr_bricks.sweep_n = sweep_n;
tr_bricks.bar_loc = bar_loc;
ii_sacc.sacc_tr = sacc_tr;
ii_sacc.sacc_bar_loc = sacc_bar_loc;
ii_cfg.tr_excl = tr_excl;
ii_cfg.bad_tr_idc = bad_tr_idc;
ii_cfg.tr_sel = tr_sel;
ii_cfg.tr_bricks = tr_bricks;
ii_cfg.rsvp_prf_params = rsvp_params;



end

