% Created by Mrugank (11/13/2022)
% Simulates eye-data in the format obtained from ii_import_edf. As of now
% it is very specific to MD_TMS_EEG/mgs_stimul task. Needs to be
% generalized, but should be easy.

clear; close all; clc;
t_array = [0.5, 2, 3/20, 4-2-3/20, 0.15, 0.7, 0.8];
p.itiDuration = [2,3];
p.nTrials = 40;
p.nBlocks = 1;
ifg_freq = 1000;
dt = 1/ifg_freq;
load('/d/DATC/datc/MD_TMS_EEG/data/phosphene_data/sub99/taskMap_sub99_day01_antitype_mirror.mat')
p.trialDuration = sum(t_array) + mean(p.itiDuration);
%sampleCount = int(p.nTrials * p.trialDuration * ifg_freq);
xC = 960; yC = 540;
ierr_mean = 0; ierr_std = 10;
ferr_mean = 0; ferr_std = 5;

for block = p.nBlocks
    correct_coords = taskMap(block).saccLocpix;
    s = 1;
    ii_data.XDAT = [];
    ii_data.TarX = [];
    ii_data.TarY = [];
    sacc_onset_arr = [];
    sacc_end_arr = [];
    feedback_onset_arr = [];
    feedback_end_arr = [];

    isacc_err_va = ierr_mean + ierr_std * randn(p.nTrials, 1);
    fsacc_err_va = ferr_mean + ferr_std * randn(p.nTrials, 1);
    for trial = 1:p.nTrials
        for epoch = 1:size(t_array, 2)+1
            if epoch < 8
                s_end = s+t_array(epoch)*ifg_freq+1;
            else
                iti_now = p.itiDuration((rand() >= 0.5) + 1);
                s_end = s+iti_now*ifg_freq+1;
            end
            ii_data.XDAT(s:s_end) = epoch;
            if epoch == 5 || epoch == 6
                ii_data.TarX(s:s_end) = correct_coords(trial, 1);
                ii_data.TarY(s:s_end) = correct_coords(trial, 2);
            else
                ii_data.TarX(s:s_end) = xC;
                ii_data.TarY(s:s_end) = yC;
            end
            if epoch == 5
                sacc_onset_arr = [sacc_onset_arr, s];
            elseif epoch == 6
                sacc_end_arr = [sacc_end_arr, s_end];
            elseif epoch == 7
                feedback_onset_arr = [feedback_onset_arr, s];
                feedback_end_arr = [feedback_end_arr, s_end];
            end
            s = s_end;
        end
    end
    t = 1:dt:s*dt-dt;
    X = xC;
    Y = yC;
    pup = 1200;
    
    for freq = 0.03:0.03:1
        X = X + 10*exp(-freq) * sin(2*pi*freq*t+10*randn());
        Y = Y + 10*exp(-freq) * sin(2*pi*freq*t+10*randn());
        pup = pup + 50*exp(-freq) * sin(2*pi*freq*t+50*randn());
    end
    
    X = X + 2*randn(1, length(t));
    Y = Y + 2*randn(1, length(t));
    pup = pup + 20*randn(1, length(t));
    
    for trial = 1:p.nTrials
        sacc_start = poissrnd(sacc_onset_arr(trial) + (sacc_end_arr(trial) - sacc_onset_arr(trial))/8);
        sacc_end = poissrnd(feedback_onset_arr(trial)-20 + (feedback_end_arr(trial) - (feedback_onset_arr(trial)-20))/8);
        X(sacc_start:sacc_end) = correct_coords(trial, 1);
        Y(sacc_start:sacc_end) = correct_coords(trial, 2);
        
        
    end
    ii_data.Pupil = pup;
    ii_data.X = X;
    ii_data.Y = Y;
    
    ii_cfg.cursel = [];
    ii_cfg.sel = zeros(s, 1);
    ii_cfg.cfg = 'p_1000hz.ifg';
    ii_cfg.vis = 'X,Y,Pupil,XDAT,TarX,TarY';
    ii_cfg.nchan = 6;
    ii_cfg.lchan = {{'X','Y','Pupil','XDAT','TarX','TarY'}'};
    ii_cfg.hz = 1000;
    ii_cfg.velocity = [];
    ii_cfg.tcursel = [];
    ii_cfg.tsel = zeros(s, 1);
    ii_cfg.tindex = 0;
    ii_cfg.saccades = [];
    ii_cfg.history = {'Simulated data'};
    ii_cfg.edf_file = ['/d/DATC/datc/MD_TMS_EEG/analysis/sim/day01/EyeData/block' num2str(block,"%02d") '/00' num2str(block,"%02d") '0000.edf'];
    ii_cfg.microsacc = [];
end



% pup = 1200;
% for freq = 2:30
%     pup = pup + 500* exp(-freq) * sin(2*pi*freq*t+500*randn());
% end
% pup = pup + 20* randn(1, length(t));
% %pup = 1200 + 500 * sin(2*pi*2/1000*(1:10000))+ 100 * sin(2*pi*8/1000*(1:10000)) + 50 * sin(2*pi*20/1000*(1:10000)) + 250* rand(1, nsamp);
% %pup = 1200 + 500 * cos(2*pi*4*t);
% figure(); plot(pup); ylim([0, max(pup)]);
    