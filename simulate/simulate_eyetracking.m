% Created by Mrugank (11/13/2022)
% Simulates eye-data in the format obtained from ii_import_edf. As of now
% it is very specific to MD_TMS_EEG/mgs_stimul task. Needs to be
% generalized, but should be easy.
function [ii_data, ii_cfg, real_error] = simulate_eyetracking(subjID, block, taskMap)
%clearvars -except block; close all; clc;
t_array = [0.5, 2, 3/20, 4-2-3/20, 0.15, 0.7, 0.8];
p.itiDuration = [2,3];
p.nTrials = 40;
p.ifg_freq = 1000;

p.trialDuration = sum(t_array) + mean(p.itiDuration);
%sampleCount = int(p.nTrials * p.trialDuration * ifg_freq);
p.xC = 960; p.yC = 540;
%ferr_std = [1, ;

correct_coords = taskMap(block).saccLocpix;

s = 1;
ii_data.XDAT = [];
ii_data.TarX = [];
ii_data.TarY = [];
sacc_arr = zeros(p.nTrials, 2);
feedback_arr = zeros(p.nTrials, 2);

for trial = 1:p.nTrials
    for epoch = 1:size(t_array, 2)+1
        if epoch < 8
            s_end = s+t_array(epoch)*p.ifg_freq+1;
        else
            iti_now = p.itiDuration((rand() >= 0.5) + 1);
            s_end = s+iti_now*p.ifg_freq+1;
        end
        ii_data.XDAT(s:s_end) = epoch;
        if epoch == 5 || epoch == 6 || epoch == 7
            ii_data.TarX(s:s_end) = correct_coords(trial, 1);
            ii_data.TarY(s:s_end) = correct_coords(trial, 2);
        else
            ii_data.TarX(s:s_end) = p.xC;
            ii_data.TarY(s:s_end) = p.yC;
        end
        if epoch == 5
            sacc_arr(trial, 1) = s;
        elseif epoch == 6
            sacc_arr(trial, 2) = s_end;
        elseif epoch == 7
            feedback_arr(trial, 1) = s;
            feedback_arr(trial, 2) = s_end;
        end
        s = s_end;
    end
end

[X, Y, pup, real_error] = generate_XYPupil(correct_coords,s,sacc_arr,feedback_arr,p);
%[X, Y, pup] = generate_blinks(X, Y, pup, s);
ii_data.Pupil = pup';
ii_data.X = X';
ii_data.Y = Y';
ii_data.TarX = ii_data.TarX';
ii_data.TarY = ii_data.TarY';
ii_data.XDAT = ii_data.XDAT';

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
ii_cfg.edf_file = ['/d/DATC/datc/MD_TMS_EEG/analysis/sim/day01/EyeData/block' num2str(block,"%02d") '/' num2str(subjID, "%02d") num2str(block,"%02d") '0000.edf'];
ii_cfg.microsacc = [];
end