function [ii_data,ii_cfg] = ii_reach_selectempty(ii_data,ii_cfg)
%II_SELECTEMPTY Ensures 'selection' is empty
%   Maintains conventions of ii_data, ii_cfg arguments
%
% Updated TCS 8/14/2017 for iEye refactor

sel = ii_cfg.reach_sel;
sel = sel*0;

cursel = [];

ii_cfg.reach_cursel = cursel;
ii_cfg.reach_sel = sel;


%ii_hideselections;
end

