function [ii_data,ii_cfg] = ii_selectreachinverse(ii_data,ii_cfg)
%II_SELECTINVERSE Select those points unselected, and de-select the points
%that were selected
% If no selections made, full timeseries will be selected

sel = ii_cfg.reach_sel;


sel = ~(sel==1);

startidx = find(diff([0; sel])== 1);
endidx   = find(diff([sel; 0])==-1);

ii_cfg.reach_sel = sel;
ii_cfg.reach_cursel = [startidx endidx];

end

