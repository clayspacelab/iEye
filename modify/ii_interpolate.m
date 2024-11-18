function [ii_data, ii_cfg] = ii_interpolate(ii_data,ii_cfg, chan_names, interp_method)
% Created by Mrugank (11/18/2024): Work in progress

% input:
%       ii_data: fields required are X, Y, and Pupil
%       ii_cfg: required filed is blinks which has sample onset and offset
%       for each blink
%       chan_names
%       interp_method: could be one of spline, pchip, linear ... (default = 'spline')

% default values
if nargin < 3
    chan_names = {'X','Y'};
end

if ~iscell(chan_names)
    chan_names = {chan_names};
end

if nargin < 4
    interp_method = 'spline';
end

sampBuffer = 50;

blinks = ii_cfg.blinks;

for cc = 1:length(chan_names)
    if ismember(chan_names{cc},fieldnames(ii_data))
        chanData = ii_data.(chan_names{cc});
        for i = 1:length(blinks)
            blnkOnset = blinks(i, 1);
            blnkOffset = blinks(i, 2);
    
            xQuery = blnkOnset-sampBuffer:blnkOffset+sampBuffer;
            chanData(xQuery) = interp1(1:length(chanData), chanData, xQuery, interp_method);
    
        
        end

        ii_data.(sprintf('%s',chan_names{cc})) = chanData;

    else
    
        error('iEye:ii_smooth:channelNotFound', 'Channel %s does not exist in ii_data',chan_names{cc})
    end
end




end