function [ii_data] = ii_checkfor_reaches (ii_data)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% need a decision point to divert saccade / reach analyses. can find this
% by looking through the ii_data.reach_velocity and determine whether any
% reaches are present. 

if max(ii_data.reach_velocity) <= 20 
    ii_data.blocktype = 'sacc';
    disp('this is likely not a reach block, skipping reach analyses')
else
    ('we have reaches!')
    ii_data.blocktype = 'reach';
end 
ii_data = ii_data; 
end

