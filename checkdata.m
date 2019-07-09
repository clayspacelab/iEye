function  [handle] = checkdata(ii_data,ii_cfg,plotnum)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
handle = figure(plotnum)
hold on; plot(ii_data.X)
plot(ii_data.Y) 

end

