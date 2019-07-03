function [ii_data,ii_cfg, rX,rY] =  upsamplereachfile_gh(no_header_reach_file,header_reach_file,ii_data,ii_cfg)
%IMPORT and UPSAMPLE to 1000K MOCAP TSV FILES
%   This function will import MOCAP files after they have been converted to
%   TSV files. Note that you MUST insert a header into the first row of the
%   TSV file with the following form:  "Frame	Time	rX	rY	rZ".  This
%   function will take the TSV file with added header as an input, and will
%   create a struct with the 5 fields.  Note that the MOCAP data is sampled at 60 Hz, while the EyeLink data is sampled at 1000 Hz. Then, it will upsample the data to
%   6000k, then downsample to 1000k, to match the Eye data. 


data = tdfread(no_header_reach_file,'\t');
% If you get errors, check to make sure you added the header row to the TSV
% file! 

filename = sprintf('%s.mat',no_header_reach_file(1:(end-4)));

% Upsample the 60 Hz X and Y channels to 6000 Hz
data6k.rX = interp((data.rX/10),100);
data6k.rY = interp((data.rY/10),100);

% Downsample the 6000 Hz X and Y channels to 1000 Hz
rX = downsample(data6k.rX,6);
rY = downsample(data6k.rY,6);

%%% TODO geh 062019

%add additional arguments which are the timing offset from header file 
%and length of edf x/y 
%make sure there are no nans in rX / rY
%trim rX/rY by length(ii_data.XDAT) to align the end

%add rX,rY too ii_data struct

%the eye data and mocap data are not the same length, and do not start at
%the same time. therefore, we first want to align them by their respective
%end indices. to do this, first take the difference in length between the
%files and indx the mocap file starting at that index:end to trim it. 
get_eyedata_size = length(ii_data.XDAT);
get_mocap_size =length(rX);

size_discrep = get_mocap_size-get_eyedata_size; %take the diff btwn our two files
trim_rX = rX((size_discrep+1):end); 
trim_rY = rY((size_discrep+1):end);

% how much do we need to adjust?
fid = fopen(header_reach_file);
s = textscan(fid,'%s %s %s %s');
tmp = s{1,4}(10,1); %location of value invariant -- always line 10 position 4 
motime = str2double(tmp)*1000; %location of value invariant

fprintf('time offset from mocap start to first fixation = %f16 milliseconds \n', round(motime))

buff = round(motime); 
rXtrim_aligned = rX(buff:(length(trim_rX)+buff-1)); %when derived originally, called this var 'rXtrim_new'
rYtrim_aligned = rY(buff:(length(trim_rY)+buff-1));

ii_data.rX = rXtrim_aligned;
ii_data.rY = rYtrim_aligned;

%store the history of the reach file 
ii_cfg.reach_file = [filename];

% Save X and Y to mat file
save(filename,'rX','rY');
