%data_compare

%TODO: load subject data files dynamically, loop over all subjects 
load('/Users/grace/GITHUB/iEye/examples/s5_1_4_iEye.mat') %eye data
load('/Users/grace/GITHUB/iEye/examples/5_1_4.mat') %upsampled mocap data

%the eye data and mocap data are not the same length, and do not start at
%the same time. therefore, we first want to align them by their respective
%end indices. to do this, first take the difference in length between the
%files and indx the mocap file starting at that index:end to trim it. 
get_eyedata_size = length(ii_data.XDAT);
get_mocap_size =length(rX);

size_discrep = get_mocap_size-get_eye_size; %take the diff btwn our two files
trim_rX = rX((size_discrep+1):end); 
trim_rY = rY((size_discrep+1):end);

%% collect the correct indices from rX by account for time onset lag
%what is the val of the first EVENT FIXATION from the .tsv file? (w/HEADER)

buff = start_time; 
rXtrim_new = rX(buff:(length(trim_rX)+buff-1)); 
rYtrim_new = rX(buff:(length(trim_rY)+buff-1));
%% plot // visual confirm 

figure (3); plot(ii_data.XDAT-10,'linewidth',2) %subtract 10 to zero align our XDATs, which start at 11
hold on; plot(ii_data.TarX,'linewidth',1.5,'linestyle','--') %target X position
hold on; plot(ii_data.TarY,'linewidth',1.5,'linestyle','--') %target Y position
plot(rXtrim_new, 'linewidth',2,'linestyle',':') %plot the correct version of rX
plot(rYtrim_new, 'linewidth',2,'linestyle',':') %plot the correct version of rY
legend({'XDAT','TarX','TarY','rX','rY'})



%% add adjusted rX/rYtrim_corr to ii_data as new fields. 

ii_data.rX = rXtrim_new;
ii_data.rY = rYtrim_new;

