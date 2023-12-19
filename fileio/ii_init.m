function ii_init()
% initializes path for ii_import_edf
% TCS for GH 8/25/2016
%
% If no preferences set, loads from CURTIS LAB default.
% NOTE: this WILL NOT work outside Curtis lab.
%
% To set a new default edf2asc path, do:
% setpref('iEye_ts','edf2asc_path','/path/to/your/binary');
%
% TODO: somehow do an automatic git pull from master at init?
[ret, hostname] = system('hostname');
if ret ~= 0
    hostname = getenv('hostname');
end
hostname = strtrim(hostname);

if ~ispref('iEye','edf2asc_path')
    if strcmp(hostname, 'syndrome') || strcmp(hostname, 'vader') || strcmp(hostname, 'zod.psych.nyu.edu')  ...
            ||strcmp(hostname, 'zod') || strcmp(hostname, 'loki.psych.nyu.edu')
        edf2asc_path = '/d/DATA/hyper/spacebin';
        %edf2asc_path = '/usr/local/bin/edf2asc';
    else
        edf2asc_path = '/Users/mrugankdake/remote/hyper/spacebin';
    end
else
    edf2asc_path = getpref('iEye','edf2asc_path');
end

path1 = getenv('PATH');
path1 = [path1 ':' edf2asc_path];
setenv('PATH', path1);
clear path1;

return