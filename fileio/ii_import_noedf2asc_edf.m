function [ii_data,ii_cfg] = ii_import_noedf2asc_edf(edf_file,ifg_file,data_file, varargin)
% This function would import EDF file without edf2asc.exe, but it requires
% edfmex, which also requires download the EyeLink Developers Kit.
% The edfmex could be found on SR research support websit after log in:
% https://www.sr-support.com/thread-54.html

% The output would be sample structure ii_data and event structure ii_cfg
% for the future use of the following process

% The required input of function are EDF file and IFG file

% The function is adapted from ii_import_edf.m, which is developed by Tommy
% Sprague, latest update on 8/17/2017:"

    % ii_import_edf(edf_file, ifg_file) saves to edf_file_iEye.mat
    %
    % ii_import_edf(edf_file, ifg_file, []) prompts user for name/location of
    % output file
    %
    % ii_import_edf(edf_file, ifg_file, data_file) reads
    % edf_file with configuration ifg_file and saves data into data_file. if
    % any arguments are empty, they are prompted via GUI.
    %
    % [ii_data, ii_cfg] = ii_import_edf(edf_file, ifg_file, data_file) returns
    % ii_data, ii_cfg structures for use by other iEye programs
    %
    % ii_import_edf(edf_file, ifg_file, data_file, 'oldstyle') saves each
    % channel separately, rather than packing into ii_data (for backwards
    % compatibility)"

% Qingqing Yang, 6/2/2022

if nargin < 1 || isempty(edf_file)
    [filename, pathname] = uigetfile('*.edf', 'Select EDF file');
    edf_file = fullfile(pathname, filename);
    if filename==0
        fprintf('User selected cancel\n');
        return
    end
    
end
if nargin < 2 || isempty(ifg_file)
    [filename_ifg, pathname] = uigetfile('*.ifg', 'Select IFG file');
    ifg_file = fullfile(pathname, filename_ifg);
    
    if filename_ifg==0
        fprintf('User selected cancel\n');
        return
    end

    
end


if nargin < 3
    
   data_file = sprintf('%s_iEye.mat',edf_file(1:(end-4)));
    
end

if isempty(data_file)
    iEye_file = sprintf('%s_iEye.mat',edf_file(1:(end-4))); %strrep(edf_file, '.edf', '.mat');
    [filename_data, pathname] = uiputfile(iEye_file, 'Create data file');
    data_file = fullfile(pathname, filename_data);    
    if filename_data==0
        fprintf('User selected cancel\n');
        return
    end    
end


%% GET CONFIG
[nchan,lchan,schan,cfg] = ii_openifg(ifg_file);
nchan = str2num(nchan);
schan = str2num(schan);
vis = lchan;
lchan = textscan(lchan,'%s','delimiter',',');
echan = nchan - 3;

%% Extract Samples
edfStruct_origin = edfmex(edf_file)
edfStruct=struct2cell(edfStruct_origin);

% For comformity, convert single to double precision data
s_num = transpose(double(edfStruct{1,1}.time));
x = transpose(double(edfStruct{1,1}.gx(2,:)));
y = transpose(double(edfStruct{1,1}.gy(2,:)));
% here the x,y should be 0, but it's 100000000, so we adjust them later
pupil = transpose(double(edfStruct{1,1}.pa(2,:)));
% pupil is 0 when it should be 0, no adjustment
size1=size(x);
size2=size(y);
if size1(1) == size2(1)
    for i =1:size1(1)
        if x(i)==100000000
            x(i)=0;
        end
        
        if y(i)==100000000
            y(i)=0;
        end
    end
else
    fprintf('User selected cancel\n');
end

M=[];
M(:,1) =s_num;
M(:,2) =x;
M(:,3) =y;
M(:,4) =pupil;

%% Extract Events
e = edfStruct{2,1};
e1 = struct2cell(e);
sizee1=size(e1);

% Get the message, has varb1 and vval
d =[]
for i = 1:sizee1(3)
    st=cell2mat(e1(35,1,i));
    st = string(st);
    if isempty(st)
        st=" ";
    end
    d=[d;st];
end
st = 0;
% Get the samp_n and MGS, Here is 'MESSAGEEVENT'
E =[];
c = e1([4,sizee1(1)],1,:);
E = transpose(squeeze(c(:,1,:)));

E = [E,d];
mline = 1;
Mess = {};
Mess = cell2mat(Mess);
%%
% GET MSG EVENTS
% Which is 'MESSAGEEVENT' events here
% token = strtok(E);

samp_n = [];
message = [];
token = E(:,2);
sam_tkn=string(E(:,1));
mess_tkn = string(E(:,3));
mline =1;
for v = 1:length(token)
%     dm = strcmpi(token(v),'MSG');
    st = string(token(v));
    dm = strcmp(string(st),'MESSAGEEVENT');
    if dm == 1
%        Mess(mline) = [Mess;st];
        nn=str2double(sam_tkn(v));
        samp_n = [samp_n;nn];
        cr = string(mess_tkn(v));        
        message = [message;string(cr)];
        mline = mline + 1;
    end
end

[varbl, vval] = strtok(message);
vval = str2double(vval);

%% SEARCH MSG EVENTS FOR VARIABLE 
% Below is the same as in the ii_import_edf.m

for i = 4:nchan
    mline = 1;
    cname = lchan{1}{i};
    MV = [];
    
    for v = 1:length(varbl)
        dm = strcmpi(varbl(v),cname); % TS: no longer case-sensitive!!!
        
        if dm == 1
            MV(mline,:) = [samp_n(v) vval(v)];
            mline = mline + 1;
        end
    end
    
    % GET INDICES & SET VALUES
    li = 1;
    ci = 1;
    cv = 0;
    M(:,(i+1)) = 0;
    
    for h = 1:length(MV)
        ci = find(M(:,1)==MV(h,1));
        if isempty(ci) == 0
            M((ci:length(M)),(i+1)) = MV(h,2);
            li = ci;
        else
            MV(h,1) = MV(h,1) - 1;
            ci = find(M(:,1)==MV(h,1));
            M((ci:length(M)),(i+1)) = MV(h,2);
            li = ci;
        end
    end
end

%% Create File Matrix
M(:,1) = [];

% also insert channel data into eyedata, to be saved out into data_file
eyedata = [];


for i = 1:nchan
    cname = lchan{1}{i};
    cvalue = M(:,i);
    eyedata.(lchan{1}{i})=M(:,i);
end

x = M(:,1);

%% Creat ii_cfg Structure
dt = datestr(now, 30);%'mmmm dd, yyyy HH:MM:SS.FFF AM');

ii_cfg.cursel = [];
ii_cfg.sel = x*0;
ii_cfg.cfg = cfg;
ii_cfg.vis = vis;
ii_cfg.nchan = nchan;
ii_cfg.lchan = lchan;
ii_cfg.hz = schan;
ii_cfg.velocity = [];
ii_cfg.tcursel = [];
ii_cfg.tsel = x*0;
ii_cfg.tindex = 0;
ii_cfg.saccades = [];
ii_cfg.history{1} = ['EDF imported ', dt];

ii_cfg.edf_file = edf_file;

ii_cfg.microsacc =[]; 

ii_data = eyedata;
%% Save File
if nargin >= 3 || (nargin<3 && nargout == 0)
    if ~isempty(varargin) && strcmpi(varargin{1},'oldstyle')
        save(data_file,'ii_cfg','edf_file','M','ii_data');
        save(data_file,'-struct','eyedata','-append');
    else % default save state - ii_cfg and ii_data encapsulate all info
        save(data_file,'ii_cfg','ii_data','edf_file');
    end
end

% maybe will help us avoid weird errors?
fclose('all');

end