function D = ICARTTmerge(files,timeInfo,alignInfo)
% function D = ICARTTmerge(files,timeInfo,alignInfo)
% tool to merge ICARTT files.
% 
% INPUTS:
% files: Can be either:
%   1) full path to directory containing files to merge. In this case, All files to be merged are in a single directory containing no other files.
%   2) cell array of file names
% timeInfo: two or three-element cell array containing info for the master time vector.
%   First cell is the data ID for the file containing the master time vector.
%   Second cell is the name of the time variable (or start time if using a window merge).
%   Third cell (optional) is name of stop time if using a window merge.
%
% alignInfo: OPTIONAL cell array indicating which variables to align by lag-covariance.
%   First cell is master alignment variable name.
%   Subsequent cells are names of variables to align.
%
% OUTPUT D is the merged dataset structure.
%
% EXAMPLE USE
% files = cd;
% timeInfo = {'FIREXAQ-MetNav','Time_Start'};
% alignInfo = {'CO_LGR_ppb','CH2O_CAMS_pptv','CH2O_ISAF'};
% D = ICARTTmerge(fldr,timeInfo,alignInfo);
%
% 20190803 GMW
% 20190821 GMW  Added ability to do a window merge (e.g. on WAS data).
% 20200311 GMW  Added option to specify file names instead of folder.

if nargin<3, alignInfo = {}; end
    
% get file names, if folder specificed
if ~iscell(files)
    fileinfo = dir(files);
    L = length(fileinfo)-2; %hack off . and ..
    files = cell(L,1);
    [files{:}] = fileinfo(3:end).name;
else
    L = length(files);
end

% read em
d = struct;
for i=1:L
    dataID = strtok(files{i},'_');
    dataID = strrep(dataID,'-','_'); %make OK for var name
    d.(dataID) = ICARTTreader(files{i});
    d.(dataID) = rmfield(d.(dataID),'header');
end

% get time variable to average to
timeDataID = strrep(timeInfo{1},'-','_');
t1name = timeInfo{2};
assert(isfield(d,timeDataID),'timeInfo dataID %s not found.',timeDataID)
assert(isfield(d.(timeDataID),t1name),'timeInfo variable name %s not found in dataID %s.',t1name,timeDataID)
D = d.(timeDataID); %initialize output structure
t1 = D.(t1name);

% get stop time for window if indicated
if length(timeInfo)==3
    t1name_stop = timeInfo{3};
    assert(isfield(d.(timeDataID),t1name_stop),'timeInfo variable name %s not found in dataID %s.',t1name_stop,timeDataID)
    t1stop = D.(t1name_stop);
else
    t1stop=[];
end

% adjust for time offsets (assuming 1 Hz)
% BinAvg routine assumes midpoint
if isempty(t1stop)
    switch(lower(t1name))
        case 'time_start'
            t1 = t1+0.5;
        case 'time_mid'
        case 'time_stop'
            t1 = t1-0.5;
    end
end

% average other variables
d = rmfield(d,timeDataID);
dnames = fieldnames(d);
for i = 1:length(dnames)
    dd = d.(dnames{i});
    ddnames = fieldnames(dd);
    t2name = ddnames{1};
    t2 = dd.(t2name);
    
    % adjust for time offsets (assuming 1 Hz)
    % assumes 1s data
    if isempty(t1stop)
        switch(lower(t2name))
            case 'time_start'
                t2 = t2+0.5;
            case 'time_mid'
            case 'time_stop'
                t2 = t2-0.5;
        end
    end
    
    % average
    for j = 2:length(ddnames)
        % avoid overwriting variables if names are the same
        n = ddnames{j};
        while isfield(D,n)
            n = [n '_'];
        end
        
        if isempty(t1stop)
            D.(n) = BinAvg(t2,dd.(ddnames{j}),t1);
        else
            D.(n) = BinAvg(t2,dd.(ddnames{j}),[t1 t1stop]);
        end
    end
    
end

% time alignment
if ~isempty(alignInfo)
    alignName = alignInfo{1};
    alignVars = alignInfo(2:end);
    assert(isfield(D,alignName),'Alignment variable %s not found in files.',alignName)
    nlags = 20; %# of lag points
    for i = 1:length(alignVars)
        if ~isfield(D,alignVars{i}), continue; end
        
        D.(alignVars{i}) = lagit(D.(alignName),D.(alignVars{i}),nlags,0);
    end
end



