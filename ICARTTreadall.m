function D = ICARTTreadall(fldr,vars)
% function D = ICARTTreadall(fldr,vars)
% Reads all ICARTT files in a folder and combines them into a single set of variables.
% For this to work, the following conditions must be met:
% 1) The target folder contains only .ict files
% 2) All .ict files must contain the same variables
%
% INPUT
% fldr is the folder containing the files
% vars is an optional cell array specifying which variables to extract.
% OUTPUT D is the data structure.
%   This will also contain several other variables:
%   "fnum" is an index denoting the file number.
%   "files" is a cell array denoting the files loaded.
%
% 20161018 GMW
% 20170209 GMW  Added provisions to fill in nans for missing variables
% 20180208 GMW  Added option to specify which variabls are output.

if nargin<2, vars=[]; end

% get files
fileinfo = dir(fldr);
L = length(fileinfo)-2; %hack off . and ..
files = cell(L,1);
[files{:}] = fileinfo(3:end).name;

% read em
for i=1:L
    
    % read file
    d = ICARTTreader(files{i});
    d = rmfield(d,'header');
    dnames = fieldnames(d);
    
    % limit to desired variables
    if ~isempty(vars)
        junk = setdiff(dnames,vars);
        d = rmfield(d,junk);
        dnames = fieldnames(d);
    end
    
    % fnum
    ld = length(d.(dnames{1}));
    d.fnum = i*ones(ld,1);
    dnames = fieldnames(d); %need to include fnum
    
    if i==1
        D = d; %initialize
    else
        Dnames = fieldnames(D);
        lD = length(D.fnum);
        
        % deal with missing fields in D
        miss = ~ismember(dnames,Dnames);
        for mname = dnames(miss)'
            D.(char(mname)) = nan(lD,1);
        end
        
        % deal with missing fields in d
        miss = ~ismember(Dnames,dnames);
        for mname = Dnames(miss)'
            d.(char(mname)) = nan(ld,1);
        end
        
        % concatenate
        D = structCat(D,d);
    end
end

D.files = files;


