function filename = ICARTTwriter_example(flightDay,data,rev,ICTdir,specialComments)
% function filename = ICARTTwriter_example(flightDay,data,rev,ICTdir,specialComments)
% Generates an ICARTT file for ISAF HCHO.
% Note that revision special comments can be entered near the top of the script as needed.
% For more info on format, see the ICARTT documentation:
% http://www-air.larc.nasa.gov/missions/etc/IcarttDataFormat.htm
%
% INPUTS:
% flightDay:     UTC beginning date for data, yyyymmdd.
%               If multiple sorties were flown in one day, this can also have a 3-character sortie tag, eg. _L1
% data:         5-column cell array containing data and information.
%               Note, it is assumed that first row is the independent variable (typically, time).
%                   1st column: names of variables to write.
%                   2nd column: actual data. All should be vectors of the same length, and missing values should be NaNs.
%                   3rd column: format strings. Default value is '%6.3f'. For more info, see help for fprintf.
%                   4th column: units string
%                   5th column: variable description string (optional)
% rev:          revision letter (for preliminary data) or number (for final data). e.g. 'RA' or 'R0'.
% ICTdir:       full path for save directory.
% specialComments: cell array of special comments for flight day (optional).
% 
% OUTPUT is the name of the file.
%
% EXAMPLE USE
%
% flightDay = '20220123';
% x = (1:100)';
% y = rand(100,1);
% ICTrev = 'R1';
% ICTdir = cd;
% data = {...
%     % var name      % data      %format     %units      %description
%     'Time_Start'          x     '%6.0f'     's'         'Time_Start, start time of measurement'
%     'Ozone'               y     '%6.3f'     'ppbv'      'Gas_O3_insitu_S_AVMR, ozone mixing ratio'
%     };
% specialComments = {
%    'Instrument failure after first hour of flight, possibly due to poor thermal management.'
%    'Fun fact: special comments can be more than one line.'
%    };
% flnm = ICARTTwriter_example(flightDay,data,ICTrev,ICTdir,specialComments);
% %
% 20160730 GMW  Modified from ICARTTwriter_example.
% 20170310 GMW  Modified to accept cell array for data, info and format
% 20170517 GMW  Modified revComments handlng from switch to cell array with string checking.
% 20170607 RAH  Added provision to save filenames with multiple sorties in one day.
% 20210812 GMW  Added specialComments input
% 20220107 GMW  Removed sampleRate input upon converting this back to an example.


%% MEASUREMENT-SPECIFIC INFO

% FILENAME INFO
dataID      = 'SOOI-NOZE-OZONE';        %DATA ID (1st part of filename)
locID       = 'WB57';                   %LOCATION ID (2nd part of filename)

% TOP OF HEADER
HeaderIntro = {...
    'PI: Wolfe, Glenn';...                                              % PI name
    'ORGANIZATION: NASA';...                                            % Organization
    'DATA SOURCE: NOZE (Nose Ozone Experiment)';...                     % Data Source
    'MISSION: SOOI (Stratospheric Ozone Olfaction Investigation)';...   % Mission name
    '1, 1';...                                                          % volume number, number of file volumes
    '1';...                                                             % time interval (see documentation)
    };

% DATA FLAGS
MISSflag = -9999;
ULODflag = -7777;
LLODflag = -8888;

% NORMAL COMMENTS
NormalComments = {...
    'PI_CONTACT_INFO: NASA GSFC, Atmospheric Chemistry and Dynamics Lab, Greenbelt, MD 20771; glenn.m.wolfe@nasa.gov';...
    'PLATFORM: NASA WB-57 (N926NA)';...
    'LOCATION: Houston, TX, USA';...
    'ASSOCIATED_DATA: Aircraft position and tongue-based temperature data available on the mission archive. Data available at DOI ###';...
    ['INSTRUMENT_INFO: Graduate student installed in starboard superpod forebody, sniffing air from an external probe.'];...
    ['DATA_INFO: Native sampling rate is 12 breaths/min. ' ...
        '10 s Gaps occur every 30 minutes so the student can do a background measurement by sniffing a cylinder of zero air. ' ...
        'Occasional larger gaps of ~30 seconds due to bathroom breaks. '...
        'Data is not time-aligned. '...
        'Actual response time (1/e flush time) is 0.5 Hz due to olfactory hysteresis.'];...
    ['UNCERTAINTY: Accuracy (systematic uncertainty) is 10%%. '...
        'Precision is +/- 100 ppbv (signal/noise = 1). '...
        'When averaging, precision reduces as (number of points)^0.5, but accuracy does not reduce. '...
        'Nominal detection limit for signal/noise = 1 is 100 ppbv at 1 Hz.'];...
    ['ULOD_FLAG: ' num2str(ULODflag)];...
    'ULOD_VALUE: N/A';...
    ['LLOD_FLAG: ' num2str(LLODflag)];...
    'LLOD_VALUE: N/A';...
    'DM_CONTACT_INFO: Ima Student, College, ima.student@college.edu';...
    'PROJECT_INFO: SOOI is the Stratospheric Ozone Olfaction Investigation. January 2022. More info at https://website.com.';...
    'STIPULATIONS_ON_USE: Please contact the PI or DM for proper use. All NASA data use policies apply.';...
    'OTHER_COMMENTS: This food tastes funny now.';...
    ['REVISION: ' rev];...
    };

% REVISION COMMENTS
% Latest revision should appear at top of list.
% Note each revision comment is in a separate cell.
revComments = {...
    'R1: This is finaler data.';...
    'R0: This is final data.';...
    'RA: This is field data.';...
    };

%% INPUT CHECKING

% revisions
assert(ischar(rev),'Input "rev" must be a string, e.g. ''RA'' or ''R0''.');
revChk = strncmp(rev,revComments,2);
assert(any(revChk),[mfilename ': rev "' rev '" not recognized. Add comments to Revision list.'])
revComments = revComments(find(revChk,1):end); %include all comments below this one
NormalComments = [NormalComments; revComments];

% data
if ~iscell(data)
    error('Input "data" must be a cell array.')
end

% break it out
names = data(:,1);
vars  = data(:,2);
forms = data(:,3);
units = data(:,4);
info  = data(:,5);

%check variables
L = length(vars{1});
numvar = length(names);
for i=1:numvar
    
    %check variable length
    d = vars{i};
    if ~isvector(d) || length(d)~=L
        error(['Data field ' names{i} ' must be same size as time.'])
    end
    vars{i} = vars{i}(:); %ensure column vector
    
    %check format
    if isempty(forms{i})
        disp(['CAUTION: Format string missing for variable ' names{i}])
        forms{i} = '%6.3f';
    end
    
    %check units
    if isempty(units{i})
        disp(['CAUTION: Units missing for variable ' names{i}])
        units{i} = 'none';
    end
    
    %check info
    if isempty(info{i})
        disp(['CAUTION: Description missing for variable ' names{i}])
        info{i} = '';
    end

end

%check directory
if ~isfolder(ICTdir)
    yn = input(['ICTdir ' ICTdir ' does not exist. Create? y/n [y]: '],'s');
    if isempty(yn) || yn=='y'
        mkdir(ICTdir)
    else
        disp(['Invalid save path ' ICTdir '. File not saved.']);
        return
    end
end

%% DATA HANDLING

%build variable names string and formatting string
% nameStr = []; fStr = [];
% for i=1:length(names)
%     nameStr = [nameStr ', ' names{i}]; %column names
%     fStr = [fStr ', ' forms{i}]; %data formatting string
% end
nameStr = strjoin(names,', ');
fStr = [strjoin(forms,', ') '\r\n'];
NormalComments = [NormalComments; nameStr];

% extract independent variable
time = vars{1};
time_name = names{1};
time_units = units{1};
time_info = info{1};
names(1)=[]; vars(1)=[]; forms(1)=[]; units(1)=[]; info(1)=[];
numvar = length(names);

% convert data to matrix
data = cell2mat(vars'); %single matrix, with one column for each variable

% replace missing data
data(isnan(data)) = MISSflag;
missStr = repmat([int2str(MISSflag) ', '],1,numvar);
missStr = missStr(1:end-2);

%scaling factor
scaleStr = repmat('1, ',1,numvar);
scaleStr = scaleStr(1:end-2);

% add info to HeaderIntro
HeaderIntro = [...
    HeaderIntro;...
    [time_name ', ' time_units ', ' time_info];...  % independent variable info
    int2str(numvar);...                             % Number of dependent variables
    scaleStr;...                                    % Scaling factors
	missStr;...                                     % Missing data indicator
    ];

%% DATES

% SORTIE TAG
if length(flightDay)>8
    sortie = flightDay(9:end);
    flightDay = flightDay(1:8);
else
    sortie = '';
end

% format dates
fltDateForm = datestr(datenum(flightDay,'yyyymmdd'),'yyyy, mm, dd'); %formatted
revDate     = datestr(now,'yyyy, mm, dd'); %revision date

% stuff into intro array
HeaderIntro = [HeaderIntro(1:5);...
    [fltDateForm ', ' revDate];...
    HeaderIntro(6:end)];

%% LINE NUMBERS
numintro  = length(HeaderIntro);               %number of beginning header lines
numspec  = length(specialComments);               %number of special comments
numnorm  = length(NormalComments);             %number of normal comments
numother = 3;                                  %other lines
numlines = numintro + numvar + numspec + numnorm + numother; %number of lines in header

%% PRINT FILE

filename    = [dataID '_' locID '_' flightDay '_' rev sortie '.ict'];
filepath    = fullfile(ICTdir,filename);
[fid,message] = fopen(filepath,'w'); %open file and overwrite if it exists
if fid==-1
    disp(message)
    return
end

fprintf(fid,[int2str(numlines) ', 1001\r\n']);      % Number of header lines, file format index

for i=1:numintro, fprintf(fid,[HeaderIntro{i} '\r\n']); end                       % Intro Comments

for i=1:numvar, fprintf(fid,[names{i} ', ' units{i} ', ' info{i} '\r\n']); end    % Dependent variables

fprintf(fid,[int2str(numspec) '\r\n']);                                     % # special comments
for i=1:numspec, fprintf(fid,[specialComments{i} '\r\n']); end                 % Special comments

fprintf(fid,[int2str(numnorm) '\r\n']);                                     % # normal comments
for i=1:numnorm, fprintf(fid,[NormalComments{i} '\r\n']); end               % Normal comments

%print data
for i=1:L
    fprintf(fid,fStr,time(i),data(i,:));
end

% close file
status = fclose(fid);
if status
    disp('Problem closing file.');
end
disp(['ICARTT file written: ' filepath]);


