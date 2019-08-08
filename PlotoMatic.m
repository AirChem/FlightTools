function varargout = PlotoMatic(varargin)
% PLOTOMATIC M-file for plotomatic.fig
% PlotoMatic is a graphical user interface for visualizing airborne in situ observations.
% Such data is 5-dimensional (x,y,z,t, and scalars), and there are many different ways to slice it.
% This GUI will simultaneously make a time series, vertical profile, map, and scalar-scalar scatter plot.
% Additional functionality is included to change data scaling and index discrete segments (legs).
% Suggestions for improvement are always welcome.
%
% EXAMPLE USE
% D = ICARTTmerge(someDirectoryWithData);
% Xnames = {''Time_Start'',''Latitude'',''Longitude'',''Altitude''};
% legTimes = [12345 67899];
% PlotoMatic(D,Xnames,legTimes)
%
% INPUTS
% D: a structure containing data vectors, e.g. from a merge file.
%       All fields must be 1-D arrays of the same length.
% Xnames: cell array of strings specifying names of variables within D for independent variables: 
%       {time,latitude,longitude,altitude}.
% legTimes: OPTIONAL 2-column matrix of start and stop times for individual legs. Default is [].
%
% PLOT OPTIONS MENU
% X1 and X2 dropdown menus: select primary and secondary variables.
%   X1 is displayed in all plots. X2 is currently only used in the scatter and map plots.
% LOG: transform the variable to a base-10 log scale.
% LAG: specify how many data points X1 or X2 should be lagged forward (positive) or backward (negative) in time. 
%   Use if variables are not properly aligned in time, as is often the case with preliminary data.
% MATH: performs difference or ratio between X1 and X2. Be sure that variables are aligned before using this.
% RESET AXES: resets all plot axes to their outer limits.
%
% LEG SELECTION MENU
% Leg Table: enter start and stop times here to flag specific segments of data.
%   Clicking the "Plot" box will plot the leg as a distinct color on all plots.
%   Note, Plot will not work if Stop>Start.
%
% Dump: prints all leg times to the command window for cutting and pasting.
%  qCan be used with plumeAverage.m or other functions.
% Add Leg: uses the data cursors on the time series plot to add a new leg to the table.
% Plot All: toggles between plotting all or none of the valid legs.
%
% COLOR BY MENU
% Found above the map, this lets you determine how the data on the map are colored.
%
% PLOT NAVIGATION
% Hot keys are coded in for quickly switching between plot tools.
% CTRL+Z: zoom
% CTRL+X: zoom out
% CTRL+A: pan
% CTRL+d: data cursors
%
% 20170929 GMW
% 20190805 GMW Initiial public offering.
% 20190807 GMW Added input checking and better readme.

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PlotoMatic_OpeningFcn, ...
                   'gui_OutputFcn',  @PlotoMatic_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before PlotoMatic is made visible.
function PlotoMatic_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

handles.output = hObject; % Choose default command line output for PlotoMatic
guidata(hObject, handles); % Update handles structure

% break out varargin and check inputs
assert(length(varargin)>=2,'Not enough input arguments.')
assert(length(varargin)<=3,'Too many input arguments.')
D = varargin{1};
Xnames = varargin{2};
assert(isstruct(D),'First input must be a structure containing 1-D arrays of the same size.')
assert(iscell(Xnames) && length(Xnames)==4,'Second input must be a 4-element cell array containing field names for time, latitude, longitude, and altitude, respectively.')
assert(all(cellfun(@ischar,Xnames)),'All cells in input Xnames must be strings.')
assert(all(isfield(D,Xnames)),'All variables in Xnames must be in input data structure.')

if length(varargin)==3
    legTimes = varargin{3};
    assert(size(legTimes,2)==2,'Input legTimes must be of size N x 2.')
else
    legTimes = [];
end

% filter data structure
c1 = structfun(@isnumeric,D);
c2 = structfun(@isvector,D);
l = structfun(@length,D);
c3 = l==mode(l);
c = c1 & c2 & c3;
if any(~c)
    Dnames = fieldnames(D);
    D = rmfield(D,Dnames(~c));
end
Dnames = fieldnames(D);

% get independent variables
time = D.(Xnames{1}); time_label = strrep(Xnames{1},'_',' ');
lat = D.(Xnames{2});
lon = D.(Xnames{3}); lon(lon>180) = lon(lon>180) - 360;
alt = D.(Xnames{4});  alt_label = strrep(Xnames{4},'_',' ');

% initialize variable selectors
set(handles.Var1,'String',Dnames);
set(handles.Var2,'String',Dnames);
x1n = Dnames{get(handles.Var1,'Value')};
x2n = Dnames{get(handles.Var2,'Value')};
x1 = D.(x1n);
x2 = D.(x2n);

set(handles.VarMath,'String',{'none','X1-X2','X1/X2'})

% initialize leg table
if ~isempty(legTimes)  
    nlegs = size(legTimes,1);
    legTable = get(handles.legTable,'Data');
    legTable(1:nlegs,1:2) = num2cell(legTimes);
    legTable(1:nlegs,3) = num2cell(true(nlegs,1));
    set(handles.legTable,'Data',legTable)
end

% add info to guidata
handles.D        = D;
handles.Dnames   = Dnames;
handles.time     = time;
handles.lat      = lat;
handles.lon      = lon;
handles.alt      = alt;
handles.alt_label = alt_label;
handles.time_label = time_label;
guidata(hObject,handles);

% update plots
UpdatePlots(handles,1)

% --- Executes on any button press...
function UpdatePlots(handles,firstCall,newX1)

if nargin<3, newX1 = 0; end
if nargin<2, firstCall = 0; end

% get variables
x1n = handles.Dnames{get(handles.Var1,'value')};
x2n = handles.Dnames{get(handles.Var2,'value')};
x1 = handles.D.(x1n);
x2 = handles.D.(x2n);
x1 = x1(:); %ensure column vectors
x2 = x2(:);
x1n = strrep(x1n,'_',' '); %for plot labels
x2n = strrep(x2n,'_',' ');

% lag if desired
lag = str2double(handles.x1lag.String);
if lag
    n = nan(abs(lag),1); %filler
    if lag<0,     x1 = [x1(-lag+1:end); n]; %shift backward
    elseif lag>0, x1 = [n; x1(1:end-lag)]; %shift forward
    end
end

lag = str2double(handles.x2lag.String);
if lag
    n = nan(abs(lag),1); %filler
    if lag<0,     x2 = [x2(-lag+1:end); n]; %shift backward
    elseif lag>0, x2 = [n; x2(1:end-lag)]; %shift forward
    end
end

% set scale
if get(handles.var1log,'value')
    x1 = log10(x1);
    x1n = ['log(' x1n ')'];
end

if get(handles.var2log,'Value')
    x2 = log10(x2);
    x2n = ['log(' x2n ')'];
end

% do math
switch handles.VarMath.Value
    case 1 % none
    case 2 % X1 - X2
        x1 = x1 - x2;
        x1n = [x1n '-' x2n];
    case 3 % X1/X2
        x1 = x1./x2;
        x1n = [x1n '/' x2n];
    otherwise
        keyboard
end

all_color = [0.8 0.8 0.8]; %color for all data

% get axes limits, cursor info
if ~firstCall
    ylim_VP = ylim(handles.VerticalProfile);
    xlim_TS = xlim(handles.TimeSeries);
    xlim_MP = xlim(handles.Map);
    ylim_MP = ylim(handles.Map);
    
    cursorInfo = getCursorInfo(datacursormode(handles.figure1));
end

% vertical profile
plot(handles.VerticalProfile,x1,handles.alt,'*','color',all_color)
xlabel(handles.VerticalProfile,x1n)
ylabel(handles.VerticalProfile,handles.alt_label)

% time series
hTS = plot(handles.TimeSeries,handles.time,x1,'.-','color',all_color);
xlabel(handles.TimeSeries,handles.time_label)
ylabel(handles.TimeSeries,x1n)
set(handles.TimeSeries,'XGrid','on','YGrid','on')

% map (currently set up for US)
if exist('IBWread','file')
    mlat = IBWread('USstates_lat.ibw'); mlat=mlat.y;
    mlon = IBWread('USstates_lon.ibw'); mlon=mlon.y;
    plot(handles.Map,mlon,mlat,'-','color',[0.5 0.5 0.5])
    hold(handles.Map,'on');
end
plot(handles.Map,handles.lon,handles.lat,'-','color',all_color)
set(handles.Map,'xlim',[min(handles.lon) max(handles.lon)])
set(handles.Map,'ylim',[min(handles.lat) max(handles.lat)])
set(handles.Map,'XGrid','on','YGrid','on')
xlabel(handles.Map,'Longitude')
ylabel(handles.Map,'Latitude')

% map colors
Nmax = 1e6; % scatter eats memory with too many points
if length(x1)<=Nmax
    i2 = true(size(x1));
else
    di = floor(length(x1)/Nmax); %increment
    i2 = 1:di:length(x1);
end

switch handles.MapColorSelector.SelectedObject.String
    case 'none'
        % do nothing
    case 'X1'
        scatter(handles.Map,handles.lon(i2),handles.lat(i2),25,x1(i2),'fill');
        c = colorbar(handles.Map,'EastOutside');
        title(c,'X1')
    case 'X2'
        scatter(handles.Map,handles.lon(i2),handles.lat(i2),25,x2(i2),'fill');
        c = colorbar(handles.Map,'EastOutside');
        title(c,'X2')
    case 'leg'
        % handled below
end

% scatter plot
plot(handles.ScatterPlot,x1,x2,'*','color',all_color);
xlabel(handles.ScatterPlot,x1n)
ylabel(handles.ScatterPlot,x2n)
set(handles.ScatterPlot,'XGrid','on','YGrid','on')

% add individual legs
legTable = get(handles.legTable,'Data');
plotflag = cell2mat(legTable(:,3));
if any(plotflag)
    legTimes = cell2mat(legTable(plotflag,1:2));
    nLegs = sum(plotflag);
    cLegs = getLegColors(50);
    cLegs = cLegs(plotflag,:); %keep same colors for each leg
    
    hold(handles.VerticalProfile,'on')
    hold(handles.TimeSeries,'on')
    hold(handles.ScatterPlot,'on')
    hold(handles.Map,'on');
    
    for i = 1:nLegs
        j = handles.time>=legTimes(i,1) & handles.time<=legTimes(i,2);
        plot(handles.VerticalProfile,x1(j),handles.alt(j),'*','color',cLegs(i,:))
        plot(handles.TimeSeries,handles.time(j),x1(j),'.','color',cLegs(i,:))
        plot(handles.ScatterPlot,x1(j),x2(j),'*','color',cLegs(i,:))
        
        if strcmp(handles.MapColorSelector.SelectedObject.String,'leg')
            plot(handles.Map,handles.lon(j),handles.lat(j),'.','color',cLegs(i,:))
        end
    end
    
    rowsPerColumn = 5;
    legend(handles.TimeSeries,[{'all'}; cellstr(num2str(find(plotflag)))],'NumColumns',ceil((nLegs+1)/rowsPerColumn))
end

hold(handles.VerticalProfile,'off')
hold(handles.TimeSeries,'off')
hold(handles.ScatterPlot,'off')
hold(handles.Map,'off');

% set axes limits
if ~firstCall
    set(handles.VerticalProfile,'ylim',ylim_VP);
    set(handles.TimeSeries,'xlim',xlim_TS);
    set(handles.Map,'xlim',xlim_MP);
    set(handles.Map,'ylim',ylim_MP);
end

% fonts
f = 14;
set(handles.VerticalProfile,'FontSize',f)
set(handles.TimeSeries,'FontSize',f)
set(handles.Map,'FontSize',f)
set(handles.ScatterPlot,'FontSize',f)

% Time Series data tips
if firstCall | newX1
    i = find(~isnan(x1+handles.time),1,'first');
    pos1 = [handles.time(i) x1(i) 0];
    i = find(~isnan(x1+handles.time),1,'last');
    pos2 = [handles.time(i) x1(i) 0];
else
    pos1 = cursorInfo(1).Position;
    pos2 = cursorInfo(2).Position;
end
cursorObj = datacursormode(handles.figure1);
dTip1 = createDatatip(cursorObj,hTS);
dTip1.Position = pos1;
dTip2 = createDatatip(cursorObj,hTS);
dTip2.Position = pos2;

% --- Executes when entered data in editable cell(s) in legTable.
function legTable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to legTable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

% check that leg times are valid if plotting
legTable = get(hObject,'Data');
plotflag = cell2mat(legTable(:,3));
if any(plotflag)
    for r = find(plotflag')
        dt = legTable{r,2} - legTable{r,1};
        if isempty(dt) || isnan(dt) || dt<=0
            legTable{r,3} = false;
        end
    end
    set(handles.legTable,'Data',legTable)
end
UpdatePlots(handles)

% --- Executes on button press in DumpButton (DUMP TO CW).
function DumpButton_Callback(hObject, eventdata, handles)

% dump indices to command window for easy cut 'n paste
legTable = get(handles.legTable,'Data');
legTimes = legTable(:,1:2);
i = ~sum(cellfun(@isempty,legTimes),2); %flag rows with valid start/stop times
legTimes = cell2mat(legTimes(i,:))';
outtie = sprintf('%1.2f, %1.2f;\n',legTimes);
outtie = [sprintf('%s\n','legTimes=[') outtie '];'];
disp(outtie)

% --- Executes on button press in PlotAllLegs.
function PlotAllLegs_Callback(hObject, eventdata, handles)

legTable = get(handles.legTable,'Data');
plotflag = cell2mat(legTable(:,3));
legTimes = legTable(:,1:2);
i = ~sum(cellfun(@isempty,legTimes),2); %flag rows with valid start/stop times
if all(plotflag(i))
    plotflag(i) = 0; % if all are checked, then uncheck all.
else
    plotflag(i) = 1; % otherwise, check all.
end
legTable(:,3) = num2cell(plotflag);
set(handles.legTable,'Data',legTable)
legTable_CellEditCallback(handles.legTable, [], handles)

% --- Executes on button press in resetAxes.
function resetAxes_Callback(hObject, eventdata, handles)
UpdatePlots(handles,1)

% --- Executes on button press in AddLeg.
function AddLeg_Callback(hObject, eventdata, handles)

% get cursor x locations
dcmObj = datacursormode(handles.figure1);
c = getCursorInfo(dcmObj);
legLoc = sort([c(1).Position(1) c(2).Position(1)]);

% add to leg Table
legTable = get(handles.legTable,'Data');
legTimes = legTable(:,1:2);
firstEmpty = find(sum(cellfun(@isempty,legTimes),2)==2,1,'first');
legTable{firstEmpty,1} = legLoc(1);
legTable{firstEmpty,2} = legLoc(2);
legTable{firstEmpty,3} = true;
set(handles.legTable,'Data',legTable)
legTable_CellEditCallback(handles.legTable, [], handles)

% --- UI Controls that cause plots to update when changed
function Var1_Callback(hObject, eventdata, handles)
UpdatePlots(handles,0,1)
function Var2_Callback(hObject, eventdata, handles)
UpdatePlots(handles)
function var1log_Callback(hObject, eventdata, handles)
UpdatePlots(handles)
function var2log_Callback(hObject, eventdata, handles)
UpdatePlots(handles)
function VarMath_Callback(hObject, eventdata, handles)
UpdatePlots(handles,0,1)
function x1lag_Callback(hObject, eventdata, handles)
UpdatePlots(handles)
function x2lag_Callback(hObject, eventdata, handles)
UpdatePlots(handles)
function MapColorSelector_SelectionChangedFcn(hObject, eventdata, handles)
UpdatePlots(handles)

% ----MENUS------------------------------------------
function ToolbarHotKeysMenu_Callback(hObject, eventdata, handles)
function HotKey_ZoomIn_Callback(hObject, eventdata, handles)
zoom on
function HotKey_ZoomOut_Callback(hObject, eventdata, handles)
zoom out
function HotKey_Pan_Callback(hObject, eventdata, handles)
pan on
function HotKey_DataCursor_Callback(hObject, eventdata, handles)
datacursormode on
function READMEMenu_Callback(hObject, eventdata, handles)
doc('PlotoMatic')

% function to set leg colors, generated from distinguishable_colors.m
function c = getLegColors(n)
c = [
         0         0    1.0000
         0    1.0000         0
    1.0000         0         0
         0         0    0.2100
    1.0000         0    0.6900
         0    0.3100         0
    1.0000    0.8300         0
         0    0.6200    1.0000
    0.6200    0.3100    0.2800
    0.2100    1.0000    0.7200
    0.4500    0.2100    0.7200
    0.0700    0.4800    0.5200
    1.0000    0.6900    1.0000
    0.4800    0.7200    0.1000
    0.8300    0.4500         0
    0.7600    0.6900    0.4100
    0.9000         0    1.0000
    0.1400    0.1000         0
    0.9300    0.0300    0.3400
    0.4800         0    0.3400
    0.3100    0.9700    1.0000
    0.4800    0.4100    0.5900
    0.3400    0.6600    0.4500
    0.4100    0.3100         0
    0.8600    1.0000         0
    0.6200         0         0
    1.0000    0.6900    0.6900
    0.7900    1.0000    0.6200
         0    0.2800    0.6200
    0.9300    0.4500    1.0000
    0.5900    0.7900    0.9700
    0.9300    0.4100    0.6200
    0.4100    1.0000    0.4500
    0.5200    0.4800    0.4100
    1.0000    0.4800    0.3800
    0.1400    0.4500    1.0000
    0.2400         0    0.1000
    0.0700    0.2100    0.2800
         0         0    0.3800
         0    0.2400    0.9300
    0.6600    0.6200         0
    0.5900    0.4800    0.8600
         0    0.6900    0.2800
    0.6900    0.9300    0.7900
    1.0000    1.0000    0.4500
    0.1700    0.3100    0.2100
    0.6900    0.1000    0.6200
    1.0000    0.7600    0.3400
    0.4100    0.5500    0.1700
    0.6200         0    0.2400
    ];

if n>length(c)
    c = repmat(c,ceil(n./length(c)),1); %wrap
end 

c = c(1:n,:);


%%%%% UNUSED JUNK FOLLOWS %%%%%
% --- Executes during object creation, after setting all properties.
function Var2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function Var1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function VarMath_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function x1lag_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function x2lag_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Outputs from this function are returned to the command line.
function varargout = PlotoMatic_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;


% THAT'S ALL FOLKS
