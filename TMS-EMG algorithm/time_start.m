function varargout = time_start(varargin)
% TIME_START MATLAB code for time_start.fig
%      TIME_START, by itself, creates a new TIME_START or raises the existing
%      singleton*.
%
%      H = TIME_START returns the handle to a new TIME_START or the handle to
%      the existing singleton*.
%
%      TIME_START('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TIME_START.M with the given input arguments.
%
%      TIME_START('Property','Value',...) creates a new TIME_START or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before time_start_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to time_start_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help time_start

%% 
% Function description: It refers to the GUI that allows the user to select
% the time instant of effective acquisition (output)

% input arguments: 
% - acqdata, a structure that contains the raw continuous EMG data
% - file, a cell array that contains the string identificative of the participant being analyzed 
%%
% Last Modified by GUIDE v2.5 19-Jul-2019 21:01:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @time_start_OpeningFcn, ...
                   'gui_OutputFcn',  @time_start_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
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


% --- Executes just before time_start is made visible.
function time_start_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to time_start (see VARARGIN)

% enables functionalities of toolbar and menubar for the present GUI
set(handles.figure1,'toolbar','figure');
set(handles.figure1,'menubar','figure');

acqdata = varargin{1}; % input structure
file = varargin{2}; 
handles.sampling_frequency = varargin{3}; 
EMGdata = acqdata.data(:,1); %gets EMG data from channel 1
x_start = acqdata.x_start; %initiate time of effective acquisition (initially is set at 0 seconds) 

% assigns to handles of type text, the name that identifies the patient being analyzed
set(handles.patient,'String',file); 

[ymin,~] = min(EMGdata); %gets the maximum amplitude value of data
[ymax,~] = max(EMGdata); %gets the minimum amplitude value of data 
ylims = [ymin ymax]; %y-axis limits 

handles.ylims=ylims;
handles.acqdata = acqdata;
handles.EMGdata = EMGdata;

plot_figure(EMGdata,handles,x_start);

% Choose default command line output 
handles.output = hObject;

% prompt to really close without saving
set(handles.figure1, 'CloseRequestFcn', @closeGUI); 

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes time_start wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = time_start_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%acqdata = handles.acqdata; %get acqdata structure, that contains the selected initiate time of 
%effective acquisition 
varargout{1} = fix(str2double(get(handles.initiate_time,'string')) * handles.sampling_frequency) + 1;
% output of GUI, relative to the selected initiate time (in samples) 

delete(handles.figure1); % exit GUI 

function closeGUI(hObject, eventdata, handles)
if isequal(get(hObject, 'waitstatus'), 'waiting')
% The GUI is still in UIWAIT, us UIRESUME
uiresume(hObject);
else
% The GUI is no longer waiting, just close it
delete(hObject);
end

%% PLOT FIGURE
function plot_figure(EMGdata,handles,x_start)

% get handles fields to plot figure
ylims=handles.ylims;
acqdata = handles.acqdata;

% buttons to display
set(handles.mouse,'visible','off');
set(handles.mouse,'visible','on');   

y=EMGdata; %EMG data from channel 1
x =(0: length(y)-1)*(1/handles.sampling_frequency); %time scale (in seconds)
axes(handles.graphic); % it gets the plot axis that was design in GUIDE
plot(x,y,'k');

title('Raw EMG data','FontSize',16, 'FontName', 'Arial'); %plot title

xlim([x_start acqdata.x_start+3]); % x-axis range 

%adds x label
xlabel('Time (s)', 'FontSize',10, 'FontName', 'Arial', 'FontWeight', 'bold');

% adds y label
ylabel('ch1(mV)','FontSize',10, 'FontName', 'Arial', 'FontWeight', 'bold');

% y limits
ylim(handles.ylims);

%adds x_start line 
x_time = acqdata.x_start;
line([x_time x_time], ylims(1, :) ,'Color',[1 0 1],'Marker','o');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%% GUI ELEMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

display('saved');
uiresume(handles.figure1); %delevering of the output, and consequent deletetion of present GUI

% --- Executes on button press in mouse
function mouse_Callback(hObject, eventdata, handles)
% hObject    handle to mouse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get handles fields
EMGdata = handles.EMGdata;
acqdata = handles.acqdata; 

axh = gca;  
axh.Toolbar = matlab.ui.controls.AxesToolbar();

%manually selects the initiate time of effective acquisition (in seconds)
x_new = ginput(1);
acqdata.x_start = x_new(1);
x_start = acqdata.x_start;

sample_start = fix(x_new(1)*handles.sampling_frequency); 
set(handles.initiate_time, 'string', x_new(1)); %assigns to handles of type text the selected initiate 
% time of effective acquisition (in seconds)

handles.acqdata = acqdata; % atualizes acqdata structure

guidata(hObject, handles); 
plot_figure(EMGdata,handles,x_start); % plots EMG data from channel 1, with the selected initiate 
%time of effective acquisition 
