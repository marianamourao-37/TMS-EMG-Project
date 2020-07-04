function varargout = specifications(varargin)
% SPECIFICATIONS MATLAB code for specifications.fig
%      SPECIFICATIONS, by itself, creates a new SPECIFICATIONS or raises the existing
%      singleton*.
%
%      H = SPECIFICATIONS returns the handle to a new SPECIFICATIONS or the handle to
%      the existing singleton*.
%
%      SPECIFICATIONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPECIFICATIONS.M with the given input arguments.
%
%      SPECIFICATIONS('Property','Value',...) creates a new SPECIFICATIONS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before specifications_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to specifications_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help specifications

%% 
% Function description: It refers to the GUI that allows the user to modify
% the value of certain parameters (input and output of the present GUI): 

% min_amp_MEP - minimum MEP's amplitude (mV) for its quantification;

% latency_MEP - MEP's latency (ms) relative to the last TMS artefact 
% (used to define MEP's window for its detection);

% quality_comb - quality factor of IIR comb filter (nondimensional);

% curve_model - Type of I/O curve Modulation (1 for statistical modulation 
% and 0 for linear regression);

% artefact_detection - select if desired to see the 2nd peaks of the 1st TMS 
%artefact in each rectified EMG segment (0 for no, 1 for yes)

% distribution_tms - select if desired to see the distribution of amplitude of the 
%detected  TMS artefacts (before and after optimization of the TMS
%detection algorithm): 0 for no, 1 for yes;

% agremment_measures - select if desired to see the agreement between automatic 
% and manual quantification (0 for no, and 1 for yes);

% comb - switch applicability of IIR comb filter (0 to not apply, 1 for
% applying);

% constant_time_tms - time constant (ms) to construct TMS templates, related to an 
%estimate of a TMS pulse duration. 
%%
% Last Modified by GUIDE v2.5 30-Oct-2019 23:12:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @specifications_OpeningFcn, ...
                   'gui_OutputFcn',  @specifications_OutputFcn, ...
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


% --- Executes just before specifications is made visible.
function specifications_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to specifications (see VARARGIN)

% state of the respective handles given the input values:

%state of the IIR comb filter handle
handles.comb = varargin{8}; 
if handles.comb == 1
    set(handles.on_comb,'value',1);
    set(handles.off_comb,'value',0);
else
    set(handles.on_comb,'value',0);
    set(handles.off_comb,'value',1);
end

set(handles.min_amp_MEP,'string',varargin{1}); % attribution of minimum MEP amplitude to respective 
% handle
set(handles.latency_MEP,'string',varargin{2}); % attribution of MEP latency to respective handle
set(handles.quality_comb,'string',varargin{3}); % attribution of quality factor to respective handle

%state of the I/O curve modulation 
handles.curve_model = varargin{4}; 
if handles.curve_model == 0 
    set(handles.regression_curve,'value',1);
    set(handles.statistical_curve,'value',0);
else 
    set(handles.regression_curve,'value',0);
    set(handles.statistical_curve,'value',1);
end

% state of various performance demonstration of the algorithm (off or on)
set(handles.artefact_detection,'value',varargin{5}); 
set(handles.distribution_tms,'value',varargin{6});  
set(handles.agremment_measures,'value',varargin{7}); 

set(handles.constant_time_tms,'string',varargin{9}); % attribution of the time duration of TMS pulse 
% to respective handle

set(handles.sampling_frequency,'string',varargin{10});

% Choose default command line output for specifications
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes specifications wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = specifications_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Deletes figure of present GUI after uiresume 
delete(hObject);

function closeGUI(hObject, eventdata, handles)
if isequal(get(hObject, 'waitstatus'), 'waiting')
% The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
% The GUI is no longer waiting, just close it
    delete(hObject);
end

% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%after pressing save, saves the values of the input parameters
if get(handles.on_comb,'value')==1
    handles.comb = 1;
else 
    handles.comb = 0;
end

handles.min_amp_MEP = str2double(get(handles.min_amp_MEP,'string'));
handles.latency_MEP = str2double(get(handles.latency_MEP,'string'));
handles.quality_comb = str2double(get(handles.quality_comb,'string'));

if get(handles.regression_curve,'value') == 1
    handles.curve_model = 0;
else
    handles.curve_model = 1;
end

handles.artefact_detection = get(handles.artefact_detection,'value');
handles.distribution_tms = get(handles.distribution_tms,'value');
handles.agremment_measures = get(handles.agremment_measures,'value');
handles.constant_time_tms = str2double(get(handles.constant_time_tms,'string'));
handles.sampling_frequency = str2double(get(handles.sampling_frequency,'string'));

choosedata(handles.min_amp_MEP,handles.latency_MEP, handles.quality_comb, handles.curve_model, ...
    handles.artefact_detection, handles.distribution_tms, handles.agremment_measures, handles.comb,...
    handles.constant_time_tms,handles.sampling_frequency); % call to next GUI 

uiresume(handles.figure1); %it proceeds to the delivering of the output, and consequent 
%deletetion of present GUI

% --- Executes on button press in default_values.
function default_values_Callback(hObject, eventdata, handles)
% hObject    handle to default_values (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%sets the input parameters and its corresponding handles to its default values 
set(handles.on_comb,'value',1);
set(handles.off_comb,'value',0);
handles.comb = 1;
set(handles.min_amp_MEP,'string',0.05);
set(handles.latency_MEP,'string',18);
set(handles.quality_comb,'string',300);
set(handles.regression_curve,'value',1);
set(handles.statistical_curve,'value',0);
handles.curve_model = 0;
set(handles.artefact_detection,'value',0);
set(handles.distribution_tms,'value',0);
set(handles.agremment_measures,'value',0);
set(handles.constant_time_tms,'string',12);
set(handles.sampling_frequency,'string',10000);

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function min_amp_MEP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to min_amp_MEP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function latency_MEP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to latency_MEP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in change_templates.
function change_templates_Callback(hObject, eventdata, handles)
% hObject    handle to change_templates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get, in seconds, the defined time constant (ms) that accounts the duration of TMS pulse to 
% construct TMS templates
constant_time_tms= 10^(-3)*str2double(get(handles.constant_time_tms,'string'));
sampling_frequency = str2double(get(handles.sampling_frequency,'string'));
templates_tms(constant_time_tms,sampling_frequency); %calls to nest GUI, that will allow to construct TMS templates

% --- Executes during object creation, after setting all properties.
function constant_time_tms_CreateFcn(hObject, eventdata, handles)
% hObject    handle to constant_time_tms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function quality_comb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to quality_comb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sampling_frequency_Callback(hObject, eventdata, handles)
% hObject    handle to sampling_frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sampling_frequency as text
%        str2double(get(hObject,'String')) returns contents of sampling_frequency as a double


% --- Executes during object creation, after setting all properties.
function sampling_frequency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sampling_frequency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function quality_comb_Callback(hObject, eventdata, handles)
% hObject    handle to quality_comb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of quality_comb as text
%        str2double(get(hObject,'String')) returns contents of quality_comb as a double

% --- Executes on button press in off_comb.
function off_comb_Callback(hObject, eventdata, handles)
% hObject    handle to off_comb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of off_comb

% --- Executes on button press in on_comb.
function on_comb_Callback(hObject, eventdata, handles)
% hObject    handle to on_comb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of on_comb

function constant_time_tms_Callback(hObject, eventdata, handles)
% hObject    handle to constant_time_tms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of constant_time_tms as text
%        str2double(get(hObject,'String')) returns contents of constant_time_tms as a double

% --- Executes on button press in agremment_measures.
function agremment_measures_Callback(hObject, eventdata, handles)
% hObject    handle to agremment_measures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of agremment_measures

% --- Executes on button press in distribution_tms.
function distribution_tms_Callback(hObject, eventdata, handles)
% hObject    handle to distribution_tms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of distribution_tms

% --- Executes on button press in artefact_detection.
function artefact_detection_Callback(hObject, eventdata, handles)
% hObject    handle to artefact_detection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of artefact_detection

% --- Executes on button press in regression_curve.
function regression_curve_Callback(hObject, eventdata, handles)
% hObject    handle to regression_curve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of regression_curve

% --- Executes on button press in statistical_curve.
function statistical_curve_Callback(hObject, eventdata, handles)
% hObject    handle to statistical_curve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of statistical_curve

function latency_MEP_Callback(hObject, eventdata, handles)
% hObject    handle to latency_MEP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of latency_MEP as text
%        str2double(get(hObject,'String')) returns contents of latency_MEP as a double

function min_amp_MEP_Callback(hObject, eventdata, handles)
% hObject    handle to min_amp_MEP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of min_amp_MEP as text
%        str2double(get(hObject,'String')) returns contents of min_amp_MEP as a double
