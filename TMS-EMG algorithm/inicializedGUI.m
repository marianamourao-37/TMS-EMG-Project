function varargout = inicializedGUI(varargin)
% INICIALIZEDGUI MATLAB code for inicializedGUI.fig
%      INICIALIZEDGUI, by itself, creates a new INICIALIZEDGUI or raises the existing
%      singleton*.
%
%      H = INICIALIZEDGUI returns the handle to a new INICIALIZEDGUI or the handle to
%      the existing singleton*.
%
%      INICIALIZEDGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INICIALIZEDGUI.M with the given input arguments.
%
%      INICIALIZEDGUI('Property','Value',...) creates a new INICIALIZEDGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before inicializedGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to inicializedGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help inicializedGUI

%% 
%Function description: It refers to the GUI that allows the user to choose 1 
%of 3 analysis modules (respective to the wanted protocol) 
%%

% Last Modified by GUIDE v2.5 30-Oct-2019 12:03:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @inicializedGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @inicializedGUI_OutputFcn, ...
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


% --- Executes just before inicializedGUI is made visible.
function inicializedGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to inicializedGUI (see VARARGIN)

% Choose default command line output for inicializedGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes inicializedGUI wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = inicializedGUI_OutputFcn(hObject, eventdata, handles) 
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

% --- Executes on button press in paired_protocol.
function paired_protocol_Callback(hObject, eventdata, handles)
% hObject    handle to paired_protocol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%global variables that connect with other functions and GUIs
global basepath
global protocol
global raw_files_i;
global segmented_files_i;
global analysed_files_i;
global visualized_files_i;

protocol = 'paired_protocol'; % variable that selects the appropriate functions to 
% process data related to this protocol

%Call user input for working folder
basepath = uigetdir(pwd,'Select Data Folder');

% files with the required extensions in the selected path
files_mat = dir(fullfile(basepath, '*.mat'));
files_acq = dir(fullfile(basepath, '*.acq'));

%for loop that analyzes if there is any .acq file without corresponding .mat 
%raw file, for the paired-pulse protocol. If there is, proceeds with the conversion
for p = 1:length(files_acq)
    file_acq = files_acq(p).name;
    if ~any(contains({files_mat(1:length(files_mat)).name},file_acq(1:end-4))) ...
            && ~contains(file_acq,'INPUTOUTPUT')
        dataconvertion(basepath,'paired_pulse',file_acq);
    end
end

files_mat = dir(fullfile(basepath, '*.mat')); %atualizes the number of .mat 
%files (eventually altered by the previous for loop)

%creation of cell arrays that identifies and categorizes the initial files in 
%the selected path 
raw_files_i = {}; %includes .mat files with raw data (continuous data)
segmented_files_i = {}; %includes .mat files with segmented data (separation of trials)
analysed_files_i = {}; %includes .mat files with analyzed data (parameters of 
%interest extracted)
visualized_files_i = {}; %includes .mat files that were edited by the user in GUI

%for loop that reads .mat files' names and categorizes them, by the analysis if 
%certain strings representative of different stages of data processing are contained 
%in the file's name
for j = 1:length(files_mat)
    file_mat = files_mat(j).name; % mat file's name
    
    %adds file's name in the end+1 position of the respective cell
    %array
    if contains(file_mat(1:end-4), 'segmented') && ...
            ~contains(file_mat(1:end-4),'INPUTOUTPUT')
        segmented_files_i{end+1} = file_mat; 
    elseif contains(file_mat(1:end-4), 'analysed') ...
            && ~contains(file_mat(1:end-4),'INPUTOUTPUT')
        analysed_files_i{end+1} = file_mat;
    elseif contains(file_mat(1:end-4), 'visualized') && ...
            ~contains(file_mat(1:end-4),'INPUTOUTPUT')
        visualized_files_i{end+1} = file_mat;
    else
        for b = 1:length(files_acq)
            file_acq = files_acq(b).name;
            if contains(file_mat(1:end-4),file_acq(1:end-4)) && ...
                    ~contains(file_mat(1:end-4),'INPUTOUTPUT')
                %if it doesn't contain any of the
                %known strings, it must be a raw file's 
                %name (which has the same name as the 
                %.acq file)
                raw_files_i{end+1} = file_mat;
            end
        end
    end
end

guidata(hObject,handles); % Update handles structure

%default values of certain variables (later, they can be modified by the
%user - in 'specification' GUI):
min_amp_MEP = 0.05; % minimum MEP's amplitude (mV) for its quantification
latency_MEP = 18; % MEP's latency (ms) relative to the last TMS artefact 
%(used to define MEP's window for its detection)
quality_comb = 300; % quality factor of IIR comb filter (nondimensional)
curve_model = 0; % Type of I/O curve Modulation (1 for statistical modulation 
% and 0 for linear regression) 
artefact_detection = 0; % select if desired to see the 2nd peaks of the 1st TMS 
%artefact in each rectified EMG segment (0 for no, 1 for yes)
distribution_tms = 0; % select if desired to see the distribution of amplitude of the 
%detected  TMS artefacts (before and after optimization of the TMS
%detection algorithm): 0 for no, 1 for yes.
agremment_measures = 0; % select if desired to see the agreement between automatic 
%and manual quantification (0 for no, and 1 for yes)
comb = 1; %switch applicability of IIR comb filter (0 to not apply, 1 to apply)
constant_time_tms = 12; % time constant (ms) to construct TMS templates, related to an 
%estimate of a TMS pulse duration 
sampling_frequency = 10000;
choosedata(min_amp_MEP,latency_MEP, quality_comb, curve_model, ...
    artefact_detection, distribution_tms, agremment_measures, comb,...
    constant_time_tms,sampling_frequency); % call to next GUI 
uiresume(handles.figure1); %it proceeds to the delivering of the output, 
%and consequent deletion of present GUI

% --- Executes on button press in curve_protocol.
function curve_protocol_Callback(hObject, eventdata, handles)
% hObject    handle to curve_protocol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%global variables that connect with other functions and GUIs
global basepath;
global protocol;
global raw_files_i;
global segmented_files_i;
global analysed_files_i;
global visualized_files_i;

% variable that selects the appropriate functions to process data related to 
%this protocol
protocol = 'curve_protocol';

%Call user input for working folder
basepath = uigetdir(pwd,'Select Data Folder');

% files with the required extensions in the selected path
files_mat = dir(fullfile(basepath, '*.mat'));
files_acq = dir(fullfile(basepath, '*.acq'));

%for loop that analyzes if there is any .acq file without corresponding
%.mat raw file, for the protocol input-output. If there is, proceeds with
%the conversion 
for p = 1:length(files_acq)
    file_acq = files_acq(p).name;
    if ~any(contains({files_mat(1:length(files_mat)).name},file_acq(1:end-4))) ...
            && contains(file_acq,'INPUTOUTPUT')
        dataconvertion(basepath,'curve_protocol',file_acq);
    end
end

files_mat = dir(fullfile(basepath, '*.mat')); %atualizes the number of .mat 
%files (eventually altered by the previous for loop)

%creation of cell arrays that identifies and categorizes the initial files in the 
%selected path 
raw_files_i = {};
visualized_files_i = {};
analysed_files_i = {};
segmented_files_i = {};

%for loop that reads .mat files' names and categorizes them, by the analysis if 
%certain strings representative of different stages of data processing are contained 
%in the file's name
for j = 1:length(files_mat)
    file_mat = files_mat(j).name; %file's name
    
    %adds file's name in the end+1 position of the respective cell
    %array
    if contains(file_mat(1:end-4), 'segmented') && ...
            contains(file_mat(1:end-4),'INPUTOUTPUT')
        segmented_files_i{end+1} = file_mat; 
    elseif contains(file_mat(1:end-4), 'analysed') && ...
            contains(file_mat(1:end-4),'INPUTOUTPUT')
        analysed_files_i{end+1} = file_mat;
    elseif contains(file_mat(1:end-4), 'visualized') && ...
            contains(file_mat(1:end-4),'INPUTOUTPUT')
        visualized_files_i{end+1} = file_mat;
    else
        for b = 1:length(files_acq)
            file_acq = files_acq(b).name;
            if contains(file_mat(1:end-4),file_acq(1:end-4)) && contains(file_mat(1:end-4),'INPUTOUTPUT')
                %if it doesn't cointain any of the
                %known strings, it must be a raw file's 
                %name (which has the same name as the 
                %.acq file)
                raw_files_i{end+1} = file_mat;
            end
        end
    end
end

guidata(hObject,handles); % Update handles structure

%default values of certain variables (later, they can be modified by the
%user - in 'specification' GUI):
min_amp_MEP = 0.05; % minimum MEP's amplitude (mV) for its quantification
latency_MEP = 18; % MEP's latency (ms) relative to the last TMS artefact 
%(used to define MEP's window for its detection)
quality_comb = 300; % quality factor of IIR comb filter (nondimensional)
curve_model = 0; % Type of I/O curve Modulation (1 for statistical modulation 
% and 0 for linear regression) 
artefact_detection = 0; % select if desired to see the 2nd peaks of the 1st TMS 
%artefact in each rectified EMG segment (0 for no, 1 for yes)
distribution_tms = 0; % select if desired to see the distribution of amplitude of the 
%detected  TMS artefacts (before and after optimization of the TMS
%detection algorithm): 0 for no, 1 for yes.
agremment_measures = 0; % select if desired to see the agreement between automatic 
%and manual quantification (0 for no, and 1 for yes)
comb = 1; % switch applicability of IIR comb filter (0 to not apply, 1 for applying)
constant_time_tms = 12; % time constant (ms) to construct TMS templates, related to an 
%estimate of a TMS pulse duration  
sampling_frequency = 10000;
choosedata(min_amp_MEP,latency_MEP, quality_comb, curve_model, ...
    artefact_detection, distribution_tms, agremment_measures, comb,...
    constant_time_tms,sampling_frequency); % call to next GUI 
uiresume(handles.figure1); %it proceeds to the delivering of the output, and 
%consequent deletion of present GUI

% --- Executes on button press in cortical_protocol.
function cortical_protocol_Callback(hObject, eventdata, handles)
% hObject    handle to cortical_protocol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%warning message that pops whenever the user pretends to automatically
%analyze Silence Period Protocol 
warndlg('Warning: It hasn''t been yet automatize the analysis of these protocol');


% --- Executes on button press in instructions_button.
function instructions_button_Callback(hObject, eventdata, handles)
% hObject    handle to instructions_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
instructions();
