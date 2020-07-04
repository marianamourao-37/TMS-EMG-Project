function varargout = choosedata(varargin)
% CHOOSEDATA MATLAB code for choosedata.fig
%      CHOOSEDATA, by itself, creates a new CHOOSEDATA or raises the existing
%      singleton*.
%
%      H = CHOOSEDATA returns the handle to a new CHOOSEDATA or the handle to
%      the existing singleton*.
%
%      CHOOSEDATA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CHOOSEDATA.M with the given input arguments.
%
%      CHOOSEDATA('Property','Value',...) creates a new CHOOSEDATA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before choosedata_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to choosedata_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help choosedata

%% 
% Function description: It refers to the GUI that allows the user to navigate 
%throughout the different stage of processing:
% 1. segmentation of the raw continuous EMG data 
% 2. extraction of wanted parameters 
% 3. visualization of the results

% It has as input variables the values of different parameters: 
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

% Last Modified by GUIDE v2.5 24-Jul-2019 22:30:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @choosedata_OpeningFcn, ...
                   'gui_OutputFcn',  @choosedata_OutputFcn, ...
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


% --- Executes just before choosedata is made visible.
function choosedata_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to choosedata (see VARARGIN)

% Values of certain variables (they can be modified by the user - in 'specification' 
%GUI):
handles.min_amp_MEP = varargin{1}; 
handles.latency_MEP = varargin{2}; 
handles.quality_comb = varargin{3}; 
handles.curve_model = varargin{4}; 
handles.artefact_detection = varargin{5}; 
handles.distribution_tms = varargin{6}; 
handles.agremment_measures = varargin{7}; 
handles.comb = varargin{8}; 
handles.constant_time_tms = varargin{9};
handles.sampling_frequency = varargin{10};
% calls global variables created in 'inicializedGUI.mat' file. It contains in 
% cell arrays the names of .mat files initially present in the selected path
global raw_files_i;
global segmented_files_i;
global analysed_files_i;
global visualized_files_i;
global basepath;

%assigns to handles of type text, the path's name selected (string). Only useful 
%for visualization and better perception of which patient you're analyzing
set(handles.path_file,'String',basepath);

% set the maximum number of permitted files' selections in each of the
% list boxes. In the stages of data processing corresponding to 
% raw data --> segmented data and segmented data --> analyzed data it is permitted 
% to select all the files. In the stages of data processing corresponding to 
% analyzed data --> visualized data and visualized data --> edited visualized data it 
% is only permitted the selection of 1 file at each time. 
set(handles.raw_list,'max',length(raw_files_i));
set(handles.segmented_list,'max',length(segmented_files_i));
set(handles.analysed_list,'max',1);
set(handles.visualized_list,'max',1);

% sets files's names to the corresponding list box (representative of a certain
% stage in data processing), for later selection
set(handles.raw_list,'String',raw_files_i);
set(handles.segmented_list,'String',segmented_files_i);
set(handles.analysed_list,'String',analysed_files_i);
set(handles.visualized_list,'String',visualized_files_i);

% Choose default command line output 
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes choosedata wait for user response (see UIRESUME)
% uiwait(handles.figure2);

% --- Outputs from this function are returned to the command line.
function varargout = choosedata_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
%delete(handles.figure2);

% --- Executes on button press in get_segmented.
function get_segmented_Callback(hObject, eventdata, handles)
% hObject    handle to get_segmented (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global basepath; %global variable related to the path's name of selected directory
global new_segmented_files; %creation of new global variable, cell array containing 
%strings which identifies the new segmented files
global protocol; %global variable created in the initial GUI ('inicializedGUI.mat'), 
% that contains a string which identifies the respective protocol 
num_mat = get(handles.raw_list,'Value'); % get the index of the selected files in the 
% listbox structure
string_mat = cellstr(get(handles.raw_list,'String')); 
files_mat = string_mat(num_mat); %get the selected files identified by their index in 
% the listbox

%for the respective protocol, calls the appropriate function that permits
%extraction of the desired parameters
if strcmp(protocol,'paired_protocol')
    data_convertion_segmentation(basepath,handles.quality_comb,handles.artefact_detection,...
    handles.distribution_tms,handles.comb,files_mat,handles.sampling_frequency)
end
if strcmp(protocol,'curve_protocol')
    data_convertion_segmentation_curve(basepath,handles.quality_comb,handles.artefact_detection,...
    handles.distribution_tms,handles.comb,files_mat,handles.sampling_frequency)
end

existed_segmentation = get(handles.segmented_list,'String'); %get the pre-existing segmented files 
% in listbox 
list_segmented = [existed_segmentation{:},new_segmented_files];
files_list_segmented = unique(list_segmented,'last');
set(handles.segmented_list,'String', files_list_segmented); %set the respective 
% listbox with the pre-existing segmented files and the new ones

disp('saved')

set(handles.segmented_list,'max',length(files_list_segmented)); %atualizes the maximum number of 
% files permitted to be selected in the respective listbox

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in get_analysed.
function get_analysed_Callback(hObject, eventdata, handles)
% hObject    handle to get_analysed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global basepath; %global variable related to the path's name of the selected directory
global new_analysed_files; %creation of new global variable, cell array containing  strings which 
% identify the new analyzed files
global protocol; %global variable created in the initial GUI ('inicializedGUI.mat'), 
% that contains a string which identifies the respective protocol 

num_mat = get(handles.segmented_list,'Value'); % get the index of the selected files in the listbox 
%structure
string_mat = cellstr(get(handles.segmented_list,'String'));
files_segmented = string_mat(num_mat); %get the selected files identified by their index in the listbox

%for the respective protocol, calls the appropriate function that permits
%extraction of the desired parameters
if strcmp(protocol,'paired_protocol')
    final_parameters(basepath,handles.min_amp_MEP, handles.latency_MEP, files_segmented,...
        handles.sampling_frequency);
end
if strcmp(protocol,'curve_protocol')
    final_parameters_curve(basepath,handles.min_amp_MEP, handles.latency_MEP,handles.curve_model,...
        files_segmented,handles.sampling_frequency);
end

existed_analysed = get(handles.analysed_list,'String'); %get the pre-existing analyzed files in listbox 

list_analysed = [existed_analysed{:},new_analysed_files];
files_list_analysed = unique(list_analysed,'last');
set(handles.analysed_list,'String', files_list_analysed); %set the respective 
% listbox with the pre-existing analyzed files and the new ones

disp('saved')

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in get_visualized.
function get_visualized_Callback(hObject, eventdata, handles)
% hObject    handle to get_visualized (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global basepath; %global variable related to the path's name of the selected directory
global new_visualized_file; %creation of new global variable, cell array containing strings 
%which identify the new visualized files
global protocol; %global variable created in the initial GUI ('inicializedGUI.mat'), 
% that contains a string which identifies the respective protocol 

num_mat = get(handles.analysed_list,'Value'); % get the index of the selected files in the listbox 
%structure
string_mat = cellstr(get(handles.analysed_list,'String'));
file_visualized = string_mat(num_mat); %get the selected files identified by their index in the listbox

%for the respective protocol, calls the appropriate GUI that permits
%visualization and correction of analyzed data
if strcmp(protocol,'paired_protocol')
    pairedpulse(basepath,file_visualized,handles.agremment_measures,handles.sampling_frequency);
elseif strcmp(protocol,'curve_protocol')
    curve(basepath,file_visualized,handles.agremment_measures,handles.sampling_frequency);
end

existed_visualized = get(handles.visualized_list,'String'); %get the pre-existing visualized files in 
%listbox 

list_visualized = [existed_visualized;new_visualized_file];
files_list_visualized = unique(list_visualized,'last');
set(handles.visualized_list,'String', files_list_visualized);  %set the respective 
% listbox with the pre-existing visualized files and the new ones

disp('saved')

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in get_edited_visualized.
function get_edited_visualized_Callback(hObject, eventdata, handles)
% hObject    handle to get_edited_visualized (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global basepath; %global variable related to the path's name of the selected directory
global new_edited_visualized_file; %creation of new global variable, cell array containing strings 
%which identify the new edited visualized files
global protocol; %global variable created in the initial GUI ('inicializedGUI.mat'), 
% that contains a string that identifies the respective protocol 

string_mat = cellstr(get(handles.visualized_list,'String'));
file_edited_visualized = string_mat(get(handles.visualized_list,'Value')); %get the selected files 
% identified by their index in the listbox

existed_visualized = get(handles.visualized_list,'String'); %get the pre-existing visualized files in 
%listbox 

%for the respective protocol, calls the appropriate GUI that permits
%visualization and correction of analyzed data
if strcmp(protocol,'paired_protocol')
    pairedpulse(basepath,file_edited_visualized,handles.agremment_measures,handles.sampling_frequency);
end
if strcmp(protocol,'curve_protocol')
    curve(basepath,file_edited_visualized,handles.agremment_measures,handles.sampling_frequency);
end

list_edited_visualized = [existed_visualized;new_edited_visualized_file];
files_list_edited_visualized = unique(list_edited_visualized,'last');
set(handles.visualized_list,'String', files_list_edited_visualized);  %sets the 
%respective listbox with the pre-existing visualized files (that weren't selected), the edited visualized 
%files (that weren't selected) and the new edited visualized files

disp('saved')

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function raw_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to raw_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function visualized_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to visualized_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function analysed_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to analysed_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function segmented_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to segmented_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes after pressing change_specifications button 
function change_specifications_Callback(hObject, eventdata, handles)
% hObject    handle to change_specifications (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% calls to next GUI, that allows changing the value of certain parameters
% (input of the present GUI)
specifications(handles.min_amp_MEP,handles.latency_MEP, handles.quality_comb, handles.curve_model, ...
    handles.artefact_detection, handles.distribution_tms, handles.agremment_measures, handles.comb,...
    handles.constant_time_tms,handles.sampling_frequency);

% --- Executes on selection change in segmented_list.
function segmented_list_Callback(hObject, eventdata, handles)
% hObject    handle to segmented_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns segmented_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from segmented_list

% --- Executes on selection change in analysed_list.
function analysed_list_Callback(hObject, eventdata, handles)
% hObject    handle to analysed_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns analysed_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from analysed_list

% --- Executes on selection change in visualized_list.
function visualized_list_Callback(hObject, eventdata, handles)
% hObject    handle to visualized_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns visualized_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from visualized_list

% --- Executes on selection change in raw_list.
function raw_list_Callback(hObject, eventdata, handles)
% hObject    handle to raw_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns raw_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from raw_list