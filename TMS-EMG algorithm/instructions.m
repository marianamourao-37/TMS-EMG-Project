function varargout = instructions(varargin)
% INSTRUCTIONS MATLAB code for instructions.fig
%      INSTRUCTIONS, by itself, creates a new INSTRUCTIONS or raises the existing
%      singleton*.
%
%      H = INSTRUCTIONS returns the handle to a new INSTRUCTIONS or the handle to
%      the existing singleton*.
%
%      INSTRUCTIONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INSTRUCTIONS.M with the given input arguments.
%
%      INSTRUCTIONS('Property','Value',...) creates a new INSTRUCTIONS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before instructions_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to instructions_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help instructions

% Last Modified by GUIDE v2.5 30-Oct-2019 12:52:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @instructions_OpeningFcn, ...
                   'gui_OutputFcn',  @instructions_OutputFcn, ...
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


% --- Executes just before instructions is made visible.
function instructions_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to instructions (see VARARGIN)

% Choose default command line output for instructions
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes instructions wait for user response (see UIRESUME)
%uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = instructions_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Deletes figure of present GUI after uiresume 
varargout{1} = handles.output;
%delete(hObject);

% --- Executes on button press in excel_files.
function excel_files_Callback(hObject, eventdata, handles)
% hObject    handle to excel_files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
excel_files();

% --- Executes on button press in statistical_results.
function statistical_results_Callback(hObject, eventdata, handles)
% hObject    handle to statistical_results (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 statistical_bland_gui();
