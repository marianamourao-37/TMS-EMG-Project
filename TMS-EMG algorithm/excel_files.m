function varargout = excel_files(varargin)
% EXCEL_FILES MATLAB code for excel_files.fig
%      EXCEL_FILES, by itself, creates a new EXCEL_FILES or raises the existing
%      singleton*.
%
%      H = EXCEL_FILES returns the handle to a new EXCEL_FILES or the handle to
%      the existing singleton*.
%
%      EXCEL_FILES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EXCEL_FILES.M with the given input arguments.
%
%      EXCEL_FILES('Property','Value',...) creates a new EXCEL_FILES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before excel_files_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to excel_files_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help excel_files

% Last Modified by GUIDE v2.5 30-Oct-2019 11:12:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @excel_files_OpeningFcn, ...
                   'gui_OutputFcn',  @excel_files_OutputFcn, ...
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


% --- Executes just before excel_files is made visible.
function excel_files_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to excel_files (see VARARGIN)

axes(handles.axes1)
matlabImage = imread('pp_excel_example.jpg');
image(matlabImage)
axis off
axis image

axes(handles.axes2)
matlabImage = imread('io_excel_example.jpg');
image(matlabImage)
axis off
axis image
% Choose default command line output for excel_files
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes excel_files wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = excel_files_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
