function varargout = pairedpulse(varargin)
% PAIREDPULSE MATLAB code for pairedpulse.fig
%      PAIREDPULSE, by itself, creates a new PAIREDPULSE or raises the existing
%      singleton*.
%
%      H = PAIREDPULSE returns the handle to a new PAIREDPULSE or the handle to
%      the existing singleton*.
%
%      PAIREDPULSE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PAIREDPULSE.M with the given input arguments.
%
%      PAIREDPULSE('Property','Value',...) creates a new PAIREDPULSE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pairedpulse_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pairedpulse_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pairedpulse

%% 
% Function description: It refers to the GUI that allows the user to visualize
%and correct the analyzed data, for the paired-pulse protocol 

%input variables:
% - pathname: string with the path's name of the selected directory (that contains the files 
%to be analyzed) 
% - file_visualized: selected file, in 'choosedata' GUI, to be visualized
% - agremment_measures: select if desired to see the agreement between automatic 
%and manual quantification (0 for no, and 1 for yes)
%%

% Last Modified by GUIDE v2.5 22-Sep-2019 13:14:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pairedpulse_OpeningFcn, ...
                   'gui_OutputFcn',  @pairedpulse_OutputFcn, ...
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
end

% --- Executes just before pairedpulse is made visible.
function pairedpulse_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pairedpulse (see VARARGIN)

% enables funcionalities of toolbar and menubar for the present GUI
set(handles.figure1,'toolbar','figure');
set(handles.figure1,'menubar','figure');

handles.agremment_measures = varargin{3}; 

handles.sampling_frequency = varargin{4};
file_visualized = varargin{2}; %cell array that contains the file's name to be visualized 
file = file_visualized{1}; %extract from cell array the file's name to be visualized 
handles.file = file; % adds to handles structure the field 'file'
handles.pathname = varargin{1}; %file's path to be visualized 

EMGdata=load([handles.pathname,'\',file]); % loads the file to be visualized
a=1; %variable that identifies the segment to be visualized

%assigns to handles of type text the name that identifies the pacient 
if contains(file,'analysed') % visualize a previous analyzed file 
    set(handles.identification,'String',file(1:strfind(file,'_analysed')-1));
else % visualize a previuos visualized file
    set(handles.identification,'String',file(1:strfind(file,'_visualized')-1));
end
        
ISI_values = unique(EMGdata.trials.ISI_sec(:,1)); %get the ISI values of paired pulse 
%sequence

% writes the desired parameters in the present GUI 
texts = sprintf(['Mean amp. 0ms (mV): %0.3f ± %0.3f\n Mean amp. ', ...
    num2str(ISI_values(2)*10^3),'ms (mV): %0.3f ± %0.3f\n Mean amp. ', ...
    num2str(ISI_values(3)*10^3),'ms (mV): %0.3f ± %0.3f\n Normalized Mean amp. ', ...
    num2str(ISI_values(2)*10^3),'ms: %0.3f ± %0.3f\n Normalized Mean amp. ',...
    num2str(ISI_values(3)*10^3),'ms: %0.3f ± %0.3f'],...
    [EMGdata.statistics.mean_mep_amplitude(1,1),...
    EMGdata.statistics.sd_mep_amplitude(1,1)],...
    [EMGdata.statistics.mean_mep_amplitude(2,1),...
    EMGdata.statistics.sd_mep_amplitude(2,1)],...
    [EMGdata.statistics.mean_mep_amplitude(3,1),...
    EMGdata.statistics.sd_mep_amplitude(3,1)],...
    [EMGdata.results.(['pp_value_ISI_',num2str(ISI_values(2)*10^3),'ms']),...
    EMGdata.results.(['std_pp_value_ISI_',num2str(ISI_values(2)*10^3),'ms'])],...
    [EMGdata.results.(['pp_value_ISI_',num2str(ISI_values(3)*10^3),'ms']),...
    EMGdata.results.(['std_pp_value_ISI_',num2str(ISI_values(3)*10^3),'ms'])]);

% assigns to handles of type text the desired parameters 
set(handles.final_parameters,'String', texts);

% calls function that plots 1st segment of EMG data 
plot_figure(EMGdata,handles,a);

handles.EMGdata = EMGdata;
handles.a = a;
xlim1 = EMGdata.trials.artloc(1,1)/handles.sampling_frequency - 0.005;
xlim2 = xlim1 + EMGdata.trials.ISI_sec(1,1) + 0.1;
len_max = xlim2 - xlim1;

for y = 1:40
    xlim1 = EMGdata.trials.artloc(y,1)/handles.sampling_frequency - 0.005; % corrects the window for 
%visualization to start 5 ms before 1st TMS artefact 
    xlim2 = xlim1 + EMGdata.trials.ISI_sec(y,1) + 0.1; %correts the window for 
%visualization to end 95 ms after 2nd TMS artefact 
    len = xlim2 - xlim1;
    
    if len < len_max
        xlim2 = len_max+xlim1;
    else 
        len_max = len;
    end
    
    data(:,1) = EMGdata.trials.EMGfilt{y,1}(round(xlim1*handles.sampling_frequency):...
    round(xlim2*handles.sampling_frequency));

    if y == 1
        save('dataset_meps.mat','data');
    else 
        save('dataset_meps.mat','data','-append');
    end
   
end

% Choose default command line output 
handles.output = hObject;

% prompt to really close without saving
set(handles.figure1, 'CloseRequestFcn', @closeGUI); 

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes visualizeEMG wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = pairedpulse_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

function closeGUI(hObject, eventdata, handles)
selection = questdlg('You may have unsaved edits. Are you sure you want to close?', ...
    'Warning', ...
    'Yes','No','Yes');
switch selection
    case 'Yes'
        close(gcf)
    case 'No'
        return
end
%      end
end
%% PLOT FIGURE
function plot_figure(EMGdata,handles,a)

% buttons to display
set(handles.ch1_TMS_art,'visible','off');
set(handles.ch1_clearTMSart,'visible','off');   
set(handles.ch1_MEP,'visible','off');
set(handles.ch1_clearMEP,'visible','off');
set(handles.ch1_MEP,'visible','on');
set(handles.ch1_clearMEP,'visible','on');
set(handles.ch1_TMS_art,'visible','on');
set(handles.ch1_clearTMSart,'visible','on');

% delete annotations relative to MEP amplitude of the previous segment
% (maintains non-normalized and normalized mean MEP amplitudes) 
delete(findall(gcf,'type','annotation'));

% set colors for markers and buttons
MEP_color = [0.5843 0.8157 0.9882];
TMS_color = [1 0 0];

y=EMGdata.trials.EMGfilt{a,1}; % gets the 'a' segment of EMG data
x = linspace(0,(length(y)/handles.sampling_frequency),length(y)); % time scale for the given segment
axes(handles.graphic); % it gets the plot axis that was designed in the present GUI
plot(x,y,'k'); %plots EMG data in time domain, with black trace 
grid on; 

xlim1 = EMGdata.trials.artloc(a,1)/handles.sampling_frequency - 0.005; % corrects the window for 
%visualization to start 5 ms before 1st TMS artefact 

xlim2 = xlim1 + EMGdata.trials.ISI_sec(a,1) + 0.1; %correts the window for 
%visualization to end 95 ms after 2nd TMS artefact 

xlim([xlim1 xlim2]); % x-axis limits 

[ymin,~] = min(EMGdata.trials.EMGfilt{a,1}); %gets the maximum amplitude value of 
%the given segment 
[ymax,~] = max(EMGdata.trials.EMGfilt{a,1}); %gets the minimum amplitude value of 
%the given segment
ylims = [ymin-0.1 ymax+0.1]; %corrects the previous values to be 0.1 above of 
%ymax or 0.1 bellow of ymin
ylim([ymin-0.1 ymax+0.1]); % y-axis limits 

% Plot Title
if EMGdata.trials.ISI_sec(a,1) == 0 
    set(handles.graphic_title,'String',sprintf(['Trial #:', num2str(a),...
        ' - Baseline']));
else
    set(handles.graphic_title,'String',sprintf(['Trial #:', num2str(a),...
        ' - ISI = ',num2str(EMGdata.trials.ISI_sec(a,1)*10^3),' ms']));
end

%adds x label
xlabel('Time (s)', 'FontSize',10, 'FontName', 'Arial', 'FontWeight', 'bold');

% adds y label
ylabel('ch1(mV)','FontSize',10, 'FontName', 'Arial', 'FontWeight', 'bold');
    
% adds MEP colored region 
if EMGdata.trials.artloc(a,1)    
    MEP_on = (EMGdata.trials.MEP_onset(a,1)-1)/handles.sampling_frequency; %lower limit of MEP
    MEP_off = (EMGdata.trials.MEP_offset(a,1)-1)/handles.sampling_frequency; %higher limit of MEP
    
    %polygon that identifies the MEP window
    patch([MEP_on MEP_on MEP_off MEP_off], [ylims(1, 1) ylims(1, 2) ylims(1, 2) ...
        ylims(1, 1)], MEP_color); 
    alpha(.7);      
    
end

%adds lines relative to TMS artefacts
if EMGdata.trials.artloc(a,1) % if it was detected the TMS artefact 
    %line that identifies the 1st TMS artefact (in seconds) 
    line([(EMGdata.trials.artloc(a,1)-1)/handles.sampling_frequency,...
        (EMGdata.trials.artloc(a,1)-1)/handles.sampling_frequency], ...
        ylims(1, :) ,'Color',TMS_color,'LineStyle', '--','Marker','o')       
    
    lim = (EMGdata.trials.artloc(a,1)-1)/handles.sampling_frequency + EMGdata.trials.ISI_sec(a,1); 
    % calculates the position of the 2nd TMS artefact (in seconds) 
    
    % line that identifies 2nd TMS artefact (in seconds)
    line([lim lim], ylims(1, :) ,'Color',TMS_color,'LineStyle', '--','Marker','o')       
end

% set MEP amplitude relative to the given EMG segment   
set(handles.MEP_amplitude_trial,'String',sprintf('MEP amp. (mV): %0.3f',...
    EMGdata.trials.MEP_amplitude(a,1)));
end

function enter_trial_Callback(hObject, eventdata, handles)
% hObject    handle to enter_trial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of enter_trial as text
%        str2double(get(hObject,'String')) returns contents of enter_trial as a double

a = str2double(get(hObject,'String')); % atualizes the 'a'th segment to be visualized 
EMGdata = handles.EMGdata; % get EMGdata structure, that contains all segments 

plot_figure(EMGdata,handles,a)

%update handles fields
handles.EMGdata = EMGdata;
handles.a = a;
handles.output = hObject;

guidata(hObject, handles);  % Update handles structure
end

% --- Executes during object creation, after setting all properties.
function enter_trial_CreateFcn(hObject, eventdata, handles)
% hObject    handle to enter_trial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,...
        'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in back_button.
function back_button_Callback(hObject, eventdata, handles)
% hObject    handle to back_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a=handles.a; % gets the index value of the present segment 
EMGdata=handles.EMGdata; % get EMGdata structure, that contains all segments 

% decrement of the segment index to be visualized  
if a>1
    a=a-1;
end

plot_figure(EMGdata,handles,a); % calls function that plots the 'a'th segment 
%of EMG data 

%update handles variables
handles.EMGdata = EMGdata;
handles.a = a;
handles.output = hObject;

guidata(hObject, handles);  % Update handles structure
end

% --- Executes on button press in forward_button.
function forward_button_Callback(hObject, eventdata, handles)
% hObject    handle to forward_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a=handles.a; % gets the index value of the present segment 
EMGdata=handles.EMGdata; % get EMGdata structure, that contains all segments 

% increment of the segment index to be visualized
if a < height(EMGdata.trials)
    a=a+1;
end

plot_figure(EMGdata,handles,a) % calls function that plots 'a'th segment of EMG data 

%update handles variables
handles.EMGdata = EMGdata;
handles.a = a;
handles.output = hObject;

guidata(hObject, handles);  % Update handles structure
end

% --- Executes on button press in export_button.
function export_button_Callback(hObject, eventdata, handles)
% hObject    handle to export_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

EMGdata = handles.EMGdata; % get EMGdata structure, that contains all segments 
ISI_values = unique(EMGdata.trials.ISI_sec(:,1)); %isi values contained in the randomized key 
file = handles.file; % name of the selected file to be visualized 

if contains(file,'analysed') % visualize a previuos analyzed file
    file = file(1:strfind(file,'_analysed')-1);
else % visualize a previuos visualized file
    file = file(1:strfind(file,'_visualized')-1);
end

filename = [handles.pathname,'\reviewdata_',file,'.xlsx']; %name of the Excel file to be created 

total_data = table();
total_data.ISI_sec = EMGdata.trials.ISI_sec(:,1); %randomized key 
total_data.mep_amp_mV = EMGdata.trials.MEP_amplitude(:,1); %MEP amplitudes 
writetable(total_data,filename,'Sheet',1,'Range','A1'); %adds information to 1st sheet 

analysed_statistics = table(); 
analysed_statistics.ISI_sec = ISI_values(:,1);
analysed_statistics.mean_mep_mV = EMGdata.statistics.mean_mep_amplitude(:,1); %mean of MEP 
%amplitudes, for each isi 
analysed_statistics.sd_mep_mV = EMGdata.statistics.sd_mep_amplitude(:,1); %standard deviation 
%of MEP amplitudes, for each isi 
writetable(analysed_statistics,filename,'Sheet',2,'Range','A1'); %adds information to 2nd sheet 

normalized_analysis = table();
normalized_analysis.amp_normalized(1,1) = EMGdata.results.(['pp_value_ISI_',...
    num2str(ISI_values(2)*10^3),'ms']);
normalized_analysis.amp_normalized(2,1) = EMGdata.results.(['pp_value_ISI_',...
    num2str(ISI_values(3)*10^3),'ms']);
normalized_analysis.sd_normalized(1,1)= EMGdata.results.(['std_pp_value_ISI_',...
    num2str(ISI_values(2)*10^3),'ms']);
normalized_analysis.sd_normalized(2,1)= EMGdata.results.(['std_pp_value_ISI_',...
    num2str(ISI_values(3)*10^3),'ms']);
writetable(normalized_analysis,filename,'Sheet',2,'Range','E2'); %adds information to 2nd sheet, 
%starting in E2 cell
end

% --- Executes on button press in ch1_clearTMSart.
function ch1_clearTMSart_Callback(hObject, eventdata, handles)
% hObject    handle to ch1_clearTMSart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a=handles.a; % gets the index value of the present segment 
EMGdata=handles.EMGdata; % get EMGdata structure, that contains all segments 

EMGdata.trials.artloc(a,1)=0; % sets to zero the position of the 1st TMS artefact 
%of the given segment 

handles.EMGdata = EMGdata; %atualizes EMGdata structure 
plot_figure(EMGdata,handles,a); % calls function that plots 'a'th segment 
%(same segment) of EMG data 
guidata(hObject, handles); % Update handles structure
end

% --- Executes on button press in ch1_TMS_art.
function ch1_TMS_art_Callback(hObject, eventdata, handles)
% hObject    handle to ch1_TMS_art (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a=handles.a; % gets the index value of the present segment 
EMGdata=handles.EMGdata; % gets EMGdata structure, that contains all segments 

axh = gca;  
axh.Toolbar = matlab.ui.controls.AxesToolbar(); 

% manually selects 1st TMS artefact 
x_new = ginput(1);
EMGdata.trials.artloc(a,1) = x_new(1)*handles.sampling_frequency; % saves in EMGdata structure the 
%position of 1st TMS artefact (in samples)

handles.EMGdata = EMGdata; %atualizes EMGdata structure 
plot_figure(EMGdata,handles,a); % calls function that plots 'a'th segment 
%(same segment) of EMG data 
guidata(hObject, handles); % Update handles structure
end

% --- Executes on button press in ch1_MEP.
function ch1_MEP_Callback(hObject, eventdata, handles)
% hObject    handle to ch1_MEP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a=handles.a; % gets the index value of the present segment 
EMGdata=handles.EMGdata; % get EMGdata structure, that contains all segments 

axh = gca;  
axh.Toolbar = matlab.ui.controls.AxesToolbar(); 

% manually selects onset and offset limits of MEP
x_new = ginput(2);
EMGdata.trials.MEP_onset(a,1) = x_new(1)*handles.sampling_frequency; % saves in EMGdata structure 
%the position of onset MEP limit in samples

EMGdata.trials.MEP_offset(a,1) = x_new(2)*handles.sampling_frequency; % saves in EMGdata structure 
%the position of offset MEP limit in samples

MEPchannel=EMGdata.trials.EMGfilt{a,1}; %gets the segment data
MEPsearchrange = MEPchannel(round(x_new(1) * handles.sampling_frequency):...
    round(x_new(2) * handles.sampling_frequency)); 
% identifies the MEP data region 

%get the maximum and minimum values in MEP data
[max_MEP_value,~] = max(MEPsearchrange);
[min_MEP_value,~] = min(MEPsearchrange);

EMGdata.trials.MEP_amplitude(a,1)=max_MEP_value-min_MEP_value; %MEP amplitude

handles.EMGdata = EMGdata; %atualizes EMGdata structure 
plot_figure(EMGdata,handles,a); % calls function that plots 'a'th segment 
%(same segment) of EMG data 
guidata(hObject, handles); % Update handles structure
end
% --- Executes on button press in ch1_clearMEP.
function ch1_clearMEP_Callback(hObject, eventdata, handles)
% hObject    handle to ch1_clearMEP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a=handles.a; % gets the index value of the present segment 
EMGdata=handles.EMGdata; % get EMGdata structure, that contains all segments 


EMGdata.trials.MEP_onset(a,1)=0; % sets to zero MEP's lower limit 
EMGdata.trials.MEP_offset(a,1)=0; % sets to zero MEP's upper limit
EMGdata.trials.MEP_amplitude(a,1)=0; % sets to zero MEP's amplitude 

handles.EMGdata = EMGdata; %atualizes EMGdata structure 
plot_figure(EMGdata,handles,a); % calls function that plots 'a'th segment 
%(same segment) of EMG data 
guidata(hObject, handles); % Update handles structure
end

% --- Executes on button press in save_button.
function save_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global new_visualized_file; %global variable created in 'choosedata.mat' GUI
global new_edited_visualized_file; %global variable created in 'choosedata.mat' GUI
           %global pathname;
           
EMGdata=handles.EMGdata; % get EMGdata structure, that contains all segments 

ISI_values = unique(EMGdata.trials.ISI_sec(:,1)); % get the ISI values applied 
%during a certain paired pulse sequence acquisition

EMGdata.statistics = table();
EMGdata.statistics.ISI_sec(:,1) = ISI_values;

%calculate the mean MEP amplitudes for each ISI value 
for u= 1:length(ISI_values)
    index = find(EMGdata.trials.ISI_sec(:,1) == ISI_values(u));
    mep_values = EMGdata.trials.MEP_amplitude(index,1); %MEP amplitudes relative to the uth isi value
    mep_values(mep_values == 0) = NaN; %doesn't account MEP amplitudes equal to 0 
    EMGdata.statistics.mean_mep_amplitude(u,1) = mean(mep_values,'omitnan'); %calculate mean
    EMGdata.statistics.sd_mep_amplitude(u,1) = std(mep_values,'omitnan')'; %calculate standard deviation
    EMGdata.statistics.cv_mep_amplitude(u,1) = ...
        EMGdata.statistics.sd_mep_amplitude(u,1)/...
        EMGdata.statistics.mean_mep_amplitude(u,1); % calculate variation coefficient
end

%calculates the normalized mean MEP amplitudes 
for z = 2:length(ISI_values)
    EMGdata.results.(['pp_value_ISI_',num2str(ISI_values(z)*10^3),'ms']) = ...
        EMGdata.statistics.mean_mep_amplitude(z,1)/EMGdata.statistics.mean_mep_amplitude(1,1);
    
    %standard deviation calculated from the formula of propagation of errors 
    EMGdata.results.(['std_pp_value_ISI_',num2str(ISI_values(z)*10^3),'ms']) = ...
        (EMGdata.results.(['pp_value_ISI_',num2str(ISI_values(z)*10^3),'ms'])*...
        ((EMGdata.statistics.sd_mep_amplitude(z,1)/EMGdata.statistics.mean_mep_amplitude(z,1))^2+...
        (EMGdata.statistics.sd_mep_amplitude(1,1)/...
        EMGdata.statistics.mean_mep_amplitude(1,1))^2))^(1/2);
end
 
%eliminates the previous parameters (non-normalized and normalized mean and std MEP amplitudes) 
set(handles.final_parameters,'String', '');

% writes the desired parameters 
texts = sprintf(['Mean amp. 0ms (mV): %0.3f ± %0.3f\n Mean amp. ', ...
    num2str(ISI_values(2)*10^3),'ms (mV): %0.3f ± %0.3f\n Mean amp. ', ...
    num2str(ISI_values(3)*10^3),'ms (mV): %0.3f ± %0.3f\n Normalized Mean amp. ', ...
    num2str(ISI_values(2)*10^3),'ms: %0.3f ± %0.3f\n Normalized Mean amp. ', ...
    num2str(ISI_values(3)*10^3),'ms: %0.3f ± %0.3f'],...
    [EMGdata.statistics.mean_mep_amplitude(1,1),...
    EMGdata.statistics.sd_mep_amplitude(1,1)],...
    [EMGdata.statistics.mean_mep_amplitude(2,1),...
    EMGdata.statistics.sd_mep_amplitude(2,1)],...
    [EMGdata.statistics.mean_mep_amplitude(3,1),...
    EMGdata.statistics.sd_mep_amplitude(3,1)],...
    [EMGdata.results.(['pp_value_ISI_',num2str(ISI_values(2)*10^3),'ms']),...
    EMGdata.results.(['std_pp_value_ISI_',num2str(ISI_values(2)*10^3),'ms'])],...
    [EMGdata.results.(['pp_value_ISI_',num2str(ISI_values(3)*10^3),'ms']),...
    EMGdata.results.(['std_pp_value_ISI_',num2str(ISI_values(3)*10^3),'ms'])]);

% assigns to handles of type text the desired parameters 
set(handles.final_parameters,'String', texts);

guidata(hObject, handles); % Update handles structure
display('saved');

% opens dialog for user to decide, after pressing save button, if it really
% wants to exit or not
selection = questdlg('Do you want to close visualization of data?', ...
      'Warning', ...
      'Yes','No','Yes');
         switch selection
         case 'Yes'
             
           file=handles.file; %gets file's name 
           
           trials = EMGdata.trials;
           results = EMGdata.results;
           statistics = EMGdata.statistics;

           if contains(file,'analysed') %if the visualized file was a 
               %previous analyzed file (1st edition) 
    
               new_visualized_file=[file(1:strfind(file,'_analysed')-1),'_visualized.mat']; % adds to global 
               % variable (cell array) the name of the new visualized file
    
               save([handles.pathname,'\',new_visualized_file],'trials','results','statistics');
               % save structures in .mat file, in the previously selected path

           else %if the visualized file was a previous visualized file (it isn't the 
               %1st edition) 
               
               new_edited_visualized_file=file; % adds to global variable 
               %(cell array) the name of the new edited visualized file
               
               save([handles.pathname,'\',new_edited_visualized_file],'trials','results',...
                   'statistics');
               % save structures in .mat file, in the previously selected path
           end
           
           % Get default command line output from handles structure, and deletes 
           %figure of present GUI after uiresume 
           varargout{1} = handles.output;
           delete(handles.figure1); 
           
         case 'No' %if user selects that it doesn't pretend to exit, it returns to GUI
           return 
         end
end

% --- Executes on button press in statistics.
function statistics_Callback(hObject, eventdata, handles)
% hObject    handle to statistics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global pathname

EMGdata=handles.EMGdata; % get EMGdata structure, that contains all segments 
algorithm_data = handles.file; % name of the selected file to be visualized 

if contains(algorithm_data,'analysed') 
    file = algorithm_data(1:strfind(algorithm_data,'_analysed')-1);
else 
    file = algorithm_data(1:strfind(algorithm_data,'_visualized')-1);
end
    
separator = strfind(algorithm_data,'_');
name_pacient = algorithm_data(1:separator(2)); %patient name is relative to the second _

% files with the required extensions in the selected path
files_mat  = dir(fullfile(handles.pathname, '*.mat'));
files_name = extractfield(files_mat,'name');

%indexes, in 'files_name' cell array, of the files relative to the wanted patient 
count_pacient = find(contains(files_name, name_pacient)~=0);

vMean = []; %mean values 
vXvalues  = []; %isi values 
vStd = []; %std values 

for j = 1:numel(count_pacient) 
    if ~strcmp(files_name{count_pacient(j)}(1:end-4),file) %file relative to the wanted patient but a 
        %different paired-pulse sequence of the one being analyzed 
        if contains(files_name{count_pacient(j)},'analysed') %file with information after being analyzed  
            data = load([handles.pathname,'\',files_name{count_pacient(j)}]); %loads file
            ISI_values = unique(data.trials.ISI_sec(:,1)); %get the ISI values of paired pulse sequence
            vMean(end+1) = data.results.(['pp_value_ISI_',num2str(ISI_values(2)*10^3),'ms']);
            vXvalues(end+1) = ISI_values(2)*10^3;
            vStd(end+1) = data.results.(['std_pp_value_ISI_',num2str(ISI_values(2)*10^3),'ms']);
            vMean(end+1) = data.results.(['pp_value_ISI_',num2str(ISI_values(3)*10^3),'ms']);
            vXvalues(end+1) = ISI_values(3)*10^3;
            vStd(end+1) = data.results.(['std_pp_value_ISI_',num2str(ISI_values(3)*10^3),'ms']);
        end
    else %file relative to the wanted patient and the same paired-pulse sequence of the one being 
        %analyzed 
        ISI_values = unique(EMGdata.trials.ISI_sec(:,1));
        vMean(end+1) = EMGdata.results.(['pp_value_ISI_',num2str(ISI_values(2)*10^3),'ms']);
        vXvalues(end+1) = ISI_values(2)*10^3;
        vStd(end+1) = EMGdata.results.(['std_pp_value_ISI_',num2str(ISI_values(2)*10^3),'ms']);
        vMean(end+1) = EMGdata.results.(['pp_value_ISI_',num2str(ISI_values(3)*10^3),'ms']);
        vXvalues(end+1) = ISI_values(3)*10^3;
        vStd(end+1) = EMGdata.results.(['std_pp_value_ISI_',num2str(ISI_values(3)*10^3),'ms']);
    end
end

[vXvalues, sortorder] = sort(vXvalues); % sort ISIs in crescent order 
vMean = vMean(sortorder); % sort in the same order as it was in vXvalues
vStd = vStd(sortorder); % sort in the same order as it was in vXvalues
figure(1);
plot(vXvalues,vMean,'bo--','MarkerSize',8); hold on;
errorbar(vXvalues,vMean,vStd,'LineStyle','none'); %errobar as being the standard deviation of 
%MEP amplitudes (calculated with propagation error formula) for each isi value 
title('MEP Amplitude in Paired-Pulse paradigm');
ylabel('Normalized Amplitude (% of baseline)');
xlabel('ISI (ms)');
xlim([0 110]);

%% Bladan-Altman analysis
if handles.agremment_measures %if desired to see the agreement between automatic 
%and manual quantification 

    % files with the required extensions in the selected path
    excel_files = dir(fullfile(handles.pathname, '*.xlsx'));
    excel_names = extractfield(excel_files,'name');

    if ~any(contains(excel_names,[file,'.xlsx']))
        %selected path doesn't contain the corresponding Excel file for making Bladan-Altman analysis
        warndlg(['Warning: Excel file needed doesn''t exist. Please load an excel file',...
            ' with name ', file,'.xlsx, that contains data analysed by hand (human_analysis).']);
    else %selected path contains the corresponding Excel file for making Bladan-Altman analysis
        human_data = xlsread([handles.pathname,'\',file]); %reads corresponding Excel file of the one 
        % being visualized  
        methods_amplitudes.alghorithm_mep = []; %automatic MEP quantification
        methods_amplitudes.human_mep = []; % manual MEP quantification
        
        ISI_values = unique(human_data(:,1));
        for i = 2:length(ISI_values) 
            index = find(human_data(:,1) == ISI_values(i)); %indexes of the ith isi value 
            %on the randomized key (excluding isi = 0) 
            
            methods_amplitudes.alghorithm_mep(:,end+1) = EMGdata.trials.MEP_amplitude(index,1);
            %mep amplitudes automatically quantified, for the ith isi value (excluding isi = 0) 
            
            methods_amplitudes.human_mep(:,end+1) = human_data(index,2);
            %mep amplitudes manually quantified, for the ith isi value (excluding isi = 0) 
        end 

        index = find(human_data(:,1) == 0); %indexes of isi = 0 on the randomized key 
        baseline_alghorithm = EMGdata.trials.MEP_amplitude(index,1); %mep amplitudes automatically 
        %quantified, for isi = 0 
        baseline_human = human_data(index,2); %mep amplitudes manually quantified, for isi = 0 

        if length(baseline_human) < size(methods_amplitudes.alghorithm_mep,1)
    
            baseline_alghorithm = [baseline_alghorithm; ...
        NaN([size(methods_amplitudes.alghorithm_mep,1)-length(baseline_human),1])];
            % making the baseline_alghorithm vector the same size as the 
            % methods_amplitudes.alghorithm_mep vector (i.e, adds 5 NaN elements, that won't be 
            % account in Bladan-Altman analysis)
            
            baseline_human = [baseline_human; ...
        NaN([size(methods_amplitudes.alghorithm_mep,1)-length(baseline_human),1])];
            % making the baseline_human vector the same size as the 
            % methods_amplitudes.human_mep vector (i.e, adds 5 NaN elements, that won't be 
            % account in Bladan-Altman analysis)
        end
        methods_amplitudes.alghorithm_mep = [baseline_alghorithm,methods_amplitudes.alghorithm_mep];
        % automatic MEP quantification 
        
        methods_amplitudes.human_mep = [baseline_human,methods_amplitudes.human_mep];
        % manual MEP quantification 

%% load data
        territories = {'ISI=0ms',['ISI=',num2str(ISI_values(2)*10^3),'ms'],...
    ['ISI=',num2str(ISI_values(3)*10^3),'ms']};

        data1= methods_amplitudes.alghorithm_mep(:,:); % automatic MEP quantification 
        data2= methods_amplitudes.human_mep(:,:); % manual MEP quantification 

        % BA plot paramters
        tit = 'MEP Amplitude Accuracy (Bland-Altman)'; % figure title
        tit1 = 'MEP Amplitude Correlation'; % figure title
        gnames = {territories}; % names of groups in data 
        label = {'Algorithm Measure','Human Measure','mV'}; % Names of data sets
        corrinfo = {'n','SSE','r2','eq','P'}; % stats to display of correlation scatter plot
        limits = 'auto'; % how to set the axes limits
        colors = 'brgmkc'; % character codes

%% Example 2 - using non-Gaussian data and one figure with 2 analyses

        [~, fig, ~] = correlationPlot(data1, data2,label,tit1,gnames,'corrInfo',...
    corrinfo,'colors',colors);

        [~, fig, statsStruct2] = BlandAltman(data1, data2,label,...
    [tit ' (using non-parametric stats)'],gnames,'corrInfo',corrinfo,...
    'axesLimits',limits,'colors',colors,'baYLimMode','square','showFitCI',' on',...
    'baStatsMode','non-parametric');
        disp('Statistical results (using non-parametric stats):')
        disp(statsStruct2)
    end
else
    warndlg(['If you desired to compare the', ...
        ' automatic quantification versus manual, please verify the settings on ''Change',...
        ' algorithm specifications'' in the previous interface (''choosedata'' - where', ...
        ' the files are organized)']);
end
end 
