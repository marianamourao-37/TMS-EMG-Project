function varargout = curve(varargin)
% CURVE MATLAB code for curve.fig
%      CURVE, by itself, creates a new CURVE or raises the existing
%      singleton*.
%
%      H = CURVE returns the handle to a new CURVE or the handle to
%      the existing singleton*.
%
%      CURVE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CURVE.M with the given input arguments.
%
%      CURVE('Property','Value',...) creates a new CURVE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before curve_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to curve_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help curve

%% 
% Function description: It refers to the GUI that allows the user to visualize
%and correct the analyzed data, for the I/O protocol 

%input variables:
% - pathname: string with the path's name of the selected directory (that contains the files 
%to be analyzed) 
% - file_visualized: selected file, in 'choosedata' GUI, to be visualized
% - agremment_measures: select if desired to see the agreement between automatic 
%and manual quantification (0 for no, and 1 for yes)
%%

% Last Modified by GUIDE v2.5 22-Sep-2019 13:14:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @curve_OpeningFcn, ...
                   'gui_OutputFcn',  @curve_OutputFcn, ...
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

% --- Executes just before curve is made visible.
function curve_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to curve (see VARARGIN)

% enables funcionalities of toolbar and menubar for the present GUI
set(handles.figure1,'toolbar','figure');
set(handles.figure1,'menubar','figure');

%global pathname;

handles.agremment_measures = varargin{3};

handles.sampling_frequency = varargin{4};
file_visualized = varargin{2}; %cell array that contains the file's name to be 
%visualized 
handles.pathname = varargin{1}; %file's path to be visualized 
file = file_visualized{1}; %extract from cell array the file's name to be visualized 
handles.file = file; % adds to handles structure the field file, that as associate 
%the file's name to be visualized

EMGdata=load([handles.pathname,'\',file]); % loads the file to be visualized
a=1; %variable that identifies the segment to be visualized

% writes the desired parameters (mean MEP amplitudes for each intensity, as well as 
%the intensity of half maximum and slope of linear region)
texts = sprintf(['Slope of linear portion (mV/stim): %0.3f\n stim intensity of half max = %0.3f\n mean amp. ', num2str(EMGdata.statistics.intensities(1,1)),':%0.3f ± %0.3f (mV)\n mean amp. ', num2str(EMGdata.statistics.intensities(2,1)),': %0.3f ± %0.3f (mV)\n mean amp. ', num2str(EMGdata.statistics.intensities(3,1)),': %0.3f ± %0.3f (mV)\n mean amp. ', num2str(EMGdata.statistics.intensities(4,1)),': %0.3f ± %0.3f (mV)\n mean amp. ', num2str(EMGdata.statistics.intensities(5,1)),': %0.3f ± %0.3f (mV)\n mean amp. ', num2str(EMGdata.statistics.intensities(6,1)),': %0.3f ± %0.3f (mV)'],...
    EMGdata.results.slope2,...
    EMGdata.results.intensity_halfmax,...
    [EMGdata.statistics.mean_mep_amplitude(1,1),...
    EMGdata.statistics.sd_mep_amplitude(1,1)],...
    [EMGdata.statistics.mean_mep_amplitude(2,1),...
    EMGdata.statistics.sd_mep_amplitude(2,1)],...
    [EMGdata.statistics.mean_mep_amplitude(3,1),...
    EMGdata.statistics.sd_mep_amplitude(3,1)],...
    [EMGdata.statistics.mean_mep_amplitude(4,1),...
    EMGdata.statistics.sd_mep_amplitude(4,1)],...
    [EMGdata.statistics.mean_mep_amplitude(5,1),...
    EMGdata.statistics.sd_mep_amplitude(5,1)],...
    [EMGdata.statistics.mean_mep_amplitude(6,1),...
    EMGdata.statistics.sd_mep_amplitude(6,1)]);

% assigns to handles of type text the desired parameters 
set(handles.curve_text,'String',texts);

%assigns to handles of type text the name that identifies the pacient 
if contains(file,'analysed') % visualize a previous analysed file 
    set(handles.identification,'String',file(1:strfind(file,'_analysed')-1));
else % visualize a previuos visualized file
    set(handles.identification,'String',file(1:strfind(file,'_visualized')-1));
end

% calls function that plots 1st EMG segment 
plot_figure(EMGdata,handles,a);

handles.EMGdata = EMGdata;
handles.a = a;

% Choose default command line output for visualizeEMG
handles.output = hObject;

% prompt to really close without saving
set(handles.figure1, 'CloseRequestFcn', @closeGUI); 

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes visualizeEMG wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = curve_OutputFcn(hObject, eventdata, handles) 
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
%(maintains mean MEP amplitudes for each intensity, as well as the intensity of 
%half maximum and slope of linear region) 
delete(findall(gcf,'type','annotation'));

% set colors for markers and buttons
MEP_color = [0.5843 0.8157 0.9882];
TMS_color = [1 0 0];

y=EMGdata.trials.EMGfilt{a,1}; % gets the 'a'th segment of EMG data
x = linspace(0,(length(y)/handles.sampling_frequency),length(y)); % time scale for the given segment
axes(handles.graphic); % it gets the plot axis that was designed in GUI
plot(x,y,'k'); %plots EMG data in time domain, with black trace 
grid on; 

xlim1 = EMGdata.trials.artloc(a,1)/handles.sampling_frequency - 0.005; % corrects the window for 
%visualization to start 5 ms before TMS artefact 

xlim2 = xlim1 + 0.1; %correts the window for visualization to end 95 ms after 
%TMS artefact 

xlim([xlim1 xlim2]); % x-axis limits 

[ymin,~] = min(EMGdata.trials.EMGfilt{a,1}); %gets the maximum amplitude value of 
%the given segment 
[ymax,~] = max(EMGdata.trials.EMGfilt{a,1}); %gets the minimum amplitude value of 
%the given segment
ylims = [ymin-0.1 ymax+0.1]; %corrects the previous values to be 0.1 above of 
%ymax or 0.1 bellow of ymin
ylim([ymin-0.1 ymax+0.1]); % y-axis limits 

% Plot Title
set(handles.graphic_title,'String',sprintf(['Trial #:', num2str(a),' - ',...
    num2str(EMGdata.trials.intensities(a,1)),'% intensity']));

%adds x label
xlabel('Time (s)', 'FontSize',10, 'FontName', 'Arial', 'FontWeight', 'bold');

% adds y label
ylabel('ch1(mV)','FontSize',10, 'FontName', 'Arial', 'FontWeight', 'bold');

% adds MEP colored region 
if EMGdata.trials.MEP_onset(a,1)    
    MEP_on = EMGdata.trials.MEP_onset(a,1)/handles.sampling_frequency; %lower limit of MEP
    MEP_off = EMGdata.trials.MEP_offset(a,1)/handles.sampling_frequency; %upper limit of MEP   

    %polygon that identifies the MEP region 
    patch([MEP_on MEP_on MEP_off MEP_off], [ylims(1, 1) ylims(1, 2) ylims(1, 2) ...
        ylims(1, 1)], MEP_color); 
    
    alpha(.7); 
end

%adds line relative to TMS artefact
if EMGdata.trials.artloc(a,1)     
    line([EMGdata.trials.artloc(a,1)/handles.sampling_frequency,...
        EMGdata.trials.artloc(a,1)/handles.sampling_frequency], ...
        ylims(1, :) ,'Color',TMS_color,'LineStyle', '--','Marker','o')       
end

% set MEP amplitude relative to the given EMG segment   
set(handles.MEP_amplitude_trial,'String',sprintf('MEP amp. (mV): %0.3f',...
    EMGdata.trials.MEP_amplitude(a,1)));
end

% --- Executes on button press in export_excel.
function export_excel_Callback(hObject, eventdata, handles)
% hObject    handle to export_excel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global pathname
EMGdata = handles.EMGdata; % get EMGdata structure, that contains all segments 
intensity_values = unique(EMGdata.trials.intensities(:,1)); %intensity values contained in 
%the randomized key
file = handles.file; % name of the selected file to be visualized 

if contains(file,'analysed') % visualize a previuos analyzed file
    file = file(1:strfind(file,'_analysed')-1);
else % visualize a previuos visualized file
    file = file(1:strfind(file,'_visualized')-1);
end

filename = [handles.pathname,'\reviewdata_',file,'.xlsx']; %name of the Excel file to be created 

total_data = table();
total_data.intensities = EMGdata.trials.intensities(:,1); %randomized key
total_data.mep_amp_mV = EMGdata.trials.MEP_amplitude(:,1); %MEP amplitudes
writetable(total_data,filename,'Sheet',1,'Range','A1'); %adds information to 1st sheet 
    
analysed_statistics = table(); 
analysed_statistics.intensities = intensity_values(:,1); 
analysed_statistics.mean_mep_mV = EMGdata.statistics.mean_mep_amplitude(:,1); %mean of MEP 
%amplitudes, for each intensity 
analysed_statistics.sd_mep_mV = EMGdata.statistics.sd_mep_amplitude(:,1); %standard deviation 
%of MEP amplitudes, for each intensity  
writetable(analysed_statistics,filename,'Sheet',2,'Range','A1'); %adds information to 2nd sheet 

curve_analysis = table();
curve_analysis.slope = EMGdata.results.slope2; %slope of the linear region of I/O curve 
curve_analysis.intensity_halfmax = EMGdata.results.intensity_halfmax; %half-maximum intensity 

writetable(curve_analysis,filename,'Sheet',2,'Range','E1'); %adds information to 2nd sheet, 
%starting in E1 cell
end

% --- Executes on button press in ch1_clearTMSart.
function ch1_clearTMSart_Callback(hObject, eventdata, handles)
% hObject    handle to ch1_clearTMSart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

a=handles.a; % gets the index value of the present segment 
EMGdata=handles.EMGdata; % get EMGdata structure, that contains all segments 

EMGdata.trials.artloc(a,1)=0; % sets to zero the position of the TMS artefact of 
%the given segment 

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
EMGdata=handles.EMGdata; % get EMGdata structure, that contains all segments 

axh = gca;  
axh.Toolbar = matlab.ui.controls.AxesToolbar(); 

% manually selects TMS artefact 
x_new = ginput(1);
EMGdata.trials.artloc(a,1) = x_new(1)*handles.sampling_frequency; % saves in EMGdata structure the 
%position of TMS artefact in samples

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
EMGdata.trials.MEP_onset(a,1) = x_new(1)*handles.sampling_frequency; % saves in EMGdata structure the 
%position of onset MEP limit in samples

EMGdata.trials.MEP_offset(a,1) = x_new(2)*handles.sampling_frequency; % saves in EMGdata structure the 
%position of offset MEP limit in samples

MEPchannel=EMGdata.trials.EMGfilt{a,1}; %it gets the 'a'th segment data
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

function enter_trial_Callback(hObject, eventdata, handles)
% hObject    handle to enter_trial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of enter_trial as text
%        str2double(get(hObject,'String')) returns contents of enter_trial as a double

a = str2double(get(hObject,'String')); % atualizes the segment to be visualized 
EMGdata = handles.EMGdata; % get EMGdata structure, that contains all segments 

plot_figure(EMGdata,handles,a) % calls function that plots 'a'th segment of EMG data 

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

%update handles fields
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

%update handles fields
handles.EMGdata = EMGdata;
handles.a = a;
handles.output = hObject;

guidata(hObject, handles);
end

% --- Executes on button press in input_output_curve.
function input_output_curve_Callback(hObject, eventdata, handles)
% hObject    handle to input_output_curve (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

EMGdata = handles.EMGdata; % get EMGdata structure, that contains all segments 

% initial values to calculate the 4 parameters in fhSigmo
mid_point = EMGdata.statistics.intensities(fix(end/2));
vBO = [nanmin(EMGdata.statistics.mean_mep_amplitude(:,1)) 0.1 ...
    mid_point max(EMGdata.statistics.mean_mep_amplitude(:,1))];

% 4 parameters sigmoid function
fhSigmo = @(p,x) abs(p(1)) + p(4)./(1+exp(p(2).*(p(3)-x)));

% non linear fit, taking the intensities and MEP values, as well as
% the initial values that estimates where the curve passes. vB
% contains the 4 parameters of fhSigmo
[vB,~,~,~,MSE] = nlinfit(EMGdata.statistics.intensities,...
    EMGdata.statistics.mean_mep_amplitude, fhSigmo, vBO);

b = fhSigmo(vB,vB(3))- EMGdata.results.slope2*vB(3); %y-intercept 

% Single point plotting
figure(1); 
plot(EMGdata.statistics.intensities(:,1),EMGdata.statistics.mean_mep_amplitude(:,1),'b*','MarkerSize',4); hold on;

vX = [EMGdata.statistics.intensities(1)-5 :  1 : EMGdata.statistics.intensities(end)+5]; 
% vector that contains the intensity range 
y = fhSigmo(vB,vX); % I/O curve for the vX intensity range 

plot(vX, y, 'r','LineWidth',2); hold on; % Sigmoid fit plotting
plot(vX, EMGdata.results.slope2.*vX + b,'--m'); % plotting of linear region 

if b<0 %y-intercept
    text = sprintf([' MSE = ',num2str(MSE),' mV \n y = ', ...
    num2str(EMGdata.results.slope2),'*x - ',num2str(abs(b)),'\n s50 = ', ...
    num2str(vB(3)),'%%RMT']);
else %y-intercept
    text = sprintf([' MSE = ',num2str(MSE),' mV \n y = ', ...
    num2str(EMGdata.results.slope2),'*x + ',num2str(b),'\n s50 = ', ...
    num2str(vB(3)),'%%RMT']);
end

annotation('textbox',[0.15 0.8 0.3 0.1],'String',text,'FitBoxToText','on'); %text in I/O curve plot, 
%with information relative to MSE, line equation and half-maximum intensity

ylim([min(y)-1 max(y)+1]);
xlim([85 145]);
xlabel('Stimulation power (%RMT)','FontSize',10); %adds x label 
ylabel('Amplitude (mV)','FontSize',10); % adds y label 
title('Input-Output Curve');
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

intensities_values = unique(EMGdata.trials.intensities(:,1));
    
EMGdata.statistics = table();
EMGdata.statistics.intensities(:,1) = intensities_values;

%calculate the mean MEP amplitudes for each intensity value 
for u= 1:length(intensities_values)
    index = find(EMGdata.trials.intensities(:,1) == intensities_values(u));
    mep_values = EMGdata.trials.MEP_amplitude(index,1); %MEP amplitudes relative to the uth 
    %intensity value
    mep_values(mep_values == 0) = NaN; %doesn't account MEP amplitudes equal to 0
    EMGdata.statistics.mean_mep_amplitude(u,1) = mean(mep_values,'omitnan'); %calculate mean
    EMGdata.statistics.sd_mep_amplitude(u,1) = std(mep_values,'omitnan'); %calculate standard deviation
    EMGdata.statistics.cv_mep_amplitude(u,1) = ...
        EMGdata.statistics.sd_mep_amplitude(u,1)/...
        EMGdata.statistics.mean_mep_amplitude(u,1); % calculate variation coefficient
end

%initial values to calculate the 4 parameters in fhSigmo
mid_point = EMGdata.statistics.intensities(fix(end/2));
vBO = [nanmin(EMGdata.statistics.mean_mep_amplitude(:,1)) 0.1 mid_point ...
    max(EMGdata.statistics.mean_mep_amplitude(:,1))];

% 4 parameters sigmoid function
fhSigmo = @(p,x) abs(p(1)) + p(4)./(1+exp(p(2).*(p(3)-x)));

% non linear fit, taking the intensities and MEP values, as well as
% the initial values that estimates where the curve passes. vB
% contains the 4 parameters of fhSigmo
[vB,~,~,~,~] = nlinfit(EMGdata.statistics.intensities,...
    EMGdata.statistics.mean_mep_amplitude, fhSigmo, vBO);

% slope of sigmoid's linear region, obtained by doing the 1st derivative of
% the sigmoid function (fhSigmo), and calculating for the mid point value of 
% the linear region (parameter in fhSigmo -> vB(3)). assumes that the linear 
% region has the same 1st derivative value). 
EMGdata.results.slope2 = vB(4)*vB(2)/4;
    
EMGdata.results.intensity_halfmax = vB(3); % structure that contains field with 
%the half maximum intensity,relative to a parameter in fhSigmo that was estimated by 
%non linar fit

handles.EMGdata = EMGdata;

set(handles.curve_text,'String', ''); %eliminates the previous parameters (mean 
%MEP amplitudes for each intensity, as well as the intensity of half maximum and 
%slope of linear region)

% writes the desired parameters (mean MEP amplitudes for each intensity, as well as 
%the intensity of half maximum and slope of linear region)
texts = sprintf(['Slope of linear portion (mV/stim): %0.3f\n stim intensity of half max = %0.3f\n mean amp. ', num2str(EMGdata.statistics.intensities(1,1)),':%0.3f ± %0.3f (mV)\n mean amp. ', num2str(EMGdata.statistics.intensities(2,1)),': %0.3f ± %0.3f (mV)\n mean amp. ', num2str(EMGdata.statistics.intensities(3,1)),': %0.3f ± %0.3f (mV)\n mean amp. ', num2str(EMGdata.statistics.intensities(4,1)),': %0.3f ± %0.3f (mV)\n mean amp. ', num2str(EMGdata.statistics.intensities(5,1)),': %0.3f ± %0.3f (mV)\n mean amp. ', num2str(EMGdata.statistics.intensities(6,1)),': %0.3f ± %0.3f (mV)'],...
    EMGdata.results.slope2,...
    EMGdata.results.intensity_halfmax,...
    [EMGdata.statistics.mean_mep_amplitude(1,1),...
    EMGdata.statistics.sd_mep_amplitude(1,1)],...
    [EMGdata.statistics.mean_mep_amplitude(2,1),...
    EMGdata.statistics.sd_mep_amplitude(2,1)],...
    [EMGdata.statistics.mean_mep_amplitude(3,1),...
    EMGdata.statistics.sd_mep_amplitude(3,1)],...
    [EMGdata.statistics.mean_mep_amplitude(4,1),...
    EMGdata.statistics.sd_mep_amplitude(4,1)],...
    [EMGdata.statistics.mean_mep_amplitude(5,1),...
    EMGdata.statistics.sd_mep_amplitude(5,1)],...
    [EMGdata.statistics.mean_mep_amplitude(6,1),...
    EMGdata.statistics.sd_mep_amplitude(6,1)]);

% assigns to handles of type text the desired parameters 
set(handles.curve_text,'String',texts);
display('saved');

guidata(hObject, handles); % Update handles structure

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
           if contains(file,'analysed') %if the visualized file was a previous 
               %analyzed file (1st edition) 
    
               new_visualized_file = [file(1:strfind(file,'_analysed')-1),'_visualized.mat']; % adds to global 
               % variable (cell array) the name of the new visualized file
    
               save([handles.pathname,'\',new_visualized_file],'trials','results','statistics');
               % save structures in .mat file, in the previously selected path  
               
           else %if the visualized file was a previous visualized file (it isn't the 1st edition) 
               
               new_edited_visualized_file=file; % adds to global variable (cell array) the name 
               % of the new edited visualized file
    
               save([handles.pathname,'\',new_edited_visualized_file], 'trials','results','statistics');
               % save structures in .mat file, in the previously selected path 
           end
   
           % Get default command line output from handles structure, and deletes figure of present GUI 
           %after uiresume
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

[intensities, sortorder] = sort(EMGdata.trials.intensities); %sort in crescent order the intensity 
%values in the randomized key 
mep_amp = EMGdata.trials.MEP_amplitude(sortorder); %sort in the same order as it was done in 
%intensities vector
mep_amp(mep_amp == 0) = NaN; %non-quantified MEPs are associate to NaN values, for not to be 
%being accounted to marginal distributions
h = figure(1);
s = scatterhist(intensities,mep_amp,'Group',intensities,'Kernel','on',...
    'Location','SouthEast','Direction','out',...
    'LineStyle',':','Legend','off','Parent',h); % plot marginal distribution 
delete(s(2));
title('MEP Dispersion for each Stimulus Intensity');
xlabel('Stimulation Intensity (%RMT)');
ylabel('MEP amplitude (mV)');

%text with the variation coefficient for each intensity value 
str_cv = sprintf('CV 90RMT = %0.3f\n CV 100RMT = %0.3f\n CV 110RMT = %0.3f\n CV 120RMT = %0.3f \n CV 130RMT = %0.3f\n CV 140RMT = %0.3f ',...
    EMGdata.statistics.cv_mep_amplitude(1,1),...
    EMGdata.statistics.cv_mep_amplitude(2,1),...
    EMGdata.statistics.cv_mep_amplitude(3,1),...
    EMGdata.statistics.cv_mep_amplitude(4,1),...
    EMGdata.statistics.cv_mep_amplitude(5,1),...
    EMGdata.statistics.cv_mep_amplitude(6,1));
 
annotation('textbox',[0.85 0.75 0.15 0.2],'String',str_cv,'FitBoxToText','on');

vMean = EMGdata.statistics.mean_mep_amplitude(:,1); %mean MEP amplitudes 
vXvalues = EMGdata.statistics.intensities(:,1);
vStd = EMGdata.statistics.sd_mep_amplitude(:,1); %standard deviation of MEP amplitudes 

figure(2);
plot(vXvalues,vMean,'bo--','MarkerSize',8); hold on;
errorbar(vXvalues,vMean,vStd,'LineStyle','none'); %errobar as being the standard deviation of 
%MEP amplitudes
title('MEP Amplitude in Input-Output Paradigm');
ylabel('MEP Amplitude');
xlabel('stimulation intensity (%RMT)');
xlim([80 150]);
    
%% Bladan-Altman analysis
if handles.agremment_measures %if desired to see the agreement between automatic 
%and manual quantification 

    if contains(algorithm_data,'analysed') % visualize a previuos analyzed file
        file = algorithm_data(1:strfind(algorithm_data,'_analysed')-1);
    else % visualize a previuos visualized file
        file = algorithm_data(1:strfind(algorithm_data,'_visualized')-1);
    end

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

        intensities = unique(human_data(:,1));
    
        methods_amplitudes.alghorithm_mep = []; %automatic MEP quantification
        methods_amplitudes.human_mep = []; % manual MEP quantification
    
        for i= 1:length(intensities)
            index = find(human_data(:,1) == intensities(i)); %indexes of the ith intensity value 
            %on the randomized key 
            methods_amplitudes.alghorithm_mep(:,end+1) = EMGdata.trials.MEP_amplitude(index,1);
            %mep amplitudes automatically quantified, for the ith intensity value 
            
            methods_amplitudes.human_mep(:,end+1) = human_data(index,2);
            %mep amplitudes manually quantified, for the ith intensity value 
        end
 
        %% load data
        territories = {'90%RMT','100%RMT','110%RMT','120%RMT','130%RMT','140%RMT'};

        data1= methods_amplitudes.alghorithm_mep(:,:); % automatic MEP quantification 
        data2= methods_amplitudes.human_mep(:,:); % manual MEP quantification

        % BA plot paramters
        tit1 = 'MEP Amplitude Correlation'; % figure title
        tit = 'MEP Amplitude Accuracy (Bland-Altman)';% figure title
        gnames = {territories}; % names of groups in data 
        label = {'Algorithm Measure','Human Measure','mV'}; % Names of data sets
        corrinfo = {'n','SSE','r2','eq','P'}; % stats to display of correlation scatter plot
        limits = 'auto'; % how to set the axes limits
        colors = 'brgmkc';      % character codes

    %% Example 2 - using non-Gaussian data and one figure with 2 analyses

        [~, fig, statsStruct2] = BlandAltman(data1, data2,label,...
    [tit ' (using non-parametric stats)'],gnames,'corrInfo',corrinfo,...
    'axesLimits',limits,'colors',colors,'baYLimMode','Square','showFitCI',...
    ' on','baStatsMode','non-parametric');
        disp('Statistical results (using non-parametric stats):')
        disp(statsStruct2)

    [~, fig, ~] = correlationPlot(data1, data2,label,tit1,gnames,'corrInfo',...
    corrinfo,'colors',colors);
    end
else
    warndlg(['If you desired to compare the', ...
        ' automatic quantification versus manual, please verify the settings on ''Change',...
        ' algorithm specifications'' in the previous interface (''choosedata'' - where', ...
        ' the files are organized)']);
end
end
