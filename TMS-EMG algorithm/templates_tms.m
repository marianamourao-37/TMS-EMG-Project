function varargout = templates_tms(varargin)
% TEMPLATES_TMS MATLAB code for templates_tms.fig
%      TEMPLATES_TMS, by itself, creates a new TEMPLATES_TMS or raises the existing
%      singleton*.
%
%      H = TEMPLATES_TMS returns the handle to a new TEMPLATES_TMS or the handle to
%      the existing singleton*.
%
%      TEMPLATES_TMS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TEMPLATES_TMS.M with the given input arguments.
%
%      TEMPLATES_TMS('Property','Value',...) creates a new TEMPLATES_TMS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before templates_tms_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to templates_tms_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help templates_tms

%% 
% Function description: It refers to the GUI that allows the user to:
% 1. modify TMS templates
% 2. see the current TMS templates 
% 3. change to the default TMS templates 

% Besides that, it compares samples of TMS artefacts from the participants
% selected to construted them, to see if they're adequate, i.e, have
% similar shape (same polarity, same duration, and they can have different
% amplitudes, as it depends on the stimulus intensity). 

% It has as input variable (constant_time_tms) the time constant (s) to construct TMS 
% templates, related to an estimate of a TMS pulse duration. 

%%

% Last Modified by GUIDE v2.5 29-Oct-2019 23:54:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @templates_tms_OpeningFcn, ...
                   'gui_OutputFcn',  @templates_tms_OutputFcn, ...
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


% --- Executes just before templates_tms is made visible.
function templates_tms_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to templates_tms (see VARARGIN)

handles.constant_time_tms = varargin{1}; % attribution of the time duration of TMS pulse 
% to respective handle
handles.sampling_frequency = varargin{2};

global basepath %global variable related to the path's name of selected directory

% counts how many files with the required extensions exists in the selected path
files_mat = dir(fullfile(basepath, '*.mat'));
files_acq = dir(fullfile(basepath, '*.acq'));

%creation of cell array that identifies .mat files with raw data
%(continuous data) from both protocols
files = {};

%for loop that reads .mat files' names from the selected directory, and identifies which 
% are relative to the category 'raw data', by comparing if they include the name of an .acq file 
% present in the selected directory, and doesn't contain, as well, strings representative of 
%different stages of data processing (segmented, analysed and visualized).
%If so, they are added in the end+1 position of the cell array
for b = 1:length(files_acq)
    file_acq = files_acq(b).name;
    for j = 1:length(files_mat)
        file_mat = files_mat(j).name; % mat file's name
        if contains(file_mat(1:end-4),file_acq(1:end-4)) && ...
                ~contains(file_mat(1:end-4), 'segmented') && ...
                ~contains(file_mat(1:end-4), 'analysed') && ...
                ~contains(file_mat(1:end-4), 'visualized')
            files{end+1} = file_mat; %adds file's name in the end+1 position of the respective cell
            %array
        end
    end
end

% sets files's names to the list box, for later selection
set(handles.participants_list,'String', files);

handles.files = files; 

% Choose default command line output for templates_tms
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes templates_tms wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = templates_tms_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in participants_list.
function participants_list_Callback(hObject, eventdata, handles)
% hObject    handle to participants_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns participants_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from participants_list

files = handles.files; % get the cell array with raw data files 

set(handles.participants_list,'max',length(files)); % set the maximum number of permitted files' 
%selections in the list box.

num_mat = get(handles.participants_list,'Value'); % get the index of the selected files in the 
%listbox structure
string_mat = cellstr(get(handles.participants_list,'String')); 
files_mat = string_mat(num_mat); %get the selected files identified by their index in the listbox

input_output_files = {}; % creation of cell array that will contain the corresponding selected 
%files relative to I/O protocol
paired_pulse_files = {}; % creation of cell array that will contain the corresponding selected 
%files relative to paired-pulse protocol

for b = 1:length(files_mat)
    file = files_mat{b}; %name of the bth .mat file contained in the cell array 
    %'files_mat'

    if ~contains(file(1:end-4),'INPUTOUTPUT')  
        paired_pulse_files{end+1} = file; 
    else
        input_output_files{end+1} = file; 
    end
end

handles.input_output_files = input_output_files;
handles.paired_pulse_files = paired_pulse_files;
handles.files_mat = files_mat; 

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function participants_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to participants_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in default_templates.
function default_templates_Callback(hObject, eventdata, handles)
% hObject    handle to default_templates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get the path that identifies the default TMS templates (contain in
%default_patterns folder)
sici_default = 'default_patterns\pattern_SICI.mat';
icf_default = 'default_patterns\pattern_ICF.mat';
lici_default = 'default_patterns\pattern_LICI.mat';
baselines_default = 'default_patterns\pattern_baselines.mat';

%copy and replace 
copyfile(sici_default);
copyfile(icf_default);
copyfile(lici_default);
copyfile(baselines_default);

% --- Executes on button press in extract_templates.
function extract_templates_Callback(hObject, eventdata, handles)
% hObject    handle to extract_templates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global basepath

input_output_files = handles.input_output_files;
paired_pulse_files = handles.paired_pulse_files;
constant_time_tms = handles.constant_time_tms;

%initialization of variables that identify the first time of a paired-pulse sequence 
%(0 -> it didn't occur, 1 -> occur)
iteration_sici = 0; 
iteration_icf = 0;
iteration_lici = 0;

%initialization of variables that will contain, for each ISI, the sum of
%the respective patterns for the participants selected (and contained in
%the cell array paired_pulse_files) 
SICI = [];
ICF = [];
LICI = [];

%initialization of variables that will count, for each ISI, the number of
%times that were extracted TMS patterns 
num_sici = 0;
num_icf = 0;
num_lici = 0;

for g = 1:length(paired_pulse_files)
    file = paired_pulse_files{g}; %name of the gth .mat file contained in the 
    %cell array 'paired_pulse_files'
    
    str = load([basepath,'\',file]); %loads the corresponding .mat file 
    
    acq_data = str.acq_data; %get the acq_data structure
    
    trials = str.trials; %get the table named trials, that will be 
    %filled with intended parameters
    
    locpulse = acq_data.data(:,2); %get data from channel 2 of the 
    %EMG acquisition (channel with trigger information)
    
    acq_data_detrend = detrend(acq_data.data(:,1)); %get data from channel 1 of the 
    %EMG acquisition (channel with raw EMG continuous data)
    
    ISI_values = unique(trials.ISI_sec(:,1)); %get the ISI values  
    
    if ~isempty(find(locpulse < 4.8,1)) %by observation, it was set the threshold of 
        %4.8, and if it's found at least one value lower than threshold, it means 
        %that contains trigger (occurrence of 1st TMS artefact of a pulse pair) 
        %information relevant for posterior analysis
        
        c = fix(constant_time_tms * handles.sampling_frequency); %constant that accounts from the start of 
        %the 2nd TMS artefact to its end (time duration of TMS pulse)
        
        if contains(file(1:end-4),'SICI') && ~iteration_sici %if file name contains the string 'SICI' and 
            %it it is the first occurrence of this paired-pulse sequence in the selected files 
            
            iteration_sici = 1; %first occurrence of SICI sequence

            %initialization of fields in SICI structure, for each ISI, that
            %will contain the sum of the respective TMS patterns. 
            for d = 1:numel(ISI_values)
                SICI.(['sum_patterns',num2str(ISI_values(d)*10^3),'ms']) = ...
                    zeros(c+(ISI_values(d)*handles.sampling_frequency) +1,1); 
                % each pattern has the length c+(ISI_values(d)*handles.sampling_frequency) +1
            end
            num_sici = zeros(numel(ISI_values),1); %initialization, for each ISI, of the number of times 
            %that were extracted TMS patterns 
        end
        
        if contains(file(1:end-4),'ICF') && ~iteration_icf %if file name contains the string 'ICF' and 
            %it is the first occurrence of this paired-pulse sequence in the selected files 
            
            iteration_icf = 1; %first occurrence of ICF sequence
            
            %initialization of fields in ICF structure, for each ISI, that
            %will contain the sum of the respective TMS patterns.
            for d = 1:numel(ISI_values)
                ICF.(['sum_patterns',num2str(ISI_values(d)*10^3),'ms']) = ...
                    zeros(c+(ISI_values(d)*handles.sampling_frequency) +1,1);
                % each pattern has the length c+(ISI_values(d)*handles.sampling_frequency) +1
            end
            num_icf = zeros(numel(ISI_values),1); %initialization, for each ISI, of the number of times 
            %that were extracted TMS patterns 
        end
        
        if contains(file(1:end-4),'LICI') && ~iteration_lici %if file name contains the string 'LICI' and 
            %it is the first occurrence of this paired-pulse sequence in the selected files 
            
            iteration_lici = 1; %first occurrence of LICI sequence
            
            %initialization of fields in LICI structure, for each ISI, that
            %will contain the sum of the respective TMS patterns.
            for d = 1:numel(ISI_values)
                LICI.(['sum_patterns',num2str(ISI_values(d)*10^3),'ms']) = ...
                    zeros(c+(ISI_values(d)*handles.sampling_frequency) +1,1);
                % each pattern has the length c+(ISI_values(d)*handles.sampling_frequency) +1
            end
            num_lici = zeros(numel(ISI_values),1); %initialization, for each ISI, of the number of times 
            %that were extracted TMS patterns 
        end
    
        acq_data.x_start = 0; % initializes the initiate time of effective acquisition
        
        x_start=time_start(acq_data,file,handles.sampling_frequency); % calls 'time_start.mat' GUI for selection of 
    %initiate time of effective acquisition
    
        acq_data.x_start = x_start; % atualizes the value in x_start field
 
        time_delta_index1 = 0; %initializes the lower limit of the segment (in samples)
   
        %for loop with the same number of iteration as the number of ISI
        %in the randomized key (i.e, 40 in this case) 
        for j = 1:height(trials)
            if j ~= length(trials.ISI_sec) % these segments contain IIPP time intervals 
                %(time interval between paired of pulses)
                
                time_delta_index2 = fix((sum(trials.IIPP_sec(1:j,1))*handles.sampling_frequency)); %upper limit of 
                %the segment (in samples)
            
                stpulse = locpulse((x_start+time_delta_index1):(x_start+time_delta_index2),1); 
                %corresponding data in channel 2, for this segment 
            
                time1 = find(stpulse < 4.8,1); 
                if ~isempty(time1) %if this segment contains trigger information 
                    delta = fix(trials.ISI_sec(j,1) * handles.sampling_frequency); % variable that accounts from the 
                    %start of 1st TMS artefact to the start of the 2nd
                
                    %pattern of the corresponding pair of TMS artefacts, for this segment 
                    trials.artefact_pattern{j,1} = ...
                    acq_data_detrend(x_start+time_delta_index1+time1:...
                    x_start+time_delta_index1+time1+delta+c);
                
                    ind = find(ISI_values == trials.ISI_sec(j,1)); %identifies the TMS pair pattern 
                    %(index) 
                    t = ISI_values(ind,1); %ISI value of TMS pair pattern
                    
                    if contains(file(1:end-4),'SICI') %if file name contains the string 'SICI'
                        
                        %sum of previous patterns for this ISI value, and the new one found
                        SICI.(['sum_patterns',num2str(t*10^3),'ms']) = ...
                            SICI.(['sum_patterns',num2str(t*10^3),'ms']) + trials.artefact_pattern{j,1};
                        
                        num_sici(ind) = num_sici(ind) + 1; % actualizes the number of patterns for 
                        %this ISI value 
                    end
                    
                    if contains(file(1:end-4),'ICF') %if file name contains the string 'ICF'
                        
                        %sum of previous patterns for this ISI value, and the new one found
                        ICF.(['sum_patterns',num2str(t*10^3),'ms']) = ...
                            ICF.(['sum_patterns',num2str(t*10^3),'ms']) + trials.artefact_pattern{j,1};
                        
                        num_icf(ind) = num_icf(ind) + 1; % actualizes the number of patterns for 
                        %this ISI value 
                    end
                    
                    if contains(file(1:end-4),'LICI') %if file name contains the string 'LICI'
                        
                        %sum of previous patterns for this ISI value, and the new one found
                        LICI.(['sum_patterns',num2str(t*10^3),'ms']) = ...
                            LICI.(['sum_patterns',num2str(t*10^3),'ms']) + trials.artefact_pattern{j,1}; 
                        
                        num_lici(ind) = num_lici(ind) + 1; % actualizes the number of patterns for 
                        %this ISI value 
                    end
                end
            
                time_delta_index1 = fix((sum(trials.IIPP_sec(1:j,1))*handles.sampling_frequency)); 
                %atualizes the lower limit for the next segment
            
            else %for the last segment, it doesn't exists IIPP time interval 
            
                stpulse = locpulse((x_start+time_delta_index1):end); %corresponding data in channel 2, 
                %for this segment 
            
                time1 = find(stpulse < 4.8,1);
            
                if ~isempty(time1) %if this segment contains trigger information
                    delta = fix(trials.ISI_sec(j,1) * handles.sampling_frequency); % variable that 
                %accounts from the start of 1st TMS artefact to the start of the 2nd

                %pattern of the corresponding pair of TMS artefacts, for this segment 
                    trials.artefact_pattern{j,1} = ...
                    acq_data_detrend(x_start+time_delta_index1+time1:...
                    x_start+time_delta_index1+time1+delta+c,1);
                
                    ind = find(ISI_values == trials.ISI_sec(j,1)); %identifies the 
                    %TMS pair pattern (index) 
                    t = ISI_values(ind,1); %ISI value of TMS pair pattern
                    
                    if contains(file(1:end-4),'SICI') %if file name contains the string 'SICI'
                        
                        %sum of previous patterns for this ISI value, and the new one found
                        SICI.(['sum_patterns',num2str(t*10^3),'ms']) = ...
                            SICI.(['sum_patterns',num2str(t*10^3),'ms']) + ...
                        trials.artefact_pattern{j,1};
                        
                        num_sici(ind) = num_sici(ind) + 1; % actualizes the number of patterns for 
                        %this ISI value
                    end
                    
                    if contains(file(1:end-4),'ICF') %if file name contains the string 'ICF'
                        
                        %sum of previous patterns for this ISI value, and the new one found
                        ICF.(['sum_patterns',num2str(t*10^3),'ms']) = ...
                            ICF.(['sum_patterns',num2str(t*10^3),'ms']) + ...
                        trials.artefact_pattern{j,1};
                        
                        num_icf(ind) = num_icf(ind) + 1; % actualizes the number of patterns for 
                        %this ISI value
                    end
                    
                    if contains(file(1:end-4),'LICI') %if file name contains the string 'LICI'
                        
                        %sum of previous patterns for this ISI value, and the new one found
                        LICI.(['sum_patterns',num2str(t*10^3),'ms']) = ...
                            LICI.(['sum_patterns',num2str(t*10^3),'ms']) + ...
                        trials.artefact_pattern{j,1}; 

                        num_lici(ind) = num_lici(ind) + 1; % actualizes the number of patterns for 
                        %this ISI value
                  
                    end
                end
            end
        end
    else
        warndlg(['File ',file(1:end-4),' doesn''t have any trigger information. It won''t',...
            ' be account for template extraction.']);
    end
end

%gets the field's names 
fields_sici = fieldnames(SICI);
fields_icf = fieldnames(ICF);
fields_lici = fieldnames(LICI);

baseline_paired = 0; % initialization of the sum baseline patterns (ISI = 0ms) in paired-pulse protocol 
num_paired = 0; % initialization of the number of baseline patterns (ISI = 0ms) in paired-pulse protocol 

%for loop that will calculate the TMS templates, for SICI sequence in
%paired-pulse protocol 
for t = 2:numel(fieldnames(SICI))
    isi = char(regexp(fields_sici{t},'\d*','Match')); %gets the ISI value contain in field's name
    if num_sici(t) ~= 0 %if for that ISI there are patterns extracted 
        results_sici.(['mean_pattern',isi,'ms']) = ...
        SICI.(['sum_patterns',isi,'ms'])/num_sici(t); %SICI TMS template for the corresponding ISI
    else
        warndlg(['It wasn''t extrated the template relative to ISI = ',isi,'ms of SICI protocol.',...
            ' It will be associate the respective default template']);
        
        %get the path that identifies the default TMS templates (contain in
        %default_patterns folder)
        sici_default = 'default_patterns\pattern_ICF.mat';
        sici_default = load(sici_default);
        sici_default = sici_default.results_sici.(['mean_pattern',isi,'ms']); 
        results_sici.(['mean_pattern',isi,'ms']) = sici_default;
    end 
end

if num_sici(1) ~= 0 %if there were extracted baselines in SICI sequence of paired-pulse protocol 
    baseline_paired = SICI.sum_patterns0ms;
    num_paired = num_sici(1);
end

%for loop that will calculate the TMS templates, for ICF sequence in
%paired-pulse protocol 
for t = 2:numel(fieldnames(ICF))
    isi = char(regexp(fields_icf{t},'\d*','Match')); %gets the ISI value contain in field's name
    if num_icf(t) ~= 0 %if for that ISI there are patterns extracted 
        results_icf.(['mean_pattern',isi,'ms']) = ...
        ICF.(['sum_patterns',isi,'ms'])/num_icf(t); %ICF TMS template for the corresponding ISI
    else 
        warndlg(['It wasn''t extrated the template relative to ISI = ',isi,'ms of ICF protocol.', ...
            ' It will be associate the respective default template']);
        
        %get the path that identifies the default TMS templates (contain in
        %default_patterns folder)
        icf_default = 'default_patterns\pattern_ICF.mat';
        icf_default = load(icf_default);
        icf_default = icf_default.results_icf.(['mean_pattern',isi,'ms']);
        results_icf.(['mean_pattern',isi,'ms']) = icf_default;
    end                 
end

if num_icf(1) ~= 0 %if there were extracted baselines in ICF sequence of paired-pulse protocol 
    if num_sici(1) ~= 0 %if there were extracted baselines in SICI sequence of paired-pulse protocol 
        baseline_paired = baseline_paired + ICF.sum_patterns0ms;
        num_paired = num_sici(1) + num_icf(1);
    else 
        baseline_paired = ICF.sum_patterns0ms;
        num_paired = num_icf(1);
    end
end

%for loop that will calculate the TMS templates, for LICI sequence in
%paired-pulse protocol 
for t = 2:numel(fieldnames(LICI))
    isi = char(regexp(fields_lici{t},'\d*','Match')); %gets the ISI value contain in field's name
    if num_lici(t) ~= 0 %if for that ISI there are patterns extracted 
        results_lici.(['mean_pattern',isi,'ms']) = ...
        LICI.(['sum_patterns',isi,'ms'])/num_lici(t); %LICI TMS template for the corresponding ISI
    else
        warndlg(['It wasn''t extrated the template relative to ISI = ',isi,'ms of LICI protocol.', ...
            ' It will be associate the respective default template']);
        
        %get the path that identifies the default TMS templates (contain in
        %default_patterns folder)
        lici_default = 'default_patterns\pattern_ICF.mat';
        lici_default = load(lici_default);
        lici_default = lici_default.results_lici.(['mean_pattern',isi,'ms']);
        results_lici.(['mean_pattern',isi,'ms']) = lici_default;
    end 
end

if num_lici(1) ~= 0 %if there were extracted baselines in LICI sequence of paired-pulse protocol 
    if num_sici(1) ~= 0 && num_icf(1) ~= 0 %if there were extracted baselines in ICF and SICI sequence 
        %of paired-pulse protocol 
        baseline_paired = baseline_paired + LICI.sum_patterns0ms;
        num_paired = num_lici(1) + num_sici(1) + num_icf(1);
    elseif num_sici(1) ~= 0 && num_icf(1) == 0 %if there were extracted baselines in SICI sequence 
        %of paired-pulse protocol 
        baseline_paired = baseline_paired + LICI.sum_patterns0ms;
        num_paired = num_lici(1) + num_sici(1);
    elseif num_sici(1) == 0 && num_icf(1) ~= 0 %if there were extracted baselines in ICF sequence 
        %of paired-pulse protocol 
        baseline_paired = baseline_paired + LICI.sum_patterns0ms;
        num_paired = num_lici(1) + num_icf(1);
    elseif num_sici(1) == 0 && num_icf(1) == 0 %if there weren't extracted baselines in ICF and SICI 
        % sequence of paired-pulse protocol 
        baseline_paired = LICI.sum_patterns0ms;
        num_paired = num_lici(1);
    end
end

if ~isempty(results_sici) %if there were extracted SICI templates from the selected files 
    save('pattern_SICI', 'results_sici'); %save SICI TMS templates extracted from the selected files 
else %get the default SICI TMS templates 
    warndlg(['It wasn''t extrated the template relative to SICI protocol. It was',...
        ' associate the default SICI template']);
    
    %get the path that identifies the default TMS templates (contain in
    %default_patterns folder)
    sici_default = 'default_patterns\pattern_SICI.mat';
    copyfile(sici_default);
end

if ~isempty(results_icf) %if there were extracted ICF templates from the selected files 
    save('pattern_ICF', 'results_icf'); %save ICF TMS templates extracted from the selected files 
else %get the default ICF TMS templates 
    warndlg(['It wasn''t extrated the template relative to ICF protocol. It was',...
        ' associate the default ICF template']);
    
    %get the path that identifies the default TMS templates (contain in
    %default_patterns folder)
    icf_default = 'default_patterns\pattern_ICF.mat';
    copyfile(icf_default);
end

if ~isempty(results_lici) %if there were extracted LICI templates from the selected files 
    save('pattern_LICI', 'results_lici'); %save LICI TMS templates extracted from the selected files 
else %get the default LICI TMS templates 
    warndlg(['It wasn''t extrated the template relative to LICI protocol. It was',...
        ' associate the default LICI template']);
    
    %get the path that identifies the default TMS templates (contain in
    %default_patterns folder)
    lici_default = 'default_patterns\pattern_LICI.mat';
    copyfile(lici_default);
end

pattern_sum = []; %initializes the sum of TMS patterns 
num = 0; %initializes the number of the extracted TMS patterns 

for g = 1:length(input_output_files)
    file = input_output_files{g}; %name of the gth .mat file contained in the 
    %cell array 'input_output_files'
    
    str = load([basepath,'\',file]); %loads the corresponding .mat file 
    
    acq_data = str.acq_data; %get the acq_data structure
    
    trials = str.trials; %get the table named trials, that will be filled with intended parameters
    
    locpulse = acq_data.data(:,2); %get data from channel 2 of the 
    %EMG acquisition (channel with trigger information)  
    
    acq_data_detrend = detrend(acq_data.data(:,1)); %get data from channel 1 of the 
    %EMG acquisition (channel with raw EMG continuous data)
    
    if ~isempty(find(locpulse < 4.8,1)) %by observation, it was set the threshold of 
        %4.8, and if it's found at least one value lower than threshold, it means 
        %that contains trigger (occurrence of 1st TMS artefact of a pulse pair) 
        %information relevant for posterior analysis
        
        c = fix(constant_time_tms * handles.sampling_frequency); %constant that accounts from the start of 
        %the 2nd TMS artefact to its end
        
        pattern_sum = zeros(c+1,1); %initializes the sum of TMS patterns, with length c+1 
        
        acq_data.x_start = 0; % initializes the initiate time of effective acquisition
        
        x_start=time_start(acq_data,file,handles.sampling_frequency); % calls 'time_start.mat' GUI for selection of 
    %initiate time of effective acquisition
    
        acq_data.x_start = x_start; % atualizes the value in x_start field
 
        time_delta_index1 = 0; %initializes the lower limit of the segment (in samples)
    
        %for loop with the same number of iteration as the number of
        %intensities values in the randomized key (i.e, 60 in this case) 
        for j = 1:height(trials)
            if j ~= height(trials) % this segments contains IIPP time intervals
                time_delta_index2 = fix((sum(trials.IIPP_sec(1:j,1))*...
                handles.sampling_frequency)); %upper limit of the segment (in samples)
            
                stpulse = locpulse((x_start+time_delta_index1):...
                (x_start+time_delta_index2),1); %corresponding data in channel 2, for this segment 
            
                time1 = find(stpulse < 4.8,1); 
                if ~isempty(time1) %if this segment contains trigger information 
                
                %pattern of the corresponding TMS artefact, for this segment 
                    trials.artefact_pattern{j,1} = ...
                    acq_data_detrend(x_start+time_delta_index1+time1:...
                    x_start+time_delta_index1+time1+c);
                
                    pattern_sum = pattern_sum + trials.artefact_pattern{j,1};
                    num = num +1;
                end
            
                time_delta_index1 = fix((sum(trials.IIPP_sec(1:j,1))*handles.sampling_frequency)); 
                %atualizes the lower limit for the next segment
            
            else %for the last segment, it doesn't exists IIPP time interval 
            
                stpulse = locpulse((x_start+time_delta_index1):end,1); %corresponding 
            %data in channel 2, for this segment 
            
                time1 = find(stpulse < 4.8,1);
            
                if ~isempty(time1) %if this segment contains trigger 
                %information
                
                    trials.locpulse(j,1)=1; % identifies that this segment contains 
                %trigger information
                
                %pattern of the correspondig pair of TMS artefacts, for this segment 
                    trials.artefact_pattern{j,1} = ...
                    acq_data_detrend(x_start+time_delta_index1+time1:...
                    x_start+time_delta_index1+time1+c);
                
                    %sum of previous TMS patterns and the new one found
                    pattern_sum = pattern_sum + trials.artefact_pattern{j,1};
                    num = num +1; %actualizes the number of TMS patterns 
                end
            end
        end
    else
        warndlg(['File ',file(1:end-4),' doesn''t have any trigger information. It won''t',...
            ' be account for template extraction.']);
    end
end

if num ~=0 || num_paired ~= 0 % if single-pulse was extracted at least from one of the protocols 
    if num~=0 && num_paired == 0 %single-pulse extracted just from I/O protocol 
        results.mean_pattern = pattern_sum/num;
    elseif num~=0 && num_paired ~= 0 %single-pulse extracted from both protocols 
        results.mean_pattern = (pattern_sum + baseline_paired)/...
        (num+num_paired);
    elseif num_paired ~= 0 && num==0 %single-pulse extracted just from paired-pulse protocol 
        results.mean_pattern = baseline_paired/num_paired;
    end
    save('pattern_baselines', 'results'); %save extracted single-pulse template 
else
    warndlg(['It wasn''t extracted the template relative to single-pulse TMS (baseline).',...
            ' It was associate the default baseline template']);
        
    %get the path that identifies the default TMS templates (contain in
    %default_patterns folder)
    baselines_default = 'default_patterns\pattern_baselines.mat';
    copyfile(baselines_default);
end


% --- Executes on button press in view_templates.
function view_templates_Callback(hObject, eventdata, handles)
% hObject    handle to view_templates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get the existing templates in the selected directory 
icf = load('pattern_ICF.mat');
lici = load('pattern_LICI.mat');
sici = load('pattern_SICI.mat');
baselines = load('pattern_baselines.mat');

%get field's names from structures 
fields_sici = fieldnames(sici.results_sici);
fields_icf = fieldnames(icf.results_icf);
fields_lici = fieldnames(lici.results_lici);

fig = 0;

%for loop that runs for each ISI in SICI sequence of paired-pulse protocol 
for t = 1:numel(fields_sici) 
    isi = char(regexp(fields_sici{t},'\d*','Match')); %gets the ISI value contain in field's name
    fig = fig +1;
    figure(fig)
    plot(sici.results_sici.(['mean_pattern',isi,'ms'])) %plots the corresponding TMS template
    title(['Template of TMS artifact for ISI = ',isi,' ms']);
    xlabel('Samples');
    ylabel('Amplitude (mV)');
end

%for loop that runs for each ISI in ICF sequence of paired-pulse protocol 
for t = 1:numel(fields_icf)
    isi = char(regexp(fields_icf{t},'\d*','Match')); %gets the ISI value contain in field's name
    fig = fig +1;
    figure(fig)
    plot(icf.results_icf.(['mean_pattern',isi,'ms'])) %plots the corresponding TMS template
    title(['Template of TMS artifact for ISI = ',isi,' ms']);
    xlabel('Samples');
    ylabel('Amplitude (mV)');
end

%for loop that runs for each ISI in LICI sequence of paired-pulse protocol 
for t = 1:numel(fields_lici)
    isi = char(regexp(fields_lici{t},'\d*','Match')); %gets the ISI value contain in field's name
    fig = fig +1;
    figure(fig)
    plot(lici.results_lici.(['mean_pattern',isi,'ms'])) %plots the corresponding TMS template
    title(['Template of TMS artifact for ISI = ',isi,' ms']);
    xlabel('Samples');
    ylabel('Amplitude (mV)');
end

fig = fig +1;
figure(fig)
plot(baselines.results.mean_pattern); %plots the baseline template
title('Template of TMS artifact for ISI = 0 ms (baseline)');
xlabel('Samples');
ylabel('Amplitude (mV)');


% --- Executes on button press in polarity.
function polarity_Callback(hObject, eventdata, handles)
% hObject    handle to polarity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global basepath 
intensities = 100; %intensity from which is pretended to get a sample of TMS pattern in I/O protocol
isi = 0; % ISI from which is pretended to get a sample of TMS pattern in paired-pulse protocol

files_mat = handles.files_mat; %selected files in list box 
constant_time_tms = handles.constant_time_tms; %time constant that accounts for time duration of 
%TMS pulse

trials_check = {}; %initializes the cell array that will contain a TMS artifact sample, for each 
%participant
identify_participants = []; % identify the participants from whom pulse characteristics will be compared

condition = 1; % variable that identifies the first occurence of intensity = 100 (I/O protocol) or 
%isi = 0 (paired-pulse protocol) in each selected file

for g = 1:length(files_mat)
    file = files_mat{g}; %name of the gth .mat file contained in the 
    %cell array 'files_mat'
    
    str = load([basepath,'\',file]); %loads the corresponding .mat file 
    
    acq_data = str.acq_data; %get the acq_data structure
    
    trials = str.trials; %get the table named trials, that will be 
    %filled with intended parameters
    
    locpulse = acq_data.data(:,2); %get data from channel 2 of the 
    %EMG acquisition (channel with trigger information)
    
    acq_data_detrend = detrend(acq_data.data(:,1));  %get data from channel 1 of the 
    %EMG acquisition (channel with raw EMG continuous data)
    
    if ~isempty(find(locpulse < 4.8,1)) %by observation, it was set the threshold of 
        %4.8, and if it's found at least one value lower than threshold, it means 
        %that contains trigger (occurrence of 1st TMS artefact of a pulse pair) 
        %information relevant for posterior analysis
        
        c = fix(constant_time_tms * handles.sampling_frequency); %constant that accounts from the start of 
        %the 2nd TMS artefact to its end
    
        acq_data.x_start = 0; % initializes the initiate time of effective acquisition
        
        x_start=time_start(acq_data,file,handles.sampling_frequency); % calls 'time_start.mat' GUI for selection of 
    %initiate time of effective acquisition
    
        acq_data.x_start = x_start; % atualizes the value in x_start field
 
        time_delta_index1 = 0; %initializes the lower limit of the segment (in samples)
    
        for j = 1:height(trials)
            if condition % this segments contains IIPP time intervals
                time_delta_index2 = fix((sum(trials.IIPP_sec(1:j,1))*...
                handles.sampling_frequency)); %upper limit of the segment (in samples)
            
                stpulse = locpulse((x_start+time_delta_index1):...
                (x_start+time_delta_index2),1); %corresponding data in channel 2, for this segment 
            
                time1 = find(stpulse < 4.8,1); 
                
                if ~isempty(time1) %if this segment contains trigger information
                    if ~contains(file(1:end-4),'INPUTOUTPUT') % if it is a file relative to paired-pulse
                        %protocol 
                        if trials.ISI_sec(j,1) == isi
                    
                            trials_check{end+1} = ...
                        acq_data_detrend(x_start+time_delta_index1+time1:...
                        x_start+time_delta_index1+time1+c); %sets in position end+1 a sample of TMS 
                        % patterns of isi = 0 ms for this participant 
                
                            identify_participants = [identify_participants,g];
                            % actualizes the participants that are going to be compared g is a number
                    
                            condition = 0;
                        end
                    elseif contains(file(1:end-4),'INPUTOUTPUT') % if it is a file relative to I/O protocol 
                        if trials.intensities(j,1) == intensities
                            
                            trials_check{end+1} = ...
                        {acq_data_detrend(x_start+time_delta_index1+time1:...
                        x_start+time_delta_index1+time1+c)}; %sets in position end+1 a sample of TMS 
                        % patterns of isi = 0 ms for this participant
                
                            identify_participants = [identify_participants,g]; % actualizes the 
                            %participants that are going to be compared g is a number
                    
                            condition = 0; %stops the analysis for this participant, since a sample of 
                            %TMS artifact has already been extracted for this participant 
                        end
                    end
                end
                time_delta_index1 = fix((sum(trials.IIPP_sec(1:j,1))*handles.sampling_frequency)); 
                %atualizes the lower limit for the next segment
            else
                break % it has already been extracted a sample of TMS artifact for this participant
            end
        end
    else
        warndlg(['File ',file(1:end-4),' doesn''t have any trigger information. It won''t',...
            ' be accounted for test TMS pulses similarity.']);
    end
    condition = 1; %for the next participant, reinitializes the analysis of the first occurrence 
    %of the pretended TMS artefact sample
end

for i = 1:length(identify_participants)
    figure(1);
    data = trials_check{1,i}; 
    plot(data'); hold on; %plot TMS artifact sample for the ith selected participant
end

figure(1);
legend(num2str(identify_participants(:))); %adds legend into the figure that has plotted the 
%TMS artifacts samples 
title('Comparation of TMS artifacts for the participants selected');
xlabel('Samples');
ylabel('Amplitude (mV)');
    
