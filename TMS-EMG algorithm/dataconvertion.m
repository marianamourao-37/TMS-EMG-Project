function dataconvertion(pathname,protocol,file_acq)
%%
%Function description: function that converts .acq files in .mat files

% - Input variables: 
%   pathname - chosen directory for analyzing the contained files
%   protocol - chosen protocol, that will analyze only the corresponding files 
%   file_acq - identified .acq files, that will be converted to .mat files
%   (output of the present function) 

trials = table(); %creation of table 
    
if strcmp(protocol,'paired_pulse') %if the input variable contains the desired string 
    acq_fullname = file_acq;
    files_xls = dir(fullfile(pathname, '*.xlsx'));
    
    %for loop for identifying the corresponding .xlsx file 
    for k = 1:length(files_xls) 
            
        xlsname = files_xls(k).name; %name of the kth .xlsx file
            
        if contains(acq_fullname(1:end-4),xlsname(1:end-5)) %if .acq file's name 
            %contains .xlsx file's name(which can be ICF,SICI and LICI)
                
            xls_fullname = xlsname; % identification of the correct .xlsx file 
            
        end
        
    end
        
    num = xlsread(fullfile(pathname, xls_fullname)); %matlab built-in function that 
    %reads the data in .xlsx file of the selected path 
        
    trials.IIPP_sec = num(:,2); %adds column in table, with the time intervals 
    %between consecutive pulse pairs (in seconds)

    trials.ISI_sec = num(:,1); %adds column in table, with the time intervals 
    % between TMS artefacts (in seconds)
    
end

if strcmp(protocol,'curve_protocol') %if the input variable contains the desired string
        
    acq_fullname = file_acq; %name of .acq file
        
    xls_fullname = 'key_aquisition.xlsx'; %for this protocol only exists one .xlsx file 
        
    list = xlsread(fullfile(pathname, xls_fullname)); %matlab built-in function that 
    % reads the data in .xlsx file of the selected path   
        
    if contains(acq_fullname, 'groupA') %if .acq file's name contains desired 
        %string name ('groupA')
            
        trials.IIPP_sec = list(:,3); % adds column in table, with time intervals 
        % between pulse pairs (in seconds)
            
        trials.intensities = list(:,1); % adds column in table, with the stimulation 
        %intensities (relative to the resting motor threshold)
        
        
    else %if .acq file's name contains desired string name ('groupB')
            
        trials.IIPP_sec = list(:,4); % adds column in table, with the intervals 
        %between pulse pairs (in seconds)
            
        trials.intensities = list(:,2); % adds column in table, with the 
        %stimulation intensities (relative to the resting motor threshold)
        
    end
    
end


acq_data = load_acq(fullfile(pathname, acq_fullname)); % calls function that converts 
%.acq file to .mat file
    
save([pathname,'\',acq_fullname(1:end-4),'_',date], 'acq_data','trials'); % save 
    %structures in .mat file, in the previously selected path. Name of the file is the 
    %same has the .acq file, adding the date at which was made the conversion
end
