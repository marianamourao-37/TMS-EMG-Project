function final_parameters_curve(pathname,min_amp_MEP, latency_MEP,curve_model,files_segmented,...
    sampling_frequency)

%% 
% Function description: analyzes the EMG segments (MEP detection and quantification) for the 
% I/O protocol. Modulates the I/O curve 

%input variables:
% - pathname: string with the path's name of the selected directory (that contains the files 
%to be analyzed) 
% - min_amp_MEP: minimum MEP's amplitude (mV) for its quantification
% - latency_MEP: MEP's latency (ms) relative to the last TMS artefact 
%(used to define MEP's window for its detection)
% - curve_model: Type of I/O curve Modulation (1 for statistical modulation 
% and 0 for linear regression)
% - files_segmented: selected files, in 'choosedata' GUI, to be analyzed 
%%

global new_analysed_files; %global variable created in 'choosedata.mat' GUI
new_analysed_files = {}; %cell array that will contain the new analysed file's names. 

%for loop that analyzes the ith .mat files that were selected in GUI
for i=1:length(files_segmented)
    
    clearvars trials results mean;
    
    file = files_segmented{i}; %name of the ith .mat file contained in the cell array 
    %files_segmented
    load([pathname,'\',file]); %loads the data from ith .mat file
    
    trials.MEP_onset(:,1) = 0; %initializes lower limit of MEP (in samples)
    
    trials.MEP_offset(:,1) = 0; % initializes upper limit of MEP (in samples) 
    
    %% trials loop
    
    intensities_values = unique(trials.intensities(:,1)); % get the intensities 
    %values applied during input-output protocol
    
    %for loop that analyzes each EMG segment
    for i = 1:height(trials)
        trials.artloc(i,1) = trials.interval_pulse{i,1}; %location of 1st TMS artefact in ith 
        %EMG segment 
        
        if trials.artloc(i,1) ~= 0 % 1st TMS artefact was detected
            
            %set lower and upper limit of MEP search window.
            MEP_onset_window = trials.artloc(i,1) + fix(latency_MEP*10^-3*sampling_frequency);
            MEP_offset_window = trials.artloc(i,1) + fix(0.04*sampling_frequency);
        
            % get data from MEP search interval
            MEPsearch = trials.EMGfilt{i,1}(MEP_onset_window:MEP_offset_window); 
        
            %get the desired points (p.e = 10) that represent the maximum
            %changes of a certain parameter (p.e mean) 
            ipoints = findchangepts(MEPsearch, 'MaxNumChanges', 10, 'Statistic', 'mean');
        
            trials.MEP_onset(i,1) = ipoints(1); %set the lower limit of MEP 
           
            trials.MEP_offset(i,1) = ipoints(end); %set the upper limit of MEP
        
            trials.MEP_onset(i,1) = fix(MEP_onset_window + trials.MEP_onset(i,1));
            %correction of the lower limit of MEP, to be relative to the start 
            %of the respective segment
        
            trials.MEP_offset(i,1) = fix(MEP_onset_window + trials.MEP_offset(i,1));
            %correction of the upper limit of MEP, to be relative to the start of 
            %the respective segment
        
            % MEP data 
            MEPsearchrange = trials.EMGfilt{i,1}(trials.MEP_onset(i,1):...
                trials.MEP_offset(i,1));
        
            %get the maximum and minimum values in MEP data
            [max_MEP_value,~] = max(MEPsearchrange);
            [min_MEP_value,~] = min(MEPsearchrange);
            
            trials.MEP_amplitude(i,1) = max_MEP_value - min_MEP_value; %MEP amplitude

            if trials.MEP_amplitude(i,1) < min_amp_MEP %if the MEP amplitude is higher than 
                %the required value 
                trials.MEP_amplitude(i,1) = 0; %it won't be accounted for posterior calculations
                trials.MEP_onset(i,1) = 0;
                trials.MEP_offset(i,1) = 0;
            end
        else %1st TMS artefact wasn't detected 
            trials.MEP_amplitude(i,1) = 0; %it won't be accounted for posterior calculations
        end
    end
 
    statistics = table();
    statistics.intensities(:,1) = intensities_values;
    
    %calculate the mean MEP amplitudes for each intensity value 
    for u= 1:length(intensities_values)
        index = find(trials.intensities(:,1) == intensities_values(u));
        mep_values = trials.MEP_amplitude(index,1); %MEP amplitudes relative to the uth intensity value
        mep_values(mep_values == 0) = NaN; %doesn't account MEP amplitudes equal to 0 
        statistics.mean_mep_amplitude(u,1) = mean(mep_values,'omitnan'); %calculate mean
        statistics.sd_mep_amplitude(u,1) = std(mep_values,'omitnan'); %calculate standard deviation
        statistics.cv_mep_amplitude(u,1) = statistics.sd_mep_amplitude(u,1)/...
            statistics.mean_mep_amplitude(u,1); % calculate variation coefficient 
    end
    
    if curve_model == 1 % statistical modulation algorithm of the I/O curve it hasn't yet been 
        %developed 
        warndlg(['It wasn''t yet implemented an alternative modulation besides de default',...
            ' (regression modulation of the I/O curve). Thus, it will be applied the',...
            ' regression model']);
    end
    
    % 4 parameters sigmoid function
    fhSigmo = @(p,x) abs(p(1)) + p(4)./(1+exp(p(2).*(p(3)-x)));
    
    % initial values to calculate the 4 parameters in fhSigmo
    mid_point = statistics.intensities(fix(end/2));
    vBO = [nanmin(statistics.mean_mep_amplitude(:,1)) 0.1 mid_point ...
        max(statistics.mean_mep_amplitude(:,1))];
    
    %nonlinear fit, taking the intensities and MEP values, as well as the initial values 
    %vBO that estimate where the curve passes. vB contains the 4 parameters of fhSigmo
    vB = nlinfit(statistics.intensities(:,1),statistics.mean_mep_amplitude(:,1), ...
        fhSigmo, vBO);

    % slope of sigmoid's linear region, obtained by doing the 1st derivative of
    % the sigmoid function (fhSigmo), and calculating for the midpoint value of 
    % the linear region (parameter in fhSigmo -> vB(3)). assumes that the linear 
    % region has the same 1st derivative value
    results.slope2 = vB(4)*vB(2)/4;
    
    results.intensity_halfmax = vB(3);  % structure that contains field with the half-maximum 
    %intensity, relative to a parameter in fhSigmo that was estimated by nonlinear fit 

    new_analysed_files{end+1}=[file(1:strfind(file,'_segmented')-1),'_analysed.mat']; 
    % adds to global variable (cell array) the name of the new analyzed file
    
    save([pathname,'\',new_analysed_files{end}], 'trials','results','statistics'); 
    % save structures in .mat file,in the previously selected path

end 
end
