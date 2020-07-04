function data_convertion_segmentation_curve(pathname,quality_comb,artefact_detection,...
    distribution_tms,comb,files_mat,sampling_frequency)
%% 
% Function description: Segments and detects the location, in each segment, of the 1st peak
% of TMS artefact, for the I/O protocol

%input variables:
% - pathname: string with the path's name of the selected directory (that contains the files 
%to be analyzed) 
% - quality_comb: quality factor of IIR comb filter (nondimensional)
% - artefact_detection: select if desired to see the 2nd peaks of the 1st TMS 
%artefact in each rectified EMG segment (0 for no, 1 for yes)
% - distribution_tms: select if desired to see the distribution of amplitude of the 
%detected  TMS artefacts (before and after optimization of the TMS
%detection algorithm): 0 for no, 1 for yes
% - comb: switch applicability of IIR comb filter (0 to not apply, 1 to apply)
% - files_mat: selected files, in 'choosedata' GUI, to be segmented 
%%
global new_segmented_files; %global variable created in 'choosedata.mat' GUI
new_segmented_files = {}; %cell array that will contain the new segmented file's names. 

%for loop that segments the rth .mat files that were selected in GUI for
%segmentation
for r=1:length(files_mat)

    file = files_mat{r}; %name of the rth .mat file contained in the cell array 
    %'files_mat'
    str = load([pathname,'\',file]); %loads the corresponding .mat file data 
    
    acq_data = str.acq_data; %get the acq_data structure, obtained by load_acq 
    %function 
    trials = str.trials; %get the table named trials, that will be filled with 
    %intended parameters
    
    acq_data_detrend = detrend(acq_data.data(:,1)); %detrend data from channel 1 of the EMG acquisition
    
    if comb % IIR comb filter 
        f0 = 50; % 1st frequency 
        q = quality_comb; %quality factor 
        bw = (f0/(sampling_frequency/2))/q; % sets bandwidth
        [b,a] = iircomb(sampling_frequency/f0,bw); %built-in function that eliminates the 50 Hz frequency and 
        %its harmonics
        [sos,g] = tf2sos(b,a); % implementation of second order sections 
        acq_data_detrend = filtfilt(sos,g,acq_data_detrend); % filters forwards and backwards to 
        %avoid signal distortion
    end
    
    pattern_baselines = load('pattern_baselines.mat'); % load the baseline 
    %pattern 
    pattern_baselines = pattern_baselines.results.mean_pattern; % get from structure 
    %the baseline pattern  

    [~,imax] = max(pattern_baselines); 
    [~,imin] = min(pattern_baselines);
    
    %input variables for interval_st_pulse function:
    distance = abs(imax-imin); %get distance between the 2 peaks in baseline pattern 
    isi = 0;
    corr_percentange = 0.2; %percentage of the maximum correlation between baseline template and 
    %EMG segment 
    st_pulse = 1; % identifies that is an estimate of the time of effective acquisition 
    trial = 0; %it won't be accounted to plots being drawn in interval_st_pulse function 
    peak_min = 0; %it hasn't an estimate of the minimum amplitude of the 2 peaks being detected 
    
    %localizes the time of effective acquisition automatically 
    [st_tms,~] = interval_st_pulse(isi,pattern_baselines,...
        acq_data_detrend(sampling_frequency:25*sampling_frequency),corr_percentange,distance,...
        st_pulse,peak_min,trial,artefact_detection,sampling_frequency);

    if st_tms == 0 % if it wasn't automatically found the time of effective acquisition
        acq_data.x_start = 0; 
        st_tms = time_start(acq_data,file,sampling_frequency); % calls GUI that manually permits the selection 
        %of the time of effective acquisition
    end
    
    time_delta_index1 = 0; %initializes the lower limit of the segment (in samples)
    constant = (0.04 * sampling_frequency);
    
    %for loop that will create segments of data that contains the TMS artefact 
    %and the corresponding MEP 
    for j = 1:height(trials)
        if j ~= height(trials) % this segments contains IIPP time intervals
            time_delta_index2 = fix((sum(trials.IIPP_sec(1:j,1))*sampling_frequency)); 
            %upper limit of the segment (in samples)
            
            % EMG data segment 
            trial_data = acq_data_detrend((st_tms+time_delta_index1-constant):...
                (st_tms+time_delta_index2-constant),1); 
            
            time_delta_index1 = fix((sum(trials.IIPP_sec(1:j,1))*sampling_frequency)); 
            %atualizes the lower limit for the next segment
            
        else %for the last segment, it doesn't exists IIPP time interval 
            
            % EMG data segment
            trial_data = acq_data_detrend((st_tms+time_delta_index1-constant):end,1);
             
        end
        
        % adds column in trials table, corresponding to jth EMG data segment 
        trials.EMGfilt{j,1} = trial_data;
        
        corr_percentange = 0.7; %percentage of the maximum correlation between baseline template and 
    %EMG segment 
        isi = 0;
        trial = j; %it will be accounted to plots being drawn in interval_st_pulse function
        peak_min = 0; %it hasn't an estimate of the minimum amplitude of the 2 peaks being detected 
        st_pulse = 0; % it doesn't correspond to an estimate of the time of effective acquisition 
        [interval_st_artefact,peak] = interval_st_pulse(isi,...
            pattern_baselines,trial_data,corr_percentange,distance,st_pulse,...
            peak_min,trial,artefact_detection,sampling_frequency);
        
        trials.interval_pulse{j,1} = interval_st_artefact; % location (in samples) of 1st peak of the 
        %detected TMS artefact 
        trials.peak(j,1) = peak; % amplitude of 2nd peak of the detected TMS artefact  

    end
    
    intensity_values = unique(trials.intensities(:,1)); %get the intensity values on the randomized 
    %key (i.e, 6 intensities: 90,100,110,120,130 and 140)
        
    % for loop that runs for each intensity value 
    for j=1:numel(intensity_values)
        index_intensities = find(trials.intensities == intensity_values(j)); %indexes in randomized key 
        %corresponding to the jth intensity value 
        peaks = trials.peak(index_intensities,1); %amplitude of 2nd TMS peak for the corresponding 
        %trials, for the jth intensity value 
        [wrong_peaks,L,U,C] = isoutlier(peaks,'ThresholdFactor',0.75); %identification of outliers, 
        %with the MAD criteria and an 0.75 threshold, with:
        %L = median(peaks) - 0.75*MAD; U = median(peaks) + 0.75*MAD; C = median(peaks)
        
        index_wrong = find(wrong_peaks ==1); % indexes of the identified outliers 
        
        if ~isempty(index_wrong) % if it were identified outliers 
            index_correct = index_intensities(find(wrong_peaks~=1)); %corresponding trials that weren't 
            %identified as outliers, for the jth intensity value 
            peak_min = 0.5*min(trials.peak(index_correct,1)); % minimum amplitude of the TMS artefacts 
            %peaks (=2) being detected
            
            %for the identified outliers, redo the identification of the 1st peak of TMS artefact 
            for l = 1:numel(index_wrong)
                n = index_intensities(index_wrong(l)); %trial 
                
                corr_percentange = 0.7; %percentage of the maximum correlation between baseline template and 
                %EMG segment 
                isi = 0;
                st_pulse = 0; % it doesn't correspond to an estimate of the time of effective acquisition 
                trial = 0; %it won't be accounted to plots being drawn in interval_st_pulse function 
                
                [interval_st_artefact,peak] = interval_st_pulse(isi,pattern_baselines,...
                    trials.EMGfilt{n,1},corr_percentange,distance,st_pulse,peak_min,...
                    trial,artefact_detection,sampling_frequency);
                
                trials.peak(n,1) = peak; % amplitude of 2nd peak of the detected TMS artefact  
                trials.interval_pulse{n,1} = interval_st_artefact; % location (in samples) of 1st peak of the 
        %detected TMS artefact 
            end
        end
        
        if distribution_tms %if desired to see the distribution of amplitude of the 
%detected  TMS artefacts (before and after optimization of the TMS
%detection algorithm)
            vector_inten = 1:length(index_intensities);
            figure(j) %for each intensity value create a figure
            subplot(2,1,1), plot(vector_inten,peaks,'b*',vector_inten(wrong_peaks),...
        peaks(wrong_peaks),'ro',vector_inten,L*ones(1,length(index_intensities)),...
        '--',vector_inten,U*ones(1,length(index_intensities)),'--',vector_inten,...
        C*ones(1,length(index_intensities)),'--');
            legend('Original Data','Outlier',['Median - 0.75*MAD =',num2str(L)],...
        ['Median + 0.75*MAD = ',num2str(U)],['Median = ',num2str(C)],...
        'Location','southeast');
            title(['Peak Amplitude of TMS Artefacts -',num2str(intensity_values(j)),'%RMT']);
            ylabel('Peak Amplitude (mV)');
            xlabel('Samples');
            ylim([-0.1 max(peaks)+0.05]);
            xlim([1 max(vector_inten)]);
        end
            
    end
    
    if distribution_tms %if desired to see the distribution of amplitude of the 
%detected  TMS artefacts (before and after optimization of the TMS
%detection algorithm)

        % for loop that runs for each intensity value 
        for j=1:numel(intensity_values)
            index_intensities = find(trials.intensities == intensity_values(j)); %indexes in randomized key 
        %corresponding to the jth intensity value 
            peaks = trials.peak(index_intensities,1); %amplitude of 2nd TMS peak for the corresponding 
        %trials, for the jth intensity value 
            [wrong_peaks,L,U,C] = isoutlier(peaks,'ThresholdFactor',0.75); %identification of outliers, 
        %with the MAD criteria and an 0.75 threshold, with:
        %L = median(peaks) - 0.75*MAD; U = median(peaks) + 0.75*MAD; C = median(peaks)
            
            vector_inten = 1:length(index_intensities);
            figure(j) %for each intensity value create a figure
            subplot(2,1,2), plot(vector_inten,peaks,'b*',vector_inten(wrong_peaks),...
        peaks(wrong_peaks),'ro',vector_inten,L*ones(1,length(index_intensities)),...
        '--',vector_inten,U*ones(1,length(index_intensities)),'--',vector_inten,...
        C*ones(1,length(index_intensities)),'--');
            legend('Corrected data','Remaning Outlier',['Median - 0.75*MAD =',num2str(L)],...
        ['Median + 0.75*MAD = ',num2str(U)],['Median = ',num2str(C)],...
        'Location','southeast');
            ylabel('Peak Amplitude (mV)');
            xlabel('Samples');
            ylim([-0.1 max(peaks)+0.05]);
            xlim([1 max(vector_inten)]);
        end
    end
    new_segmented_files{end+1}=[file(1:end-4),'_segmented.mat']; % adds to global 
    %variable (cell array) the name of the new segmented file
    
    save([pathname,'\',new_segmented_files{end}], 'trials'); % save structures in 
    %.mat file, in the previously selected path
    
end
end