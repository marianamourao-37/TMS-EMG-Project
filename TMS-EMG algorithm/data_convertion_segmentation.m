function data_convertion_segmentation(pathname,quality_comb,artefact_detection,...
    distribution_tms,comb,files_mat,sampling_frequency)

%% 
% Function description: Segments and detects the location, in each segment, of the 1st peak
% of TMS artefact, for the paired-pulse protocol

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
    
    str = load([pathname,'\',file]); %loads the corresponding .mat file 
    
    acq_data = str.acq_data; %get the acq_data structure, obtained by load_acq 
    %function 
    trials = str.trials; %get the table named trials, that will be filled with 
    %intended parameters
   
    acq_data_detrend = detrend(acq_data.data(:,1)); %detrend data from channel 1 
    %of the EMG acquisition
    
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
    
    st_pattern = trials.ISI_sec(1,1); 
    isi = fix(st_pattern*10^3); %isi (ms) of the first trial identified in the randomized key 

    if isi == 0 %get the baseline template 
        mean_patterns = load('pattern_baselines.mat');
        mean_patterns = mean_patterns.results.mean_pattern;
    else 
        if contains(file,'SICI') %get the corresponding isi template 
            mean_patterns = load('pattern_SICI.mat');
            mean_patterns = mean_patterns.results_sici.(['mean_pattern',num2str(isi),'ms']);
        elseif contains(file,'ICF') %get the corresponding isi template 
            mean_patterns = load('pattern_ICF.mat');
            mean_patterns = mean_patterns.results_icf.(['mean_pattern',num2str(isi),'ms']);
        elseif contains(file,'LICI') %get the corresponding isi template 
            mean_patterns = load('pattern_LICI.mat');
            mean_patterns = mean_patterns.results_lici.(['mean_pattern',num2str(isi),'ms']);
        end
    end
            
    pattern_baseline = load('pattern_baselines.mat'); % load the baseline 
    %pattern 
    pattern_baseline = pattern_baseline.results.mean_pattern; % get from structure 
    %the baseline pattern  
    [~,imax] = max(pattern_baseline);
    [~,imin] = min(pattern_baseline);
    
    %input variables for interval_st_pulse function:
    distance = abs(imax-imin); %get distance between the 2 peaks in baseline pattern
    corr_percentange = 0.4; %percentage of the maximum correlation between EMG segment and the 
    %respective TMS template 
    peak_min = 0; %it hasn't an estimate of the minimum amplitude of the 2 peaks being detected 
    st_pulse = 1; % identifies that is an estimate of the time of effective acquisition
    trial = 0; %it won't be accounted to plots being drawn in interval_st_pulse function 

    %localizes the time of effective acquisition automatically 
    [st_tms,~] = interval_st_pulse(isi,mean_patterns,...
        acq_data_detrend(sampling_frequency:25*sampling_frequency),corr_percentange,distance,st_pulse,...
        peak_min,trial,artefact_detection,sampling_frequency);
    
    if st_tms == 0 %if it wasn't automatically found the time of effective acquisition
        acq_data.x_start = 0;
        st_tms = time_start(acq_data,file,sampling_frequency); % calls GUI that manually permits the selection 
        %of the time of effective acquisition
    end

    time_delta_index1 = 0; %initializes the lower limit of the segment (in samples)
    
    constant = (0.04 * sampling_frequency);
    
    %for loop that will create segments of data that contains 1st and 2nd TMS
    %artefacts and the corresponding MEP 
    for j = 1:height(trials)
        isi = fix(trials.ISI_sec(j,1)*sampling_frequency);
        if j ~= length(trials.ISI_sec) % this segments contains IIPP time intervals 
            time_delta_index2 = fix((sum(trials.IIPP_sec(1:j,1))+...
                sum(trials.ISI_sec(1:j,1)))*sampling_frequency);
            
            % EMG data segment 
            trial_data = acq_data_detrend((st_tms+time_delta_index1-constant):...
                (st_tms+time_delta_index2-constant),1);
       
            time_delta_index1 = fix((sum(trials.IIPP_sec(1:j,1))+...
                sum(trials.ISI_sec(1:j,1)))*sampling_frequency); %atualizes the lower 
            %limit for the next segment
            
        else %for the last segment, it doesn't exists IIPP time interval 
            trial_data = acq_data_detrend((st_tms+time_delta_index1-constant):end,1);
             
        end
        
        % adds column in trials table, corresponding to jth EMG data segment 
        trials.EMGfilt{j,1} = trial_data;
        
        st_pulse =  0; % it doesn't correspond to an estimate of the time of effective acquisition
        trial = j; %it will be accounted to plots being drawn in interval_st_pulse function
        corr_percentange = 0.8; %percentage of the maximum correlation between baseline template and 
        %EMG segment 
        peak_min = 0; %it hasn't an estimate of the minimum amplitude of the 2 peaks being detected 
        [interval_st_artefact,peak] = interval_st_pulse(isi,mean_patterns,...
            trial_data,corr_percentange,distance,st_pulse,peak_min,trial,artefact_detection,...
            sampling_frequency);
        
        
        trials.interval_pulse{j,1} = interval_st_artefact; % location (in samples) of 1st peak of the 
        %detected TMS artefact 
        trials.peak(j,1) = peak; % amplitude of 2nd peak of the detected TMS artefact 
        
    end
    
    baseline = find(trials.ISI_sec == 0); %indexes in randomized key corresponding to isi = 0 
    other = find(trials.ISI_sec ~=0); %indexes in randomized key not corresponding to isi = 0 
    
    p_b = trials.peak(baseline,1); %amplitude of 2nd TMS peak for the corresponding 
        %trials, for isi = 0 (supposed to have the same value) 
    p_o = trials.peak(other,1); %amplitude of 2nd TMS peak for the corresponding 
        %trials, for isi ~= 0 (supposed to have the same value) 
    
    [wrong_peak_baseline,L_b,U_b,C_b] = isoutlier(p_b,'ThresholdFactor',0.75); %identification of 
    %outliers for isi = 0, with the MAD criteria and an 0.75 threshold
    [wrong_peak_other,L_o,U_o,C_o] = isoutlier(p_o,'ThresholdFactor',0.75); %identification of 
    %outliers for isi ~= 0, with the MAD criteria and an 0.75 threshold
    
    % L = median(peaks) - 0.75*MAD; U = median(peaks) + 0.75*MAD; C = median(peaks)
    
    vector_b = 1:length(baseline);
    vector_o = 1:length(other);
    
    %%
    if distribution_tms %if desired to see the distribution of amplitude of the 
%detected  TMS artefacts (before and after optimization of the TMS
%detection algorithm)
        figure(1)
        subplot(2,1,1), plot(vector_b,p_b,'b*',vector_b(wrong_peak_baseline),...
        p_b(wrong_peak_baseline),'ro',vector_b,L_b*ones(1,length(baseline)),...
        '--',vector_b,U_b*ones(1,length(baseline)),'--',vector_b,...
        C_b*ones(1,length(baseline)),'--');
        legend('Original Data','Outlier',['Median - 0.75*MAD =',num2str(L_b)],...
        ['Median + 0.75*MAD = ',num2str(U_b)],['Median = ',num2str(C_b)],...
        'Location','southeast')
        title('Peak Amplitude of TMS Artefacts - Baseline');
        ylabel('Peak Amplitude (mV)');
        xlabel('Samples');
        ylim([-0.1 max(p_b)+0.05]);
        xlim([1 max(vector_b)]);
        subplot(2,1,2), plot(vector_o,p_o,'b*',vector_o(wrong_peak_other),...
        p_o(wrong_peak_other),'ro',vector_o,L_o*ones(1,length(other)),'--',...
        vector_o,U_o*ones(1,length(other)),'--',vector_o,...
        C_o*ones(1,length(other)),'--');
        legend('Original Data','Outlier',['Median - 0.75*MAD =',num2str(L_o)],...
        ['Median + 0.75*MAD = ',num2str(U_o)],['Median = ',num2str(C_o)],...
        'Location','southeast')
        title('Peak Amplitude of 1st TMS Artefacts - ISIs');
        ylabel('Peak Amplitude (mV)');
        xlabel('Samples');
        ylim([-0.1 max(p_o)+0.05]);
        xlim([1 max(vector_o)]);
    end
    %%
    
    index_baseline_wrong = find(wrong_peak_baseline ==1); % indexes of the identified outliers, for isi=0
    index_other_wrong = find(wrong_peak_other ==1); % indexes of the identified outliers, for isi~=0
    
    if ~isempty(index_baseline_wrong) % if it were identified outliers 
        index_baseline_correct = baseline(find(wrong_peak_baseline ~=1)); %corresponding trials that 
        %weren't identified as outliers, for the isi = 0  
        index_other_correct = other(find(wrong_peak_other ~=1)); %corresponding trials that 
        %weren't identified as outliers, for the isi ~= 0 
        peak_min_baseline = 0.5*min(trials.peak(index_baseline_correct,1)); % minimum amplitude of the 
        %TMS artefacts peaks (=2) being detected, for isi = 0
        peak_min_other = 0.5*min(trials.peak(index_other_correct,1)); % minimum amplitude of the TMS 
        %artefacts peaks (=2) being detected, for isi ~=0
        index = [index_baseline_wrong;index_other_wrong];
        
        %for the identified outliers, redo the identification of the 1st peak of TMS artefact 
        for c = 1:numel(index)
            if c<=length(index_baseline_wrong) %isi = 0 
                n = baseline(index(c)); %trial/Segment 
                isi = 0;
                corr_percentange = 0.8; %percentage of the maximum correlation between EMG segment and 
                %respective TMS template
                st_pulse = 0; % it doesn't correspond to an estimate of the time of effective acquisition 
                trial = 0; %it won't be accounted to plots being drawn in interval_st_pulse function
                [interval_st_artefact,peak] = interval_st_pulse(isi,mean_patterns,...
                    trials.EMGfilt{n,1},corr_percentange,distance,st_pulse,...
                    peak_min_baseline,trial,artefact_detection,sampling_frequency);
            else %isi ~=0
                n = other(index(c)); %trial/Segment 
                isi = trials.ISI_sec(n,1);
                corr_percentange = 0.8; %percentage of the maximum correlation between EMG segment and 
                %respective TMS template
                st_pulse = 0; % it doesn't correspond to an estimate of the time of effective acquisition
                trial = 0; %it won't be accounted to plots being drawn in interval_st_pulse function
                [interval_st_artefact,peak] = interval_st_pulse(isi,...
                    mean_patterns,trials.EMGfilt{n,1},corr_percentange,distance,st_pulse,...
                    peak_min_other,trial,artefact_detection,sampling_frequency);
            end
            trials.peak(n,1) = peak; % amplitude of 2nd peak of the detected TMS artefact
            trials.interval_pulse{n,1} = interval_st_artefact; % location (in samples) of 1st peak of the 
        %detected TMS artefact 
        end
    end
    %%
    if distribution_tms %if desired to see the distribution of amplitude of the 
%detected  TMS artefacts (before and after optimization of the TMS
%detection algorithm)
        figure(2) 
    
        p_b = trials.peak(baseline,1); %amplitude of 2nd TMS peak for the corresponding 
        %trials, for isi = 0 (supposed to have the same value)
        p_o = trials.peak(other,1); %amplitude of 2nd TMS peak for the corresponding 
        %trials, for isi ~= 0 (supposed to have the same value)
    
        [wrong_peak_baseline,L_b,U_b,C_b] = isoutlier(trials.peak(baseline,1),'ThresholdFactor',0.75);
        %identification of outliers for isi = 0, with the MAD criteria and an 0.75 threshold
        [wrong_peak_other,L_o,U_o,C_o] = isoutlier(trials.peak(other,1),'ThresholdFactor',0.75);
        %identification of outliers for isi ~= 0, with the MAD criteria and an 0.75 threshold
        
        %L = median(peaks) - 0.75*MAD; U = median(peaks) + 0.75*MAD; C = median(peaks)
        
        subplot(2,1,1), plot(vector_b,p_b,'b*',vector_b(wrong_peak_baseline),...
        p_b(wrong_peak_baseline),'ro',vector_b,L_b*ones(1,length(baseline)),...
        '--',vector_b,U_b*ones(1,length(baseline)),'--',vector_b,...
        C_b*ones(1,length(baseline)),'--');
        title('Peak Amplitude of TMS Artefacts - Baseline');
        ylabel('Peak Amplitude (mV)');
        xlabel('Samples');
        ylim([-0.1 max(p_b)+0.05]);
        xlim([1 max(vector_b)]);
        legend('Corrected data','Remaning Outlier',['Median - 0.75*MAD =',num2str(L_b)],...
        ['Median + 0.75*MAD = ',num2str(U_b)],['Median = ',num2str(C_b)],'Location',...
        'southeast');
        subplot(2,1,2), plot(vector_o,p_o,'b*',vector_o(wrong_peak_other),...
        p_o(wrong_peak_other),'ro',vector_o,L_o*ones(1,length(other)),'--',...
        vector_o,U_o*ones(1,length(other)),'--',vector_o,...
        C_o*ones(1,length(other)),'--');
        title('Peak Amplitude of 1st TMS Artefacts - ISIs');
        ylabel('Peak Amplitude (mV)');
        xlabel('Samples');
        ylim([-0.1 max(p_o)+0.05]);
        xlim([1 max(vector_o)]);
        legend('Corrected data','Remaning Outlier',['Median - 0.75*MAD =',num2str(L_o)],...
        ['Median + 0.75*MAD = ',num2str(U_o)],['Median = ',num2str(C_o)],...
        'Location','southeast');
    end
    
    new_segmented_files{end+1}=[file(1:end-4),'_segmented.mat']; % adds to global variable 
    %(cell array) the name of the new ith segmented file
    
    save([pathname,'\',new_segmented_files{end}], 'trials'); % save structures in .mat 
    %file, in the previously selected path. 
    
end
end