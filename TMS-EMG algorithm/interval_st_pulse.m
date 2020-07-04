function [interval_st_artefact,peak] = interval_st_pulse(isi,mean_pattern,...
    trial_data,corr_percentange,distance,st_pulse,peak_min,trial,artefact_detection,...
    sampling_frequency)
%% 
% Function description: Detects the two peaks of the 1st TMS artefact in
% the EMG segment (trial_data)

%input variables:
% - isi: isi value of the respetive EMG segment (= 0 for the I/O protocol) 
% - mean_pattern: TMS template relative to the isi value (always the baseline for the I/O protocol)
% - trial_data: EMG segment (contains TMS artefact and the subsequent MEP)
% - corr_percentange: percentage of the maximum correlation between respective template 
%(mean_pattern) and EMG segment 
% - distance: distance between the 2 peaks in baseline pattern (estimate of
% the distance between the two peaks to be detected) 
% - st_pulse: identifies if it is an estimate of the time of effective
% acquisition (= 1 if it is, = 0 if it isn't) 
% - peak_min: estimate of the minimum amplitude of the 2 peaks being
% detected (= 0 if there isn't a pre-estimate)
% - trial: identifies the trials/segments that will be accounted to plots being 
%drawn (= 0 doesn't accounts, ~= 0 accounts)
% - artefact_detection: select if desired to see the 2nd peaks of the 1st TMS 
%artefact in each rectified EMG segment (0 for no, 1 for yes)

%output variables:
% - interval_st_artefact: location (in samples, and relative to the respective segment's
% length) of the 1st peak of the TMS artefact 
% - peak: amplitude of the 2nd peak of the TMS artefact, being used for the
% detection of outliers 
%%

[correlation,lags] = xcorr(trial_data,mean_pattern); % correlation between respective template 
%(mean_pattern) and EMG segment (trial_data)

i_maxcorr = find(abs(correlation) > corr_percentange*max(abs(correlation)),1,'first');
% index of the first point in the EMG segment being higher than corr_percentage times the maximum 
%of the correlation 

constant1 = (0.03*sampling_frequency); 
constant2 = (0.05*sampling_frequency); %constant to define, relative to the i_maxcorr, the upper limit of 
%the interval where is going to be searched the TMS artefact
inter =constant1+isi; %constant to define, relative to the i_maxcorr, the lower limit of 
%the interval where is going to be searched the TMS artefact 

st_artefact = ceil(lags(i_maxcorr)-inter); %lower limit of 
%the interval where is going to be searched the TMS artefact 
st_artefact_up = lags(i_maxcorr)+constant2; %upper limit of 
%the interval where is going to be searched the TMS artefact

%while the limits of the interval where is going to be searched the TMS
%artefact aren't correct
while (st_artefact<= 0) || (st_artefact_up >length(trial_data))
    if st_pulse == 1 %if it corresponds to an estimate of the time of effective acquisition 
        corr_percentange = corr_percentange + 0.05; %increment the percentage 
        if corr_percentange<=1 %percentage has to be contained in  0<corr_percentange<=1
            i_maxcorr = find(correlation > corr_percentange*max(correlation),1,'first');
            % index of the first point in the EMG segment being higher than corr_percentage times 
            %the maximum of the correlation 

            st_artefact = ceil(lags(i_maxcorr)-inter); %lower limit of 
            %the interval where is going to be searched the TMS artefact
            st_artefact_up = lags(i_maxcorr)+constant2; %upper limit of 
            %the interval where is going to be searched the TMS artefact
        else 
            st_artefact = 1; %lower limit of 
            %the interval where is going to be searched the TMS artefact =
            %lower limit of the EMG segment trial_data
            st_artefact_up = length(trial_data); %upper limit of 
            %the interval where is going to be searched the TMS artefact =
            %upper limit of the EMG segment trial_data
            break 
        end
    else
        st_artefact = 1; %lower limit of 
            %the interval where is going to be searched the TMS artefact =
            %lower limit of the EMG segment trial_data
        st_artefact_up = lags(i_maxcorr)+constant2; %upper limit of 
            %the interval where is going to be searched the TMS artefact
        if (st_artefact_up >length(trial_data))
            st_artefact_up = length(trial_data); %upper limit of 
            %the interval where is going to be searched the TMS artefact =
            %upper limit of the EMG segment trial_data
        end
    end
end

if peak_min == 0 % there isn't a pre-estimate of the minimum amplitude of the 2 peaks being
% detected
    minimum = 3.5*rms(abs(trial_data(1:st_artefact))); %estimate of the noise magnitude, 
    %being 3.5* RMS (root mean square) of the rectified EMG signal before the lower limit 
    %of the interval where is going to be searched the TMS artefact
else 
    minimum = peak_min; % there is a pre-estimate of the minimum amplitude of the 2 peaks being
% detected (after the detection of outliers) 
end

x = st_artefact:st_artefact_up; % interval where is going to be searched the TMS artefact
y = abs(trial_data(x)); %rectified EMG signal contained in x  

[p,LOCS] = findpeaks(y,'NPeaks',2,...
    'WidthReference','halfheight','MinPeakHeight',minimum,'MinPeakDistance',...
    0.5*distance); %search the two peaks in y, where p are its amplitude and LOCS its locations 
%(in samples, and relative to the length of y) 

condition = 1; 
iteration = 0; 
while condition %while it wasn't detected two peaks that satisfies the given characteristics 
    if numel(LOCS) == 2 %if it are detected 2 peaks 
        if abs(LOCS(1)-LOCS(2))>distance-(0.5*distance) && ...
                abs(LOCS(1)-LOCS(2))< 2.5*distance && ...
                round(p(1),4)>=round(max(y(LOCS(1):LOCS(1)+5)),4) && ...
                iteration ~=2
            %if the amplitude of the 1st peak is the maximum relative to
            %its next 5 neighbors, and if the distance between 1st and 2nd
            %peak is contained in a given interval 
            if round(p(2),4)>=round(max(y(LOCS(2):LOCS(2)+5)),4)
                %if the amplitude of the 2nd peak is the maximum relative to
                %its next 5 neighbors
                condition = 0; %there were found the two peaks that satisfy the given 
                %characteristics  
                
                if artefact_detection %if desired to see the 2nd peaks of the 1st TMS 
                    %artefact in each rectified EMG segment
                    if trial ~=0 %trials/segments that will be accounted to plots being 
                        %drawn 
                        figure(trial+2)
                        plot(y); hold on; %plot EMG signal contained in interval x 
                        plot([LOCS(1) LOCS(1)],[min(y),max(y)],'-b');  hold on; %1st peak detected
                        plot([LOCS(2) LOCS(2)],[min(y),max(y)],'-r'); hold on; %2nd peak detected
                        plot(LOCS(2),p(2),'bo');
                        legend('EMG segment rectified','1st Peak detected','2nd peak detected');
                        title(['Detection of the 1st TMS artefact in the', num2str(trial),...
                            'th rectified EMG segment - Before Optimization Algorithm']);
                        xlabel('Samples');
                        ylabel('mV');
                    end
                end
            else %if the amplitude of the 2nd peak ins't the maximum relative to
                %its next 5 neighbors
                iteration = iteration +1;
                [p,LOCS] = findpeaks(y,...
                    'NPeaks',2,'WidthReference','halfheight','MinPeakHeight',...
                    p(2),'MinPeakDistance',0.5*distance); %search the two peaks in y, establishing 
                %the minimum amplitude as the previous amplitude of the 2nd peak detected 
            end
            
        else
            iteration = 0;
            if st_artefact+LOCS(1)+1 < st_artefact_up %if the lower limit (after incrementation 
                %of one sample relative to the location of the 1st peak) of the interval being 
                %searched is lower than the upper limit 
                st_artefact = st_artefact + LOCS(1) + 1; %incrementation of one sample relative 
                %to the location of the 1st peak
                x = st_artefact:st_artefact_up; % interval where is going to be searched the 
                %TMS artefact
                y = abs(trial_data(x)); %rectified EMG signal contained in x  
                [p,LOCS] = findpeaks(y,...
                    'NPeaks',2,'WidthReference','halfheight','MinPeakHeight',...
                    minimum,'MinPeakDistance',0.5*distance); %search the two peaks in y, where 
                %p are its amplitude and LOCS its locations (in samples, and relative to the 
                %length of y) 
                
            else %if the lower limit (after incrementation of one sample relative to the location 
                %of the 1st peak) of the interval being searched is higher than the upper limit, 
                %the entire interval has already been searched  
                break
            end
        end
    else %if not detected 2 peaks, it never will be, even in posterior iterations 
        break 
    end

end

if numel(LOCS) == 1 || isempty(LOCS) % only detected 1 or 0 peaks
    % function output identificative of the failure in detect TMS artefact 
    interval_st_artefact = 0; 
    peak = 0;
else
    if st_pulse ==1 %if it corresponds to an estimate of the time of effective acquisition
        interval_st_artefact = st_artefact+LOCS(1)+sampling_frequency-1; % location (in samples, and relative 
        % to the respective segment's length) of the 1st peak of the TMS artefact
    else %if it doesn't corresponds to an estimate of the time of effective acquisition
        interval_st_artefact = st_artefact+LOCS(1)-1; % location (in samples, and relative 
        % to the respective segment's length) of the 1st peak of the TMS artefact
    end
    peak = p(2); %amplitude of the 2nd peak of the TMS artefact, being used for the
% detection of outliers 
end
end