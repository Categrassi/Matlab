%[metadata,edf_data] = edfread3('L:\LovbeskyttetMapper01\Sleep Microstructure\PSG Files\temp edf files\CJN.edf'); %super slow with two otputs
clear all

[metadata] = edfread3('L:\LovbeskyttetMapper01\Sleep Microstructure\PSG Files\edf KU\11927.edf');%only meta data
[edf_data,edf_header] = lab_read_edf2('L:\LovbeskyttetMapper01\Sleep Microstructure\PSG Files\edf KU\11927.edf'); % cell structure with data - no context
[edf_struct, manual_annotation] = edfread('L:\LovbeskyttetMapper01\Sleep Microstructure\PSG Files\edf KU\11927.edf');

data_clean = ('patientID');
data = rmfield(metadata,data_clean);

%% start code

t2b_removed_h = 5.5; % time to be removed in h
t2b_removed_s = t2b_removed_h*60*60;

scoring = edf_data {22};
scoring_cut = scoring(t2b_removed_s:end);
edf_dur = length(scoring_cut);
edf_length = length(edf_data{1});

edf_C3 = edf_data {18};
edf_dur_C3 = size(edf_struct, 1);
edf_dur_C3_cut = edf_dur_C3-t2b_removed_s;

sampling_freq = edf_header.hdr.numbersperrecord(18)/edf_header.hdr.duration; %this is sampling frequency of C3 channel
C3_timetrace = (0:(edf_length-1))/sampling_freq;
C3_fs = length(C3_timetrace)/edf_dur_C3; % sampling freq of C3
t2b_removed_samples = t2b_removed_s*C3_fs; %samples to be removed
%edf_cut_dur_C3 = (edf_dur_C3-t2b_removed_s)*C3_fs;
C3_timetrace_cut = C3_timetrace(t2b_removed_samples:end);
edf_C3_cut = edf_C3(t2b_removed_samples:end);


sampling_frequency_hypno = edf_header.hdr.numbersperrecord(22)/edf_header.hdr.duration; %this is sampling frequency of the hypnogram 
hypnogram_timetrace = 0:(edf_dur_C3-1)/sampling_frequency_hypno;
hypnogram_timetrace_cut = hypnogram_timetrace(t2b_removed_s:end);

figure;
a = subplot(2,1,1);
plot(hypnogram_timetrace_cut, scoring_cut);
title('hipnogram');
    xlabel('time');
    ylabel('sleep phase');
b = subplot(2,1,2);
plot(C3_timetrace_cut, edf_C3_cut);
title('C3');
    xlabel('time')
    ylabel('freq (Hz)');
linkaxes([a,b],'x');

%% here we create binary vectors for each sleep state (wake, REM, NREM, N1, N2 and N3)

sleepstates = edf_struct{:, 22};
sleepstates_cut = sleepstates(t2b_removed_s:end);

wakeVector = sleepstates_cut == 0;
REMVector = sleepstates_cut == 5;
NREMVector = ismember(sleepstates_cut, [1, 2, 3]);
N1Vector = sleepstates_cut == 1;
N2Vector = sleepstates_cut == 2;
N3Vector = sleepstates_cut == 3;

disp('Wake Vector:');
disp(wakeVector);

disp('REM Vector:');
disp(REMVector);

disp('N1 Vector:');
disp(N1Vector);

disp('N2 Vector:');
disp(N2Vector);

disp('N3 Vector:');
disp(N3Vector);

disp('NREM Vector:');
disp(NREMVector);

[wake_onset, wake_offset] = binary_to_OnOff(wakeVector');
wake_periods = [wake_onset, wake_offset]+t2b_removed_s; % we + seconds so we shift the values since we are cutting the trace 
[sws_onset, sws_offset] = binary_to_OnOff(NREMVector');
sws_periods = [sws_onset, sws_offset]+t2b_removed_s;
[REM_onset, REM_offset] = binary_to_OnOff(REMVector');
REM_periods = [REM_onset, REM_offset]+t2b_removed_s;
[N1_onset, N1_offset] = binary_to_OnOff(N1Vector');
N1_periods = [N1_onset, N1_offset]+t2b_removed_s;
[N2_onset, N2_offset] = binary_to_OnOff(N2Vector');
N2_periods = [N2_onset, N2_offset]+t2b_removed_s;
[N3_onset, N3_offset] = binary_to_OnOff(N3Vector');
N3_periods = [N3_onset, N3_offset]+t2b_removed_s;

wake_duration = wake_offset-wake_onset;
duration_sws = sws_offset-sws_onset;
REM_duration = REM_offset-REM_onset;
N1_duration = N1_offset-N1_onset;
N2_duration = N2_offset-N2_onset;
N3_duration = N3_offset-N3_onset;

C3_timetrace_cut_ds = downsample(C3_timetrace_cut, 100);
edf_C3_cut_ds = downsample(edf_C3_cut, 100);

hypnogram_timetrace_ds = downsample(hypnogram_timetrace,50);
hypnogram_timetrace_cut_ds = downsample(hypnogram_timetrace_cut,50);
wakeVector_ds = downsample(wakeVector,50);
NREMVector_ds = downsample(NREMVector,50);
REMVector_ds = downsample(REMVector,50);
N1Vector_ds = downsample(N1Vector, 50);
N2Vector_ds = downsample(N2Vector, 50);
N3Vector_ds = downsample(N3Vector, 50);

fig = figure; % graph of EEG C3 with NREM phase all in the same color 
    plot_sleep(C3_timetrace_cut_ds, edf_C3_cut_ds, hypnogram_timetrace_cut_ds, wakeVector_ds', REMVector_ds', NREMVector_ds');
    xlabel('time (s)');
    ylabel('EEG (V)');

fig = figure;  % graph of EEG C3  with different colors for each of the NREM phases 
    plot_sleep(C3_timetrace_cut_ds, edf_C3_cut_ds, hypnogram_timetrace_cut_ds, wakeVector_ds', REMVector_ds', N1Vector_ds', N2Vector_ds', N3Vector_ds');
    xlabel('time (s)');
    ylabel('EEG (V)');
    legend('wake', 'REM', 'N1','N2','N3')

%% EEG power spectrum analysis

% here you plot a power spectrogram (heatmap) for your EEG data and band
% power traces for the defined EEG frequency bands (incl sigma trace).
 

Data_EEG = edf_C3_cut; % should be a vector containing EEG data to perform analysis on
power_bands = {[1, 4], [4, 10], [10, 15], [15, 30], [30, 50], [50, 75]}; % define delta, theta, sigma, beta, and gamma low and high (it's up t 120 Hz)repsectively

total_power_band = [0, 75];      
frw = 0:0.2:75;
frq = sampling_freq; % sampling frequnecy of EEG data
window = 5; %sec. 1 for 30 sec


[transition_spectrogram, F, T] = spectrogram(Data_EEG,round(frq*window),[],frw,frq,'yaxis');
mean_spectrogram = log(abs(transition_spectrogram));
time_spectrogram_zero = T;
time_offset = t2b_removed_s; % Offset value in milliseconds
T_adjusted = T+time_offset;
time_spectrogram_adjusted = time_spectrogram_zero + time_offset;
filtered_mean_spectrogram = imgaussfilt(mean_spectrogram, 4);



figure;
a = subplot(3, 1, 1);
    imagesc(time_spectrogram_adjusted, F, filtered_mean_spectrogram); %plot the log spectrum
    set(gca,'YDir', 'normal'); % flip the Y Axis so lower frequencies are at the bottom
    ylim([0, 75]);
    xlim([-8, -4])
    colormap(gca, 'parula');
    hold on
    for band_i = 1:length(power_bands)
        plot([-295, -295], power_bands{band_i}, 'LineWidth', 5)
    end
    title('EEG power');
    ylabel('freq (Hz)');
b = subplot(3, 1, 2);
    band_power_collector = T;
    for band_i = 1:length(power_bands)
        power_band = power_bands{band_i};
        power_trace = mean(mean_spectrogram(find(F==power_band(1)):find(F==power_band(2)), :), 1);
        normalized_power_trace = power_trace;
        band_power_collector = [band_power_collector; normalized_power_trace];
        plot(time_spectrogram_adjusted, normalized_power_trace)
        hold on
    end
    legend({'delta','theta','sigma','beta', 'gamma low', 'gamma high'});
c = subplot(3, 1, 3);
    sigma = band_power_collector(4,:); % sigma power trace
    plot(time_spectrogram_adjusted, sigma)
    title ('sigma')
    xlabel('time (s)');
linkaxes([a,b,c],'x');

%% sigma power trace analysis

% Smoothing the trace in sigma power graph with filtfilt function
% Find peaks and troughs of the sigma trace (used to estimate amplitude of NE oscillations = peak – trough + estimate spindle density (max value))

 
window_size = 5; %smoothing sigma
b = (1/window_size) * ones(1, window_size);
a = 1;

smoothed_sigma = filtfilt(b, a, sigma);
diff_sigma = diff(smoothed_sigma);
inverted_diff_sigma = -diff_sigma;
inverted_diff_sigma_NaN = [inverted_diff_sigma, NaN];

% Plot smoothed sigma, find peaks and troughs
figure;
plot_sleep(time_spectrogram_adjusted, smoothed_sigma, hypnogram_timetrace_cut_ds, wakeVector_ds', REMVector_ds', NREMVector_ds');
hold on;
[peaks, peak_locations] = findpeaks(smoothed_sigma);
scatter(time_spectrogram_adjusted(peak_locations), peaks, 'yellow', 'filled');
inverted_smoothed_sigma = -smoothed_sigma;
[troughs,trlocs] = findpeaks(inverted_smoothed_sigma);
scatter(time_spectrogram_adjusted(trlocs), -troughs, 'g', 'filled');
hold off



window_size_inverted_diff = 5; %smoothing inverted diff sigma
b_inverted_diff = (1/window_size_inverted_diff) * ones(1, window_size_inverted_diff);
a_inverted_diff = 1;
smoothed_inverted_diff_sigma = filtfilt(b_inverted_diff, a_inverted_diff, inverted_diff_sigma); 

smoothed_inverted_diff_sigma_NaN = [smoothed_inverted_diff_sigma, NaN];
[troughs,trlocs_diff] = findpeaks(smoothed_inverted_diff_sigma_NaN);

%calculating the mean amplitude between peaks and troughs within NREM periods 

NREM_peaks = peaks(ismember(peak_locations, sws_periods));
NREM_troughs = troughs(ismember(trlocs, sws_periods)); %this is giving an empty vector ????

mean_amplitude_NREM = mean(NREM_peaks - NREM_troughs);


figure
a = subplot(3, 1, 1);
plot(time_spectrogram_adjusted, smoothed_inverted_diff_sigma_NaN)
    title ('smoothed inverted diff sigma')
    xlabel('time (s)');
b = subplot(3, 1, 2);
plot(time_spectrogram_adjusted, smoothed_sigma)
    title ('smoothed sigma')
    xlabel('time (s)');
c = subplot(3, 1, 3);
plot(time_spectrogram_adjusted, inverted_diff_sigma_NaN) %just to compare with the not smoothed
    title ('inverted diff sigma')
    xlabel('time (s)');
linkaxes([a,b,c],'x');
%% %% select only NREM sleep in the graph and mean amplitude 

peak_times = trlocs_diff;

peak_timestamps = T(peak_times); % time stamps for all peaks
[~, ~, peak_timestamps_NREM] = IsInInterval(peak_timestamps, sws_periods);% select timestamps for peaks only during NREM
peak_threshold = 0.05;
sigma_fs = length(T_adjusted)/(T_adjusted(end)-T_adjusted(1));
trough_timestamps = T_adjusted(peak_times); % time stamps for all peaks
[~, ~, trough_timestamps_NREM] = IsInInterval(trough_timestamps, sws_periods); % select timestamps for peaks only during NREM

figure
plot_sleep(time_spectrogram_adjusted, smoothed_sigma, hypnogram_timetrace_cut_ds, wakeVector_ds', REMVector_ds', NREMVector_ds');
hold on
plot(trough_timestamps_NREM, smoothed_sigma(round((trough_timestamps_NREM-t2b_removed_s)*sigma_fs)), '*');
title ('sigma')
xlabel('time (s)');


figure
a = subplot(2, 1, 1);
    plot_sleep(C3_timetrace_cut_ds, edf_C3_cut_ds, hypnogram_timetrace_cut_ds, wakeVector_ds', REMVector_ds', NREMVector_ds');
    xlabel('time (s)');
    ylabel('EEG (V)');
b = subplot(2, 1, 2);
    plot_sleep(time_spectrogram_adjusted, smoothed_sigma, hypnogram_timetrace_cut_ds, wakeVector_ds', REMVector_ds', NREMVector_ds');
    hold on
    plot(trough_timestamps_NREM, smoothed_sigma(round((trough_timestamps_NREM-t2b_removed_s)*sigma_fs)), '*');
    title ('sigma')
    xlabel('time (s)');
linkaxes([a,b],'x');

%% Extracting other values from edf data

RespThorax = edf_data{11}; %Extracting Respiration thorax data 
ECG = edf_data {7}; %Extracting ECG data 
TIBV = edf_data {9}; %Extracting TIBV muscle data 
SAO2 = edf_data {15}; %Extracting oxigen saturation data 

RespThorax_cut = RespThorax(t2b_removed_samples:end);
RespThorax_cut_ds = downsample(RespThorax_cut, 100);
ECG_cut = ECG(t2b_removed_samples:end);
ECG_cut_ds = downsample(ECG_cut, 100);
TIBV_cut = TIBV(t2b_removed_samples:end);
TIBV_cut_ds = downsample(TIBV_cut,100);
SAO2_cut = SAO2(t2b_removed_samples:end);
SAO2_cut_ds = downsample(SAO2_cut, 100);
scoring_cut = scoring(t2b_removed_s:end);
scoring_cut_ds = downsample(scoring_cut, 50);

fig = figure;
a = subplot(6,1,1);
plot(hypnogram_timetrace_cut, scoring_cut);
    title('hipnogram');
    xlabel('time');
    ylabel('sleep phase');
b = subplot(6,1,2);
plot_sleep(C3_timetrace_cut_ds, edf_C3_cut_ds, hypnogram_timetrace_cut_ds, wakeVector_ds', REMVector_ds', N1Vector_ds', N2Vector_ds', N3Vector_ds');
    xlabel('time (s)');
    ylabel('EEG (V)');
c = subplot(6,1,3);
plot_sleep(C3_timetrace_cut_ds, RespThorax_cut_ds, hypnogram_timetrace_cut_ds, wakeVector_ds', REMVector_ds', N1Vector_ds', N2Vector_ds', N3Vector_ds');
xlabel('Time');
ylabel('Respiration');
d = subplot(6,1,4);
plot_sleep(C3_timetrace_cut_ds, ECG_cut_ds, hypnogram_timetrace_cut_ds, wakeVector_ds', REMVector_ds', N1Vector_ds', N2Vector_ds', N3Vector_ds');
xlabel('Time');
ylabel('ECG');
e = subplot(6,1,5);
plot_sleep(C3_timetrace_cut_ds, TIBV_cut_ds, hypnogram_timetrace_cut_ds, wakeVector_ds', REMVector_ds', N1Vector_ds', N2Vector_ds', N3Vector_ds');
xlabel('Time');
ylabel('EMG');
f = subplot(6,1,6);
plot_sleep(C3_timetrace_cut_ds, SAO2_cut_ds, hypnogram_timetrace_cut_ds, wakeVector_ds', REMVector_ds', N1Vector_ds', N2Vector_ds', N3Vector_ds');
xlabel('Time');
ylabel('SAO2');
linkaxes([a, b, c, d, e, f],'x');
 
figure;
plot_sleep(C3_timetrace_cut_ds, ECG_cut_ds, hypnogram_timetrace_cut_ds, wakeVector_ds', REMVector_ds', N1Vector_ds', N2Vector_ds', N3Vector_ds');
xlabel('Time');
ylabel('ECG');

figure;
plot_sleep(C3_timetrace_cut_ds, RespThorax_cut_ds, hypnogram_timetrace_cut_ds, wakeVector_ds', REMVector_ds', N1Vector_ds', N2Vector_ds', N3Vector_ds');
xlabel('Time');
ylabel('Respiration');

figure;
plot_sleep(C3_timetrace_cut_ds, SAO2_cut_ds, hypnogram_timetrace_cut_ds, wakeVector_ds', REMVector_ds', N1Vector_ds', N2Vector_ds', N3Vector_ds');
xlabel('Time');
ylabel('SAO2');

figure;
plot_sleep(C3_timetrace_cut_ds, TIBV_cut_ds, hypnogram_timetrace_cut_ds, wakeVector_ds', REMVector_ds', N1Vector_ds', N2Vector_ds', N3Vector_ds');
xlabel('Time');
ylabel('EMG');

%% Analysis of sleep structure
% In this section we are analysing how much time is spent in each phase, to
% be able to compare different patients.

%calculate time spent in each phase 

total_dur_wake = sum(wake_duration); %when you sum the duration of all your bouts you get the total time spent in this state (in seconds)
total_dur_NREM = sum(duration_sws);
total_dur_REM = sum(REM_duration);
total_dur_NREM1 = sum(N1_duration);
total_dur_NREM2 = sum(N2_duration);
total_dur_NREM3 = sum(N3_duration);

% Calculate percentages

calculate_percentage = @(binary_vector) (sum(binary_vector) / numel(binary_vector)) * 100;

percentage_wake = calculate_percentage(wakeVector_ds);
percentage_NREM = calculate_percentage(NREMVector_ds);
percentage_NREM1 = calculate_percentage(N1Vector_ds);
percentage_NREM2 = calculate_percentage(N2Vector_ds);
percentage_NREM3 = calculate_percentage(N3Vector_ds);
percentage_REM = calculate_percentage(REMVector_ds);

% Display percentages
disp(['Percentage of time spent awake: ', num2str(percentage_wake), '%']);
disp(['Percentage of time spent in NREM: ', num2str(percentage_NREM), '%']);
disp(['Percentage of time spent in NREM1: ', num2str(percentage_NREM1), '%']);
disp(['Percentage of time spent in NREM2: ', num2str(percentage_NREM2), '%']);
disp(['Percentage of time spent in NREM3: ', num2str(percentage_NREM3), '%']);
disp(['Percentage of time spent in REM: ', num2str(percentage_REM), '%']);

% number of bouts/h
total_dur_s = length(scoring_cut);
total_dur_h = total_dur_s/60/60;
freq_wake_h = length(wake_duration)/total_dur_h;
freq_REM_h = length(REM_duration)/total_dur_h;
freq_NREM_h = length(duration_sws)/total_dur_h;
freq_NREM1_h = length(N1_duration)/total_dur_h;
freq_NREM2_h = length(N2_duration)/total_dur_h;
freq_NREM3_h = length(N3_duration)/total_dur_h;

%calculation of peaks/time to get a crude measure of the oscillation
%frequency (frequency/phase).
peak_timestamps = T(peak_times); % time stamps for all peaks
[~, ~, peak_timestamps_NREM] = IsInInterval(peak_timestamps, sws_periods);
[~, ~, peak_timestamps_N1] = IsInInterval(peak_timestamps, N1_periods);
[~, ~, peak_timestamps_N2] = IsInInterval(peak_timestamps, N2_periods);
[~, ~, peak_timestamps_N3] = IsInInterval(peak_timestamps, N3_periods);

freq_peaks_NREM = length(peak_timestamps_NREM)/total_dur_NREM; % peacks frequency in NREM period
freq_peaks_N1 = length(peak_timestamps_N1)/total_dur_NREM1;
freq_peaks_N2 = length(peak_timestamps_N2)/total_dur_NREM2;
freq_peaks_N3 = length(peak_timestamps_N3)/total_dur_NREM3;

%% Power spectral density (EEG)

% this will result in a plot of power across EEG frequencies (very common way to look at EEG sleep data)


% power spectral densities

t1 = sws_periods(:,1); 
t2 = sws_periods(:,2);

tsamp1 = floor((t1-t2b_removed_s)*sampling_freq); %eeg start time
tsamp2 = floor((t2-t2b_removed_s)*sampling_freq); %eeg end time

NREM_data = cell(1, numel(tsamp1));

PXX = [];

%NREM_data_collect = [];

for i=1:numel(tsamp1)
    NREM_data{i} = edf_C3_cut(tsamp1(i):tsamp2(i)); % "EEG_rawtrace_cut_filter" is the EEG_trace - so you should put your C3 trace here
    %NREM_data_cut = edf_C3_cut(tsamp1(i):tsamp2(i));
    %NREM_data_collect = [NREM_data_collect NREM_data_cut];
    [pxx, f] = pwelch(NREM_data{i}, [], [],[0:0.2:100], sampling_freq);
    logpxx = 10*log10(pxx);
    FX{i} = f;
    PXX(:,i) = logpxx;
    PXX(:,i) = pxx;
end 

mean_PXX = mean(PXX,2);
prism_psd = mean_PXX(f<45);
prism_freq = f(f<45);

% power spectral density plot

figure
plot(prism_freq,prism_psd)
mean_sigma_power_density = mean(mean_PXX(f>7 & f<15)); % look at sigma power across NREM states
mean_delta_power_density = mean(mean_PXX(f>1 & f<4));
mean_theta_power_density = mean(mean_PXX(f>4 & f<7));
mean_beta_power_density = mean(mean_PXX(f>15 & f<30));

prism_band_collect = [mean_delta_power_density mean_theta_power_density mean_sigma_power_density mean_beta_power_density]';

%% PSD on sigma trace

% This is the power spectral denity of your sigma trace. So it's similar to
% the analysis to get EEG PSD, since you will get a plot of how much power
% sigma oscillates a different frequencies (though these should range a lot lower than the EEG frequencies)

analysis_state_periods = N3_periods; %insert here the different periods, so N1_periods, sws_periods, N2_periods and N3_periods
period = analysis_state_periods-t2b_removed_s;
sigma_analysis = smoothed_sigma; %or just sigma 

% power spectral densities

t1 = analysis_state_periods(:,1);
t2 = analysis_state_periods(:,2);

tsamp1 = floor((t1-t2b_removed_s)*sigma_fs); %eeg start time
tsamp2 = floor((t2-t2b_removed_s)*sigma_fs); %eeg end time

NREM_data = cell(1, numel(tsamp1));

PXX = [];
PXXlog = [];
PXX_pk_f = [];
PXX_pk = [];

NREM_data_collect = [];
period_duration = [];

for i=1:numel(tsamp1)                                                          % This loop will run through all the bouts (one bout per iteration) of your selected state
    period_length_i = tsamp2(i)-tsamp1(i);
    if period_length_i < 120*sigma_fs % periods shorter than 120 s are excluded from analysis
        continue
    end
    if tsamp2(i) > length(sigma_analysis) % if last period ends after trace
       tsamp2(i) = length(sigma_analysis);
    end
    period_duration = [period_duration period_length_i/sigma_fs];
    NREM_data{i} = sigma_analysis(tsamp1(i):tsamp2(i));                                  
    timetrace_i = time_spectrogram_adjusted(tsamp1(i):tsamp2(i));
    %center around 0 with polyfit
    [p,s,mu] = polyfit((1:numel(NREM_data{i}))',NREM_data{i},5);
    f_y = polyval(p,(1:numel(NREM_data{i}))',[],mu);
    detrend_data = NREM_data{i} - f_y';        % Detrend data

    [pxx, f] = pwelch(detrend_data, [], [],[0:0.002:0.1], sigma_fs); %
    logpxx = 10*log10(pxx);
    FX{i} = f;
    [pxx_pk_psd, max_idx] = max(pxx);
    PXX_pk = [PXX_pk pxx_pk_psd];
    pxx_pk_f = f(max_idx);
    PXX_pk_f = [PXX_pk_f pxx_pk_f];
    PXX = [PXX pxx'];

    figure                                                                       % <<<< Consider commenting this figure out, since it will make a figure for each iteration
    set(gcf, 'Position',  [100, 300, 1500, 250])                                 % You could set i = [some random number] and run the code inside the loop to get an idea of what the code does
    a = subplot(1,2,1);                                                          % the plot shows you how it detrends the data and show you the PSD plot for each bout
        plot(timetrace_i,NREM_data{i});
        hold on
        plot(timetrace_i,detrend_data);
        legend({'raw','fitted'})
        hold off
    b = subplot(1,2,2);
        plot(f,pxx);
end

mean_PXX_sigma = mean(PXX,2);
weighted_mean_PXX_sigma = sum(period_duration.*PXX,2)/sum(period_duration); %period duration is used as weights



PXX_sigma_f_weighted_mean = sum(PXX_pk_f.*period_duration)/sum(period_duration);
PXX_sigma_pk_weighted_mean = sum(PXX_pk.*period_duration)/sum(period_duration);

[PXX_sigma_pk_mean, PXX_sigma_pk_idx] = max(weighted_mean_PXX_sigma);
PXX_sigma_f_mean = f(PXX_sigma_pk_idx);

%prism_psd_sigma = weighted_mean_PXX_sigma/mean(weighted_mean_PXX_sigma); % normalised power

prism_psd_sigma = weighted_mean_PXX_sigma;
prism_freq_sigma = f;

% power spectral density plot

figure                                                                          % This figure will be a mean PSD for all of your bouts
    plot(prism_freq_sigma, prism_psd_sigma)
    xlabel('frequency (Hz)');
    ylabel('PSD');
    period_fs = period*sigma_fs;
 %% 90-10 percentile of each NREM bout

period_fs = (analysis_state_periods-t2b_removed_s)*sigma_fs;

period_ampl = [];
trace_sigma_stitch = [];

for i = 1:length(period_fs)
   period_i = period_fs(i,:);
   tracesigma = sigma_analysis(period_i(1):period_i(2));                                      
   perc90 = prctile(tracesigma,90);
   perc10 = prctile(tracesigma,10);
   perc_amplitude = perc90-perc10;
   period_ampl = [period_ampl perc_amplitude];
   trace_sigma_stitch = [trace_sigma_stitch tracesigma];

 % figure
 %  plot(tracesigma)
 %  hold on
 %  yline(perc90)
 % yline(perc10)
end

perc90_full = prctile(trace_sigma_stitch,90);
perc10_full = prctile(trace_sigma_stitch,10);
perc_amplitude_full = perc90_full-perc10_full;

figure
    plot(trace_sigma_stitch)
    hold on
    yline(perc90_full)
    yline(perc10_full)

% weighted average based on period duration (short NREM periods weigh less than long ones)
period_dur = period(:,2)-period(:,1);
prism_mean_amplitude_weighted = sum(period_ampl*period_dur)/sum(period_dur);

%% EMG root mean square

% perhaps filter EMG to remove noise first

rms_wndw = 500; % sliding windows used for envelope rms
rms_EMG = envelope(TIBV_cut_ds,rms_wndw,'rms');

%% epoc plots of various signals

epoc_start = 30; % number of seconds prior to detected event
epoc_end = 30; % number of seconds after to detected event

epoc_time = (1:1:ceil((epoc_start+epoc_end)*(sampling_freq/100)))/(sampling_freq/100)-epoc_start;
 
%ds_factor = 30;
analysis_window = 5; %sec. 1 for 30 sec
epoc_event_times = trough_timestamps_NREM;
epoc_event_times_cut = epoc_event_times - t2b_removed_s; 

EMG_collector = [];
for i = 1:length(epoc_event_times_cut)
   pk_i = epoc_event_times_cut(i);
     if pk_i > C3_timetrace_cut_ds(end) - epoc_end  % skip if beh_onset is too close to end of recording (i.e. tail is too short) %se il mio punto e' a meno di 30 secondi dalla fine del timetrace allora non ho i 30 secondi dopo per l'analisi
     continue
    end
    pk_epoc_i = rms_EMG((pk_i - epoc_start)*(sampling_freq/100) : (pk_i + epoc_end)*(sampling_freq/100)); %here says in the raw trace take from my detected point, 30 seconds before, to 30 seconds after the detected point 
    EMG_collector = [EMG_collector; pk_epoc_i]; %this is filling the collector, so that when the code repetes, it accumulates and it doesnt take only the last one (racoon doesn't leave all the cookies that he took until now)
end

mean_EMGpk_epocs = mean(EMG_collector); %here u make the mean, because the EMG collector would make all different traces, so you take the mean one to simplify it 

figure
plot(epoc_time, mean_EMGpk_epocs)
title('NREM EMG peak');


ECG_collector = [];
for i = 1:length(epoc_event_times_cut)
    pk_i = epoc_event_times_cut(i);
    if pk_i > C3_timetrace_cut_ds(end) - epoc_end
        continue 
    end
    pk_epoc_i = ECG((pk_i - epoc_start)*(sampling_freq/100) : (pk_i + epoc_end)*(sampling_freq/100));
    ECG_collector = [ECG_collector; pk_epoc_i];
end 

mean_ECGpk_epocs = mean(ECG_collector,1);

figure
plot(epoc_time, mean_ECGpk_epocs)
title('NREM ECG peak');

SAO2_collector = [];
for i = 1:length(epoc_event_times_cut)
    pk_i = epoc_event_times_cut(i);
    if pk_i > C3_timetrace_cut_ds(end) - epoc_end
        continue 
    end
    pk_epoc_i = SAO2((pk_i - epoc_start)*(sampling_freq/100) : (pk_i + epoc_end)*(sampling_freq/100));
    SAO2_collector = [SAO2_collector; pk_epoc_i];
end 

mean_SAO2pk_epocs = mean(SAO2_collector,1);

figure
plot(epoc_time, mean_SAO2pk_epocs)
title('NREM SAO2 peak');

Resp_collector = [];
for i = 1:length(epoc_event_times_cut)
    pk_i = epoc_event_times_cut(i);
    if pk_i > C3_timetrace(end) - epoc_end
        continue 
    end
    pk_epoc_i = RespThorax((pk_i - epoc_start)*(sampling_freq/100) : (pk_i + epoc_end)*(sampling_freq/100));
    Resp_collector = [Resp_collector; pk_epoc_i];
end 

mean_Resppk_epocs = mean(Resp_collector,1);

figure
plot(epoc_time, mean_Resppk_epocs)
title('NREM RespThorax peak');

C3_collector = [];
for i = 1:length(epoc_event_times_cut)
    pk_i = epoc_event_times_cut(i);
    if pk_i > C3_timetrace(end) - epoc_end
        continue 
    end
    pk_epoc_i = edf_C3((pk_i - epoc_start)*(sampling_freq/100) : (pk_i + epoc_end)*(sampling_freq/100));
    C3_collector = [C3_collector; pk_epoc_i];
end 

mean_C3_epocs = mean(C3_collector,1);

figure
plot(epoc_time, mean_C3_epocs)
title('NREM C3 peak');

    delta_epocs = [];
    theta_epocs = [];
    sigma_epocs = [];
    beta_epocs = [];
    gamma_low_epocs = [];
    gamma_high_epocs = [];
    epoc_spectrogram_collector = [];
for i = 1:length(epoc_event_times_cut)
    pk_i = epoc_event_times_cut(i);
     if pk_i > C3_timetrace(end) - epoc_end
         continue
    end
    pk_epoc_i = edf_C3_cut((pk_i - epoc_start)*(sampling_freq/100) : (pk_i + epoc_end)*(sampling_freq/100));

   [epoc_spectrogram, F, T] = spectrogram(pk_epoc_i,round(frq*analysis_window),[],frw,frq,'yaxis'); % F = frequenciy vector, T=time vector
        log_epoc_spectrogram = log(abs(epoc_spectrogram));
        
    delta_power_trace = mean(log_epoc_spectrogram(find(F==power_bands{1}(1)):find(F==power_bands{1}(2)), :), 1);
    theta_power_trace = mean(log_epoc_spectrogram(find(F==power_bands{2}(1)):find(F==power_bands{2}(2)), :), 1);
    sigma_power_trace = mean(log_epoc_spectrogram(find(F==power_bands{3}(1)):find(F==power_bands{3}(2)), :), 1);
    beta_power_trace = mean(log_epoc_spectrogram(find(F==power_bands{4}(1)):find(F==power_bands{4}(2)), :), 1);
    gamma_low_power_trace = mean(log_epoc_spectrogram(find(F==power_bands{5}(1)):find(F==power_bands{5}(2)), :), 1); 
    gamma_high_power_trace = mean(log_epoc_spectrogram(find(F==power_bands{6}(1)):find(F==power_bands{6}(2)), :), 1);
        
    delta_epocs = [delta_epocs delta_power_trace'];
    theta_epocs = [theta_epocs theta_power_trace'];
    sigma_epocs = [sigma_epocs  sigma_power_trace'];
    beta_epocs = [beta_epocs beta_power_trace'];
    gamma_low_epocs = [gamma_low_epocs gamma_low_power_trace'];
    gamma_high_epocs = [gamma_high_epocs gamma_high_power_trace'];
    epoc_spectrogram_collector = cat(3, epoc_spectrogram_collector, epoc_spectrogram);
end 

log_spectrogram = log(abs(epoc_spectrogram_collector));
mean_spectrogram = nanmean(log_spectrogram, 3);
norm_time = abs(T_adjusted-(epoc_start/3*2)); 
norm_sampling = find(norm_time == min(norm_time));
normalization_factor = mean(mean(mean_spectrogram(find(F==total_power_band(1)):find(F==total_power_band(2)), 1:norm_sampling)));
time_spectrogram = T_adjusted-epoc_start;


figure()
a = subplot(2, 1, 1);
    filtered_mean_spectrogram = imgaussfilt(mean_spectrogram, 4);
    imagesc(time_spectrogram, F, filtered_mean_spectrogram); %plot the log spectrum
    set(gca,'YDir', 'normal'); % flip the Y Axis so lower frequencies are at the bottom
    ylim([0, 30]);
    caxis([-6.5, -4.7])
    title('spectrogram');
    colormap(gca, 'parula');
    hold on
    for band_i = 1:length(power_bands)
        plot([-295, -295], power_bands{band_i}, 'LineWidth', 5)
    end
b = subplot(2, 1, 2);
    band_power_collector = [T];
    for band_i = 1:length(power_bands)
        power_band = power_bands{band_i};
        band_power_trace = mean(mean_spectrogram(find(F==power_band(1)):find(F==power_band(2)), :), 1);
        normalized_band_power = band_power_trace/-normalization_factor+2;
        band_power_collector = [band_power_collector; normalized_band_power]; % matrix of timevector and aligned frequency band power traces
        plot(time_spectrogram, normalized_band_power)
        hold on
    end
     title('power traces');
linkaxes([a,b],'x');

%% Epoc trace extraction

 
analysis_window = 5; %sec. 1 for 30 sec

% all these "collectors" will contain individual traces for all the detected transitions. You need on for every parameter you will extract.
EMG_collector = [];
Resp_collector = [];
C3_collector = [];
ECG_collector = [];
SAO2_collector = [];
delta_epocs = [];
theta_epocs = [];
sigma_epocs = [];
beta_epocs = [];
gamma_low_epocs = [];
gamma_high_epocs = [];

epoc_spectrogram_collector = [];

for i = 1:length(epoc_event_times_cut)
    pk_i = epoc_event_times_cut(i);
    if pk_i > C3_timetrace(end) - epoc_end
        continue 
    pk_epoc_i = edf_C3((pk_i - epoc_start)*sampling_freq : (pk_i + epoc_end)*sampling_freq); 
        
        % % create EEG epoc trace
        % epoc_idx = round(epoc_time*frq);
        % epoc_before_idx = round(epoc_idx-frq*epoc_start);
        % epoc_after_idx = round(epoc_idx+frq*epoc_end);       
        % eeg_epoc_trace = edf_C3(:, epoc_before_idx:epoc_after_idx);       
        % rms_emg_trace = rms_EMG(:, epoc_before_idx:epoc_after_idx);
        % 
        % powerspectrogrtam on EEG epoc trace
        [epoc_spectrogram, F, T] = spectrogram(pk_epoc_i,round(frq*analysis_window),[],frw,frq,'yaxis'); % F = frequenciy vector, T=time vector
        log_epoc_spectrogram = log(abs(epoc_spectrogram));
        
        delta_power_trace = mean(log_epoc_spectrogram(find(F==power_bands{1}(1)):find(F==power_bands{1}(2)), :), 1);
        theta_power_trace = mean(log_epoc_spectrogram(find(F==power_bands{2}(1)):find(F==power_bands{2}(2)), :), 1);
        sigma_power_trace = mean(log_epoc_spectrogram(find(F==power_bands{3}(1)):find(F==power_bands{3}(2)), :), 1);
        beta_power_trace = mean(log_epoc_spectrogram(find(F==power_bands{4}(1)):find(F==power_bands{4}(2)), :), 1);
        gamma_low_power_trace = mean(log_epoc_spectrogram(find(F==power_bands{5}(1)):find(F==power_bands{5}(2)), :), 1);
        gamma_high_power_trace = mean(log_epoc_spectrogram(find(F==power_bands{6}(1)):find(F==power_bands{6}(2)), :), 1);
        
        delta_epocs = [delta_epocs delta_power_trace'];
        theta_epocs = [theta_epocs theta_power_trace'];
        sigma_epocs = [sigma_epocs  sigma_power_trace'];
        beta_epocs = [beta_epocs beta_power_trace'];
        gamma_low_epocs = [gamma_low_epocs gamma_low_power_trace'];
        gamma_high_epocs = [gamma_high_epocs gamma_high_power_trace'];
        
        epoc_spectrogram_collector = cat(3, epoc_spectrogram_collector, epoc_spectrogram);
        rms_EMG_epocs = [rms_EMG_epocs rms_emg_trace'];
    else continue
    end
end

log_spectrogram = log(abs(epoc_spectrogram_collector));
mean_spectrogram = nanmean(log_spectrogram, 3);
norm_time = abs((T)-(epoc_start/3*2)); 
norm_sampling = find(norm_time == min(norm_time));
normalization_factor = mean(mean(mean_spectrogram(find(F==total_power_band(1)):find(F==total_power_band(2)), 1:norm_sampling)));

% normalization of EEG band power traces to avoid negative values
normalized_delta_power_epocs = delta_epocs/-normalization_factor+2;
normalized_theta_power_epocs = theta_epocs/-normalization_factor+2;
normalized_sigma_power_epocs = sigma_epocs/-normalization_factor+2;
normalized_beta_power_epocs = beta_epocs/-normalization_factor+2;
normalized_gamma_lo_power_epocs = gamma_low_epocs/-normalization_factor+2;
normalized_gamma_hi_power_epocs = gamma_high_epocs/-normalization_factor+2;

time_spectrogram_zero = T-epoc_start; % time vector matching band power traces
time_EMG = (1:1:length(rms_EMG))/sampling_freq-epoc_start; % time vector matching EMG (rms) trace
EEG_bands_fs = length(T_adjusted)/(T_adjusted(end)-T_adjusted(1));

% figure
% plot(time_spectrogram_zero,mean(normalized_delta_power_epocs,2))
% hold on
% plot(time_spectrogram_zero, mean(normalized_theta_power_epocs,2))
% plot(time_spectrogram_zero, mean(normalized_sigma_power_epocs,2))
% plot(time_spectrogram_zero, mean(normalized_beta_power_epocs,2))
% plot(time_spectrogram_zero, mean(normalized_gamma_lo_power_epocs,2))
% plot(time_spectrogram_zero, mean(normalized_gamma_hi_power_epocs,2))

% figure of mean traces and mean power                                      <<< Look at this figure to decide which time windows to use for you pre/post analysis.
figure()
a = subplot(3, 1, 1);
    filtered_mean_spectrogram = imgaussfilt(mean_spectrogram, 4);
    imagesc(time_spectrogram_adjusted, F, filtered_mean_spectrogram); %plot the log spectrum
    set(gca,'YDir', 'normal'); % flip the Y Axis so lower frequencies are at the bottom
    ylim([0, 30]);
    caxis([-6.5, -4.7])
    title('spectrogram');
    colormap(gca, 'parula');
    hold on
    for band_i = 1:length(power_bands)
        plot([-295, -295], power_bands{band_i}, 'LineWidth', 5)
    end
b = subplot(3, 1, 2);
    band_power_collector = [T];
    for band_i = 1:length(power_bands)
        power_band = power_bands{band_i};
        band_power_trace = mean(mean_spectrogram(find(F==power_band(1)):find(F==power_band(2)), :), 1);
        normalized_band_power = band_power_trace/-normalization_factor+2;
        band_power_collector = [band_power_collector; normalized_band_power]; % matrix of timevector and aligned frequency band power traces
        plot(time_spectrogram_zero, normalized_band_power)
        hold on
    end
     title('power traces');
c = subplot(3, 1, 3);
    plot(time_EMG,mean(rms_EMG_epocs,2))
    title('EMG rms');
linkaxes([a,b,c],'x');


%% my code 
delta = band_power_collector(2,:);
Theta = band_power_collector(3,:);
sigma = band_power_collector(4,:);
beta =  band_power_collector(5,:);
gamma = band_power_collector(6,:);

EEG_bands_fs = length(T_adjusted)/(T_adjusted(end)-T_adjusted(1));

Delta_collector = [];
Theta_collector = [];
sigma_collector = [];
beta_collector = [];
gamma_collector = [];

for i = 1:length(epoc_event_times_cut)
    pk_i = epoc_event_times_cut(i);
    if pk_i > C3_timetrace_cut_ds(end)-epoc_end % skip if beh_onset is too close to end of recording (i.e. tail is too short)
     continue
    end

 
pk_delta_i = delta((pk_i - epoc_start)*EEG_bands_fs:(pk_i + epoc_end)*EEG_bands_fs);
pk_Theta_i = Theta((pk_i - epoc_start)*EEG_bands_fs:(pk_i + epoc_end)*EEG_bands_fs);
pk_sigma_i = sigma((pk_i - epoc_start)*EEG_bands_fs:(pk_i + epoc_end)*EEG_bands_fs);
pk_beta_i = beta((pk_i - epoc_start)*EEG_bands_fs:(pk_i + epoc_end)*EEG_bands_fs);
pk_gamma_i = gamma((pk_i - epoc_start)*EEG_bands_fs:(pk_i + epoc_end)*EEG_bands_fs);


Delta_collector = [Delta_collector; pk_delta_i];
Theta_collector = [Theta_collector; pk_Theta_i];
sigma_collector = [sigma_collector; pk_sigma_i];
beta_collector = [beta_collector; pk_beta_i];
gamma_collector = [gamma_collector; pk_gamma_i];
end
 

epoc_time = (1:1:ceil((epoc_start+epoc_end)*EEG_bands_fs))/EEG_bands_fs-epoc_start;
mean_Delta_pk_epocs = mean(Delta_collector,1);
mean_Theta_pk_epocs = mean(Theta_collector,1);
mean_sigma_pk_epocs = mean(sigma_collector,1);
mean_beta_pk_epocs = mean(beta_collector,1);
mean_gamma_pk_epocs = mean(gamma_collector,1);
 

figure
plot(epoc_time, mean_Delta_pk_epocs)
hold on
plot(epoc_time, mean_Theta_pk_epocs)
plot(epoc_time, mean_sigma_pk_epocs)
plot(epoc_time, mean_beta_pk_epocs)
plot(epoc_time, mean_gamma_pk_epocs)
hold off
title('NREM bands peak');
legend({'delta','Theta','sigma','beta', 'gamma'});


 %% Before/after levels (to be used for cluster analysis)
% Each epoc should end up with before after level for NE, all EEG bands and EMG rms

epoc_MA_time0 = epoc_start; % time zero of detected transition within epoc traces
pre_start = 30; % s before detected event for baseline calculation (start)       <<<< Here you adjust pre/post windows based on what makes sense looking at the figure above.
pre_end = 0; % s before detected event for baseline calculation (end)                 The pre/post windows should be the same across all patients of course
post_start = 0; % s after detected event for effect calculation (start)
post_end = 30; % s after detected event for effect calculation (end)

% Here you prepare empty vectors that will store the changes (pre-post) for each transition of the various parameters.
SAO2_diff_prepost = NaN(length(epoc_event_times_cut),1);
Respthorax_diff_prepost = NaN(length(epoc_event_times_cut),1);
C3_diff_prepost = NaN(length(epoc_event_times_cut),1);
EMG_diff_prepost = NaN(length(epoc_event_times_cut),1);
ECG_diff_prepost = NaN(length(epoc_event_times_cut),1);
delta_diff_prepost = NaN(length(epoc_event_times_cut),1);
theta_diff_prepost = NaN(length(epoc_event_times_cut),1);
sigma_diff_prepost = NaN(length(epoc_event_times_cut),1);
beta_diff_prepost = NaN(length(epoc_event_times_cut),1);
gamma_diff_prepost = NaN(length(epoc_event_times_cut),1);

% levels pre and post transitions are determined, and the change in signal is claculated for all variabels

% NB! Decide for each parameter what output measure is fitting (e.g. percentile, mean, max, min, median, area under curve (AUC)).
% Mean/median/AUC are probably more safe choices, if the change can both be positive and negative, since it's not biasing toward a peak/trough, while
% percentil/max/min are used if you are looking for specifically for a peak or trough follwoing the transition
for epoc_i = 1:length(epoc_event_times_cut)
    ECG_epoc_i = ECG_collector(:,epoc_i);
    ECG_baseline_pre_event_i = prctile(ECG_epoc_i(round((epoc_start-pre_start)*(sampling_freq/100):(epoc_start-pre_end)*(sampling_freq/100))),5);
    ECG_baseline_post_on_i = prctile(ECG_epoc_i(round((epoc_start+post_start)*(sampling_freq/100):(epoc_start+post_end)*(sampling_freq/100))),5);
    ECG_diff_prepost_i = ECG_baseline_post_on_i-ECG_baseline_pre_event_i;
    
    Respthorax_epoc_i = Resp_collector(:,epoc_i);
    Respthorax_baseline_pre_event_i = prctile(Respthorax_epoc_i(round((epoc_start-pre_start)*(sampling_freq/100):(epoc_start-pre_end)*(sampling_freq/100))),5);
    Respthorax_baseline_post_on_i = prctile(Respthorax_epoc_i(round((epoc_start+post_start)*(sampling_freq/100):(epoc_start+post_end)*(sampling_freq/100))),5);
    Respthorax_diff_prepost_i = Respthorax_baseline_post_on_i-Respthorax_baseline_pre_event_i;
    
    SAO2_epoc_i = SAO2_collector(:,epoc_i);
    SAO2_baseline_pre_event_i = prctile(SAO2_epoc_i(round((epoc_start-pre_start)*(sampling_freq/100):(epoc_start-pre_end)*(sampling_freq/100))),5);
    SAO2_baseline_post_on_i = prctile(SAO2_epoc_i(round((epoc_start+post_start)*(sampling_freq/100):(epoc_start+post_end)*(sampling_freq/100))),5);
    SAO2_diff_prepost_i = SAO2_baseline_post_on_i-SAO2_baseline_pre_event_i;
    
    EMG_epoc_i = rms_EMG_epocs(:,epoc_i);
    EMG_baseline_pre_event_i = prctile(EMG_epoc_i(round((epoc_start-pre_start)*(sampling_freq/100):(epoc_start-pre_end)*(sampling_freq/100))),5);
    EMG_baseline_post_on_i = prctile(EMG_epoc_i(round((epoc_start+post_start)*(sampling_freq/100):(epoc_start+post_end)*(sampling_freq/100))),5);
    EMG_diff_prepost_i = EMG_baseline_post_on_i-EMG_baseline_pre_event_i;

    C3_epoc_i = epoc_spectrogram_collector(:,epoc_i);
    C3_baseline_pre_event_i = prctile(C3_epoc_i(round((epoc_start-pre_start)*(sampling_freq/100):(epoc_start-pre_end)*(sampling_freq/100))),5);
    C3_baseline_post_on_i = prctile(C3_epoc_i(round((epoc_start+post_start)*(sampling_freq/100):(epoc_start+post_end)*(sampling_freq/100))),5);
    C3_diff_prepost_i = C3_baseline_post_on_i-C3_baseline_pre_event_i;
    
    delta_epoc_i = Delta_collector(:,epoc_i);
    delta_baseline_pre_event_i = prctile(delta_epoc_i(round((epoc_start-pre_start)*(sigma_fs/100):(epoc_start-pre_end)*(sigma_fs/100))),5);
    delta_baseline_post_on_i = prctile(delta_epoc_i(round((epoc_start+post_start)*(sigma_fs/100):(epoc_start+post_end)*(sigma_fs/100))),5);
    delta_diff_prepost_i = delta_baseline_post_on_i-delta_baseline_pre_event_i;
    
    Theta_epoc_i = Theta_collector(:,epoc_i);
    Theta_baseline_pre_event_i = prctile(Theta_epoc_i(round((epoc_start-pre_start)*(sigma_fs/100):(epoc_start-pre_end)*(sigma_fs/100))),5);
    Theta_baseline_post_on_i = prctile(Theta_epoc_i(round((epoc_start+post_start)*(sigma_fs/100):(epoc_start+post_end)*(sigma_fs/100))),5);
    Theta_diff_prepost_i = Theta_baseline_post_on_i-Theta_baseline_pre_event_i;
    
    Sigma_epoc_i = sigma_collector(:,epoc_i);
    Sigma_baseline_pre_event_i = prctile(Sigma_epoc_i(round((epoc_start-pre_start)*(sigma_fs/100):(epoc_start-pre_end)*(sigma_fs/100))),5);
    Sigma_baseline_post_on_i = prctile(Sigma_epoc_i(round((epoc_start+post_start)*(sigma_fs/100):(epoc_start+post_end)*(sigma_fs/100))),5);
    Sigma_diff_prepost_i = Sigma_baseline_post_on_i-Sigma_baseline_pre_event_i;
    
    beta_epoc_i = beta_collector(:,epoc_i);
    beta_baseline_pre_event_i = prctile(beta_epoc_i(round((epoc_start-pre_start)*(sigma_fs/100):(epoc_start-pre_end)*(sigma_fs/100))),5);
    beta_baseline_post_on_i = prctile(beta_epoc_i(round((epoc_start+post_start)*(sigma_fs/100):(epoc_start+post_end)*(sigma_fs/100))),5);
    beta_diff_prepost_i = beta_baseline_post_on_i-beta_baseline_pre_event_i;
    
    gamma_epoc_i = gamma_collector(:,epoc_i);
    gamma_baseline_pre_MA_i = prctile(gamma_epoc_i(round((epoc_start-pre_start)*(sigma_fs/100):(epoc_start-pre_end)*(sigma_fs/100))),5);
    gamma_baseline_post_on_i = prctile(gamma_epoc_i(round((epoc_start-pre_start)*(sigma_fs/100):(epoc_start-pre_end)*(sigma_fs/100))),5);
    gam_lo_diff_prepost_i = gamma_baseline_post_on_i-gamma_baseline_pre_MA_i;
    
    
    SAO2_diff_prepost(epoc_i,1) = SAO2_diff_prepost_i;
    Respthorax_diff_prepost(epoc_i,1) = Respthorax_diff_prepost_i;
    ECG_diff_prepost(epoc_i,1) = ECG_diff_prepost_i;
    EMG_diff_prepost(epoc_i,1) = EMG_diff_prepost_i;
    C3_diff_prepost(epoc_i,1) = C3_diff_prepost_i;
    delta_diff_prepost(epoc_i,1) = delta_diff_prepost_i;
    theta_diff_prepost(epoc_i,1) = theta_diff_prepost_i;
    sigma_diff_prepost(epoc_i,1) = sigma_diff_prepost_i;
    beta_diff_prepost(epoc_i,1) = beta_diff_prepost_i;
    gamma_diff_prepost(epoc_i,1) = gam_lo_diff_prepost_i;
    
end


%% write table containing MA variables
% you can try to see if this work. Should create an excel file with all the
% variables

mat2prism_xlsxfilepath = 'H:\MATLAB\mat2prism_export_MAvariables.xlsx'; % select folder location

rec_ID = folder_name(1);     % <<< Put in patient ID. This will be the row name in excel file
rec_ID(1:length(NE_diff_prepost))= rec_ID;
rec_ID=rec_ID';

prism_table_3 = table(rec_ID, NE_diff_prepost, As_diff_prepost, LC_diff_prepost, EMG_diff_prepost, delta_diff_prepost, theta_diff_prepost, sigma_diff_prepost, beta_diff_prepost, gam_lo_diff_prepost, gam_hi_diff_prepost); %here sheet 1 of the table is created

if ~isfile(mat2prism_xlsxfilepath) % no existing file
    writetable(prism_table_3,mat2prism_xlsxfilepath,'WriteRowNames',true)
    prism_table_4 = table(scored_MA_idx','VariableNames',rec_ID(1)); %here sheet 2 of the table is created
    writetable(prism_table_4,mat2prism_xlsxfilepath, 'sheet', 'scored_MA_idx')
else                                                                                            % <<<< The next patients should be added as new columns/rows in your excel sheet here
    existing_table_3 = readtable(mat2prism_xlsxfilepath,'sheet',1);
    IDs_in_table = existing_table_3.rec_ID;
    if ~any(IDs_in_table==rec_ID(1)) % this only runs if the ID is not already in he table (i.e. to prevent duplicates)
        writetable(prism_table_3,mat2prism_xlsxfilepath,'WriteMode','Append', 'WriteVariableNames',false,'WriteRowNames',true)
        prism_table_4 = readtable(mat2prism_xlsxfilepath, 'sheet', 'scored_MA_idx');

        % Add a new column to the end of the array (when vectors are not the same length, NaNs are put in the end of columns to match size)
        if height(prism_table_4) >= length(scored_MA_idx)
            newcol = NaN(height(prism_table_4),1);
            newcol(1:length(scored_MA_idx)) = scored_MA_idx;
            prism_table_4_update = [prism_table_4 table(newcol, 'VariableNames', rec_ID(1))]; %the existing table is concatenated with added column to create new table
        else    % add NaNs (rows) to table
            add_rows_No = length(scored_MA_idx)-height(prism_table_4);
            newrows = [prism_table_4; array2table(nan(add_rows_No,width(prism_table_4)),'variablenames',prism_table_4.Properties.VariableNames)];
            prism_table_4_update = [newrows table(scored_MA_idx', 'VariableNames', rec_ID(1))];
        end

        writetable(prism_table_4_update,mat2prism_xlsxfilepath, 'sheet', 'scored_MA_idx')
    end
end

%% PCA or cluster analysis (single patient)
% Here is the PCA or cluster analysis. You can run this section once you
% have vectors with all your pre-post changes (one per detected transition)
% for a single patient, or you can run the next section when you have all
% the pre-post measures from more patients

%MA_pca_matrix = [NE_diff_prepost As_diff_prepost LC_diff_prepost EMG_diff_prepost delta_diff_prepost theta_diff_prepost sigma_diff_prepost beta_diff_prepost gam_lo_diff_prepost gam_hi_diff_prepost ];
MA_pca_matrix = [As_diff_prepost LC_diff_prepost EMG_diff_prepost delta_diff_prepost theta_diff_prepost sigma_diff_prepost beta_diff_prepost gam_lo_diff_prepost gam_hi_diff_prepost ];


%[pca_coeff, pca_scores, ~, ~, pca_explained] = pca(MA_pca_matrix, 'Algorithm','eig');
[pca_coeff, pca_scores, ~, ~, pca_explained] = pca(MA_pca_matrix);

figure
    scatter3(pca_scores(:,1),pca_scores(:,2),pca_scores(:,3))
    hold on
    scatter3(pca_scores(scored_MA_idx,1),pca_scores(scored_MA_idx,2),pca_scores(scored_MA_idx,3), 'r', 'filled') % scored MA within detection threshold
    %axis equal
    xlabel('1st Principal Component')
    ylabel('2nd Principal Component')
    zlabel('3nd Principal Component')

figure
scatter(NE_diff_prepost,sigma_diff_prepost) % correlation between NE and other outputs

figure
plot(time_FP, mean(epoc_FP_NE_collector(:,delta_diff_prepost>0),2))
hold on
plot(time_FP, mean(epoc_FP_NE_collector(:,delta_diff_prepost<0),2))


%% Cluster analysis on all recordings
%PCA or Means Shift clustering or Gaussian Mixture Model

mat2prism_xlsxfilepath = 'H:\MATLAB\mat2prism_export_MAvariables.xlsx';         % << This should be the same file location as the one where you save the excel sheet
MA_table = readtable(mat2prism_xlsxfilepath,'sheet',1,'TextType','string');
MA_scored_table = readtable(mat2prism_xlsxfilepath,'sheet',2);

rec_IDs = unique(MA_table.rec_ID, 'stable');

for i = 1:length(rec_IDs)
    number_of_MA(i) = length(MA_table{MA_table.rec_ID==rec_IDs(i),:});
    MA_table_IDstart_idx(i) = find(MA_table.rec_ID==rec_IDs(i),1);
end

table_idx = table2array(MA_scored_table)+MA_table_IDstart_idx-1;
vector_idx = reshape(table_idx,[],1);
vector_idx = vector_idx(~isnan(vector_idx)); % remove NaNs

%how to extract rows based on row variable (e.g. extracting only a specific recording)
%MA_table{MA_table.rec_ID=="EEG_FP_triple_237",:};

%PCA
MA_pca_matrix = [ MA_table.EMG_diff_prepost MA_table.delta_diff_prepost MA_table.theta_diff_prepost MA_table.sigma_diff_prepost MA_table.beta_diff_prepost MA_table.gam_lo_diff_prepost MA_table.gam_hi_diff_prepost ];

%[pca_coeff, pca_scores, ~, ~, pca_explained] = pca(MA_pca_matrix, 'Algorithm','eig');
[pca_coeff, pca_scores, ~, ~, pca_explained] = pca(MA_pca_matrix);

figure
    scatter(pca_scores(:,1),pca_scores(:,2))
    hold on
    scatter(pca_scores(vector_idx,1),pca_scores(vector_idx,2), 'r', 'filled') % scored MA within detection threshold
    %axis equal
    xlabel('1st Principal Component')
    ylabel('2nd Principal Component')

figure
    scatter3(pca_scores(:,1),pca_scores(:,2),pca_scores(:,3))
    hold on
    scatter3(pca_scores(vector_idx,1),pca_scores(vector_idx,2),pca_scores(vector_idx,3), 'r', 'filled') % scored MA within detection threshold
    %axis equal
    xlabel('1st Principal Component')
    ylabel('2nd Principal Component')
    zlabel('3nd Principal Component')

