% install BioSig toolbox and add to path
% type eeglab in command window

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Prepare acc data - the slower wave %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify the path to the EDF file
file_path = '/Users/simonknogler/Desktop/MS1/High/C/12.edf';

% Load the EDF file
[s, h] = sload(file_path);

% The accelerometer channel for the Z-axis is the 4th channel
Accel_Z = s(:,4);

% Convert ADC values to g-values
adc_resolution = 1023;
mid_point = adc_resolution / 2;
accel_range = 2;
Accel_Z_g = ((Accel_Z - mid_point) / mid_point) * accel_range;

% Normalize the magnitude of acceleration of the Z-axis
Magnitude_Accel_Z_normalized = (Accel_Z_g - mean(Accel_Z_g)) / std(Accel_Z_g);

% Specify the sample rate and cutoff frequency for the low-pass filter
sample_rate = 1000; % Hz
cutoff_frequency = 5; % Hz

% Trim the Z-axis signal to 60 seconds worth of data
samples_to_keep = 60 * sample_rate; % 60 seconds times the sample rate
Accel_Z_g = Accel_Z_g(1:samples_to_keep);
Magnitude_Accel_Z_normalized = Magnitude_Accel_Z_normalized(1:samples_to_keep);

% Design the low-pass filter
[b, a] = butter(4, cutoff_frequency / (sample_rate / 2), 'low');

% Apply the filter to the trimmed normalized magnitude of the acceleration for the Z-axis
Magnitude_Accel_Z_filtered = filtfilt(b, a, Magnitude_Accel_Z_normalized);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Prepare ECG data - the faster wave %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the EDF file
[s, h] = sload('/Users/simonknogler/Desktop/MS1/High/C/10.edf');

% SignalLabel1_RAW corresponds to the first signal (ECG).
ECG_Values = s(:,1);

% Display the first few values
disp('First few ECG Values:');
disp(ECG_Values(1:10));

% Trim the signal at 60 seconds
sample_rate = 1000; % Sampling rate is 1000 Hz
time_to_trim = 60; % Duration in seconds to keep the signal
samples_to_keep = sample_rate * time_to_trim; % 60 seconds worth of samples
ECG_Values = ECG_Values(1:samples_to_keep);

%% 1st Filter the data to make it smoother and easier to detect R-peaks 
% The following lines run with plugin HEPLAB, installed in EEGLAB Toolbox
         
lcf = 3; % Low cut-off in Hz (HEPLAB manual recommendations)
hcf = 30; % High cut-off in Hz (HEPLAB recommendations)
ECG_original = ECG_Values; % ECG_Values is a vector with subjects' ECG data, i.e., one long vector
ECG_Values = (heplab_ecg_filt(ECG_Values, sample_rate, lcf, hcf));
            
% Finding R peaks via findpeaks function in Matlab
M = max(ECG_Values);
valuePeaks = M./4;
PeakDistance = 450;
[PKSecg, Sub_ecgPoints] = findpeaks(ECG_Values, 'MinPeakHeight', valuePeaks, 'MinPeakDistance', PeakDistance);
            
% Convert R peak indices to times (in seconds)
R_peak_times = Sub_ecgPoints / sample_rate;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Visualize n:m phase synchronization %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Phase extraction for the filtered acceleration data (slower wave)
phi_footsteps = angle(hilbert(Magnitude_Accel_Z_filtered));

% Convert R_peak_times to sample indices if they're in seconds
R_peak_indices = round(R_peak_times * sample_rate);

% Calculate the synchrogram: phase of footsteps at the peak times of the ECG
sync = phi_footsteps(R_peak_indices);

% Visualization 
figSig = figure('Name','Synchrogram');

% Plot the two signals (Magnitude_Accel_Z_filtered and ECG_Values)
subplot(211)
plot(Magnitude_Accel_Z_filtered,'b'),hold on,plot(ECG_Values,'r')
title('Acceleration and ECG signal')
xlabel('Sample Number');

% Plot the synchrogram
subplot(212)
plot(sync,'o')
title('Synchrogram')
xlabel('Sample Number');
ylabel('Phase');

% Adjust the aspect ratio of the synchrogram to flatten the slopes of the diagonals
pbaspect([5 1 1]); % This will stretch the plot horizontally

sgtitle('Visualization of Phase Synchronization');


