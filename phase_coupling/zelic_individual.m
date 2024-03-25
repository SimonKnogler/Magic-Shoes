% Specify the path to the EDF file
file_path = '/Users/simonknogler/Desktop/MS1/Low/C/31.edf';

% Load the EDF file
[s, h] = sload(file_path);

% Assume that the accelerometer channel for the Z-axis is the 4th channel
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

% Design a low-pass filter
[b, a] = butter(4, cutoff_frequency / (sample_rate / 2), 'low');

% Apply the filter to the trimmed normalized magnitude of the acceleration for the Z-axis
Magnitude_Accel_Z_filtered = filtfilt(b, a, Magnitude_Accel_Z_normalized);

% Visualize the unfiltered and filtered signals for the trimmed data
figure;
subplot(2,1,1);
plot(Magnitude_Accel_Z_normalized);
title('Unfiltered Normalized Magnitude of Accelerometer (Z-axis, 60 seconds)');

subplot(2,1,2);
plot(Magnitude_Accel_Z_filtered);
title('Filtered Normalized Magnitude of Accelerometer (Z-axis, 60 seconds)');

% Define initial parameters for peak detection
peak_prominence_factor = 1; % Factor to multiply with std deviation for peak prominence
min_distance_samples = 500; % Minimum distance between peaks in samples

% Find peaks in the filtered Z-axis signal
[peaks, locs] = findpeaks(Magnitude_Accel_Z_filtered, 'MinPeakProminence', peak_prominence_factor * std(Magnitude_Accel_Z_filtered), 'MinPeakDistance', min_distance_samples);

% Visual inspection of the findpeaks result
figure;
plot(Magnitude_Accel_Z_filtered);
hold on;
plot(locs, peaks, 'r*', 'MarkerSize', 10);
hold off;
title('Peaks in Filtered Normalized Magnitude of Acceleration (Z-axis, 60 seconds)');
xlabel('Sample Number');
ylabel('Acceleration (g)');
legend('Filtered Signal', 'Detected Peaks');
