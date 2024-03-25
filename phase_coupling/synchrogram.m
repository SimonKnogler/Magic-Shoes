% Determine the new common sampling rate (choose one of the original rates or another value)
new_sample_rate = max([length(ECG_Values), length(Magnitude_Accel_filtered)]);

% Resample ECG_Values to match the length of Magnitude_Accel_filtered
ECG_Values_resampled = resample(ECG_Values, length(Magnitude_Accel_filtered), length(ECG_Values));

% Now ECG_Values_resampled and Magnitude_Accel_filtered have the same number of samples.
% You can now proceed with your original synchrogram method.

% Find the R-peaks in the resampled ECG signal
[pks, RPeakIndices] = findpeaks(ECG_Values_resampled, 'MinPeakHeight', mean(ECG_Values_resampled) + std(ECG_Values_resampled), 'MinPeakDistance', 200);

% Calculate the phase of the slower signal using the Hilbert transform
phi_magnitude = angle(hilbert(Magnitude_Accel_filtered));

% Extract the phase information at the R-peak indices
% Ensure that RPeakIndices does not exceed the length of phi_magnitude
RPeakIndices = RPeakIndices(RPeakIndices <= length(phi_magnitude));
sync = phi_magnitude(RPeakIndices);

% Plot the synchrogram
figSig = figure('Name', 'Synchrogram');
subplot(211)
plot(Magnitude_Accel_filtered, 'b')
title('Filtered Accelerometer Data')

subplot(212)
plot(sync, 'o')
title('Synchrogram')
xlabel('R-Peak Indices')
ylabel('Phase of Accelerometer Signal at ECG R-Peaks')
