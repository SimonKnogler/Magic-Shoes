% Define the folder where the EDF files are located
folderPath = '/Users/simonknogler/Desktop/MS1/Low/L';

% Get a list of all EDF files in the folder
edfFiles = dir(fullfile(folderPath, '*.edf'));

% Initialize the dataframe to store step event results
stepEventsTable = table();

% Loop over each file
for iFile = 1:length(edfFiles)
    % Construct the full path to the EDF file
    filePath = fullfile(folderPath, edfFiles(iFile).name);

    % Load the EDF file
    [s, h] = sload(filePath);

    % Cut at 60 seconds
    num_samples_to_take = min(size(s,1), 60000);
    Accel_X = s(1:num_samples_to_take, 2);
    Accel_Y = s(1:num_samples_to_take, 3);
    Accel_Z = s(1:num_samples_to_take, 4);

    % Convert ADC values to g-values
    adc_resolution = 1023;
    mid_point = adc_resolution / 2;
    accel_range = 2;
    Accel_X_g = ((Accel_X - mid_point) / mid_point) * accel_range;
    Accel_Y_g = ((Accel_Y - mid_point) / mid_point) * accel_range;
    Accel_Z_g = ((Accel_Z - mid_point) / mid_point) * accel_range;

    % Calculate the magnitude of the acceleration vector
    Accel_Magnitude = sqrt(Accel_X_g.^2 + Accel_Y_g.^2 + Accel_Z_g.^2);

    % Normalize the magnitude of acceleration
    Magnitude_Accel_normalized = (Accel_Magnitude - mean(Accel_Magnitude)) / std(Accel_Magnitude);

    % Trim the normalized magnitude of acceleration to 60 seconds worth of data
    samples_to_keep = 60 * 1000; % Using 1000 Hz as the sample rate
    Magnitude_Accel_normalized = Magnitude_Accel_normalized(1:samples_to_keep);

    % Specify the sample rate and cutoff frequency for the low-pass filter
    sample_rate = 1000; % Hz
    cutoff_frequency = 5; % Hz

    % Design the low-pass filter
    [b, a] = butter(4, cutoff_frequency / (sample_rate / 2), 'low');

    % Apply the filter to the trimmed normalized magnitude of acceleration
    Accel_Magnitude_filtered = filtfilt(b, a, Magnitude_Accel_normalized);

    % Initial peak detection parameters
    peak_prominence_factor = 1.5; % Factor to multiply with std deviation for peak prominence
    min_distance_samples = 800; % Minimum distance between peaks in samples

    % Loop for interactive adjustment of parameters
    exit_loop = false;
    while ~exit_loop
        % Detect step events with the detectStepEvents function
        step_events = detectStepEvents(Accel_Magnitude_filtered, sample_rate, peak_prominence_factor * std(Accel_Magnitude_filtered), min_distance_samples);

        % Plot the filtered signal with the detected events
        figure('Name','Step Event Detection','NumberTitle','off');
        plot(Accel_Magnitude_filtered, 'b'); % Plot the filtered magnitude of acceleration in blue
        hold on;
        plot(step_events, Accel_Magnitude_filtered(step_events), 'ro', 'MarkerFaceColor', 'r'); % Plot step events as red circles
        title('Detected Step Events');
        xlabel('Samples');
        ylabel('Acceleration Magnitude');
        legend('Filtered Acceleration Magnitude', 'Step Events');
        hold off;

        % Prompt the user for input
        user_input = input('Approve results? (y)es, (n)o, or e(x)it: ', 's');

        % Handle user input
        switch user_input
            case 'y'
                % User approved the results, process and append data to dataframe
                exit_loop = true; % Exit the loop
            case 'n'
                % User disapproved, allow parameter adjustments
                peak_prominence_factor = input('Enter new peak prominence factor: ');
                min_distance_samples = input('Enter new minimum distance between peaks (in samples): ');
            case 'x'
                % User chose to exit
                disp('Exiting...');
                exit_loop = true;
            otherwise
                % Handle invalid input
                disp('Invalid input.');
        end
    end

    % If the user approved the results
    if strcmp(user_input, 'y')
        % Create a table with the results for this participant
        participantTable = table(step_events / sample_rate, repmat({'Step'}, length(step_events), 1),...
                                 'VariableNames', {'Time', 'EventType'});
        
        % Add participant ID
        participant_id = erase(edfFiles(iFile).name, '.edf'); % Extract participant ID from filename
        participantTable.ParticipantID = repmat({participant_id}, height(participantTable), 1);
        
        % Append to the main stepEventsTable
        stepEventsTable = [stepEventsTable; participantTable];
    end

    % If the user chooses to exit
    if strcmp(user_input, 'x')
        break; % Break out of the file loop
    end
end


% After all files are processed, or the user exits
% Save the combined dataframe to a CSV file
writetable(stepEventsTable, fullfile(folderPath, 'StepEventResults.csv'));
disp('All data processed and saved to StepEventResults.csv');

% Auxiliary function to detect step events
function step_events = detectStepEvents(Accel_Magnitude_filtered, sample_rate, peak_prominence, min_distance_samples)
    [~, step_events] = findpeaks(Accel_Magnitude_filtered, 'MinPeakProminence', peak_prominence, 'MinPeakDistance', min_distance_samples);
end