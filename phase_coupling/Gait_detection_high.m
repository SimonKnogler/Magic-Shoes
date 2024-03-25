% Define the folder where the EDF files are located
folderPath = '/Users/simonknogler/Desktop/MS1/Low/C';

% Get a list of all EDF files in the folder
edfFiles = dir(fullfile(folderPath, '*.edf'));

% Initialize the dataframe to store gait event results
gaitEventsTable = table();

% Loop over each file
for iFile = 1:length(edfFiles)
    % Construct the full path to the EDF file
    filePath = fullfile(folderPath, edfFiles(iFile).name);

    % Load the EDF file
    [s, h] = sload(filePath);

    % Cut at 60 seconds
    num_samples_to_take = min(size(s,1), 60000);
    Accel_Z = s(1:num_samples_to_take, 4);

    % Convert ADC values to g-values
    adc_resolution = 1023;
    mid_point = adc_resolution / 2;
    accel_range = 2;
    Accel_Z_g = ((Accel_Z - mid_point) / mid_point) * accel_range;

    % Normalize the magnitude of acceleration of the Z-axis
    Magnitude_Accel_Z_normalized = (Accel_Z_g - mean(Accel_Z_g)) / std(Accel_Z_g);

    % Trim the Z-axis signal to 60 seconds worth of data
    samples_to_keep = 60 * 1000; % Using 1000 Hz as the sample rate
    Accel_Z_g = Accel_Z_g(1:samples_to_keep);
    Magnitude_Accel_Z_normalized = Magnitude_Accel_Z_normalized(1:samples_to_keep);

    % Specify the sample rate and cutoff frequency for the low-pass filter
    sample_rate = 1000; % Hz
    cutoff_frequency = 5; % Hz

    % Design the low-pass filter
    [b, a] = butter(4, cutoff_frequency / (sample_rate / 2), 'low');

    % Apply the filter to the trimmed normalized magnitude of the acceleration for the Z-axis
    Accel_Z_filtered = filtfilt(b, a, Magnitude_Accel_Z_normalized);

    % Initial filtering and peak detection parameters
    peak_prominence_factor = 2; % Factor to multiply with std deviation for peak prominence
    min_distance_samples = 1500; % Minimum distance between peaks in samples
    window_size = 2000; % Size of the sliding window in samples
    overlap_step = 100; % Overlap step size in samples
    cadence_threshold_seconds = 0.5; % Cadence threshold in seconds

    % Loop for interactive adjustment of parameters
    exit_loop = false;
    while ~exit_loop
        % Detect gait events with the detectGaitEvents function
        [filtered_IC, filtered_TO] = detectGaitEvents(Accel_Z_filtered, sample_rate, peak_prominence_factor * std(Accel_Z_filtered), min_distance_samples, window_size, overlap_step, cadence_threshold_seconds);

        % Plot the filtered signal with the detected events
        figure('Name','Gait Event Detection','NumberTitle','off');
        plot(Accel_Z_filtered, 'b'); % Plot the filtered Z-axis acceleration in blue
        hold on;
        plot(filtered_IC, Accel_Z_filtered(filtered_IC), 'ro', 'MarkerFaceColor', 'r'); % Plot IC as red circles
        plot(filtered_TO, Accel_Z_filtered(filtered_TO), 'go', 'MarkerFaceColor', 'g'); % Plot TO as green circles
        title('Detected Gait Events');
        xlabel('Samples');
        ylabel('Acceleration (Z-axis)');
        legend('Filtered Acceleration', 'Initial Contact (IC)', 'Toe-Off (TO)');
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
        participantTable = table(filtered_IC / sample_rate, repmat({'IC'}, length(filtered_IC), 1),...
                                 'VariableNames', {'Time', 'EventType'});
        toeOffsTable = table(filtered_TO / sample_rate, repmat({'TO'}, length(filtered_TO), 1),...
                             'VariableNames', {'Time', 'EventType'});
        % Combine heel strikes and toe-offs
        participantTable = [participantTable; toeOffsTable];
        % Sort by time
        participantTable = sortrows(participantTable, 'Time');
        
        % Add participant ID
        participant_id = erase(edfFiles(iFile).name, '.edf'); % Extract participant ID from filename
        participantTable.ParticipantID = repmat({participant_id}, height(participantTable), 1);
        
        % Append to the main gaitEventsTable
        gaitEventsTable = [gaitEventsTable; participantTable];
    end

    % If the user chooses to exit
    if strcmp(user_input, 'x')
        break; % Break out of the file loop
    end
end


% After all files are processed, or the user exits
% Save the combined dataframe to a CSV file
writetable(gaitEventsTable, fullfile(folderPath, 'GaitEventResults.csv'));
disp('All data processed and saved to GaitEventResults.csv');

% Auxiliary functions to detect gait events
function [filtered_IC, filtered_TO] = detectGaitEvents(Accel_Z_filtered, sample_rate, peak_prominence, min_distance_samples, window_size, overlap_step, cadence_threshold_seconds)
    n = length(Accel_Z_filtered);
    window_start = 1;
    potential_IC = [];
    potential_TO = [];
    while window_start + window_size <= n
        window_end = window_start + window_size - 1;
        window_data = Accel_Z_filtered(window_start:window_end);
        [window_IC_peaks, window_IC_locs] = findpeaks(window_data, 'MinPeakProminence', peak_prominence, 'MinPeakDistance', min_distance_samples);
        [window_TO_peaks, window_TO_locs] = findpeaks(-window_data, 'MinPeakProminence', peak_prominence, 'MinPeakDistance', min_distance_samples);
        window_IC_locs = window_IC_locs + window_start - 1;
        window_TO_locs = window_TO_locs + window_start - 1;
        potential_IC = unique([potential_IC; window_IC_locs]);
        potential_TO = unique([potential_TO; window_TO_locs]);
        window_start = window_start + overlap_step;
    end
    filtered_IC = filterICEvents(potential_IC, sample_rate, cadence_threshold_seconds);
    filtered_TO = filterTOEvents(potential_TO, filtered_IC, sample_rate, cadence_threshold_seconds);
end

function filtered_IC = filterICEvents(potential_IC, sample_rate, cadence_threshold_seconds)
    cadence_samples = cadence_threshold_seconds * sample_rate;
    filtered_IC = [];
    for i = 1:length(potential_IC)
        if i == 1 || (potential_IC(i) - potential_IC(i-1)) > cadence_samples
            filtered_IC = [filtered_IC; potential_IC(i)];
        end
    end
end

function filtered_TO = filterTOEvents(potential_TO, filtered_IC, sample_rate, cadence_threshold_seconds)
    cadence_samples = cadence_threshold_seconds * sample_rate;
    filtered_TO = [];
    for i = 1:length(filtered_IC)
        TO_after_IC = potential_TO(potential_TO > filtered_IC(i));
        if ~isempty(TO_after_IC)
            filtered_TO = [filtered_TO; TO_after_IC(1)];
        end
    end
end
