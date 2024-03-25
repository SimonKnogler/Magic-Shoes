% Assuming ecgEventsTable_high and gaitEventsTable_high are loaded elsewhere

% Convert the tables to arrays for processing for the second group
E_high = table2array(ecgEventsTable_high); % n x 3, time of peaks in second row
G_hi = table2array(gaitEventsTable_high(:, [1, 3])); % Assuming we only need columns 1 and 3

% Initialize arrays to store values for participants 31 to 59 for the second group
all_meanT_ratios_high = zeros(1, 29);
all_couplings_24_high = zeros(1, 29); % For (2,4) parameters
all_couplings_35_high = zeros(1, 29); % For (3,5) parameters
all_couplings_45_high = zeros(1, 29); % For (4,5) parameters

% Loop over participants 31 to 59 for the second group
for part = 31:59
    % Heart beats
    E_S16_high = E_high(E_high(:,3) == part, 2); % Select rows for participant and second column

    % Gait events
    G_S16_high = G_hi(G_hi(:,2) == part, 1); % Select rows for participant and first column

    % Compute mean period ratio for the second group
    meanT_ratio_high = mean(diff(E_S16_high)) / mean(diff(G_S16_high));
    all_meanT_ratios_high(part-30) = meanT_ratio_high; % Store mean period ratio

    % Perform IS analysis with (2,4) parameters for the second group
    ISres_24_high = ISanalysis(E_S16_high, G_S16_high, 2, 4);
    all_couplings_24_high(part-30) = ISres_24_high.RM; % Store (2,4) coupling strength

    % Perform IS analysis with (3,5) parameters for the second group
    ISres_35_high = ISanalysis(E_S16_high, G_S16_high, 3, 5);
    all_couplings_35_high(part-30) = ISres_35_high.RM; % Store (3,5) coupling strength

    % Perform IS analysis with (4,5) parameters for the second group
    ISres_45_high = ISanalysis(E_S16_high, G_S16_high, 4, 5);
    all_couplings_45_high(part-30) = ISres_45_high.RM; % Store (4,5) coupling strength
end

% Plotting the results to compare the robustness of the method for both groups
figure;

% Plot coupling strengths for different parameters for the first group
plot(31:59, all_couplings_24, '-og', 'DisplayName', 'Group Low (2,4)');
hold on;  % Hold the plot for superimposing the next curves
plot(31:59, all_couplings_35, '-ob', 'DisplayName', 'Group Low (3,5)');
plot(31:59, all_couplings_45, '-or', 'DisplayName', 'Group Low (4,5)');

% Plot coupling strengths for different parameters for the second group
plot(31:59, all_couplings_24_high, '--dg', 'DisplayName', 'Group High (2,4)');
plot(31:59, all_couplings_35_high, '--db', 'DisplayName', 'Group High (3,5)');
plot(31:59, all_couplings_45_high, '--dr', 'DisplayName', 'Group High (4,5)');

hold off;  % Release the plot
xlabel('Participant');
ylabel('N:M Coupling');
title('N:M Coupling for Participants 31 to 59 with Different Parameters');
legend('show');  % Show legend with display names
grid on;

% Adjust figure window size if needed
set(gcf, 'Position', [100, 100, 1049, 895]);

% Save the figure as a file (optional)
saveas(gcf, 'participant_synchronization_analysis_comparison_both_groups.png');
