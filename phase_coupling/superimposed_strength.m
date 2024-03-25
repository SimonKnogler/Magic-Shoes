% ecgEventsTable_low and gaitEventsTable_low have to be loaded 

% Convert the tables to arrays for processing
E_low = table2array(ecgEventsTable_low); % n x 3, time of peaks in second row
G_lo = table2array(gaitEventsTable_low(:, [1, 3])); % Assuming we only need columns 1 and 3

% Initialize arrays to store values for participants 31 to 59
all_meanT_ratios_low = zeros(1, 29);
all_couplings_24 = zeros(1, 29); % For (2,4) parameters
all_couplings_35 = zeros(1, 29); % For (3,5) parameters
all_couplings_45 = zeros(1, 29); % For (4,5) parameters

% Loop over participants 31 to 59
for part = 31:59
    % Heart beats
    E_S16 = E_low(E_low(:,3) == part, 2); % Select rows for participant and second column

    % Gait events
    G_S16 = G_lo(G_lo(:,2) == part, 1); % Select rows for participant and first column

    % Compute mean period ratio
    meanT_ratio = mean(diff(E_S16)) / mean(diff(G_S16));
    all_meanT_ratios_low(part-30) = meanT_ratio; % Store mean period ratio

    % Perform IS analysis with (2,4) parameters
    ISres_24 = ISanalysis(E_S16, G_S16, 2, 4);
    all_couplings_24(part-30) = ISres_24.RM; % Store (2,4) coupling strength

    % Perform IS analysis with (3,5) parameters
    ISres_35 = ISanalysis(E_S16, G_S16, 3, 5);
    all_couplings_35(part-30) = ISres_35.RM; % Store (3,5) coupling strength

    % Perform IS analysis with (4,5) parameters
    ISres_45 = ISanalysis(E_S16, G_S16, 4, 5);
    all_couplings_45(part-30) = ISres_45.RM; % Store (4,5) coupling strength
end

% Plotting the results to compare the robustness of the method
figure;

% Plot coupling strengths for different parameters
plot(31:59, all_couplings_24, '-og', 'DisplayName', '(2,4)');
hold on;  % Hold the plot for superimposing the next curves
plot(31:59, all_couplings_35, '-ob', 'DisplayName', '(3,5)');
plot(31:59, all_couplings_45, '-or', 'DisplayName', '(4,5)');
hold off;  % Release the plot
xlabel('Participant');
ylabel('N:M Coupling');
title('N:M Coupling for low DP with Different Parameters');
legend('show');  % Show legend with display names
grid on;

% Adjust figure window size if needed
set(gcf, 'Position', [100, 100, 1049, 895]);

% Save the figure as a file (optional)
saveas(gcf, 'participant_synchronization_analysis_comparison.png');





