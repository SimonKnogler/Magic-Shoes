% Assuming ecgEventsTable_low and gaitEventsTable_low are loaded elsewhere

% Convert the tables to arrays for processing
E_low = table2array(ecgEventsTable_low); % n x 3, time of peaks in second row
G_low = table2array(gaitEventsTable_low(:, [1, 3])); % Assuming we only need columns 1 and 3

% Initialize arrays to store values for participants 31 to 59
all_meanT_ratios_low_24 = zeros(1, 29); % For (2,4) parameters
all_meanT_ratios_low_35 = zeros(1, 29); % For (3,5) parameters
all_meanT_ratios_low_45 = zeros(1, 29); % For (4,5) parameters
all_couplings_24_low = zeros(1, 29); % For (2,4) parameters
all_couplings_35_low = zeros(1, 29); % For (3,5) parameters
all_couplings_45_low = zeros(1, 29); % For (4,5) parameters
all_bestRatios_24_low = zeros(1, 29); % For (2,4) parameters
all_bestRatios_35_low = zeros(1, 29); % For (3,5) parameters
all_bestRatios_45_low = zeros(1, 29); % For (4,5) parameters
all_bestRatios_45_low = zeros(1, 29); % For (4,5) parameters 

% Loop over participants 31 to 59
for part = 31:59
    idx = part - 30; % Adjust index to start from 1 for array storage
    % Heart beats
    E_S16 = E_low(E_low(:,3) == part, 2); % Select rows for participant and second column

    % Gait events
    G_S16 = G_low(G_low(:,2) == part, 1); % Select rows for participant and first column

    % Perform IS analysis with (2,4) parameters
    ISres_24 = ISanalysis(E_S16, G_S16, 2, 4);
    all_meanT_ratios_low_24(idx) = mean(diff(E_S16)) / mean(diff(G_S16)); % Store (2,4) mean period ratio
    all_couplings_24_low(idx) = ISres_24.RM; % Store (2,4) coupling strength
    all_bestRatios_24_low(idx) = ISres_24.CIstats(1,1); % Store (2,4) best ratio

    % Perform IS analysis with (3,5) parameters
    ISres_35 = ISanalysis(E_S16, G_S16, 3, 5);
    all_meanT_ratios_low_35(idx) = mean(diff(E_S16)) / mean(diff(G_S16)); % Store (3,5) mean period ratio
    all_couplings_35_low(idx) = ISres_35.RM; % Store (3,5) coupling strength
    all_bestRatios_35_low(idx) = ISres_35.CIstats(1,1); % Store (3,5) best ratio

    % Perform IS analysis with (4,5) parameters
    ISres_45 = ISanalysis(E_S16, G_S16, 4, 5);
    all_meanT_ratios_low_45(idx) = mean(diff(E_S16)) / mean(diff(G_S16)); % Store (4,5) mean period ratio
    all_couplings_45_low(idx) = ISres_45.RM; % Store (4,5) coupling strength
    all_bestRatios_45_low(idx) = ISres_45.CIstats(1,1); % Store (4,5) best ratio
end

% Plotting the results to compare the robustness of the method for participants 31 to 59
figure;

% Plot Mean Period Ratios for different parameters
subplot(3,1,1);
plot(31:59, all_meanT_ratios_low_24, '-og', 'DisplayName', '(2,4)');
hold on;  % Hold the plot for superimposing the next curves
plot(31:59, all_meanT_ratios_low_35, '-ob', 'DisplayName', '(3,5)');
plot(31:59, all_meanT_ratios_low_45, '-or', 'DisplayName', '(4,5)');
hold off;  % Release the plot
xlabel('Participant');
ylabel('Mean Period Ratio');
title('Mean Period Ratio for Low DP Group with Different Parameters')
legend('show');  % Show legend with display names
grid on;

% Plot Mean Period Ratios for different parameters
subplot(3,1,2);
plot(31:59, all_couplings_24_low, '-og', 'DisplayName', '(2,4)');
hold on;  % Hold the plot for superimposing the next curves
plot(31:59, all_couplings_35_low, '-ob', 'DisplayName', '(3,5)');
plot(31:59, all_couplings_45_low, '-or', 'DisplayName', '(4,5)');
hold off;  % Release the plot
xlabel('Participant');
ylabel('N:M Coupling');
title('N:M Coupling for Low DP Group with Different Parameters')
legend('show');  % Show legend with display names
grid on;

% Plot Mean Period Ratios for different parameters
subplot(3,1,3);
plot(31:59, all_bestRatios_24_low, '-og', 'DisplayName', '(2,4)');
hold on;  % Hold the plot for superimposing the next curves
plot(31:59, all_bestRatios_35_low, '-ob', 'DisplayName', '(3,5)');
plot(31:59, all_bestRatios_45_low, '-or', 'DisplayName', '(4,5)');
hold off;  % Release the plot
xlabel('Participant');
ylabel('N:M Best Ratio');
title('N:M Best Ratio for Low DP Group with Different Parameters')
legend('show');  % Show legend with display names
grid on;
