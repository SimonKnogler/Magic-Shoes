% Assuming ecgEventsTable_high and gaitEventsTable_high are loaded elsewhere

% Convert the tables to arrays for processing
E_high = table2array(ecgEventsTable_high); % n x 3, time of peaks in second row
G_high = table2array(gaitEventsTable_high(:, [1, 3])); % Assuming we only need columns 1 and 3

% Initialize arrays to store values for participants 1 to 30
all_meanT_ratios_high_24 = zeros(1, 30); % For (2,4) parameters
all_meanT_ratios_high_35 = zeros(1, 30); % For (3,5) parameters
all_meanT_ratios_high_45 = zeros(1, 30); % For (4,5) parameters
all_couplings_24_high = zeros(1, 30); % For (2,4) parameters
all_couplings_35_high = zeros(1, 30); % For (3,5) parameters
all_couplings_45_high = zeros(1, 30); % For (4,5) parameters
all_bestRatios_24_high = zeros(1, 30); % For (2,4) parameters
all_bestRatios_35_high = zeros(1, 30); % For (3,5) parameters
all_bestRatios_45_high = zeros(1, 30); % For (4,5) parameters

% Loop over participants 1 to 30
for part = 1:30
    % Heart beats
    E_S16 = E_high(E_high(:,3) == part, 2); % Select rows for participant and second column

    % Gait events
    G_S16 = G_high(G_high(:,2) == part, 1); % Select rows for participant and first column

    % Perform IS analysis with (2,4) parameters
    ISres_24 = ISanalysis(E_S16, G_S16, 2, 4);
    all_meanT_ratios_high_24(part) = mean(diff(E_S16)) / mean(diff(G_S16)); % Store (2,4) mean period ratio
    all_couplings_24_high(part) = ISres_24.RM; % Store (2,4) coupling strength
    all_bestRatios_24_high(part) = ISres_24.CIstats(1,1); % Store (2,4) best ratio

    % Perform IS analysis with (3,5) parameters
    ISres_35 = ISanalysis(E_S16, G_S16, 3, 5);
    all_meanT_ratios_high_35(part) = mean(diff(E_S16)) / mean(diff(G_S16)); % Store (3,5) mean period ratio
    all_couplings_35_high(part) = ISres_35.RM; % Store (3,5) coupling strength
    all_bestRatios_35_high(part) = ISres_35.CIstats(1,1); % Store (3,5) best ratio

    % Perform IS analysis with (4,5) parameters
    ISres_45 = ISanalysis(E_S16, G_S16, 4, 5);
    all_meanT_ratios_high_45(part) = mean(diff(E_S16)) / mean(diff(G_S16)); % Store (4,5) mean period ratio
    all_couplings_45_high(part) = ISres_45.RM; % Store (4,5) coupling strength
    all_bestRatios_45_high(part) = ISres_45.CIstats(1,1); % Store (4,5) best ratio
end

% Plotting

% Plotting the results to compare the robustness of the method
figure;

% Plot Mean Period Ratios for different parameters
subplot(3,1,1);
plot(1:30, all_meanT_ratios_high_24, '-og', 'DisplayName', '(2,4)');
hold on;  % Hold the plot for superimposing the next curves
plot(1:30, all_meanT_ratios_high_35, '-ob', 'DisplayName', '(3,5)');
plot(1:30, all_meanT_ratios_high_45, '-or', 'DisplayName', '(4,5)');
hold off;  % Release the plot
xlabel('Participant');
ylabel('Mean Period Ratio');
title('Mean Period Ratio for High DP Group with Different Parameters');
legend('show');  % Show legend with display names
grid on;

% Plot N:M Couplings for different parameters
subplot(3,1,2);
plot(1:30, all_couplings_24_high, '-og', 'DisplayName', '(2,4)');
hold on;  % Hold the plot for superimposing the next curves
plot(1:30, all_couplings_35_high, '-ob', 'DisplayName', '(3,5)');
plot(1:30, all_couplings_45_high, '-or', 'DisplayName', '(4,5)');
hold off;  % Release the plot
xlabel('Participant');
ylabel('N:M Coupling');
title('N:M Coupling for High DP Group with Different Parameters');
legend('show');  % Show legend with display names
grid on;

% Plot N:M Couplings for different parameters
subplot(3,1,3);
plot(1:30, all_bestRatios_24_high, '-og', 'DisplayName', '(2,4)');
hold on;  % Hold the plot for superimposing the next curves
plot(1:30, all_bestRatios_35_high, '-ob', 'DisplayName', '(3,5)');
plot(1:30, all_bestRatios_45_high, '-or', 'DisplayName', '(4,5)');
hold off;  % Release the plot
xlabel('Participant');
ylabel('N:M Best Ratio');
title('N:M Best Ratio for High DP Group with Different Parameters');
legend('show');  % Show legend with display names
grid on;

