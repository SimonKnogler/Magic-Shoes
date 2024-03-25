% Prepare the figure
figure;

% Plot Mean Period Ratios for both groups
subplot(3,1,1); % This creates a subplot in a 3x1 grid, at position 1
hold on; % This command allows us to plot multiple datasets on the same subplot
plot(all_meanT_ratios_high, '-or'); % Plot the "high" group data in red with circle markers
plot(all_meanT_ratios_low, '-ob'); % Plot the "low" group data in blue with circle markers
xlabel('Data Point'); % Adjusting xlabel since we are disregarding the participant number
ylabel('Mean Period Ratio');
legend('High DP Group', 'Low DP Group'); % Add a legend to distinguish between groups
grid on; % Enable grid for better readability

% Plot N:M Couplings for both groups
subplot(3,1,2); % This creates a subplot in a 3x1 grid, at position 2
hold on; % Allows plotting multiple datasets on the same subplot
plot(all_couplings_high, '-or'); % Plot the "high" group data in red
plot(all_couplings_low, '-ob'); % Plot the "low" group data in blue
xlabel('Data Point'); % Again, adjusting xlabel
ylabel('N:M Coupling');
legend('High DP Group', 'Low DP Group'); % Add a legend
grid on; % Enable grid

% Plot N:M Best Ratios for both groups
subplot(3,1,3); % This creates a subplot in a 3x1 grid, at position 3
hold on; % Allows multiple datasets on the same subplot
plot(all_bestRatios_high, '-or'); % Plot "high" group data in red
plot(all_bestRatios_low, '-ob'); % Plot "low" group data in blue
xlabel('Data Point'); % Adjusting xlabel
ylabel('N:M Best Ratio');
legend('High DP Group', 'Low DP Group'); % Add a legend
grid on; % Enable grid

% Final adjustments
hold off; % No more plots to add
