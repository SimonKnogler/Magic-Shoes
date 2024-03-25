% Assume the data has been processed and the following variables are filled
% with the data for the two groups:
% - all_meanT_ratios for group 1-30
% - all_couplings for group 1-30
% - all_bestRatios for group 1-30
% - group2_meanT_ratios for group 31-59
% - group2_couplings for group 31-59
% - group2_bestRatios for group 31-59

% Now, plotting the results for both groups superimposed in one figure

% Plot Mean Period Ratios for both groups
figure;
subplot(3,1,1);
plot(1:30, all_meanT_ratios, '-ob', 'DisplayName', 'Group 1-30'); % Blue for group 1-30
hold on; % Hold on to plot on the same figure
plot(31:59, group2_meanT_ratios, '-or', 'DisplayName', 'Group 31-59'); % Red for group 31-59
xlabel('Participant');
ylabel('Mean Period Ratio');
title('Mean Period Ratio for Two Groups');
legend show; % Show legend
grid on;

% Plot N:M Coupling for both groups
subplot(3,1,2);
plot(1:30, all_couplings, '-ob', 'DisplayName', 'Group 1-30');
hold on; % Hold on to plot on the same figure
plot(31:59, group2_couplings, '-or', 'DisplayName', 'Group 31-59');
xlabel('Participant');
ylabel('N:M Coupling');
title('N:M Coupling for Two Groups');
legend show; % Show legend
grid on;

% Plot N:M Best Ratio for both groups
subplot(3,1,3);
plot(1:30, all_bestRatios, '-ob', 'DisplayName', 'Group 1-30');
hold on; % Hold on to plot on the same figure
plot(31:59, group2_bestRatios, '-or', 'DisplayName', 'Group 31-59');
xlabel('Participant');
ylabel('N:M Best Ratio');
title('N:M Best Ratio for Two Groups');
legend show; % Show legend
grid on;

% Adjust the x-axis limits to accommodate both groups
xlim(subplot(3,1,1), [1, 59]);
xlim(subplot(3,1,2), [1, 59]);
xlim(subplot(3,1,3), [1, 59]);
