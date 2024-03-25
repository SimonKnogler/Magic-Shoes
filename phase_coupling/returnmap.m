%% Intent to estimate the "generalized" (polyrhythmic) synchronization
%% heart - steps, from Simon time series of steps cycles (times) and heart max (ECG max times)
%% 1) ratio between mean periods, 2) use of Zelic's approach to get the ratio at which sync is best captured
%% and an estimate of the coupling (inspired by Mc Dermott et al. 2003),
%% sum of distances in return map of actual points to line of perfect synchrony.

% McDermott, W. J., Van Emmerik, R. E., & Hamill, J. (2003). Running
% training and adaptive strategies of locomotor-respiratory coordination.
% European Journal of Applied Physiology, 89, 435–444.

% Zelic, G., Varoqui, D., Kim, J. et al. A flexible and accurate method to
% estimate the mode and stability of spontaneous coordinated behaviors:
% The index-of-stability (IS) analysis. Behav Res 50, 182–194 (2018). https://doi.org/10.3758/s13428-017-0861-2

% Load the data
E_high = ecgEventsTable_high;
G_high = gaitEventsTable_high;

% Change from table to array
E_high = table2array(E_high); % n x 3, time of peaks in second row
x = G_high; % 'cause table2array pb to concatenate directly GaitPeaks
x(:,2) = [];
G_high = table2array(x); % n x 2, time of step in first row
clear x

% Specify the participant number
part = 16; % Change this to the desired participant number

% Heart beats
E2 = E_high(:,3) == part;
E3 = E_high.*E2;
E_S16 = E3(E3(:,2)~=0, 2); % Filter and select second column
clear E2 E3

% Gait events
G2 = G_high(:,2) == part;
G3 = G_high.*G2;
G_S16 = G3(G3(:,1)~=0, 1); % Filter and select first column
clear G2 G3

% Compute mean period ratio
meanT_ratio = mean(diff(E_S16)) / mean(diff(G_S16));

% Perform IS analysis
ISres = ISanalysis(E_S16, G_S16, 3, 0);
bestRatio = ISres.CIstats(1,1);
coupling = ISres.RM;
dwellTime = ISres.RM; % Dwell Time measure

% Extract the return map data
RP = ISres.ReturnMap.RP;
Dn = ISres.ReturnMap.Dn;

% Plot the return map and Euclidean distances
fig1 = figure('Name', strcat('participant-', num2str(part), '_return_map_and_distances'));
subplot(1, 2, 1)
plot(RP(1:end-1), RP(2:end), '.');
xlabel('RP(n)');
ylabel('RP(n+1)');
title(sprintf('Return Map - Participant %d', part));
xlim([0 2*pi]);
ylim([0 2*pi]);
line([0 2*pi], [0 2*pi], 'Color', 'r', 'LineStyle', '--');

subplot(1, 2, 2)
plot(Dn, '-b');
ylabel('Euclidean Distance');
xlabel('Index');
title(sprintf('Euclidean Distance - Participant %d', part));

% Save the figure as a PNG file
filename1 = sprintf('participant-%d_return_map_and_distances.png', part);
print(fig1, '-dpng', '-r300', filename1);

% Plot the return map and Euclidean distances
figSig = figure('Name', strcat('participant-', num2str(part)), 'Units', 'normalized', 'OuterPosition', [0 0 1 1]); % Fullscreen figure
set(figSig, 'DefaultTextFontSize', 16); % Increase the default font size

subplot(3, 1, 1)
plot(diff(E_S16), 'o'), ylabel('Heart periods', 'FontSize', 18), ylim([0 2])

subplot(3, 1, 2)
plot(diff(G_S16), 'o'), ylabel('Gait periods', 'FontSize', 18), ylim([0 2])

subplot(3, 4, 9)
bar(meanT_ratio), ylabel('Mean Period Ratio', 'FontSize', 16), ylim([0 1])

subplot(3, 4, 10)
bar(coupling), ylabel('N:M Coupling', 'FontSize', 16), ylim([0 1]) % Fixed y-axis limits

subplot(3, 4, 11)
bar(bestRatio), ylabel('N:M Best Ratio', 'FontSize', 16)

subplot(3, 4, 12)
bar(dwellTime), ylabel('Dwell Time', 'FontSize', 16) % Plot Dwell Time measure

subplot(3, 4, 12)
histogram(diff(E_S16), 16), hold on, histogram(diff(G_S16), 16)
title('Hist Heart/Gait Periods', 'FontSize', 16)

% Save the figure as a PNG file
filename = sprintf('participant-%d_all_plots.png', part);
print(figSig, '-dpng', '-r300', filename);