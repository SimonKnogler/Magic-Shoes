%% intent to estimate the "generalized" (polyrhythmic) synchronization
%% heart - steps, from Simon time series of steps cycles (times) and heart max (ECG max times)
%% 1) ratio between mean periods, 2) use of Zelic's approach to get the ratio at which sync is best captured
%% and an estimate of the coupling (inspired by Mc Dermott et al. 2003),
%% sum of distances in return map of actual points to line of perfect synchrony.
% McDermott, W. J., Van Emmerik, R. E., & Hamill, J. (2003). Running
% training and adaptive strategies of locomotor-respiratory coordination.
% % European Journal of Applied Physiology, 89, 435–444. Zelic, G., Varoqui, D., Kim, J. et al. A flexible and accurate method to
% estimate the mode and stability of spontaneous coordinated behaviors: 
% The index-of-stability (IS) analysis. Behav Res 50, 182–194 (2018). https://doi.org/10.3758/s13428-017-0861-2


% follow example of High DP for Low DP group further down - currently with
% no proper output 

%%%%%%%%%%%%%%%
%%% High DP %%%
%%%%%%%%%%%%%%%

% load the data
E_high = ecgEventsTable_high
G_high = gaitEventsTable_high

% change from table to array
E_high = table2array(E_high); % n x 3, time of peaks in second row
x = G_high;% 'cause table2array pb to concatenate directly GaitPeaks
x(:,2) = [];
G_high = table2array(x);% n x 2, time of step in first row
clear x

% Initialize arrays to store values for all participants
all_meanT_ratios_high = zeros(1, 30);
all_couplings_high = zeros(1, 30);
all_bestRatios_high = zeros(1, 30);

% Loop over all participants
for part = 1:30
    % heart beats
    E2 = E_high(:,3) == part;
    E3 = E_high.*E2;
    E_S16 = E3(E3(:,2)~=0, 2);  % Filter and select second column
    clear E2 E3

    % gait events
    G2 = G_high(:,2) == part;
    G3 = G_high.*G2;
    G_S16 = G3(G3(:,1)~=0, 1);  % Filter and select first column
    clear G2 G3

    % Compute mean period ratio
    meanT_ratio = mean(diff(E_S16)) / mean(diff(G_S16));
    all_meanT_ratios_high(part) = meanT_ratio;

    % Perform IS analysis
    ISres = ISanalysis(E_S16, G_S16, 2, 0);
    bestRatio = ISres.CIstats(1,1);
    all_bestRatios_high(part) = bestRatio;
    coupling = ISres.RM;
    all_couplings_high(part) = coupling;
end

% Plotting the results
figure;
subplot(3,1,1);
plot(all_meanT_ratios_high, '-o');
xlabel('Participant');
ylabel('Mean Period Ratio');
grid on;

subplot(3,1,2);
plot(all_couplings_high, '-o');
xlabel('Participant');
ylabel('N:M Coupling');
grid on;

subplot(3,1,3);
plot(all_bestRatios_high, '-o');
xlabel('Participant');
ylabel('N:M Best Ratio');
grid on;


%%%%%%%%%%%%%%%
%%% Low DP %%%
%%%%%%%%%%%%%%%

% load the data
E_low = ecgEventsTable_low
G_low = gaitEventsTable_low

% change from table to array
E_low = table2array(E_low); % n x 3, time of peaks in second row
x = G_low;% 'cause table2array pb to concatenate directly GaitPeaks
x(:,2) = [];
G_low = table2array(x);% n x 2, time of step in first row
clear x

% Initialize arrays to store values for all participants
all_meanT_ratios_high = zeros(1, 30);
all_couplings_high = zeros(1, 30);
all_bestRatios_high = zeros(1, 30);

% Loop over all participants
for part = 1:30
    % heart beats
    E2 = E_low(:,3) == part;
    E3 = E_low.*E2;
    E_S16 = E3(E3(:,2)~=0, 2);  % Filter and select second column
    clear E2 E3

    % gait events
    G2 = G_low(:,2) == part;
    G3 = G_low.*G2; 
    G_S16 = G3(G3(:,1)~=0, 1);  % Filter and select first column
    clear G2 G3

    % Compute mean period ratio
    meanT_ratio = mean(diff(E_S16)) / mean(diff(G_S16));
    all_meanT_ratios_high(part) = meanT_ratio;

    % Perform IS analysis
    ISres = ISanalysis(E_S16, G_S16, 2, 0);
    bestRatio = ISres.CIstats(1,1);
    all_bestRatios_high(part) = bestRatio;
    coupling = ISres.RM;
    all_couplings_high(part) = coupling;
end

% Plotting the results
figure;
subplot(3,1,1);
plot(all_meanT_ratios_high, '-o');
xlabel('Participant');
ylabel('Mean Period Ratio');
grid on;

subplot(3,1,2);
plot(all_couplings_high, '-o');
xlabel('Participant');
ylabel('N:M Coupling');
grid on;

subplot(3,1,3);
plot(all_bestRatios_high, '-o');
xlabel('Participant');
ylabel('N:M Best Ratio');
grid on;
