%% How to visualize n:m phase synchronization
% %% using the synchrogram method from Rosenblum et al.
% Schäfer, C., Rosenblum, M. G., Abel, H. H., & Kurths, J. (1999). 
% Synchronization in the human cardiorespiratory system. Physical Review E, 60(1), 857.

clear all; close all

omega_x = 1;omega_y = 5;
t = 0:0.01:100;% time vector
x = 1*sin(omega_x*t); y = 1*sin(omega_y*t);% 2 waves

phi_x = angle(hilbert(x));% wave of the slowest wave (x)
[pval, pdates] = findpeaks(y);% peaks values and dates of peaks of the fastest wave (y)
sync = phi_x(pdates);% values of phase of the slowest at dates of peaks of the fastest

figSig = figure('Name','Synchrogram');
subplot(211)
plot(x,'r'),hold on,plot(y,'b')
subplot(212)
plot(sync,'o')