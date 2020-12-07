%% Init
clear; clc; clf;
mic.y.data = importdata("calibration.mat");

%% 1
mic.y.err = zeros(size(mic.y.data));

for i=1:7
   mic.y.err(:,i) = detrend(mic.y.data(:,i)); 
end

mic.y.errAvg = mean(mic.y.err,2);
mic.bias = mean(mic.y.err,1);
mic.var = var(mic.y.err,1);


figure(1)
hold on; grid on;
title('Histogram of Measurement Error of Average Mic')
xlabel('Error Size (bin width 4*10^{-7})')
ylabel('Prevelance of Error')
[N, l] = hist(mic.y.errAvg,20);
Wb=l(2)-l(1); % Bin width
Ny = length(mic.y.errAvg); % Nr of samples
bar(l, N/(Ny*Wb));
% clear i l N Ny Wb;

%% 2
experiment = importdata("experiment.mat");
experiment.mic_locations = 100.*experiment.mic_locations;

% substracting the mic bias from measurements
calibrated = zeros(size(experiment.y));
for i=1:length(experiment.y(:,1))
    calibrated(i,:) = experiment.y(i,:)-mic.bias;
end

figure(2)
hold on; grid on;
for i = 1:7
    plot(experiment.mic_locations(i,1),experiment.mic_locations(i,2), 'ro')
end

% theta = [x y t_pulse]
th_hat0 = [10 60 0];
th_hat = [th_hat0 ; zeros(length(experiment.y(:,1))-1,3)];

%ftheta = (%?? ,experiment.mic_locations);

% min eps|sig --> min [y-f(theta)], because y = f(theta)+eps
% With f(theta) equals tau+d/c 
% y does not incluse bias(calibrated data)

%% Functions
function [th_hat, diagP] = nls(yk,stds,th_hat0,maxiter,mic_locations)
    th_hat = [th_hat0 ; zeros(length(experiment.y(:,1))-1,3)];

end

function dF = Jacobian(theta,mic_locations)
    c = 343; % speed of sound in [m/s]
    
end

function ftheta = f(theta,mic_locations)
    c = 343; % speed of sound in [m/s]
    %d_m(theta)  =((P_r,m|x-theta(1))^2+(P_r,m|y-theta(2))^2)^1/2
    %ftheta = t_k +_d(theta)/c
    
end