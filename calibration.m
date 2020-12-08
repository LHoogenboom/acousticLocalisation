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

%% 2 a
experiment = importdata("experiment.mat");
experiment.mic_locations = 100.*experiment.mic_locations;

% substracting the mic bias from measurements
calibrated = zeros(size(experiment.y));
for i=1:length(experiment.y(:,1))
    calibrated(i,:) = experiment.y(i,:)-mic.bias;
end

%% 2 b

% theta = [x y tau_k]
th_hat0 = [10 60 0];

%ftheta = (%?? ,experiment.mic_locations);

% min eps|sig --> min [y-f(theta)], because y = f(theta)+eps
% With f(theta) equals tau+d/c 
% y does not incluse bias(calibrated data)

%% plotting 2

figure(2)
hold on; grid on;
for i = 1:7
    plot(experiment.mic_locations(i,1),experiment.mic_locations(i,2), 'ro')
end

%%
test = Jacobian(th_hat0, experiment.mic_locations)

%% Functions
function [th_hat, diagP] = nls(yk,stds,th_hat0,maxiter,mic_locations)
    th_hat = [th_hat0 ; zeros(length(experiment.y(:,1))-1,3)];

end

function dF = Jacobian(theta,mic_locations)
    c = 343; % speed of sound in [m/s]
    dF = zeros(length(mic_locations(:,1)),length(theta(1,:)));
    for m=1:7
        dF(m,:) = [(theta(1)-mic_locations(m,1))/(c*d(theta,mic_locations(m,:)))  (theta(2)-mic_locations(m,2))/(c*d(theta,mic_locations(m,:)))  1];
    end
end

function ftheta = f(th_hat,mic_locations)
    % size th_hat = 3 by 7
    % size ftheta = 7 by 1
    ftheta = zeros(length(mic_locations(:,1)),1);
    size(ftheta)
    c = 343; % speed of sound in [m/s]
    
    for m=1:7
        % f(theta)  = tau_k    + (1/c) * d
        %d_m(theta) = theta(3) + (1/c) * (theta-p_m)^T * (theta-p_m)
        %ftheta(m) = th_hat(3) + (1/c) * sqrt((th_hat(1:2)-mic_locations(m,:)) * (th_hat(1:2)-mic_locations(m,:))')
        ftheta(m) = th_hat(3) + (1/c) * d(th_hat,mic_locations(m,:))
    end
end

function dist = d(th_hat,micloc)
    dist = sqrt((th_hat(1:2)-micloc) * (th_hat(1:2)-micloc)')
end