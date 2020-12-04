%% Init
clear; clc;
ycal = importdata("calibration.mat");
exp = importdata("experiment.mat");


%% 1a
err = zeros(length(ycal(:,1)),1);

% least squared to make trend line for each mic. and that equals tau
% so calibration measurements-trend is microphone error.
errors = zeros(length(ycal(:,1)),7);

for i=1:7
   errors(:,i) = detrend(ycal(:,i));
end

avgMicErr = mean(errors,2);
bias = mean(avgMicErr);
var = var(avgMicErr);

[N, l] = hist(avgMicErr,20);
Wb=l(2)-l(1); % Bin width
Ny = length(avgMicErr); % Nr of samples
bar(l, N/(Ny*Wb));

figure(1)
hold on; grid on;
plot(t,ycal(:,1))
%% 2.1


%% Functions

function [trend] = lls(input)
    % ad hoc trend function
end

function [th_hat, diagP] = nls(yk,stds,th_hat0,maxiter,mic_locations)

end

function dF = Jacobian(theta,mic_locations)
    c = 343; % speed of sound in [m/s]
    
end

function ftheta = f(theta,mic_locations)
    c = 343; % speed of sound in [m/s]
    
end

