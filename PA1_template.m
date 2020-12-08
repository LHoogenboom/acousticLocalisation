clear; clf; clc;
%% A1
load('calibration.mat')
%% A2
load('experiment.mat')
maxiter = 100;
for k = 1:size(y,1) 
    % initial estimate of theta
    th_hat0 = [10 60 0];
    % non linear LS estimate 
    [th_hat(k,:),diagP(k,:)] = nls(yk,stds,th_hat0,maxiter,mic_locations);
end

%% A3

%% A4

%% Functions
function [th_hat, diagP] = nls(yk,stds,th_hat0,maxiter,mic_locations)

end

function dF = Jacobian(theta,mic_locations)
    c = 343; % speed of sound in [m/s]
    
end

function ftheta = f(theta,mic_locations)
    c = 343; % speed of sound in [m/s]
    
end
