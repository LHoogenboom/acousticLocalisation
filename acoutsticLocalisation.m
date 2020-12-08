%% Init
clc; clear; clf;
y.data = importdata("calibration.mat");
exp = importdata("experiment.mat");

%% 1
y.err = zeros(size(y.data));

y.avg = mean(y.data,2);

for m=1:7
    y.err(:,m) = y.data(:,m)-y.avg;
    mic.bias(m) = mean(y.err(:,m));
end
mic.var = var(y.err);

figure(1)
grid on;
hist(y.err(:,1))
title('Histogram of Measurement Error of Mic 1')  
set(gcf,'color','w')

clear y % cal data is now redundant

%% 2

for m=1:7
    y.cal(:,m) = exp.y(:,m)-mic.bias(m);
end

th_hat0 = [.1 .6 0];
maxiter = 100;
stds = diag(sqrt(mic.var));

th_hat =  th_hat0 ;zeros(length(y.cal(:,1))-1,3);
for i=1:length(y.cal(:,1))-1
    [th_hat(i+1,:), diagP(i+1,:)] = nls(y.cal(i,:),stds,th_hat(i,:),maxiter,exp.mic_locations);
end

plotresults(th_hat(:,1:2)',diagP,exp.mic_locations')

%% 3



%% Functions

function [th_hat, diagP] = nls(y,stds,th_hat0,maxiter,miclocs)
    i=0;
    th = th_hat0;
    th_hat = th_hat0;
    sig = pinv(stds*stds');
    while i<maxiter
        J = Jacobian(th_hat,miclocs);
        ftheta = f(th_hat,miclocs);
        e = y-ftheta';
        P = (J'*sig*J);
        dth = P^-1*J'*sig*e';
        th_hat = th_hat + dth';
        i = i+1;
    end
    th = [th_hat0 ;th_hat];
    diagP =diag(P^-1);
end

function dist = d(th_hat,micloc)
    dist = sqrt((th_hat(1:2)-micloc) * (th_hat(1:2)-micloc)');
end

function ftheta = f(th_hat,mic_locations)
    ftheta = zeros(length(mic_locations(:,1)),1);
    c = 343; % speed of sound in [m/s]
    for m=1:7
        ftheta(m) = (th_hat(3) + (1/c) * d(th_hat,mic_locations(m,:)));
    end
end

function dF = Jacobian(theta,miclocs)
    c = 343; % speed of sound in [m/s]
    dF = zeros(length(miclocs(:,1)),length(theta(1,:)));
    for m=1:7
        dF(m,:) = [(theta(1)-miclocs(m,1))/(c*d(theta,miclocs(m,:)))  (theta(2)-miclocs(m,2))/(c*d(theta,miclocs(m,:)))  1];
    end
end