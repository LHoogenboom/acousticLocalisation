%% Init
clear; clc;
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