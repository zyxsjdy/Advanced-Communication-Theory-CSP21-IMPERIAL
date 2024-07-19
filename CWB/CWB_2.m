% TimZ, MSc, 2021, Imperial College.
% 4/1/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RSS Localisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% 4 received signals of 4 receivers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% r_m = estimated location obtained by RSS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
clc;
clear;
close all;

%% Define Parameters
r = [0  60 100 60;   % Locations of Rxs, (x,y)
     0 -88   9 92];  % z is 0 hence ignored
Fc = 2.4 * 10^9;     % Carrier frequency
Tcs = 5 * 10^-9;     % Symbol duration
N = 4;               % Number of Rx
Pn = 5;              % Noise power (dB)
c = 3 * 10^8;        % Propagation speed
a = 2;               % Path Loss exponent
SNR = 20;            % SNR (dB)
Ts = 5 * 10^-9;      % Sampling period
l = c/Fc;            % Wave length

%% RSS Localisation
% Association Stage
P_Tx = 10^(150/10) * 10^(-3);  % Power of received signal, 150dBm
d = zeros(1,N);
for i = 1:N
    x_RSS = cell2mat(struct2cell(load(['t2_Rx' num2str(i) '.mat'])));  % Load the reveived signal by reveiver i
    P_Rx = x_RSS * x_RSS' / size(x_RSS,2);  % Power of received signal
    d(i) = sqrt(P_Tx/P_Rx) * l/(4*pi);  % Distances/ranges, find using Friis transmission equation
end

% Metric Fusion Stage
H = r(:,2:end)';
b = [sum(r(:,2).^2) - d(2)^2 + d(1)^2; 
     sum(r(:,3).^2) - d(3)^2 + d(1)^2; 
     sum(r(:,4).^2) - d(4)^2 + d(1)^2] / 2; 
r_m = H \ b;

%% Show the result (RSS)
disp('----------------------------------------------------------------------------');
disp(['The estimated location of Tx (RSS) is (' num2str(r_m(1)) ',' num2str(r_m(2)) ').']);
disp('----------------------------------------------------------------------------');
figure;
plot(r(1,:),r(2,:),'.b','MarkerSize',30);
hold on;
plot(r_m(1),r_m(2),'.r','MarkerSize',30);
draw_circle(r(:,1),d(1));
draw_circle(r(:,2),d(2));
draw_circle(r(:,3),d(3));
draw_circle(r(:,4),d(4));
axis equal
grid on;
xlabel('x');
ylabel('y');
legend('Rx','Tx');
title('Received Signal Strength Localisation');

tilefigs([0 0.4 0.4 1.0]);