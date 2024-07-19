% TimZ, MSc, 2021, Imperial College.
% 2/1/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TOA and TDOA Localisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% 4 received signals of 4 receivers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% r_m_TOA = estimated location obtained by TOA
% r_m_TDOA_est = estimated location obtained by TDOA
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

%% TOA Localisation

% Association Stage
t0 = 20 * Ts;  % Transmission time instant of the Tx signal
TOA = zeros(1,N);
d = zeros(1,N);
for i = 1:N
    x_Time = cell2mat(struct2cell(load(['t1_Rx' num2str(i) '.mat'])));  % Load the reveived signal by reveiver i
    TOA(i) = toa(x_Time, SNR, Pn) * Ts;  % Determine the time of arrival for each Rx
    d(i) = (TOA(i) - t0) * c;  % Distances/ranges
end

% Metric Fusion Stage
H_TOA = r(:,2:end)';
b_TOA = [sum(r(:,2).^2) - d(2)^2 + d(1)^2; 
         sum(r(:,3).^2) - d(3)^2 + d(1)^2; 
         sum(r(:,4).^2) - d(4)^2 + d(1)^2] / 2; 
r_m_TOA = H_TOA \ b_TOA;

%% Show the result (TOA)
disp('----------------------------------------------------------------------------');
disp(['The estimated location of Tx (TOA) is (' num2str(r_m_TOA(1)) ',' num2str(r_m_TOA(2)) ').']);
figure;
plot(r(1,:), r(2,:),'.b','MarkerSize',30);
hold on;
plot(r_m_TOA(1), r_m_TOA(2),'.r','MarkerSize',30);
draw_circle(r(:,1), d(1));
draw_circle(r(:,2), d(2));
draw_circle(r(:,3), d(3));
draw_circle(r(:,4), d(4));
axis equal
grid on;
xlabel('x');
ylabel('y');
legend('Rx','Tx');
title('Time of Arrival Localisation');

%% TDOA Localisation
% Association Stage
d_21 = (TOA(2) - TOA(1)) * c;  % Distances/ranges
d_31 = (TOA(3) - TOA(1)) * c;
d_41 = (TOA(4) - TOA(1)) * c;

% Metric Fusion Stage
% Solve the equation, where d_1_TDOA is the positive root
syms d_1_TDOA
H_TDOA = r(:,2:end)';
b_TDOA_est = [sum(r(:,2).^2) - d_21^2 - 2*d_21*d_1_TDOA;
              sum(r(:,3).^2) - d_31^2 - 2*d_31*d_1_TDOA;
              sum(r(:,4).^2) - d_41^2 - 2*d_41*d_1_TDOA] / 2;
r_m_TDOA = ((H_TDOA' * H_TDOA) \ H_TDOA') * b_TDOA_est;

equation = d_1_TDOA ^ 2 == r_m_TDOA' * r_m_TDOA;
root = double(solve(equation,d_1_TDOA));  % Find the root
d_1_est = root(root>0);  % Positive root is the estimated distance

b_TDOA_est = [sum(r(:,2).^2) - d_21^2 - 2*d_21*d_1_est;
              sum(r(:,3).^2) - d_31^2 - 2*d_31*d_1_est;
              sum(r(:,4).^2) - d_41^2 - 2*d_41*d_1_est] / 2;
r_m_TDOA_est = H_TDOA \ b_TDOA_est;

%% Show the result (TDOA)
disp('----------------------------------------------------------------------------');
disp(['The estimated location of Tx (TDOA) is (' num2str(r_m_TDOA_est(1)) ',' num2str(r_m_TDOA_est(2)) ').']);
disp('----------------------------------------------------------------------------');
figure;
plot(r(1,:),r(2,:),'.b','MarkerSize',30);
hold on;
plot(r_m_TDOA_est(1), r_m_TDOA_est(2),'.r','MarkerSize',30);
draw_circle(r(:,1), d_1_est);
draw_hyperbola(r(:,1), r(:,2), d_21);
draw_hyperbola(r(:,1), r(:,3), d_31);
draw_hyperbola(r(:,1), r(:,4), d_41);
axis equal
grid on;
xlabel('x');
ylabel('y');
legend('Rx','Tx');
title('Time Difference of Arrival Localisation');

tilefigs([0 0.4 0.7 1.0]);