% TimZ, MSc, 2021, Imperial College.
% 10/1/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DOA Localisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% 4 received signals of 4 receivers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% r_m = estimated location obtained by DOA
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
array = [0.1250,       0, 0;
         0.0625,  0.1083, 0;
        -0.0625,  0.1083, 0;
        -0.1250,       0, 0;
        -0.0625, -0.1083, 0;
         0.0625, -0.1083, 0] / (l/2); % transfer from meter to half wavelength

%% DOA Localisation
% Association Stage
% Estimate DOA
DOA = zeros(N,1);
for i = 1:N
    DOA(i) = doa(cell2mat(struct2cell(load(['t3_Xmatrix_' num2str(i) '_DFarray.mat']))), array, Pn);
end
% Find the angles and distances, graphic description is given in the report
y_1 = DOA(2) - DOA(1);
y_3 = DOA(4) - DOA(3);
d_12 = r(:,2) - r(:,1);
d_34 = r(:,4) - r(:,3);
phi_12 = atan(d_12(2) / d_12(1)) / pi * 180 + 180;
phi_34 = atan(d_34(2) / d_34(1)) / pi * 180 + 180;
bet_1 = DOA(1) + 180 - phi_12;
bet_2 = phi_12 - DOA(2);
bet_3 = DOA(3) - phi_34;
bet_4 = phi_34 - DOA(4) + 180;
d(1) = norm(d_12) / sin(y_1/180*pi) * sin(bet_2/180*pi);
d(2) = norm(d_12) / sin(y_1/180*pi) * sin(bet_1/180*pi);
d(3) = norm(d_34) / sin(y_3/180*pi) * sin(bet_4/180*pi);
d(4) = norm(d_34) / sin(y_3/180*pi) * sin(bet_3/180*pi);

% Metric fusion stage (the alternative one)
b = zeros(2*N,1);
for i = 1:N
    b(2*i-1:2*i) = [r(1,i) + d(i) * cos(DOA(i)*pi/180),  r(2,i) + d(i) * sin(DOA(i)*pi/180)];
end
H = kron(ones(N,1),eye(2));
r_m = H \ b;

%% Show the result (RSS)
disp('----------------------------------------------------------------------------');
disp(['The estimated location of Tx (DOA) is (' num2str(r_m(1)) ',' num2str(r_m(2)) ').']);
disp('----------------------------------------------------------------------------');
figure;
plot(r(1,:),r(2,:),'.b','MarkerSize',30);
hold on;
plot(r_m(1),r_m(2),'.r','MarkerSize',30);
draw_line(r(:,1),DOA(1)*pi/180);
draw_line(r(:,2),DOA(2)*pi/180);
draw_line(r(:,3),DOA(3)*pi/180);
draw_line(r(:,4),DOA(4)*pi/180);
axis([-100 200 -150 150])
grid on;
xlabel('x');
ylabel('y');
legend('Rx','Tx');
title('Direction of Arrival Localisation');

tilefigs([0 0.4 0.4 1.0]);