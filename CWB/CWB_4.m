% TimZ, MSc, 2021, Imperial College.
% 14/1/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LAA Localisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% 4 received signals of 4 receivers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% r_m = estimated location obtained by LAA
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

%% LAA Localisation
% Association Stage
lambda = zeros(1,N);
for i = 1:N
    % 4 measurements for 4 different reference points
    x_LAA = cell2mat(struct2cell(load(['t4_Xmatrix_LAA_' num2str(i) '.mat'])));  % Load the reveived signal by reveiver i
    cov = x_LAA * x_LAA' / size(x_LAA,2);  % Calculate the covariance matrix
    eig_val = sort(eig(cov),'descend');  % Find the eigen values and sort
    y = eig_val(1);   % The largest eigen value
    lambda(i) = y - Pn;
end
K = ((lambda(2:end) ./ lambda(1)) .^ (1/(2*a)))';

% Metric fusion stage (the alternative one)
H = [2*(ones(N-1,1) * r(:,1)' - r(:,2:end)'), ones(N-1,1) - K.^2];
b = norm(r(:,1)).^2 * ones(N-1,1) - (r(1,2:end)').^2 - (r(2,2:end)').^2;
r_m = H \ b;

% Calculate the centre and radius of reference circle
r_1 = 1 / (1 - K(1)^2) * r(:,2) - (K(1)^2) / (1 - K(1)^2) * r(:,1);
d_1 = abs(K(1) / (1 - K(1)^2)) * norm(r(:,2) - r(:,1));
r_2 = 1 / (1 - K(2)^2) * r(:,3) - (K(2)^2) / (1 - K(2)^2) * r(:,1);
d_2 = abs(K(2) / (1 - K(2)^2)) * norm(r(:,3) - r(:,1));
r_3 = 1 / (1 - K(3)^2) * r(:,4) - (K(3)^2) / (1 - K(3)^2) * r(:,1);
d_3 = abs(K(3) / (1 - K(3)^2)) * norm(r(:,4) - r(:,1));

%% Show the result (LAA)
disp('----------------------------------------------------------------------------');
disp(['The estimated location of Tx (LAA) is (' num2str(r_m(1)) ',' num2str(r_m(2)) ').']);
disp('----------------------------------------------------------------------------');
figure;
plot(r(1,:),r(2,:),'.b','MarkerSize',30);
hold on;
plot(r_m(1), r_m(2),'.r','MarkerSize',30);
draw_circle(r_1,d_1);
draw_circle(r_2,d_2);
draw_circle(r_3,d_3);
plot(r_1(1),r_1(2),'.k','MarkerSize',10);
plot(r_2(1),r_2(2),'.k','MarkerSize',10);
plot(r_3(1),r_3(2),'.k','MarkerSize',10);
text(r_1(1)+10,r_1(2)+10,'r_1');
text(r_2(1)+10,r_2(2)+10,'r_2');
text(r_3(1)+10,r_3(2)+10,'r_3');
axis equal
grid on;
xlabel('x');
ylabel('y');
legend('Rx','Tx');
title('Large Aperture Array Localisation');

figure;
plot(r(1,:),r(2,:),'.b','MarkerSize',30);
hold on;
plot(r_m(1), r_m(2),'.r','MarkerSize',30);
draw_circle(r_1,d_1);
draw_circle(r_2,d_2);
draw_circle(r_3,d_3);
plot(r_1(1),r_1(2),'.k','MarkerSize',10);
plot(r_2(1),r_2(2),'.k','MarkerSize',10);
plot(r_3(1),r_3(2),'.k','MarkerSize',10);
axis equal
axis([-100 200 -150 150])
grid on;
xlabel('x');
ylabel('y');
legend('Rx','Tx');
title('Large Aperture Array Localisation');

tilefigs([0 0.4 0.7 1.0]);