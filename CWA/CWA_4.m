% TimZ, MSc, 2021, Imperial College.
% 11/1/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover the text in the given data file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% Data file 'yz6121_FastFading.mat'
% Beta_1 (3x1 Complex) = fading coefficients
% phase_shift (Integer) = delay used to generate the gold-sequence
% phi_mod (3x1 Complex) = angle phi for the QPSK constellation
% Xmatrix (5x7471 Complex) = text signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% delay_est = estimated delays of each path of the desired signal
% DOA_est = estimated azimuth and elevation of each path of the desired signal
% text_bit_est = estimated text bit stream
% text = estimated original text
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
clc;
clear;
close all;

%% Define Parameters
% Personal Setting
load('yz6121_FastFading.mat');


% Modulation
myphi = phi_mod/180 * pi ;  % Angle phi for the QPSK constellation, phi_mod is given in degrees

% Channel
text_sym_rec = Xmatrix.';  % Received text symbol stream
[~,n_ant] = size(text_sym_rec);  % Number of antennas
array = zeros(n_ant,3);  % A uniform circular array at point T_1_hat, (5 antennas and x/y/z positions)
phase_1st = pi/6;  % position of the first antenna
for i = 1:n_ant  % The unit of positions is half wavelength
    array(i,:) = [cos(phase_1st + 2*pi/5 * (i-1)), sin(phase_1st + 2*pi/5 * (i-1)), 0];
end
fading_coeff = Beta_1;  % The fading coefficient in each path
n_path = length(fading_coeff);  % Number of the paths
n_char = 60;  % Number of characters in the desired text-message

%% PN-codes
% primitive polynomials
primpoly_1 = [1 0 0 1 0 1]';
primpoly_2 = [1 0 1 1 1 1]';
d = phase_shift;  % The delay, used to generate the gold-sequence

% generate 2 m-sequences containing 0s & 1s and generate the gold sequence
MSeq_1 = fMSeqGen(primpoly_1);
MSeq_2 = fMSeqGen(primpoly_2);
GoldSeq = fGoldSeq(MSeq_1,MSeq_2,d);

%% Estimation
P_noise = 0.3;  % P_noise is about 0.3 in task 3 for 0 dB SNR
[delay_est, DOA_est] = fChannelEstimation_3_4(text_sym_rec, GoldSeq, n_path, array, P_noise);  % Estimate the delay

%% Demodulation
% Calculate the weight and retrun the weighted discretised signal
weighted_discretised_signal = spatiotemporal(text_sym_rec, GoldSeq, n_path, array, delay_est, DOA_est, fading_coeff);
text_bit_est = fDSQPSKDemodulator_4(weighted_discretised_signal, myphi);  % Demodulate the received signal

%% Show the recovered text
disp('----------------------------------------------------------------------------');
text = bit2int((reshape(text_bit_est, length(text_bit_est)/n_char, n_char)), 8);  % Convert it from binary stream to decimal stream
disp(['Recovered text = "' char(text) '"']);  
disp(['Estimated delay = [' num2str(delay_est') ']']);
disp(['Estimated DOA = [' num2str(reshape(DOA_est', 1, numel(DOA_est))) ']']);
disp('----------------------------------------------------------------------------');