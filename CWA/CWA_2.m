% TimZ, MSc, 2021, Imperial College.
% 30/12/2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover the desired input photo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% 3 photos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% delay_est = estimated delays of each path of the desired signal
% photo_bit_est = estimated photo bit stream
% recovered photo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
clc;
clear;
close all;

%% Define Parameters
% Personal Setting
X = 26;  % alphabetical order of the 1st letter of my surname (Z -> z -> 26)
Y = 25;  % alphabetical order of the 1st letter of my formal firstname (Y -> y -> 25)

% Photos
% The maximum number of bits in a photo
% 160 represents maximum width
% 112 represents maximum height
% 3 represents the R,G,B components of the photos
% each colour component has 256 values [0-255], hence use 8 bits
b_max = 160 * 112 * 3 * 8;  

% Modulation
myphi = (X + 2*Y) * pi / 180; % angle phi for the QPSK constellation

% Channel
n_signal = 3;  % number of transmitted signals
n_desired = 1;   % index of the desired signal
array = [0 0 0];  % Position of the input of the receiver, point T_1_hat
paths = [3; 1; 1];   % Each signal only has 1 path in Task 1
delay = [mod(X+Y,4); 4+mod(X+Y,5); 9+mod(X+Y,6); 8; 13];  % The relative delay in each path
fading_coeff = [0.8; 0.4*exp(-40*1i/180*pi); 0.8*exp(80*1i/180*pi); 0.5; 0.2];    % The fading coefficient in each path
DOA = [30 0; 45 0; 20 0; 80 0; 150 0]; % Direction of arrival (azimuth, elevation), has no effect in Task 2

% Noise
n_SNR = 2;  % number of different noises
SNR_dB = [0 40];
SNR = 10.^(SNR_dB/10);

%% PN-codes
% primitive polynomials
primpoly_1 = [1 0 0 1 1]';
primpoly_2 = [1 1 0 0 1]';
M_Length = 2 ^ (length(primpoly_1)-1) - 1; % length of the generated m-sequences
d_min = 1 + mod((X+Y),12);  % The minimum delay, used to generate the gold-sequence

% generate 2 m-sequences containing 0s & 1s
MSeq_1 = fMSeqGen(primpoly_1);
MSeq_2 = fMSeqGen(primpoly_2);

% generate (n_signal) gold-sequences containing 1s & -1s, forming a (M_Length * n_signal) matrix
GoldSeq = zeros(M_Length,n_signal);
temp = d_min; % delay start from the minimum delay
while(1)
    GoldSeq(:,1) = fGoldSeq(MSeq_1, MSeq_2, temp);  % generate the first gold-sequence with delay (temp)
    if sum(GoldSeq(:,1)) == -1  % check whether the first gold-sequence is "balanced"
        for i = 2:n_signal  % if true, generate the remaining gold-sequences with delays (temp+1) & (temp+2) 
            GoldSeq(:,i) = fGoldSeq(MSeq_1, MSeq_2, temp + i - 1);
        end
        break
    else
        temp = temp + 1;
    end
end

%% Transmitter
% Initialize matrices to contain the information of n_signal photos
w_infor = zeros(1,n_signal);  % width of 3 photos
h_infor = zeros(1,n_signal);  % height of 3 photos
b_infor = zeros(1,n_signal);  % number of bits in each photo, (w_infor * h_infor * z_rgb * b_colour)
photo_bit = zeros(b_max,n_signal);  % the bit streams of 3 photos
photo_sym = zeros(b_max*M_Length/2,n_signal);  % the transmitted symbols from the transmitter
for i = 1:n_signal
    % Read the image and get the required information
    [photo_bit(:,i), w_infor(i), h_infor(i)] = fImageSource([num2str(i),'.png'], b_max);
    b_infor(i) = w_infor(i) * h_infor(i) * 3 * 8;

    % Modulate the bit streams of 3 photos
    photo_sym(:, i) = fDSQPSKModulator(photo_bit(:,i), GoldSeq(:,i), myphi);
end
% Show the original photos
fImageSink(photo_bit, b_infor, w_infor, h_infor);
title('Original Photos')

%% Channel & Receiver
disp('----------------------------------------------------------------------------');
disp(['Desired source is Source ' num2str(n_desired)]);
for i = 1:n_SNR
    % Passing through the channel
    [photo_sym_rec, ~] = fChannel(photo_sym, paths, array, SNR(i), delay, DOA, fading_coeff);

    % Arriving at the receiver, receive symbol streams
    delay_est = fChannelEstimation_1_2(photo_sym_rec(:,n_desired), GoldSeq, paths);  % Estimate the delay
    photo_bit_est = fDSQPSKDemodulator_1_2_3(photo_sym_rec(:,n_desired), GoldSeq, paths, delay_est, fading_coeff, myphi);  % Demodulate the received signal with estimated delay 

    % calculate bit error rate of the desired signal
    BER = sum(mod(photo_bit_est(:, n_desired) + photo_bit(:, n_desired), 2)) / b_infor(1);
    
    disp('----------------------------------------------------------------------------');
    disp(['SNR = ' num2str(SNR_dB(i)) ' dB']);
    disp(['Actual delay for each path = [' num2str(delay(1:paths(1))') ']']);
    disp(['Estimated delay for each path = [' num2str(delay_est') ']']);
    disp(['Bit Error Rate = ' num2str(BER*100) '%']);

    % Show the recovered photos
    fImageSink(photo_bit_est(:, n_desired), b_infor, w_infor, h_infor);
    title(['SNR = ' num2str(SNR_dB(i)) ' dB'])
end
disp('----------------------------------------------------------------------------');

tilefigs([0 0.4 1.0 1.0]);