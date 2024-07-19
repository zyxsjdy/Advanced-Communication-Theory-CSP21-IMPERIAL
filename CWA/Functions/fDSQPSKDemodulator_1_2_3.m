% TimZ, MSc, 2021, Imperial College.
% 5/1/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform demodulation of the received data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% photo_sym_rec (Fx1 Complex) = R channel symbol chips received
% GoldSeq (Wx1 Integers) = W bits of 1's and 0's representing the gold
% sequence of the desired signal to be used in the demodulation process
% paths (Mx1 Integers) = Number of paths for each source in the system.
% For example, if 3 sources with 1, 3 and 2 paths respectively then
% paths = [1;3;2]
% delay_est = Vector of estimates of the delays of each path of the
% desired signal
% fading_coeff (Cx1 Integers) = Fading Coefficient for each path in the
% system starting with source 1
% myphi (Integer) = Angle index in degrees of the QPSK constellation points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% photo_bit_est (Px1 Integers) = P demodulated bits of 1's and 0's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [photo_bit_est] = fDSQPSKDemodulator_1_2_3(photo_sym_rec, GoldSeq, paths, delay_est, fading_coeff, myphi)

[d_max, ~] = size(GoldSeq);  % Get the maximum relative delay
n_sym = (length(photo_sym_rec)-d_max) / d_max;  % The number of symbols in the photo

% For the desired source 1
photo_sym_path = zeros(n_sym,paths(1));  % Initialize the symbol streams in all paths
for i = 1:paths(1)
    photo_sym_syn = circshift(photo_sym_rec, -delay_est(i));  % Synchronize the input symbols
    photo_sym_path(:,i) = reshape(photo_sym_syn(1:length(photo_sym_rec)-d_max), d_max, n_sym).' * GoldSeq(:,1);  % Despread the symbols
end
QPSK = sum(photo_sym_path .* fading_coeff(1:paths(1))', 2);  % Obtain the symbol stream for the desired source

% Quantization
photo_bit2 = zeros(2,n_sym);  % The matrix used to store the decoded bits, each column is used for one QPSK symbol
coordinate = [myphi; myphi+pi/2; myphi-pi; myphi-pi/2];
for i = 1:n_sym
    [~,Index] = min(abs(angle(QPSK(i)) - coordinate));  % Find the closest constellation point
    switch Index
        case 1
            photo_bit2(:,i) = [0;0];
        case 2
            photo_bit2(:,i) = [0;1];
        case 3
            photo_bit2(:,i) = [1;1];
        case 4
            photo_bit2(:,i) = [1;0];
        otherwise
            error('Error occurred')
    end
end
photo_bit_est = reshape(photo_bit2,[],1);  % Reshape to a bit stream
end