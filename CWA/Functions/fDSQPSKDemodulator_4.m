% TimZ, MSc, 2021, Imperial College.
% 11/1/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform demodulation of the received data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% weighted_discretised_signal (Fx1 Complex) = F weighted discretised symbols
% myphi (Integer) = Angle index in degrees of the QPSK constellation points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% photo_bit_est (Px1 Integers) = P demodulated bits of 1's and 0's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [photo_bit_est] = fDSQPSKDemodulator_4(weighted_discretised_signal, myphi)

% Quantization
n_sym = length(weighted_discretised_signal);  % The number of symbols in the photo
photo_bit2 = zeros(2,n_sym);  % The matrix used to store the decoded bits, each column is used for one QPSK symbol
coordinate = [myphi; myphi+pi/2; myphi-pi; myphi-pi/2];
for i = 1:n_sym
    [~,Index] = min(abs(angle(weighted_discretised_signal(i)) - coordinate));  % Find the closest constellation point
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