% TimZ, MSc, 2021, Imperial College.
% 3/1/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform DS-QPSK Modulation on a vector of bits using a gold sequence
% with channel symbols set by a phase phi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% photo_bit (Px1 Integers) = P demodulated bits of 1's and 0's
% GoldSeq (Wx1 Integers) = W bits of 1's and 0's representing the gold
% sequence of the desired signal to be used in the demodulation process
% myphi (Integer) = Angle index in degrees of the QPSK constellation points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% photo_sym (Rx1 Complex) = Signals being transmitted in the channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [photo_sym] = fDSQPSKModulator(photo_bit, GoldSeq, myphi)

% length(bitsIn) is a even number (160 * 112 * 3 * 8), hence can be reshaped to a (2*x) matrix
% Each column is used for one QPSK symbol
photo_bit2 = reshape(photo_bit,2,[]);  

% QPSK
QPSK = zeros(length(photo_bit2),1);  % Initialize the QPSK result matrix
for i = 1:length(photo_bit2)
    switch string(photo_bit2(:,i))
        case ["0";"0"]
            QPSK(i) = sqrt(2) * (cos(myphi) + 1i*sin(myphi));
        case ["0";"1"]
            QPSK(i) = sqrt(2) * (cos(myphi + pi/2) + 1i*sin(myphi + pi/2));
        case ["1";"1"]
            QPSK(i) = sqrt(2) * (cos(myphi - pi) + 1i*sin(myphi - pi));
        case ["1";"0"]
            QPSK(i) = sqrt(2) * (cos(myphi - pi/2) + 1i*sin(myphi - pi/2));
        otherwise
            error('Error occurred')
    end
end

% Encode the signal and reshape it to a (R x 1) matrix
photo_sym = reshape(GoldSeq * QPSK.',[],1);  
end