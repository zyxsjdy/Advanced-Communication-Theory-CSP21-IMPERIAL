% TimZ, MSc, 2021, Imperial College.
% 31/12/2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate the arrival time of signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% x_Time (NxL Complex) = received signal
% SNR (Integer) = Signal to Noise Ratio (dB)
% Pn (Integer) = power of noise (dB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% TOA (Integer) = estimated time of arrival
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [TOA] = toa(x_Time, SNR, Pn)

x_power = abs(x_Time);  % Get the power of received signal
L = length(x_power);

% Quantization, find the maximum between the present and previous values
x_quan = x_power(1);
for i = 2:L
    x_quan(i) = max(x_power(i), x_power(i-1));
end

d_x_quan = [0 diff(x_quan)];  % Apply the derivative

% Determine the time of arrival
threshold = (SNR+Pn)/2;  % The level between the power of signal and noise
%threshold = Pn*2;  % The level of twice noise power
for i = 1:L
    if d_x_quan(i) > threshold
        TOA = i;
        break
    elseif i == L
        error('Fail to determine the time of arrival, please decrease the threshold.')
    end
end
end