% TimZ, MSc, 2021, Imperial College.
% 4/1/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Models the channel effects in the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% photo_sym (LxR Complex) = Signals being transmitted in the channel
% paths (Mx1 Integers) = Number of paths for each source in the system
% For example, if 3 sources with 1, 3 and 2 paths respectively then
% paths = [1;3;2]
% array = Array locations in half unit wavelength. If no array then should
% be [0,0,0]
% SNR = Signal to Noise Ratio in dB
% delay (Cx1 Integers) = Delay for each path in the system starting with
% source 1
% DOA = Direction of Arrival for each source in the system in the form
% [Azimuth, Elevation]
% fading_coeff (Cx1 Integers) = Fading Coefficient for each path in the
% system starting with source 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% photo_sym_rec (FxN Complex) = F channel symbol chips received from each
% antenna
% P_noise (Integer) = The power of the noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [photo_sym_rec, P_noise] = fChannel(photo_sym, paths, array, SNR, delay, DOA, fading_coeff)
                                    
[L,R] = size(photo_sym);  % L is the length of input signal, R is the number of signals
L = L + 15;          % Expand the input streams by 15, 
photo_sym(L,R) = 0;  % since the maximum relative delay < 15 (tau mod 15)

% Dealing with the multi-path effect for each source
photo_sym_path = zeros(L,sum(paths));  % Initialize the symbol streams in all paths
p = 1;  % A counter used to count the path
for i = 1:R
    for n = 1:paths(i)
        photo_sym_path(:,p) = fading_coeff(p) * circshift(photo_sym(:, i), delay(p));  % signal i pathing through path p
        p = p + 1;
    end

    % The desired signal is 1, hence use it to determine the power of the noise
    if i == 1
        photo_sym_desired = (spv(array,DOA(1:paths(1),:)) * photo_sym_path(:,1:paths(1)).').';
    end
end

% Noise
P_Signal = mean(sum(abs(photo_sym_desired).^2) / length(photo_sym_desired));  % Calculate the power of the desired signal
P_noise = P_Signal / SNR;  % Calculate the power of the noise
noise = sqrt(P_noise/2) * (randn(size(photo_sym_desired)) + 1i*randn(size(photo_sym_desired)));  % Generate the noise

% Received symbol streams
SPV = spv(array,DOA);  % Calculate the gain of elements on all directions, the function "spv" is given in AM1 experiment
photo_sym_rec = photo_sym_path * SPV.' + noise;
end