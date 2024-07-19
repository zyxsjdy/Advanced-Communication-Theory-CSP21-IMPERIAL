% TimZ, MSc, 2021, Imperial College.
% 11/1/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spatiotemporal Beamformer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% text_sym_rec (Fx1 Complex) = R channel symbol chips received
% GoldSeq (Wx1 Integers) = W bits of 1's and 0's representing the gold
% sequence of the desired source used in the modulation process
% n_path (Integer) = Number of path
% For example, if 3 sources with 1, 3 and 2 paths respectively then
% paths = [1;3;2]
% array = Array locations in half unit wavelength. If no array then should
% be [0,0,0]
% delay_est = Vector of estimates of the delays of each path of the
% desired signal
% DOA_est = Estimates of the azimuth and elevation of each path of the
% desired signal
% fading_coeff (Cx1 Integers) = Fading Coefficient for each path in the
% system starting with source 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% weighted_discretised_signal (Fx1 Complex) = F weighted discretised symbols
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [weighted_discretised_signal] = spatiotemporal(text_sym_rec, GoldSeq, n_path, array, delay_est, DOA_est, fading_coeff)

[d_max, n_signal] = size(GoldSeq);  % Get the maximum relative delay and the number of sources

% Transform the received signal to discretised signal
n_ant = length(array);  % oversampling factor
discretised_signal = zeros(2*n_ant*d_max, (length(text_sym_rec)-d_max)/d_max);  % Initialize the discretised signal matrix, M_Length is the code period
for n = 1:n_ant
    sym_ant_i = text_sym_rec(:,n);
    % Odd columns
    discretised_signal(2*(n-1)*d_max + 1:2*n*d_max, 1:2:end-1) = reshape(sym_ant_i(1 : length(sym_ant_i)-d_max), 2*d_max, (length(sym_ant_i)-d_max) / (2*d_max));
    % Even columns
    discretised_signal(2*(n-1)*d_max + 1:2*n*d_max, 2:2:end) = reshape(sym_ant_i(d_max+1 : length(sym_ant_i)), 2*d_max, (length(sym_ant_i)-d_max) / (2*d_max));
end

% STAR
J = [zeros(1,2*d_max-1) 0; eye(2*d_max-1) zeros(2*d_max-1,1)];  % Shifting matrix J
c = [GoldSeq;zeros(d_max,n_signal)];  % Extend the PN code to fit the shifting matrix
p = 1;
for n = 1: n_signal
    h = zeros(2 * n_ant * d_max, n_path(n));  % Initialize the basic STAR manifold
    for i = 1: n_path(n)
        SPV = spv(array, DOA_est(p,:));  % Array manifold vector
        h(:, i) = kron(SPV, J^delay_est(p) * c(:,n));  % Basic STAR manifold
        p = p + 1;
    end
end
weight = h * fading_coeff;  % Weights of the spatiotemporal beamformer
weighted_discretised_signal = (weight' * discretised_signal).';
end