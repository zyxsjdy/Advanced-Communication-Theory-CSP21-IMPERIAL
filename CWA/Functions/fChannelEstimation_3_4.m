% TimZ, MSc, 2021, Imperial College.
% 10/1/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs channel estimation for the desired source using the received signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% photo_text_sym_rec (Fx1 Complex) = R channel symbol chips received
% GoldSeq (Wx1 Integers) = W bits of 1's and 0's representing the gold
% sequence of the desired source used in the modulation process
% paths (Mx1 Integers) = Number of paths for each source in the system
% For example, if 3 sources with 1, 3 and 2 paths respectively then
% paths = [1;3;2]
% array = Array locations in half unit wavelength. If no array then should
% be [0,0,0]
% P_noise (Integer) = The power of the noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% delay_est = Vector of estimates of the delays of each path of the
% desired signal
% DOA_est = Estimates of the azimuth and elevation of each path of the
% desired signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [delay_est, DOA_est] = fChannelEstimation_3_4(photo_text_sym_rec, GoldSeq, paths, array, P_noise)

[d_max,n_signal] = size(GoldSeq);  % Get the maximum relative delay and the number of sources
n_path = sum(paths);
DOA_est = zeros(n_path,2);  % Initialize the matrix for estimated DOA 
delay_est = zeros(n_path,1);  % Initialize the matrix for estimated delay

% Transform the received signal to discretised signal
% d_max is the code period
q = length(array);  % oversampling factor
discretised_signal = zeros(2*q*d_max, (length(photo_text_sym_rec)-d_max)/d_max);  % Initialize the discretised signal matrix
for i = 1:q
    sym_ant_i = photo_text_sym_rec(:, i);
    % Odd columns
    discretised_signal(2*(i-1)*d_max + 1:2*i*d_max, 1:2:(end-1)) = reshape(sym_ant_i(1:(length(sym_ant_i)-d_max)), 2*d_max, (length(sym_ant_i)-d_max) / (2*d_max));
    % Even columns
    discretised_signal(2*(i-1)*d_max + 1:2*i*d_max, 2:2:end) = reshape(sym_ant_i(d_max+1:length(sym_ant_i)), 2*d_max, (length(sym_ant_i)-d_max) / (2*d_max));
end

% Detection
COV = discretised_signal * discretised_signal' / length(discretised_signal);  % Calculate the covariance matrix
[eig_vec, eig_val] = eig(COV);
eig_vec_source = eig_vec(:, abs(diag(eig_val)) > 3 * P_noise);  % find the eigen vectors of sources
Pn = fpoc(eig_vec_source);  % Projection operator

% STAR
J = [zeros(1,2*d_max-1) 0; eye(2*d_max-1) zeros(2*d_max-1,1)];  % Shifting matrix J
c = [GoldSeq; zeros(d_max,n_signal)];  % Extend the PN code to fit the shifting matrix
azimuth = 0:359;  % Using different DOAs to find the ones with highest values
elevation = 0;
cost_fun = zeros(length(azimuth),d_max);  % Initialize the cost function
for n = 1:n_signal
    for a = azimuth
        SPV = spv(array, [a elevation]);  % Array manifold vector for the ath azimuth
        for i = 1:d_max
            h = kron(SPV, J^i * c(:,n));  % Basic STAR manifold
            cost_fun(a+1,i) = 1 ./ (h' * Pn * h);  % Calculate the cost function
        end
    end
    [cost_fun_max, I_DOA] = max(abs(cost_fun));  % Return the maximum values in each column (for each azi_i) and the row indices
    [~, I_delay] = maxk(cost_fun_max, paths(n_signal));  % Find the indices of paths(n_signal) maximum values
    if n == 1
        delay_est(1:paths(n),1) = sort(I_delay);  % Find the delay
        DOA_est(1:paths(n),1) = I_DOA(sort(I_delay)) - 1;  % Find the DOA
    else  % Dealing with multipath effect
        temp = sum(paths(1:n-1));
        delay_est((1:paths(n))+temp, 1) = sort(I_delay);
        DOA_est((1:paths(n))+temp, 1) = I_DOA(sort(I_delay)) - 1;
    end
end
end