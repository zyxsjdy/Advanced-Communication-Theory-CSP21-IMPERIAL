% TimZ, MSc, 2021, Imperial College.
% 10/1/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DOA estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% x_DOA (NxL Complex) = received signal
% array = Array locations in half unit wavelength
% Pn (Integer) = power of noise (dB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% DOA_est = Estimates of the azimuth of each path of the received signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [DOA_est] = doa(x_DOA, array, Pn)

COV = x_DOA * x_DOA' / size(x_DOA,2);  % Calculate the covariance matrix
[eig_vec, eig_val] = eig(COV);
eig_vec_source = eig_vec(:, abs(diag(eig_val)) > 3 * Pn);  % find the eigen vectors of sources
n_source = size(eig_vec_source,2);  % The number of sources
P = fpoc(eig_vec_source);  % Projection operator

azimuth = 0:359;  % Using different DOAs to find the ones with highest values
elevation = 0;
DOA_est = zeros(n_source,2);  % Initialize the matrix for estimated DOA 
cost_fun = zeros(length(azimuth),1);  % Initialize the cost function
for a = azimuth
    SPV = spv(array, [a elevation]);  % Array manifold vector for the ath azimuth
    cost_fun(a+1) = 1 ./ (SPV' * P * SPV);  % Calculate the cost function
end
[~, I_DOA] = maxk(abs(cost_fun), n_source);  % Find the indices of n_source maximum values
DOA_est(1:n_source,1) = I_DOA - 1;
DOA_est = DOA_est(:,1);  % Only take the azimuth, since elevation is 0
end