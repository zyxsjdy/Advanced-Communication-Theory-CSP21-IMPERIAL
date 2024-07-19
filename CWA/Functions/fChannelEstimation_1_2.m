% TimZ, MSc, 2021, Imperial College.
% 29/12/2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs channel estimation for the desired source using the received signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% photo_sym_rec (Fx1 Complex) = R channel symbol chips received
% GoldSeq (Wx1 Integers) = W bits of 1's and 0's representing the gold
% sequence of the desired source used in the modulation process
% paths (Mx1 Integers) = Number of paths for each source in the system
% For example, if 3 sources with 1, 3 and 2 paths respectively then
% paths = [1;3;2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% delay_est = Vector of estimates of the delays of each path of the
% desired signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [delay_est] = fChannelEstimation_1_2(photo_sym_rec, GoldSeq, paths)

[d_max,~] = size(GoldSeq);  % Get the maximum relative delay

% Calculate the correlation
% N is the possible delay, use a large number to increase the accuracy, e.g. N = 10000, possible delay = 1:N
N = 10000;
corr = zeros(N,1);
for d = 1:N
    corr(d) = abs(photo_sym_rec(d:(d+d_max-1)).' * GoldSeq(:,1));
end

[~, I_delay] = sort(corr,'descend');  % Order the values and get the indexs (possible delays)
I_delay_u = unique(mod(I_delay,d_max),'stable');  % Since there are many delays (N = 10000), consider unique relative delays
delay_est = mod(sort(I_delay_u(1:paths(1)))-1,d_max);  % Consider the first (paths(1)) delays as the estimated delays
end