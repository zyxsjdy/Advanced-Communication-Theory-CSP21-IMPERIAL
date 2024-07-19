%% Initialization
clc;
clear;
close all;
load 'Ex04_Array_Signal_Snapshots'
%% Main
array = [-2 0 0; -1 sqrt(3) 0; 1 sqrt(3) 0; 2 0 0; 1 -sqrt(3) 0; -1 -sqrt(3) 0];
L_a = length(array);
L = 1000;
signal = x;
signal_cov= signal * signal' / length(signal(1, :));
[Eig_vector, Eig_value] = eig(signal_cov);
snapshots = randn(1,L) + 1i * randn(1,L);
for k = 1:L
    signal_sample(:,(k-1) * L_a + 1:k * L_a) = Eig_vector * sqrt(Eig_value) * snapshots(k);
end
signal_sample_cov = (signal_sample * signal_sample') / L;

sourceNo_AIC = AIC(L, signal_sample_cov)
sourceNo_MDL = MDL(L, signal_sample_cov)













