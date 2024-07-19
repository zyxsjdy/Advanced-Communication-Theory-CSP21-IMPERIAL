% TimZ, MSc, 2021, Imperial College.
% 6/1/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Super-Resolution Beamformer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% photo_sym_rec (Fx1 Complex) = R channel symbol chips received
% array = Array locations in half unit wavelength. If no array then should
% be [0,0,0]
% DOA_est = Estimates of the azimuth and elevation of each path of the
% desired signal
% DOA_Target = The DOA of the desired source
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% weighted_photo_sym (Fx1 Complex) = F weighted photo symbols
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [weighted_photo_sym] = super_resolution(photo_sym_rec, array, DOA_est, DOA_target)

SPV_interf = spv(array, setdiff(DOA_est, DOA_target, 'rows'));   % manifold vector of the interferences
SPV_target = spv(array, DOA_target);  % manifold vector of the desired signal

weight = fpoc(SPV_interf) * SPV_target;  % Weight of the super-resolution beamformer
weighted_photo_sym = photo_sym_rec * conj(weight);  % Weight the received signals form different antennas
end