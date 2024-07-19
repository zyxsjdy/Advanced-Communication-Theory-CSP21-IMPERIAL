% TimZ, MSc, 2021, Imperial College.
% 18/12/2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Takes polynomial weights and produces an M-Sequence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% primpoly (Px1 Integers) = Polynomial coefficients. For example, if the
% polynomial is D^5+D^3+D^1+1 then the coeffs vector will be [1;0;1;0;1;1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% MSeq (Wx1 Integers) = W bits of 1's and 0's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MSeq] = fMSeqGen(primpoly)

P = length(primpoly); % Degree of the primitive polynomial 
W = 2 ^ (P-1) - 1;    % Length of the generated m-sequence
MSeq = zeros(W,1);   % Initialize the m-sequence
shift_reg = ones(1,P-1);  % Initialize the shift register

% Generate the m-sequence
for q = 1:W
    new_shift = mod(shift_reg * primpoly(2:P),2);
    shift_reg(2:P-1) = shift_reg(1:P-2);
    shift_reg(1) = new_shift;
    MSeq(q) = shift_reg(P-1);
end

% Make sure the generated m-sequence is right
if shift_reg ~= ones(1,P-1)
    error('Error occurred')
end
end