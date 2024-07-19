% TimZ, MSc, 2021, Imperial College.
% 18/12/2021

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Takes two M-Sequences of the same length and produces a gold sequence by
% adding a delay and performing modulo 2 addition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% MSeq_1 (Wx1 Integer) = First M-Sequence
% MSeq_2 (Wx1 Integer) = Second M-Sequence
% d (Integer) = Number of chips to shift second M-Sequence to the right
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% GoldSeq (Wx1 Integer) = W bits of 1's and -1's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [GoldSeq] = fGoldSeq(MSeq_1, MSeq_2, d)

MSeq_2_delay = circshift(MSeq_2, d);  % Apply the circular shift to the second m-sequence
goldSeq_01 = mod(MSeq_1+MSeq_2_delay,2);  % Generate the 0/1 version gold sequence
GoldSeq = 1-2*goldSeq_01;  % Get the 1/-1 version gold sequence
end