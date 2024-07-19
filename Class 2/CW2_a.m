%% Initialization
clc;
clear;
close all;

%% Set basic parameters

c = [1 0 1 1];
alpha_01 = fMSeqGen(c);
alpha_11 = 1-2*alpha_01;

L = length(c);
Nc = 2^(L-1)-1;

k = -10:1:10;
R = zeros(1,length(k));
for q = 1:length(k)
    for p = 1:Nc
        if mod(p+k(q),Nc) == 0
            nk = Nc;
        else
            nk = mod(p+k(q),Nc);
        end
        R(q) = R(q) + alpha_11(p) * alpha_11(nk);
    end
end
plot(k,R)






























