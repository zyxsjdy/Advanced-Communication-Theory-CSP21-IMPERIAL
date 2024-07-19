%% Initialization
clc;
clear;
close all;

%% Set basic parameters

c1 = [1 0 0 1 0 1];
alpha1_01 = fMSeqGen(c1);
alpha1_11 = 1-2*alpha1_01;
c2 = [1 1 1 1 0 1];
alpha2_01 = fMSeqGen(c2);
alpha2_11 = 1-2*alpha2_01;

L = length(c1);
Nc = 2^(L-1)-1;

b = zeros(Nc,Nc+2);
b(:,Nc+1) = alpha1_11';
b(:,Nc+2) = alpha2_11';
for q = 1:Nc
    delay = q;
    b(:,q) = (fGoldSeqGen(alpha1_01,alpha2_01,delay))';
end

check = sum(b,1);
flag = 0;
k1 = 0;
k2 = 0;
for q = 1:Nc+2
    if check(q) == -1
        if flag ==0
            k1 = q;
            flag = 1;
        elseif flag == 1
            k2 = q;
            break
        end
    end
end

k = -30:1:30;
Rbb = zeros(1,length(k));
for q = 1:length(k)
    for n = 1:Nc
        if mod(n+k(q),Nc) == 0
            nk = Nc;
        else
            nk = mod(n+k(q),Nc);
        end
        Rbb(q) = Rbb(q) + b(n,k1) * b(nk,k1);
    end
end
figure(1)
plot(k,Rbb)

Rbd = zeros(1,length(k));
for q = 1:length(k)
    for n = 1:Nc
        if mod(n+k(q),Nc) == 0
            nk = Nc;
        else
            nk = mod(n+k(q),Nc);
        end
        Rbd(q) = Rbd(q) + b(n,k1) * b(nk,k2);
    end
end
figure(2)
plot(k,Rbd)



% plot all in b
% for i = 1:Nc+2
%     k = -100:1:100;
%     R = zeros(1,length(k));
%     for q = 1:length(k)
%         for p = 1:Nc
%             if mod(p+k(q),Nc) == 0
%                 nk = Nc;
%             else
%                 nk = mod(p+k(q),Nc);
%             end
%             R(q) = R(q) + b(p,i) * b(nk,i);
%         end
%     end
%     figure(i)
%     plot(k,R)
% end

