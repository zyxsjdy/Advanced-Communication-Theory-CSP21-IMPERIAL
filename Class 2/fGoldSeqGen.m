function GoldSeq_11 = fGoldSeqGen(a1,a2,k)

L = length(a2);
a2_delay = a2;
a2_delay(1:k) = a2(L-k+1:L);
a2_delay(k+1:L) = a2(1:L-k);

GoldSeq_01 = mod(a1+a2_delay,2);
GoldSeq_11 = 1-2*GoldSeq_01;
end

