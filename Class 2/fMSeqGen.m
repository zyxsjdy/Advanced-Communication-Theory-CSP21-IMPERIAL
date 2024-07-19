function alpha_01 = fMSeqGen(c)

L = length(c);
shift_register = ones(1,L-1);
m_seq = zeros(1,2^(L-1)-1);

for q = 1:2^(L-1)-1
    new_shift = mod(c(2:L) * shift_register',2);
    shift_register(2:L-1) = shift_register(1:L-2);
    shift_register(1) = new_shift;
    m_seq(q) = shift_register(L-1);
end

if shift_register ~= ones(1,L-1)
    error('Error occurred')
end

alpha_01 = m_seq;
end

