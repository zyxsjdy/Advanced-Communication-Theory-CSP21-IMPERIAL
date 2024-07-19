function M=fpoc(A)
%complement projection operator of A
M=A/(A'*A)*A';
M=eye(size(M))-M;

