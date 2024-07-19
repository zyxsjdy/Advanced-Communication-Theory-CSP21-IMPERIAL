function N = MDL(L, signal_sample_cov)

[~, Eig_value] = eig(signal_sample_cov);
Eig_value = sort(abs(diag(Eig_value)));
L_r = length(Eig_value);

MDL_k = (-L) * (log(flip(cumprod(Eig_value))) + (L_r:-1:1)' .* (log((L_r:-1:1)') - log(flip(cumsum(Eig_value))))) ...
    + 1/2 * log(L_r) * (0:L_r-1)' .* (2 * L_r:-1:L_r + 1)'; 

[~, Row_min] = min(MDL_k);
N = Row_min - 1;
end