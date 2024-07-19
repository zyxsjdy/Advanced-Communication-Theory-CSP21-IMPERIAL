function N = AIC(L, signal_sample_cov)

[~, Eig_value] = eig(signal_sample_cov);
Eig_value = sort(abs(diag(Eig_value)));
L_r = length(Eig_value);

AIC_k = (-2) * L * (log(flip(cumprod(Eig_value))) + (L_r:-1:1)' .* (log((L_r:-1:1)') - log(flip(cumsum(Eig_value))))) ...
    + 2 * (0:L_r-1)' .* (2 * L_r:-1:L_r + 1)';

[~, Row_min] = min(AIC_k);
N = Row_min - 1;
end