function [x_t_rec, t_k, x_k] = rec_tls(s_m, K)
%REC_TLS Reconstruct a Dirac stream of K from its moments, by TLS. 
%   Total Least Squares (TLS) solution to the reconstruction procedure.

% Find the annihilating filter by TLS
h = tls_find_h(s_m, K);
% Get the t_ks as the roots of the coefficients
t_k = af_find_t_k(h);
% Find the a_ks by solving the Vandermonde system
x_k = af_find_a_k(s_m, t_k);

% % % % % Reconstruct the signal % % % % % 
x_t_rec = zeros(1, 2048);

% We need to protect against negative or complex t_k. 
if (t_k(1) > 0 && imag(t_k(1)) == 0),
    x_t_rec(round(t_k(1)*2048)) = x_k(1);
end
if (t_k(2) > 0 && imag(t_k(2)) == 0),
    x_t_rec(round(t_k(2)*2048)) = x_k(2);
end
x_t_rec = x_t_rec(1:2048);

end

