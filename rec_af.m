function [x_t_rec, t_k, x_k] = rec_af(s_m, K)
%REC_AF Reconstruct a Dirac stream of K from its moments, by AF. 
%   Simple annihilating filter approach. The original signal is assumed to 
%   have length 2048.

% Find the annihilating filter
h = af_find_h(s_m, K);
% Get the t_ks as the roots of the coefficients
t_k = af_find_t_k(h);
% Find the a_ks by solving the Vandermonde system
x_k = af_find_a_k(s_m, t_k);

% % % % % Reconstruct the signal % % % % % 
x_t_rec = zeros(1, 2048);

% We need to protect against passing noisy moments, which may lead to
% either negative or complex t_k. 
if (t_k(1) > 0 && imag(t_k(1)) == 0),
    x_t_rec(round(t_k(1)*2048)) = x_k(1);
end
if (t_k(2) > 0 && imag(t_k(2)) == 0),
    x_t_rec(round(t_k(2)*2048)) = x_k(2);
end
x_t_rec = x_t_rec(1:2048);

end

