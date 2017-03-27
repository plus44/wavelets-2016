function [t_k] = af_find_t_k(h)
%AF_FIND_T_K Find the t_ks of a provided annihilating filter.
%   Find the roots of the h[m] filter coefficients - these are the t_ks of
%   the original tau signal.
t_k = roots(h);

end

