function [a_k] = af_find_a_k(tau, t_k)
%AF_FIND_A_K Finds the a_k coefficients of a signal, given tau and t_k
%   Solves the Vandermonde system that finds the a_k values of a signal
%   through the use of an annihilating filter

% K is simply the length of t_k
K = length(t_k);

% Arrange the t_ks into a Vandermonde matrix
v = vander(t_k);
V = flipud(v');
% Find the a_ks by taking the inverse of the Vandermonde matrix
% NOTE: This can be done with the analytical formula for the Vandermonde 
% inverse in order to save computations and reduce complexity.
V_inv = V^-1;
% Take K samples of tau
p = tau(1:K)';

% Solve the Vandermonde system
a_k = V_inv*p;

end

