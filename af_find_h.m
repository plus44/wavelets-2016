function [h] = af_find_h(tau, K)
%AF_FIND_H Finds the annihilating filter for a given tau and K value.
%   For a given tau[m] = sum_{k=0}^{K-1} a_k t_k^m, find the filter
%   coefficients, h[m], that when convolved with tau[m] produce zero output
%
%   The length of tau[m] must be at least 2K.

% Length of the input signal, tau[m]
N = length(tau); 

% Arrange the tau vector in a convolution/Toeplitz matrix
A = zeros(N-K, K);
for i = 1:(N-K),
    A(i,:) = fliplr(tau(i:(K-1)+i));
end

% The result of the matrix multiplication in the Yule-Walker system
b = -tau(K+1:N)';

% Solve Ax = b, for x
x = A\b;
h = [1; x];

end

