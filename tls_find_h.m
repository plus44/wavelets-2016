function [h] = tls_find_h(tau, K)
%TLS_FIND_H Finds the annihilating filter given a signal tau, by TLS.
%   For a given tau[m] = sum_{k=0}^{K-1} a_k t_k^m, find the filter
%   coefficients, h[m], that when convolved with tau[m] produce zero 
%   output. The method used to find this filter is the Total Least Squares 
%   (TLS) method.
%
%   The length of tau[m] must be at least 2K. h will be the last column of
%   V, at the end of the algorithm.

% % % % % Set up the Toeplitz matrix S % % % % %
len_tau = length(tau);
len_h = K + 1;
S = zeros(len_tau - len_h + 1, len_h);
for i = 1:length(S(:,1)),
    S(i,:) = fliplr(tau(i:i + len_h - 1));
end

% % % % % Find the SVD of S % % % % % 
[~, ~, V] = svd(S);

% % % % % Extract h as the last col of V % % % % % 
h = V(:, length(V(1,:)));

end

