function [s_m_noisy] = gen_noisy_moments(x_t, daub, var)
%GEN_NOISY_MOMENTS Generates moments + noise for an input stream of Diracs
%   Given an input stream of Diracs, 'x_t', of size 2048, we sample it 
%   using the specified Daubechies scaling function, 'daub', at a sampling
%   period of T = 64. We compute the weights of the sampling kernel and use
%   that to generate the moments of the sampled signal. Zero-mean Gaussian
%   noise of variance, 'var', is added to the moments.

% % % % % Find the polynomial reproducing degree % % % % % 

N = 1; % The number of vanishing moments

% If there are p vanishing moments, polynomials of degree p-1 can be
% reproduced.
switch daub
    case 'db1'
        N = 1;
    case 'db2'
        N = 2;
    case 'db3'
        N = 3;
    case 'db4'
        N = 4;
    case 'db5';
        N = 5;
    case 'db6'
        N = 6;
    case 'db7'
        N = 7;
    case 'db8'
        N = 8;
    case 'db9'
        N = 9;
    case 'db10'
        N = 10;
    case 'db11'
        N = 11;
    case 'db12'
        N = 12;
    otherwise
        return
end

phi = zeros(1, 2048);
[phi_T, ~, ~] = wavefun(daub, 6);
phi_n = zeros(2048, 32); % Shifted versions of the wavelet
c_mn = zeros(N, 32); % c_{m,n} weights
phi(1:length(phi_T)) = phi_T; 

t = [0:1:2048-1]./2048;

% Find the c_mn weights for reproducing polynomials with this wavelet
for m = 0:N-1,
    for n = 0:32-1,
        % Calculate and store the shifted wavelets
        phi_n(:, n+1) = circshift(phi', 64*n);
        c_mn(m+1, n+1) = (1/64) * t.^m * phi_n(:, n+1);
    end
end

% % % % % Sample x(t) % % % % % 
y_n = zeros(1, 32);
for n = 0:32-1,
    y_n(n+1) = x_t*phi_n(:,n+1);
end

% % % % % Find moments % % % % % 
s_m = zeros(1, N);
for m = 0:N-1,
   s_m(m+1) = c_mn(m+1, :) * y_n'; 
end

eps = sqrt(var) * randn(1, N);
s_m_noisy = s_m + eps;

end

