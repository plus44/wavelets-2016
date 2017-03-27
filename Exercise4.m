close all;
clear all;

% % % % % Set up the signal % % % % % 
x_k_original = [-1.100, 0.369]; 
t_k_original = [0.22, 0.88]; % Normalised to 1
x_t = zeros(1, 2048);
for i = 1:length(x_k_original),
   x_t(round(t_k_original*2048)) = x_k_original; 
end

% % % % % Exercise 1 % % % % % 
phi = zeros(1, 2048); % Original wavelet at t = 0
phi_n = zeros(2048, 32); % Shifted versions of the wavelet
c_mn = zeros(4, 32); % c_{m,n} weights
[phi_T, psi_T, xval] = wavefun('db4', 6); % 2^6 = 64
phi(1:length(phi_T)) = phi_T;

t = [0:1:2048-1]./2048;

% Find the c_mn coefficients for reproducing polynomials with this wavelet
for m = 0:3,
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

% % % % % Find s_m % % % % % 
s_m = zeros(1, 4);
for m = 0:3,
   s_m(m+1) = c_mn(m+1, :) * y_n'; 
end

% % % % % Exercise 3 % % % % % 
K = 2;

% Find the annihilating filter
h = af_find_h(s_m, K);
% Get the t_ks as the roots of the coefficients
t_k = af_find_t_k(h);
% Find the a_ks by solving the Vandermonde system
x_k = af_find_a_k(s_m, t_k);

% % % % % Reconstruct the signal % % % % % 
x_t_rec = zeros(1, 2048);
x_t_rec(round(t_k(1)*2048)) = x_k(1);
x_t_rec(round(t_k(2)*2048)) = x_k(2);

% % % % % Plot everything % % % % %
subplot(2, 2, 1); % Original signal
plot(t, x_t, 'k-');
title('Original signal, x(t)');
xlabel('t/s');
subplot(2, 2, 2); % Sampling kernel
plot(xval, phi_T);
title('Sampling kernel, \phi(t)');
subplot(2, 2, 3); % Sampled signal
stem(linspace(0, 1, 32), y_n, 'b-');
title('Sampled signal, y_n');
xlabel('t/s');
subplot(2, 2, 4); % Reconstructed signal
plot(t, x_t_rec, 'b--');
title('Reconstructed signal, x_{rec}(t)');
xlabel('t/s');
saveas(gcf, 'figures/Exercise4_2.png');

