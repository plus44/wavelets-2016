close all;
clear all;

y_n_s = open('provided/samples.mat');
y_n = y_n_s.y_sampled;
clear y_n_s;

% % % % % Find the c_mn coefficients % % % % % 
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
subplot(3, 1, 1); % Sampling kernel
plot(xval, phi_T);
title('Sampling kernel (dB4, T=64)');
xlabel('t/s');
subplot(3, 1, 2); % Sampled signal
plot(linspace(0, 1, 32), y_n, 'b-');
title('Sampled signal, y_n');
xlabel('t/s');
subplot(3, 1, 3); % Reconstructed signal
plot(t, x_t_rec, 'b--');
title('Reconstructed signal, x_{rec}(t)');
xlabel('t/s');
saveas(gcf, 'figures/Exercise5_1.png');

