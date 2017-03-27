close all;
clear all;

% % % % % Set up the signal % % % % % 
K = 2;
x_k_original = [-1.25, 1.679]; 
t_k_original = [0.3, 0.8]; % Normalised to 1
x_t = zeros(1, 2048);
for i = 1:length(x_k_original),
   x_t(round(t_k_original*2048)) = x_k_original; 
end

t = [0:1:2048-1]./2048;

% % % % % Generate noisy moments % % % % % 
s_m = gen_noisy_moments(x_t, 'db8', 0.01);

% % % % % Noiseless method % % % % % 
[x_t_rec_af, t_k_af, x_k_af] = rec_af(s_m, K);

% % % % % TLS method % % % % % 
[x_t_rec_tls, t_k_tls, x_k_tls] = rec_tls(s_m, K);

% % % % % Cadzow's method % % % % % 
[x_t_rec_cadzow, t_k_cadzow, x_k_cadzow] = rec_cadzow(s_m, K, 5);

% % % % % Plot everything % % % % %
subplot(4, 1, 1); % Original signal
stem(t, x_t, 'k-');
title('Original signal, x(t)');
xlabel('t/s');
subplot(4, 1, 2); % Reconstructed signal annihilating filter
stem(t, x_t_rec_af, 'b--');
title('Reconstructed signal by annihilating filter, x_{af}(t)');
xlabel('t/s');
subplot(4, 1, 3); % Reconstructed signal TLS
stem(t, x_t_rec_tls, 'b--');
title('Reconstructed signal by TLS, x_{tls}(t)');
xlabel('t/s');
subplot(4, 1, 4); % Reconstructed signal Cadzow
stem(t, x_t_rec_cadzow, 'b--');
title('Reconstructed signal by Cadzow, x_{cad}(t)');
xlabel('t/s');
saveas(gcf, 'figures/Exercise7_1.png');

