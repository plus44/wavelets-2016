close all;
clear all;

% % % % % Set up the signal % % % % % 
x_k_original = [0.192, 1.682]; 
t_k_original = [0.329, 0.851]; % Normalised to 1
x_t = zeros(1, 2048);
for i = 1:length(x_k_original),
   x_t(round(t_k_original*2048)) = x_k_original; 
end

% % % % % Generate noisy moments % % % % % 
s_m = gen_noisy_moments(x_t, 'db4', 0.1);