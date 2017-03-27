close all;
clear all;

tau_s = open('provided/tau.mat');
tau = tau_s.tau;
clear tau_s;
% K = 2 for this exercise
K = 2;

% Find the annihilating filter
h = af_find_h(tau, K);
% Get the t_ks as the roots of the coefficients
t_k = af_find_t_k(h);
% Find the a_ks by solving the Vandermonde system
a_k = af_find_a_k(tau, t_k);

% Reconstruct tau
for i = 1:length(tau)
   tau_rec(i) = a_k'*(t_k.^(i-1)); 
end

figure;
tau_stem = stem(tau, 'b-', 'LineWidth', 4);
tau_labels = {num2str(tau(1)); num2str(tau(2)); num2str(tau(3)); ...
    num2str(tau(4))};
Xd = get(tau_stem, 'XData');
Yd = get(tau_stem, 'YData');
for n = 1:length(Xd)
   text(Xd(n) - 0.2 - 0.05*n, Yd(n), tau_labels{n}, ...
       'HorizontalAlignment','center', ...
       'Color', 'blue'); 
end

hold on;
tau_rec_stem = stem(tau_rec, 'r--', 'LineWidth', 2);
tau_rec_labels = {num2str(tau_rec(1)); num2str(tau_rec(2)); ...
    num2str(tau_rec(3)); num2str(tau_rec(4))};
Xd = get(tau_rec_stem, 'XData');
Yd = get(tau_rec_stem, 'YData');
for n = 1:length(Xd)
   text(Xd(n) - 0.2 - 0.05*n, Yd(n) - 500, tau_rec_labels{n}, ...
       'HorizontalAlignment', 'center', ...
       'Color', 'red'); 
end

lgd = legend([tau_stem tau_rec_stem], 'Original', 'Reconstructed');
lgd.FontSize = 12;
lgd.Location = 'northwest';
saveas(gcf, 'figures/Exercise3_1.png');


