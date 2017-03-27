close all;
clear all;

phi = zeros(1, 2048); % Original wavelet at t = 0
phi_n = zeros(2048, 32); % Shifted versions of the wavelet
c_mn = zeros(4, 32); % c_{m,n} weights
[phi_T, psi_T, xval] = wavefun('db4', 6);
phi(1:length(phi_T)) = phi_T;
t = (0:1:2048-1)./2048;

% Find the c_mn coefficients for reproducing polynomials with this wavelet
for m = 0:3
    for n = 0:32-1
        % Calculate and store the shifted wavelets
        phi_n(:, n+1) = circshift(phi', 64*n);
        c_mn(m+1, n+1) = (1/64) * t.^m * phi_n(:, n+1);
    end
end

figure;
x = zeros(1, 2048);
x(1025) = 1.5;
for m = 0:3,
    subplot(2, 2, m+1);
    stem(0.5, x(1025), 'k');
    ylim([-0.5, 2]);
    xlabel('t/s');
    hold on;
    
    % Initialise the sum of phis weighted by c_mn
    phi_c_mn_tot = zeros(2048, 1);
    for n = 0:32-1,
        phi_c_mn = phi_n(:, n+1) .* c_mn(m+1, n+1);
        phi_c_mn_tot = phi_c_mn_tot + phi_c_mn;
        plot(t, phi_c_mn, 'r--');
    end
    
    plot(t, phi_c_mn_tot, 'b-');
end
saveas(gcf, 'figures/Exercise1_1.png');

figure;
for m = 0:3,
   subplot(2, 2, m+1);
   xlabel('t/s');
   ylim([-0.5 0.5]);
   hold on;
   
   % Initialise the sum of phis weighted by c_mn
   phi_c_mn_tot = zeros(2048, 1);
   for n = 0:32-1,
       phi_c_mn = phi_n(:, n+1) .* c_mn(m+1, n+1);
       phi_c_mn_tot = phi_c_mn_tot + phi_c_mn;
   end
   
   polynomial = t.^m;
   error = phi_c_mn_tot - polynomial';
   plot(t, error, 'r-');
end
saveas(gcf, 'figures/Exercise1_2.png');
