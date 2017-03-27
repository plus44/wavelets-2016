function [x] = gen_b_spline(M, n)
%GEN_B_SPLINE Generate a B-spline of order m, N samples, centre N/2.
%   Repeated convolution

N = n/(2^M);
p = nextpow2(N);
N = 2^p;

x = zeros(1, 7*N);
x(3*N+1:4*N) = ones(1, N);

if M == 0,
    x = x(3*N+1: 4*N);
    return
end

for m = 1:M,
   x = conv(x, x)./(N);
  
%    j = 1;
%    for i = 1:length(x),
%        if x(i) ~= 0,
%           new_x(j) = x(i); 
%           j = j + 1;
%        end
%    end
%    
%    x = new_x;
end

idx_lo = ceil((length(x) - n)/2) + 1;
idx_hi = ceil((length(x) + n)/2);
x = x(idx_lo:idx_hi);

end

