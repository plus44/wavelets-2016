function [h] = cadzow_find_h(tau, K, iter)
%CADZOW_FIND_H Use Cadzow's algo to denoise tau, iter times and find h.
%   1. Set up the Toeplitz convolution matrix of tau, S
%   2. Perform SVD on S = U * L * V'
%   3. Keep the K highest entries on the diagonal of L, set the others to 0
%   4. Let S_new = U * L_new * V'.
%   5. S_new is not Toeplitz. To make it so, average along the diagonals.
%   6. Iterate steps 2-5 'iter' times on S_new.
% 
%   Once the iterations are complete, h will be the last column of V.

% % % % % 1. Set up the Toeplitz matrix S % % % % %
len_tau = length(tau);
len_h = K + 1;
S = zeros(len_tau - len_h + 1, len_h);
for i = 1:length(S(:,1)),
    S(i,:) = fliplr(tau(i:i + len_h - 1));
end

% c = [tau zeros(1, len_h - 1)];
% r = [tau(1) zeros(1, len_h - 1)];
% S = toeplitz(c, r); % Construct a Toeplitz matrix, S, s.t. SH = 0

for k = 1:iter,
    
    % % % % % 2. Perform the SVD of S % % % % % 
    [U, L, V] = svd(S);

    % % % % % 3. Keep the K highest entries along the diagonal % % % % % 
    L_sorted = sort(diag(L), 'descend');
    L_sorted = L_sorted(1:K); % Take note of the K highest entries
    L_new = L;

    % Assume that L is always a tall matrix (true as long as N >= 2K)
    for i = 1:length(L_new(1,:)),
        % If this value along the diagonal is part of the K largest values
       if ismember(L_new(i, i), L_sorted),
           % Remove it from the list of K largest elements
           L_sorted = L_sorted(L_sorted ~= L_new(i, i));
       else,
           % Otherwise clear this element.
           L_new(i, i) = 0;
       end
    end

    % % % % % 4. Reconstruct S_new from L_new % % % % % 
    S_new = U * L_new * V';

    % % % % % 5. Convert S_new to Toeplitz % % % % % 
    % Average along the longest diagonals
    n_diag = length(S_new(:,1)) - length(S_new(1,:)) + 1;
    for i = 0:n_diag - 1,
        % Shift it in place for this iteration
        S_new = circshift(S_new, -i, 1);

        % Average this diagonal
        avg = sum(diag(S_new))/len_h;
        for j = 1:len_h,
           S_new(j, j) = avg; 
        end

        % Shift it back to original
        S_new = circshift(S_new, i, 1);
    end

    % Average in the top and bottom edge cases
    for i = len_h:-1:2,
        height_S = length(S_new(:,1));

        % Get submatrices from the top right and bottom left corners
        mat_top = S_new(1:len_h-i+1, i:len_h);
        mat_bot = S_new(height_S-len_h+i:height_S, 1:len_h-i+1);

        % Find the averages along the diagonals of the submatrices
        avg_top = trace(mat_top)/length(mat_top(1,:));
        avg_bot = trace(mat_bot)/length(mat_bot(1,:));

        % The submatrices will always be square by the way we obtain them.
        for j = 1:length(mat_top(1,:)),
            % Set the diagonal elements to the average
            mat_top(j, j) = avg_top;
            mat_bot(j, j) = avg_bot;
        end

        % Assign the new, averaged matrices to S_new.
        S_new(1:len_h-i+1, i:len_h) = mat_top;
        S_new(height_S-len_h+i:height_S, 1:len_h-i+1) = mat_bot;
    end

    % % % % % 6. Assign the output to the input and iterate % % % % % 
    S = S_new;
end

% % % % % Perform one final SVD after the final iteration % % % % % 
[~, ~, V] = svd(S);

% % % % % Extract h as the last col of V, as in TLS % % % % % 
h = V(:, length(V(1,:)));

end

