function [L, Y, eig_values, expvar] = eofcore(X, n_eof)
%EOFCORE Leading EOF modes of a space-by-time matrix.
%
%   [L, Y, EIG_VALUES, EXPVAR] = EOFCORE(X, N_EOF)
%
%   Input
%     X       (n_locations x n_timesteps) data matrix. Rows with any NaN
%             should already have been removed. Time-mean is subtracted
%             here.
%     N_EOF   number of leading modes to return.
%
%   Output
%     L            (n_locations x n_eof) spatial modes (EOFs).
%                  Columns are orthonormal: L(:,a)'*L(:,b) = delta_ab.
%     Y            (n_eof x n_timesteps) principal components.
%                  Reconstruction: X_anom ≈ L * Y.
%     EIG_VALUES   (n_eof x 1) eigenvalues of the covariance X*X'/T,
%                  i.e. variance of each mode (T = n_timesteps).
%     EXPVAR       (n_eof x 1) percent of *total* (not truncated)
%                  variance explained by each mode.
%
%   Algorithm
%     Demean in time, then solve the smaller of the two covariance
%     eigenproblems (space-space vs time-time) following Björnsson &
%     Venegas (1997). If all modes are requested, a full eig is used
%     because eigs cannot return a complete spectrum.
%
%     X = L Y
%     L(:,k) is the k-th spatial EOF; Y(k,:) is the k-th PC.
%
%   Sign convention
%     The largest-magnitude entry of each EOF column is forced positive
%     so results are reproducible across MATLAB / Python / C++.
%
%   See also EOF, MEOF, EIG, EIGS, SVD.
%
%   Ref: Björnsson, H. and Venegas, S.A., 1997. A manual for EOF and SVD
%   analyses of climatic data. CCGCR Report, 97(1), pp.112-134.
%
%   Zelun Wu, Xiamen University & University of Delaware, 15 May 2020.
%   Total-variance expvar, full-eig fallback, sign convention: 2026.

if nargin < 2
    n_eof = min(size(X));
end

X = X - repmat(mean(X, 2), [1, size(X, 2)]);
[n_locations, n_timesteps] = size(X);
n_eof = min([n_eof, n_locations, n_timesteps]);
total_var = sum(sum(X .* X)) / n_timesteps;

if n_locations <= n_timesteps
    S = X * X';
    [L, eig_values] = eof_eigensolve(S, n_eof);
    eig_values = eig_values / n_timesteps;
    Y = (X' * L)';
else
    S = X' * X;
    [R, eig_values] = eof_eigensolve(S, n_eof);
    eig_values = eig_values / n_timesteps;
    % L = X * R * Lambda^{-1/2} maps time-eigenvectors back to space.
    scale = 1 ./ sqrt(max(eig_values * n_timesteps, realmin));
    L = X * (R * diag(scale));
    L = real(L);
    Y = L' * X;
end

% Orthonormal columns can flip sign; pin each mode.
for k = 1:size(L, 2)
    [~, idx] = max(abs(L(:, k)));
    if L(idx, k) < 0
        L(:, k) = -L(:, k);
        Y(k, :) = -Y(k, :);
    end
end

if total_var > 0
    expvar = eig_values / total_var * 100;
else
    expvar = zeros(size(eig_values));
end
end

function [V, eig_values] = eof_eigensolve(S, n_eof)
% Return the n_eof leading eigenpairs of symmetric S, descending.
n = size(S, 1);
% eigs cannot request a complete spectrum (k must be < n).
if n_eof >= n - 1 || n <= 2
    [V, D] = eig((S + S') / 2);
    [eig_values, order] = sort(diag(D), 'descend');
    V = V(:, order);
    V = V(:, 1:n_eof);
    eig_values = eig_values(1:n_eof);
else
    [V, D] = eigs((S + S') / 2, n_eof, 'lm');
    [eig_values, order] = sort(diag(D), 'descend');
    V = V(:, order);
end
eig_values = max(real(eig_values), 0);
end
