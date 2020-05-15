function [L, Y, eig_values, expvar] = eofcore(X, n_eof)
% [L, Y, eig_values] = EOFCORE(X, n_eof)
% Version 1.0
% Calculate the first n_eof modes of data field X(N_locations,N_timesteps)
% Output:
%     L: eigenmatrix. L(:,a)*L(:,b)' = 0 when a ~= b
%     Y: priciple component, coordinates values of the new coordinate
%     eig_values: eigenvalues 
%     expvar: percent of variance explained by each mode
% Agorithm:
%     X = LY
%     where X is the original signals measured in N_locations with N_timesteps, L is the bases of new coordinates, Y is the magnitude of the vector in the new corrdinates
%     In EOF, L(:,n_th) is the n_th mode eof_maps (spatial component), 
%     Y(n_th,:) is the n_th principle components time series.
% 
% See also EIG, EIGS, SVD, SVDS
% 
% Ref: Björnsson, H. and Venegas, S.A., 1997. A manual for EOF and SVD analyses of climatic data. CCGCR Report, 97(1), pp.112-134.
% 
% Author:
%     Zelun Wu,
%     Ph.D. student of Physical Oceanography,
%     Xiamen University & University of Delaware
%     zelunwu@stu.xmu.edu.cn, zelunwu@udel.edu
%     15th May, 2020

%% 
X = X - mean(X,2); % Remove mean in time dimension
[N_locations, N_timesteps] = size(X);

if N_locations < N_timesteps
    S = X * X'; %covariance matrix (N_locations*N_locations)
    [L,V] = eigs(S,n_eof);
    [eig_values, I] = sort(diag(V),'descend');
    L = L(:,I);
    eig_values = eig_values / N_timesteps;
    Y = (X'*L)';
else % time-space conversion
    S = X' * X; %covariance matrix (N_timesteps*N_timesteps)
    [L,V] = eigs(S,n_eof);
    [eig_values,I] = sort(diag(V),'descend');
    L = L(:,I);
    V(logical(V)) = eig_values;
    % convert back
    std_matrix = 1./(sqrt(V));
    std_matrix(std_matrix == Inf) = 0;
    L = real(X*L*std_matrix);
    eig_values = eig_values / N_timesteps;
    Y = L'*X;
end
% Percent of variance explained by each mode
expvar = eig_values/sum(eig_values) * 100;
end