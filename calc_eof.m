function [eofs, pcs, lambda, expvar] = calc_eof(X,n_eof,method)
% Compute the Nth first EOFs of matrix X(TIME,MAP).
% Inputs:
%     X: data matrix, X(time,space);
%     n_eof: the first n_eof
%     method: 'svd' or 'eig'. If it is not defined, the default method would be eigs.
%     
%  Outputs:
%     eofs : eofs(N,MAP)
%     pcs : principle components, pcs(n,:)
%     lambda : eigenvalue
%     expvar: The fraction of total variance explained by each EOF mode.
%
%
% See also EIG, EIGS, SVD, SVDS
%
% Ref: H. Bjornson and S.A. Venegas: "A manual for EOF and SVD - 
%      Analyses of climatic Data" 1997
%---------------------------------------------------------------%
% Author:
%	Zelun Wu
%	zelunwu@stu.xmu.edu.cn, zelunwu@udel.edu
%	Xiamen University, University of Delaware
%	6th, January, 2020

narginchk(1,inf);

%% Input checking
if nargin == 1
    n_eof = 5;
elseif nargin == 2
    method = 'eigs';
end

[n,p] = size(X); % n is the time demention, p is the space dimention
%% error checking.
assert(n_eof<=n, 'Input error: The number n cannot exceed the number of time steps in your data matrix A.  Time steps are inferred as the third dimension of A.') 
%---------------------------------------------------------------%
%% EOF
% Remove time mean
X = double(detrend(X,'constant'));

% Eigen analysis
switch method
    case {'eigs','eig','e'}
        % Covariance Matrix
        if n >= p
            C = X' * X;
        else
            C = X * X';
        end
%         C = C/n;
        [eigen_matrix,lambda_matrix] = eigs(C,n_eof);
        if n < p
            eigen_matrix = X' * eigen_matrix;
            %   sq = (sqrt(diag(L))+eps)';
            sq = (sqrt(diag(lambda_matrix)))';
            sq = sq(ones(1,p),:);
            eigen_matrix = eigen_matrix ./ sq;
        end
        % get pcs
        
        if n >= p
            pcs = (X * eigen_matrix)';
        else
            pcs =  eigen_matrix' * X';
        end
        
        eofs = eigen_matrix';
        
        % Amount of variance explained a 0.1 pres et en %
        lambda = diag(lambda_matrix)';
        sum_lambda = sum(lambda);
        expvar = lambda./sum_lambda * 100;
    case {'svd','svds','s'}
        if n>p
            [U,S,V] = svds(X',n_eof);
    %         [U,S,V] = svd(X');
            eofs = U';
            pcs = S*V';
            lambda = diag(S.^2/n)';
            expvar = lambda/sum(lambda)*100;
        else
            [U,S,V] = svds(X,n_eof);
    %         [U,S,V] = svd(X');
            eofs = V';
            pcs = S*U';
            lambda = diag(S.^2/n)';
            expvar = lambda/sum(lambda)*100;
        end
end
end