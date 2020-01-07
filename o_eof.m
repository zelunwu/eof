function [eof_maps,pcs,lambda,expvar] = o_eof(data,lat,varargin)
% [eof_map,pc,lambda,expvar] = o_eof(data,lat,n_eof,option)
% Compute the Nth first EOFs of matrix data(lon,lat,time).
% Inputs:
%     data: data matrix, data(lon,lat,time).
%     n_eof: the first n_eof
%     option: can be 'method':{'svd','eig'}, 'unit':{'unit'} or
%     'test':{'North';'Monte'}
%     method: 'svd' or 'eig'. If it is not defined, the default method would be eigs.
%     
%  Outputs:
%     eofs : eofs(N,MAP)
%     pcs : principle components, pcs(n,:)
%     lambda : eigenvalue
%     expvar: The fraction of total variance explained by each EOF mode.
%
% Ref: H. Bjornson and S.A. Venegas: "A manual for EOF and SVD - 
%      Analyses of climatic Data" 1997
%---------------------------------------------------------------%
% Author:
%	Zelun Wu
%	zelunwu@stu.xmu.edu.cn, zelunwu@udel.edu
%	Xiamen University, University of Delaware
%	6th, January, 2020
%---------------------------------------------------------------%


%% error check
narginchk(1,inf);
size_data = size(data);
if ndims(data) > 1 
    p_space_dim = prod(size_data(1:end-1)); % space dimention number
    n_time_dim = size_data(end); % time dimention number
else
    error('Input data should be 2 or 3 dimentional.');
end
%% nan values check
if nargin>1
	if isscalar(varargin{1})
        n_eof = varargin{1}; 
    else
        n_eof = 5;
    end
    assert(n_eof<=n_time_dim, 'Input error: The number n cannot exceed the number of time steps in your data matrix A.  Time steps are inferred as the third dimension of A.') 
    % Has the user defined a mask? 
    tmp = strcmpi(varargin,'mask'); 
    if any(tmp) 
        mask = varargin{find(tmp)+1}; 
        assert(isequal(size(mask),[size(A,1) size(A,2)])==1,'Input error: Mask dimensions must match the first two dimensions of the data matrix A.') 
        assert(islogical(mask)==1,'Input error: mask must be logical.') 
    else
        mask = ~any(isnan(data),3); 
    end
   
	% Has the user want to do sinificance test?
    tmp = strcmpi(varargin,'test'); 
    logic_test = any(tmp);
	if logic_test
        test_method = varargin{find(tmp)+1};
    end
    % Has the user define unit?
	tmp = strcmpi(varargin,'unit'); 
	if any(tmp)
        unit = varargin{find(tmp)+1};
    end
    
    % Has the user define eigenvalue calculation method?
	tmp = strcmpi(varargin,'method'); 
	if any(tmp)
        method = varargin{find(tmp)+1};
    else
        method = 'eigs';
    end
end

%% Reshape
lat = repmat(reshape(lat,[1,size(data,2)]),[size(data,1),1,size(data,3)]);

data = data.*cosd(lat); % weight
data = double(reshape(data,[p_space_dim,n_time_dim])'); % reshape to n_time_dim row and p_space_dim column.
data = data(:,mask(:)); % remove nan
[eofs, pcs, lambda, expvar] = calc_eof(data,n_eof,method);

eof_maps = nan(n_eof,numel(mask));   
eof_maps(:,mask(:)) = eofs; 
eof_maps = reshape(eof_maps,[n_eof,size(mask,1),size(mask,2)]); 
eof_maps = permute(eof_maps,[2 3 1]); 

%% Flip signs to provide consistent results:
ind = sign(pcs(:,1))<0; 
eof_maps(:,:,ind) = eof_maps(:,:,ind)*(-1); 
pcs(ind,:) = pcs(ind,:)*-1; 

if logic_test
    switch test_method 
        case {'North','north','n'} % North 1982
            L = lambda;
            L_e = L.*sqrt(2/n_time_dim); % Error
            figure;
            x = [1:n_eof];
            errorbar(x,L,L_e,'LineWidth',2);
            xlim([0,n_eof+1]);
            xticks([1:n_eof]);
            xlabel('Modes');
            ylabel(['Eigenvalue (',unit,'^2)']);
            set(gca,'fontsize',15);
        case {'monte','Monte','m','mentecarlo','Mente Carlo','MenteCarlo','mc'}
            ;
    end
end
end
