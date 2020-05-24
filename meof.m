function [eof_maps1, eof_maps2, pcs, expvar, eig_values] = meof(data1, data2, varargin)
% [eof_maps1, eof_maps2, pcs, expvar, eig_values] = meof(data1, data2, varargin)
% Version 1.0
% Calculate the Multivariable EOF of a two 3d data1 and data2(lon,lat,time)
%%   Syntax
%       [eof_maps, pcs, expvar, eig_values] = eof(data1, data2) 
%         
%       [eof_maps, pcs, expvar, eig_values] = eof(data1, data2, lat1, lat2) 
%           Weighted data with cosd(lat), length(lat*) should be equal to size(data*,2)
% 
%       [eof_maps, pcs, expvar, eig_values] = eof(data1, data2, n_eof) 
%           Only calculate the first n_eof mode. Notice that the expvar will be calculated only with the first n_eof mode. 
%           E.g., if n_eof = 1, then expvar = 100%
%                 
%       [eof_maps, pcs, expvar, eig_values] = eof(data1, data2, lat1, lat2, n_eof) 
%           Weighted data with cosd(lat), then only calculate the first n_th mode
%             
%       [eof_maps, pcs, expvar, eig_values] = eof(data1, data2, lat1, lat2, n_eof, 'std') 
%           Standardize the eof_maps and pcs with sqrt(eigenvalues)
%% Author:
%	Zelun Wu,
%   Ph.D. student of Physical Oceanography,
%	Xiamen University & University of Delaware
%	zelunwu@stu.xmu.edu.cn, zelunwu@udel.edu
%	15th May, 2020

  
%% Error checks 
narginchk(1,inf) 
[N_lon1, N_lat1, N_time1] = size(data1);
[N_lon2, N_lat2, N_time2] = size(data2);
assert(N_time1 == N_time2,'Timesteps of data1 and data3 should be the same, i.e. size(data1,3) = size(data2,3)');
N_time = N_time1;
%% Input parsing

% Set default parameters: 
[data1_2d, in_nonan_locations1] = reshape3dto2d(data1);
[data2_2d, in_nonan_locations2] = reshape3dto2d(data2);

N_loc = size(data1_2d,1)+size(data2_2d,1);
n_eof = min([N_loc, N_time1]);
flag_std = 0;
if nargin>2
    
    % Specified n_eof
    for i = 1:length(varargin)
        if isscalar(varargin{i})
            n_eof = varargin{i};
            assert(n_eof<=min([N_loc, N_time]),'Input error: n_eof cannot exceed the number of location points and the numbers of timesteps of data');
        end
    end
    
    if nargin > 3 & ismatrix(varargin{1}) & ismatrix(varargin{2}) & length(varargin{1}) == N_lat1 & length(varargin{2}) == N_lat2
        lat1 = varargin{1}; 
        assert(isempty(find(lat1>90|lat1<-90)),'Latitude should be at the range of [-90,90]');
        lat2 = varargin{2};
        assert(isempty(find(lat2>90|lat2<-90)),'Latitude should be at the range of [-90,90]');
        % weighted with latitude
        [lat1_mesh, ~] = meshgrid(lat1,[1:N_lon1]);
        lat1_weight = cosd(repmat(lat1_mesh,[1,1,N_time]));
        data1 = data1.*lat1_weight;
        [lat2_mesh, ~] = meshgrid(lat2,[1:N_lon2]);
        lat2_weight = cosd(repmat(lat2_mesh,[1,1,N_time]));
        data2 = data2.*lat2_weight;
    end

    if isscalar(varargin{1})
        n_eof = varargin{1}; 
        assert(n_eof<=min([N_loc, N_time]),'Input error: n_eof cannot exceed the number of location points and the numbers of timesteps of data');
    end

    % standardize
    tmp = strcmpi(varargin,'std'); 
    if any(tmp) 
        flag_std = 1;
	end
end

%%

[data1_2d, in_nonan_location1] = reshape3dto2d(data1);
[data2_2d, in_nonan_location2] = reshape3dto2d(data2);
[N_loc1, N_time] = size(data1_2d);

data_2d = cat(1,data1_2d,data2_2d);

[eof_maps, pcs, eig_values, expvar] = eofcore(data_2d, n_eof);

if flag_std
    for n = 1:n_eof
        std_n = sqrt(abs(eig_values(n)));
        eof_maps(:,n) = eof_maps(:,n)*std_n;
        pcs(n,:) = pcs(n,:)./std_n;
    end
end

eof_maps1 = eof_maps([1:N_loc1],:);
eof_maps1 = reshape2dto3d(eof_maps1, [N_lon1, N_lat1, n_eof], in_nonan_locations1);

eof_maps2 = eof_maps([N_loc1+1:end],:);
eof_maps2 = reshape2dto3d(eof_maps2, [N_lon2, N_lat2, n_eof], in_nonan_locations2);
end