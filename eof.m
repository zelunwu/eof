function [eof_maps, pcs, expvar, eig_values] = eof(data, varargin)
% [eof_maps, pcs, expvar, eig_values] = EOF(data, varargin)
% Version 1.0
% Calculate the EOF of a 3d data(lon,lat,time)
%%   Syntax
%       [eof_maps, pcs, expvar, eig_values] = eof(data) 
%         
%       [eof_maps, pcs, expvar, eig_values] = eof(data,lat) 
%           Weighted data with cosd(lat), length(lat) should be equal to size(data,2)
% 
%       [eof_maps, pcs, expvar, eig_values] = eof(data, n_eof) 
%           Only calculate the first n_eof mode. Notice that the expvar will be calculated only with the first n_eof mode. 
%           E.g., if n_eof = 1, then expvar = 100%
%                 
%       [eof_maps, pcs, expvar, eig_values] = eof(data, lat, n_eof) 
%           Weighted data with cosd(lat), then only calculate the first n_th mode
%             
%       [eof_maps, pcs, expvar, eig_values] = eof(data, lat, n_eof, 'std') 
%           Standardize the eof_maps and pcs with sqrt(eigenvalues)
%% Author:
%	Zelun Wu,
%   Ph.D. student of Physical Oceanography,
%	Xiamen University & University of Delaware
%	zelunwu@stu.xmu.edu.cn, zelunwu@udel.edu
%	15th May, 2020

  
%% Error checks 
narginchk(1,inf) 
%% Input parsing

% Set default parameters: 
[N_lon, N_lat, N_time] = size(data);
N_loc = N_lon*N_lat;
n_eof = min([N_loc, N_time]);
flag_std = 0;
if nargin>1
    % Specified n_eof
    for i = 1:length(varargin)
        if isscalar(varargin{i})& (~isstr(varargin{1}))
            n_eof = varargin{i};
            assert(n_eof<=min([N_loc, N_time]),'Input error: n_eof cannot exceed the number of location points and the numbers of timesteps of data');
        end
    end

    if ismatrix(varargin{1})& (~isstr(varargin{1}))
        lat = varargin{1}; 
        assert(isempty(find(lat>90|lat<-90)),'Latitude should be at the range of [-90,90]');
        % weighted with latitude
        [lat_mesh, ~] = meshgrid(lat,[1:N_lon]);
        lat_weight = cosd(repmat(lat_mesh,[1,1,N_time]));
        data = data.*lat_weight;
    end

    % standardize
    tmp = strcmpi(varargin,'std'); 
    if any(tmp) 
        flag_std = 1;
	end
end


%%
[data_2d, in_nonan_locations] = reshape3dto2d(data);
[eof_maps, pcs, eig_values, expvar] = eofcore(data_2d, n_eof);
eof_maps = reshape2dto3d(eof_maps, [N_lon, N_lat, n_eof], in_nonan_locations);

if flag_std
    for n = 1:n_eof
        std_n = sqrt(abs(eig_values(n)));
        eof_maps(:,:,n) = eof_maps(:,:,n)*std_n;
        pcs(n,:) = pcs(n,:)./std_n;
    end
end
end