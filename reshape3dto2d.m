function [data_2d, in_nonan_locations] = reshape3dto2d(data)
% [data_2d, in_nonan_locations] = RESHAPE3DTO2D(data)
% Version 1.0
% Subroutine of EOF toolbox.
% Reshape 3d data(lon,lat,time) into a 2d data(lon*lat,time) and remove the
% nan values.
%% Author:
%	Zelun Wu,
%   Ph.D. student of Physical Oceanography,
%	Xiamen University & University of Delaware
%	zelunwu@stu.xmu.edu.cn, zelunwu@udel.edu
%	15th May, 2020

size_3d = size(data);
data = reshape(data, [size_3d(1)*size_3d(2),size_3d(3)]);
[in_nonan_locations, in_nonan_timesteps] = find(~isnan(data));
in_nonan_locations = unique(in_nonan_locations);
data_2d = data(in_nonan_locations,:);
end