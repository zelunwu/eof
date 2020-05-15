function data_3d = reshape2dto3d(data_2d, size_3d, in_nonan_locations)
% data_3d = RESHAPE2DTO3D(data_2d, size_3d, in_nonan_locations)
% Version 1.0
% Subroutine of EOF toolbox.
% Reshape 2d data(lon*lat,time) into a 3d data(lon*lat,time) and insert the
% nan locations.
%% Author:
%	Zelun Wu,
%   Ph.D. student of Physical Oceanography,
%	Xiamen University & University of Delaware
%	zelunwu@stu.xmu.edu.cn, zelunwu@udel.edu
%	15th May, 2020

data_3d = nan(size_3d);
data_3d = reshape(data_3d, [size_3d(1)*size_3d(2),size_3d(3)]);
data_3d(in_nonan_locations,:) = data_2d;
data_3d = reshape(data_3d,size_3d);
end