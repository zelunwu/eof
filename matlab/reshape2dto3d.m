function data_3d = reshape2dto3d(data_2d, size_3d, in_nonan_locations)
%RESHAPE2DTO3D Scatter a (space, mode) matrix back onto a (lon, lat, mode) grid.
%
%   DATA_3D = RESHAPE2DTO3D(DATA_2D, SIZE_3D, IN_NONAN_LOCATIONS)
%
%   SIZE_3D is [n_lon, n_lat, n_mode]. Locations not listed in
%   IN_NONAN_LOCATIONS are filled with NaN (e.g. land mask).
%
%   See also RESHAPE3DTO2D, EOF.
%
%   Zelun Wu, 15 May 2020.

data_3d = nan(size_3d);
flat = reshape(data_3d, [size_3d(1) * size_3d(2), size_3d(3)]);
flat(in_nonan_locations, :) = data_2d;
data_3d = reshape(flat, size_3d);
end
