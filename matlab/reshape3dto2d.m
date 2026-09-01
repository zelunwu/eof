function [data_2d, in_nonan_locations] = reshape3dto2d(data)
%RESHAPE3DTO2D Flatten (lon, lat, time) to (space, time) and drop bad rows.
%
%   [DATA_2D, IN_NONAN_LOCATIONS] = RESHAPE3DTO2D(DATA)
%
%   DATA is (n_lon x n_lat x n_time). Flattening uses MATLAB column-major
%   order, so space index = lon + (lat-1)*n_lon.
%
%   A spatial location is kept only if it is finite at *every* time step.
%   (v1 kept any cell with at least one valid value, which left NaNs in
%   the matrix and broke the eigenproblem.)
%
%   DATA_2D is (n_valid x n_time). IN_NONAN_LOCATIONS indexes the kept
%   rows of the flattened (n_lon*n_lat) space axis.
%
%   See also RESHAPE2DTO3D, EOF.
%
%   Zelun Wu, 15 May 2020; NaN policy clarified 2026.

size_3d = size(data);
if numel(size_3d) < 3
    size_3d(3) = 1;
end
n_space = size_3d(1) * size_3d(2);
data = reshape(data, [n_space, size_3d(3)]);
in_nonan_locations = find(all(isfinite(data), 2));
data_2d = data(in_nonan_locations, :);
end
