function [eof_maps, pcs, expvar, eig_values] = eof(data, varargin)
%EOF Empirical Orthogonal Functions of a (lon, lat, time) field.
%
%   [EOF_MAPS, PCS, EXPVAR, EIG_VALUES] = EOF(DATA, ...)
%
%   DATA is an n_lon x n_lat x n_time array (NaNs mark unused cells).
%   The time mean is removed. Spatial modes are returned on the same
%   horizontal grid; missing cells remain NaN.
%
%   Spatial (area / latitude) weighting
%   -----------------------------------
%   Equal-area grids do not need extra weights. On a regular lat-lon
%   grid, polar cells are smaller, so they must not count as much as
%   equatorial cells in the covariance. Pass LAT to apply the standard
%   weight sqrt(max(cosd(lat),0)) (NCAR Climate Data Guide).
%
%   After the eigenproblem, spatial maps are *unweighted* (divided by
%   W) so EOF_MAPS is in the physical units of DATA. Set 'unweight'
%   false to keep maps in weighted space (orthonormal columns).
%
%   Syntax
%     eof(data)
%     eof(data, lat)
%     eof(data, n_eof)
%     eof(data, lat, n_eof)
%     eof(data, lat, n_eof, 'std')
%     eof(data, 'lat', lat, 'n_eof', k)
%     eof(data, 'weights', W)                 % W multiplied as-is
%     eof(data, 'lat', lat, 'weighting', 'sqrtcos')  % default
%     eof(data, 'lat', lat, 'weighting', 'cos')      % legacy v1
%     eof(data, 'lat', lat, 'unweight', false)
%
%   Name/value options
%     'lat'        latitude vector, length n_lat, degrees
%     'weights'    (n_lat,), (n_lon,n_lat), or (n_lon*n_lat,) weights
%                  multiplied onto DATA before the eigenproblem
%     'weighting'  'sqrtcos' (default) | 'cos' | 'none'  (used with lat)
%     'n_eof'      number of modes (default: min(space, time))
%     'unweight'   divide maps by W after analysis (default true)
%     'std'        scale maps by sqrt(lambda) and PCs by 1/sqrt(lambda)
%
%   Outputs
%     EOF_MAPS     n_lon x n_lat x n_eof spatial patterns
%     PCS          n_eof x n_time principal components
%     EXPVAR       n_eof x 1 percent of total (weighted) variance
%     EIG_VALUES   n_eof x 1 covariance eigenvalues (variance per mode)
%
%   Reconstruction of anomalies (with default unweight true):
%     data_anom ≈ (eof_maps .* W) contracted with pcs, or equivalently
%     maps_weighted * pcs, where maps_weighted = eof_maps .* W.
%
%   See also MEOF, EOFCORE, EOF_MAKE_WEIGHTS.
%
%   Zelun Wu, Xiamen University & University of Delaware, 15 May 2020.
%   Area weighting (sqrt(cos lat)), unweighting, name/value API: 2026.

narginchk(1, inf)
[n_lon, n_lat, n_time] = size(data);
n_loc = n_lon * n_lat;
opt = eof_parse_options(varargin, n_lon, n_lat, n_time, n_loc);

W = eof_make_weights(n_lon, n_lat, opt.lat, opt.weights, opt.weighting);
W3 = repmat(W, [1, 1, n_time]);
data_w = data .* W3;

[data_2d, in_nonan_locations] = reshape3dto2d(data_w);
n_valid = size(data_2d, 1);
n_eof = min(opt.n_eof, min(n_valid, n_time));

[eof_2d, pcs, eig_values, expvar] = eofcore(data_2d, n_eof);

if opt.unweight
    w_flat = W(:);
    w_valid = w_flat(in_nonan_locations);
    w_valid(w_valid == 0) = Inf;  % pole weights: leave map at 0
    eof_2d = eof_2d ./ repmat(w_valid, [1, n_eof]);
    eof_2d(~isfinite(eof_2d)) = 0;
end

if opt.flag_std
    for n = 1:n_eof
        std_n = sqrt(abs(eig_values(n)));
        if std_n > 0
            eof_2d(:, n) = eof_2d(:, n) * std_n;
            pcs(n, :) = pcs(n, :) / std_n;
        end
    end
end

eof_maps = reshape2dto3d(eof_2d, [n_lon, n_lat, n_eof], in_nonan_locations);
end
