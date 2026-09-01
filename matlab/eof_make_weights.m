function W = eof_make_weights(n_lon, n_lat, lat, weights, method)
%EOF_MAKE_WEIGHTS Build a (n_lon x n_lat) spatial weight map for EOF.
%
%   W = EOF_MAKE_WEIGHTS(N_LON, N_LAT, LAT, WEIGHTS, METHOD)
%
%   Precedence: WEIGHTS (applied as-is) > LAT (converted by METHOD) > ones.
%
%   METHOD when LAT is given (default 'sqrtcos'):
%     'sqrtcos'  W(i,j) = sqrt(max(cosd(lat(j)), 0))
%                Standard area weighting: grid-cell area ~ cos(lat), and
%                the covariance inner product uses sqrt(area) on the field
%                (NCAR Climate Data Guide; Baldwin et al., 2009).
%     'cos'      W(i,j) = max(cosd(lat(j)), 0)
%                Legacy v1.x behaviour of this toolbox (generally too
%                strong; kept for bit-reproducible comparisons).
%     'none'     W(:,:) = 1
%
%   WEIGHTS accepted shapes:
%     (n_lat,)            broadcast along longitude
%     (n_lon, n_lat)      one weight per grid cell
%     (n_lon*n_lat, 1)    flattened in MATLAB column-major order
%
%   See also EOF, MEOF.
%
%   Zelun Wu, 2020; spatial-weight rewrite 2026.

if nargin < 5 || isempty(method)
    method = 'sqrtcos';
end
if nargin < 4
    weights = [];
end
if nargin < 3
    lat = [];
end

if ~isempty(weights)
    W = reshape_weights(weights, n_lon, n_lat);
    return
end

W = ones(n_lon, n_lat);
if isempty(lat)
    return
end

lat = lat(:);
assert(numel(lat) == n_lat, ...
    'eof:BadLat', ...
    'length(lat) must equal size(data,2) (the latitude dimension).');
assert(isempty(find(lat > 90 | lat < -90, 1)), ...
    'eof:BadLat', 'Latitude must lie in [-90, 90].');

method = lower(method);
switch method
    case 'sqrtcos'
        wlat = sqrt(max(cosd(lat), 0));
    case 'cos'
        wlat = max(cosd(lat), 0);
    case 'none'
        wlat = ones(size(lat));
    otherwise
        error('eof:BadWeighting', ...
            'Unknown weighting method ''%s''. Use sqrtcos, cos, or none.', method);
end
W = repmat(wlat', [n_lon, 1]);
end

function W = reshape_weights(weights, n_lon, n_lat)
n_loc = n_lon * n_lat;
sz = size(weights);
if isvector(weights) && numel(weights) == n_lat
    W = repmat(weights(:)', [n_lon, 1]);
elseif isequal(sz(1:min(2, ndims(weights))), [n_lon, n_lat]) && numel(weights) == n_loc
    W = reshape(weights, [n_lon, n_lat]);
elseif isvector(weights) && numel(weights) == n_loc
    W = reshape(weights, [n_lon, n_lat]);
else
    error('eof:BadWeights', ...
        ['weights must be (n_lat,), (n_lon, n_lat), or (n_lon*n_lat,).', ...
         ' Got size %s.'], mat2str(sz));
end
end
