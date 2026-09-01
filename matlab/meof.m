function [eof_maps1, eof_maps2, pcs, expvar, eig_values] = meof(data1, data2, varargin)
%MEOF Multivariate (combined) EOF of two (lon, lat, time) fields.
%
%   [MAPS1, MAPS2, PCS, EXPVAR, EIG_VALUES] = MEOF(DATA1, DATA2, ...)
%
%   The two fields are flattened, optionally area-weighted and/or
%   variance-normalized, concatenated along space, and passed to EOFCORE.
%   They must share the time length: size(DATA1,3) == size(DATA2,3).
%   Horizontal grids may differ.
%
%   Variable (dimension) weighting
%   ------------------------------
%   SST and SLP, for example, have incompatible units and variances.
%   Pass 'normalize', true (default) to divide each field by its own
%   scalar standard deviation before concatenation so neither field
%   dominates by unit choice. Spatial maps are scaled back afterwards.
%
%   Spatial weighting is the same as EOF: pass lat1/lat2 or weights.
%
%   Syntax
%     meof(data1, data2)
%     meof(data1, data2, lat1, lat2)
%     meof(data1, data2, n_eof)
%     meof(data1, data2, lat1, lat2, n_eof)
%     meof(data1, data2, lat1, lat2, n_eof, 'std')
%     meof(data1, data2, 'lat1', lat1, 'lat2', lat2, 'n_eof', k, ...
%          'normalize', true, 'weighting', 'sqrtcos')
%
%   Name/value options
%     'lat1','lat2'        latitude vectors
%     'weights1','weights2'  explicit spatial weights (see EOF)
%     'weighting'          'sqrtcos' | 'cos' | 'none'
%     'n_eof'              number of modes
%     'normalize'          equalize field variances (default true)
%     'unweight'           unweight spatial maps (default true)
%     'std'                sqrt(lambda) scaling of maps/PCs
%
%   See also EOF, EOFCORE.
%
%   Zelun Wu, 15 May 2020; area + variable weighting 2026.

narginchk(2, inf)
[n_lon1, n_lat1, n_time1] = size(data1);
[n_lon2, n_lat2, n_time2] = size(data2);
assert(n_time1 == n_time2, 'eof:TimeMismatch', ...
    'size(data1,3) must equal size(data2,3).');
n_time = n_time1;

[opt, extra] = parse_meof_options(varargin, n_lon1, n_lat1, n_lon2, n_lat2, n_time);

W1 = eof_make_weights(n_lon1, n_lat1, extra.lat1, extra.weights1, opt.weighting);
W2 = eof_make_weights(n_lon2, n_lat2, extra.lat2, extra.weights2, opt.weighting);

data1w = data1 .* repmat(W1, [1, 1, n_time]);
data2w = data2 .* repmat(W2, [1, 1, n_time]);

[d1, loc1] = reshape3dto2d(data1w);
[d2, loc2] = reshape3dto2d(data2w);

s1 = 1;
s2 = 1;
if extra.normalize
    s1 = field_std(d1);
    s2 = field_std(d2);
    if s1 > 0
        d1 = d1 / s1;
    end
    if s2 > 0
        d2 = d2 / s2;
    end
end

n_loc = size(d1, 1) + size(d2, 1);
n_eof = min(opt.n_eof, min(n_loc, n_time));
[maps, pcs, eig_values, expvar] = eofcore([d1; d2], n_eof);

n1 = size(d1, 1);
maps1 = maps(1:n1, :);
maps2 = maps(n1+1:end, :);

if extra.normalize
    maps1 = maps1 * s1;
    maps2 = maps2 * s2;
end

if opt.unweight
    maps1 = unweight_rows(maps1, W1, loc1);
    maps2 = unweight_rows(maps2, W2, loc2);
end

if opt.flag_std
    for n = 1:n_eof
        std_n = sqrt(abs(eig_values(n)));
        if std_n > 0
            maps1(:, n) = maps1(:, n) * std_n;
            maps2(:, n) = maps2(:, n) * std_n;
            pcs(n, :) = pcs(n, :) / std_n;
        end
    end
end

eof_maps1 = reshape2dto3d(maps1, [n_lon1, n_lat1, n_eof], loc1);
eof_maps2 = reshape2dto3d(maps2, [n_lon2, n_lat2, n_eof], loc2);
end

function s = field_std(d)
d = d - repmat(mean(d, 2), [1, size(d, 2)]);
v = d(:);
s = sqrt(sum(v .* v) / max(numel(v) - 1, 1));
end

function maps = unweight_rows(maps, W, loc)
w = W(:);
w = w(loc);
w(w == 0) = Inf;
maps = maps ./ repmat(w, [1, size(maps, 2)]);
maps(~isfinite(maps)) = 0;
end

function [opt, extra] = parse_meof_options(args, n_lon1, n_lat1, n_lon2, n_lat2, n_time)
% MEOF has two latitude vectors, so positional parsing differs from EOF.
extra.lat1 = [];
extra.lat2 = [];
extra.weights1 = [];
extra.weights2 = [];
extra.normalize = true;

opt.n_eof = inf;
opt.weighting = 'sqrtcos';
opt.unweight = true;
opt.flag_std = false;

i = 1;
numeric_pos = {};
while i <= numel(args)
    a = args{i};
    if ischar(a) || (exist('isstring', 'builtin') && isstring(a))
        key = char(a);
        if strcmpi(key, 'std')
            if i < numel(args) && (islogical(args{i+1}) || ...
                    (isnumeric(args{i+1}) && isscalar(args{i+1})))
                opt.flag_std = logical(args{i+1});
                i = i + 2;
            else
                opt.flag_std = true;
                i = i + 1;
            end
        else
            if i >= numel(args)
                error('eof:MissingValue', 'Option ''%s'' requires a value.', key);
            end
            val = args{i+1};
            switch lower(key)
                case 'lat1'
                    extra.lat1 = val;
                case 'lat2'
                    extra.lat2 = val;
                case {'lat'}
                    extra.lat1 = val;
                    extra.lat2 = val;
                case 'weights1'
                    extra.weights1 = val;
                case 'weights2'
                    extra.weights2 = val;
                case {'n_eof', 'neof', 'n'}
                    opt.n_eof = val;
                case {'weighting', 'weight_method'}
                    opt.weighting = val;
                case 'unweight'
                    opt.unweight = logical(val);
                case {'normalize', 'normalise', 'equalize'}
                    extra.normalize = logical(val);
                otherwise
                    error('eof:UnknownOption', 'Unknown option ''%s''.', key);
            end
            i = i + 2;
        end
    elseif isnumeric(a)
        numeric_pos{end+1} = a; %#ok<AGROW>
        i = i + 1;
    else
        error('eof:BadArgument', 'Unrecognized argument at position %d.', i);
    end
end

% Positional leftovers: (lat1, lat2), (n_eof), or (lat1, lat2, n_eof)
if numel(numeric_pos) == 1
    v = numeric_pos{1};
    if isscalar(v)
        opt.n_eof = v;
    else
        extra.lat1 = v;
        extra.lat2 = v;
    end
elseif numel(numeric_pos) >= 2
    extra.lat1 = numeric_pos{1};
    extra.lat2 = numeric_pos{2};
    if numel(numeric_pos) >= 3 && isscalar(numeric_pos{3})
        opt.n_eof = numeric_pos{3};
    end
end

n_loc_guess = n_lon1 * n_lat1 + n_lon2 * n_lat2;
if ~isfinite(opt.n_eof)
    opt.n_eof = min([n_loc_guess, n_time]);
end
assert(opt.n_eof >= 1, 'eof:BadNEof', 'n_eof must be >= 1.');
end
