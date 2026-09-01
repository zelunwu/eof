function opt = eof_parse_options(args, n_lon, n_lat, n_time, n_loc)
%EOF_PARSE_OPTIONS Parse varargin for EOF / MEOF.
%
%   Accepts the original positional API and name/value pairs:
%
%     eof(data)
%     eof(data, lat)
%     eof(data, n_eof)
%     eof(data, lat, n_eof)
%     eof(data, lat, n_eof, 'std')
%     eof(data, 'lat', lat, 'n_eof', k, 'weights', W, ...
%              'weighting', 'sqrtcos', 'unweight', true, 'std', true)
%
%   Fields of OPT:
%     n_eof      number of modes (default min(n_loc, n_time))
%     lat        latitude vector or []
%     weights    user weight map or []
%     weighting  'sqrtcos' | 'cos' | 'none'
%     unweight   if true, divide spatial maps by W after the eigenproblem
%                so maps are in the physical units of the input field
%     flag_std   if true, map *= sqrt(lambda), pc /= sqrt(lambda)
%
%   See also EOF, MEOF.
%
%   Zelun Wu, 2020; option parser 2026.

opt.n_eof = min([n_loc, n_time]);
opt.lat = [];
opt.weights = [];
opt.weighting = 'sqrtcos';
opt.unweight = true;
opt.flag_std = false;

i = 1;
while i <= numel(args)
    a = args{i};
    if ischar(a) || (exist('isstring', 'builtin') && isstring(a))
        key = char(a);
        if strcmpi(key, 'std')
            % Bare flag (legacy) or name/value.
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
                case {'lat', 'latitude'}
                    opt.lat = val;
                case {'weights', 'weight'}
                    opt.weights = val;
                case {'n_eof', 'neof', 'n'}
                    opt.n_eof = val;
                case {'weighting', 'weight_method', 'lat_weight'}
                    opt.weighting = val;
                case 'unweight'
                    opt.unweight = logical(val);
                otherwise
                    error('eof:UnknownOption', 'Unknown option ''%s''.', key);
            end
            i = i + 2;
        end
    elseif isnumeric(a) && isscalar(a)
        opt.n_eof = a;
        i = i + 1;
    elseif isnumeric(a) && ~isscalar(a)
        opt.lat = a;
        i = i + 1;
    else
        error('eof:BadArgument', 'Unrecognized argument at position %d.', i);
    end
end

assert(opt.n_eof >= 1, 'eof:BadNEof', 'n_eof must be >= 1.');
assert(opt.n_eof <= min([n_loc, n_time]), 'eof:BadNEof', ...
    'n_eof cannot exceed min(n_locations, n_timesteps).');
end
