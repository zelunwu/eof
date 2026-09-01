function results = test_eof
%TEST_EOF Unit tests for the MATLAB EOF toolbox.
%   Run from the repository root:
%     addpath(pwd); results = test_eof
%   Each field of RESULTS is true on success.

this_dir = fileparts(mfilename('fullpath'));
matlab_dir = fileparts(this_dir);
addpath(matlab_dir);

results = struct();
results.rank2_recovers_modes = test_rank2();
results.expvar_sums_near_100 = test_expvar();
results.nan_mask_preserved = test_nan_mask();
results.area_weight_prefers_equator = test_area_weight();
results.meof_normalize_balances_fields = test_meof_normalize();
results.parse_name_value = test_parse();

names = fieldnames(results);
failed = {};
for i = 1:numel(names)
    if ~results.(names{i})
        failed{end+1} = names{i}; %#ok<AGROW>
    end
end
if isempty(failed)
    fprintf('All %d MATLAB EOF tests passed.\n', numel(names));
else
    error('eof:TestFailure', 'Failed: %s', strjoin(failed, ', '));
end
end

function ok = test_rank2()
n_lon = 8; n_lat = 6; n_time = 40;
[x, y, t] = ndgrid(1:n_lon, 1:n_lat, 1:n_time);
p1 = sin(2*pi*x/n_lon) .* cos(2*pi*y/n_lat);
p2 = cos(2*pi*x/n_lon) .* sin(2*pi*y/n_lat);
a1 = sin(2*pi*(0:n_time-1)/n_time);
a2 = cos(4*pi*(0:n_time-1)/n_time);
data = zeros(n_lon, n_lat, n_time);
for k = 1:n_time
    data(:,:,k) = p1 * a1(k) + 0.3 * p2 * a2(k);
end
[maps, pcs, expvar] = eof(data, 2);
anom = data - mean(data, 3);
reco = zeros(size(data));
for k = 1:n_time
    reco(:,:,k) = maps(:,:,1)*pcs(1,k) + maps(:,:,2)*pcs(2,k);
end
rel = norm(reco(:) - anom(:)) / (norm(anom(:)) + eps);
ok = rel < 1e-8 && expvar(1) > expvar(2);
end

function ok = test_expvar()
n_lon = 5; n_lat = 4; n_time = 20;
rng(0);
data = randn(n_lon, n_lat, n_time);
[~, ~, expvar] = eof(data);
ok = abs(sum(expvar) - 100) < 1e-6;
end

function ok = test_nan_mask()
n_lon = 5; n_lat = 4; n_time = 12;
data = randn(n_lon, n_lat, n_time);
data(1,1,:) = NaN;
[maps, ~] = eof(data, 3);
ok = all(isnan(maps(1,1,:))) && all(isfinite(maps(2,2,:)));
end

function ok = test_area_weight()
% Identical-amplitude signals at equator vs high latitude.
% sqrt(cos lat) must give the equatorial mode more variance.
lat = [-75; 0; 75];
n_lon = 4; n_time = 80;
data = zeros(n_lon, 3, n_time);
t = 0:n_time-1;
for i = 1:n_lon
    data(i, 2, :) = sin(2*pi*t/20);          % equator
    data(i, 3, :) = sin(2*pi*t/13 + 0.4);    % high lat, incommensurate
end
[~, ~, exp_unw] = eof(data, 2);
[~, ~, exp_w] = eof(data, 'lat', lat, 'n_eof', 2, 'weighting', 'sqrtcos');
% Without weights the two spatial peaks have similar power.
% With weights, mode 1 must explain more than the unweighted mode 1
% relative to mode 2, i.e. equator/high-lat ratio increases.
ok = (exp_w(1)/exp_w(2)) > (exp_unw(1)/exp_unw(2)) * 1.05;
end

function ok = test_meof_normalize()
n_lon = 4; n_lat = 3; n_time = 40;
t = 0:n_time-1;
a = sin(2*pi*t/10);
b = cos(2*pi*t/10);
d1 = zeros(n_lon, n_lat, n_time);
d2 = zeros(n_lon, n_lat, n_time);
for k = 1:n_time
    d1(:,:,k) = a(k);
    d2(:,:,k) = 100 * b(k);
end
[~, ~, ~, exp_raw] = meof(d1, d2, 'normalize', false, 'n_eof', 2);
[~, ~, ~, exp_norm] = meof(d1, d2, 'normalize', true, 'n_eof', 2);
ok = (exp_raw(1)/exp_raw(2) > 10) && (abs(exp_norm(1)-exp_norm(2)) < 5);
end

function ok = test_parse()
n_lon = 3; n_lat = 3; n_time = 10;
data = randn(n_lon, n_lat, n_time);
lat = [-10, 0, 10];
[m1, p1] = eof(data, lat, 2);
[m2, p2] = eof(data, 'lat', lat, 'n_eof', 2);
ok = max(abs(m1(:)-m2(:))) < 1e-10 && max(abs(p1(:)-p2(:))) < 1e-10;
end
