# Empirical Orthogonal Function Toolbox

EOF / MEOF (combined EOF) for oceanic and atmospheric fields, with **spatial
area weighting** and matching implementations in **MATLAB**, **Python**, and
**C++**.

On a regular latitude–longitude grid, polar cells cover less area than
equatorial cells. If every cell is given equal weight, high latitudes
dominate the covariance. The standard fix (NCAR Climate Data Guide;
Baldwin et al., 2009) is to multiply the field by

```text
W = sqrt(max(cos(latitude), 0))
```

before the eigenproblem, then (by default) divide the spatial maps by `W`
so they are returned in the physical units of the input.

v1 of this toolbox either skipped that step or multiplied by `cos(lat)`
without unweighting the maps. v2 makes weighting a first-class argument
and uses `sqrt(cos(lat))` as the default when latitudes are supplied.

## Layout

```text
matlab/     .m functions (File Exchange)
python/     package eof_toolbox  (classes Eof, MultivariateEof, …)
cpp/        include/eof  public headers
            src/         class implementations (Eigen)
```

C++ keeps both `include/` and `src/`: declarations stay in the public
header, method bodies live in a compiled library so dependents do not
recompile the eigenproblem.

## MATLAB

```matlab
addpath /path/to/eof/matlab
```

The original positional API still works.

```matlab
[maps, pcs, expvar, eig] = eof(data);
[maps, pcs, expvar, eig] = eof(data, lat);
[maps, pcs, expvar, eig] = eof(data, lat, 8);
[maps, pcs, expvar, eig] = eof(data, lat, 8, 'std');

[maps, pcs, expvar, eig] = eof(data, 'lat', lat, 'n_eof', 8, ...
    'weighting', 'sqrtcos', 'unweight', true);
[maps, pcs, expvar, eig] = eof(data, 'weights', sqrt(area));

[m1, m2, pcs, expvar, eig] = meof(sst, slp, lat_sst, lat_slp, 5);
[m1, m2, pcs, expvar, eig] = meof(sst, slp, 'lat1', lat_sst, 'lat2', lat_slp, ...
    'normalize', true, 'n_eof', 5);
```

`data` is `n_lon × n_lat × n_time`. Cells that are NaN at any time are
dropped (land mask). `expvar` is the percent of **total** (weighted)
variance, including modes you did not request.

| `weighting` | weight              | when to use                         |
|-------------|---------------------|-------------------------------------|
| `sqrtcos`   | `sqrt(cosd(lat))`   | default; correct area inner product |
| `cos`       | `cosd(lat)`         | bit-compare with toolbox v1         |
| `none`      | `1`                 | ignore latitudes                    |

```matlab
cd matlab/tests
results = test_eof
```

## Python

Primary API is object-oriented (`Eof`, `MultivariateEof`, `PackedEof`,
`SpatialWeights`). `eof()` / `meof()` are thin wrappers with the MATLAB
return signature.

```bash
cd python
pip install -e ".[test]"
pytest -q
```

```python
from eof_toolbox import Eof, MultivariateEof

solver = Eof(data, lat=lat, n_eof=8)          # (n_lon, n_lat, n_time)
maps = solver.maps()                          # physical units (unweighted)
pcs = solver.pcs()
expvar = solver.expvar
anom = solver.reconstruct()

me = MultivariateEof(sst, slp, lat1=lat_sst, lat2=lat_slp, n_eof=5)
maps_sst, maps_slp = me.maps()
```

## C++

Compiled library: `include/eof/eof.hpp` + `src/eof.cpp`. Requires
[Eigen](https://eigen.tuxfamily.org/) 3.3+ (fetched by CMake if missing).

```bash
cmake -S cpp -B cpp/build
cmake --build cpp/build
ctest --test-dir cpp/build
```

```cpp
#include <eof/eof.hpp>

eof::PackedEof packed(X, n_eof);   // X: n_loc x n_time

eof::GridOptions opt;
opt.n_eof = 8;
eof::GridEof solver(data, n_lon, n_lat, n_time, lat, n_lat, opt);
auto maps = solver.maps();         // physical units
auto pcs = solver.pcs();

eof::SpatialWeights W(n_lon, n_lat);
eof::MultivariateEof me(sst, n_lon, n_lat, slp, n_lon, n_lat, n_time, W, W);
```

Grid arrays are column-major (lon, then lat, then time), matching MATLAB.

## Reconstruction

With default unweighting,

```text
anomaly(x, t) ≈ Σ_k  maps(x, k) * W(x) * pcs(k, t)
```

## References

- Björnsson, H. and Venegas, S.A., 1997. A manual for EOF and SVD analyses
  of climatic data. CCGCR Report, 97(1).
- Baldwin, M.P., et al., 2009. On the weighting of climate data.
- NCAR Climate Data Guide, Empirical Orthogonal Function analysis.

## Author

**Zelun Wu**  
zelunwu@stu.xmu.edu.cn, zelunwu@udel.edu  
College of Ocean and Earth Science, Xiamen University  
College of Earth, Ocean and Environment, University of Delaware
