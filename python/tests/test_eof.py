from __future__ import annotations

import numpy as np
import pytest

from eof_toolbox import Eof, MultivariateEof, PackedEof, SpatialWeights, eof, eofcore, latitude_weights, meof
from eof_toolbox.weights import make_weights


def _rank2_field(n_lon=8, n_lat=6, n_time=40):
    x = np.arange(n_lon)[:, None, None]
    y = np.arange(n_lat)[None, :, None]
    t = np.arange(n_time)[None, None, :]
    p1 = np.sin(2 * np.pi * x / n_lon) * np.cos(2 * np.pi * y / n_lat)
    p2 = np.cos(2 * np.pi * x / n_lon) * np.sin(2 * np.pi * y / n_lat)
    a1 = np.sin(2 * np.pi * t / n_time)
    a2 = np.cos(4 * np.pi * t / n_time)
    return p1 * a1 + 0.3 * p2 * a2


def test_core_orthonormal_and_reconstruction():
    rng = np.random.default_rng(0)
    X = rng.standard_normal((12, 30))
    L, Y, eig, expvar = eofcore(X, n_eof=5)
    gram = L.T @ L
    assert np.allclose(gram, np.eye(5), atol=1e-8)
    Xa = X - X.mean(axis=1, keepdims=True)
    reco = L @ Y
    # 5 modes cannot reconstruct all of a random matrix; residual must
    # still drop vs the mean-only baseline.
    assert np.linalg.norm(Xa - reco) < np.linalg.norm(Xa)
    assert expvar[0] >= expvar[-1]
    assert eig[0] >= eig[-1]


def test_full_rank_expvar_sums_to_100():
    rng = np.random.default_rng(1)
    data = rng.standard_normal((5, 4, 20))
    _, _, expvar, _ = eof(data)
    assert expvar.shape[0] == min(20, 20)
    assert abs(expvar.sum() - 100.0) < 1e-6


def test_rank2_grid_reconstruction():
    data = _rank2_field()
    maps, pcs, expvar, _ = eof(data, n_eof=2)
    anom = data - data.mean(axis=2, keepdims=True)
    reco = np.einsum("ijm,mt->ijt", maps, pcs)
    rel = np.linalg.norm(reco - anom) / np.linalg.norm(anom)
    assert rel < 1e-8
    assert expvar[0] > expvar[1]


def test_nan_mask():
    rng = np.random.default_rng(2)
    data = rng.standard_normal((5, 4, 12))
    data[1, 1, :] = np.nan
    maps, *_ = eof(data, n_eof=3)
    assert np.all(np.isnan(maps[1, 1, :]))
    assert np.all(np.isfinite(maps[0, 0, :]))


def test_latitude_weights_sqrtcos():
    lat = np.array([-60.0, 0.0, 60.0])
    W = latitude_weights(lat, n_lon=2, method="sqrtcos")
    assert W.shape == (2, 3)
    assert np.allclose(W[0], W[1])
    assert W[0, 1] > W[0, 0]
    assert np.allclose(W[0, 1], 1.0)


def test_area_weight_prefers_equator():
    lat = np.array([-75.0, 0.0, 75.0])
    n_lon, n_time = 4, 80
    t = np.arange(n_time)
    data = np.zeros((n_lon, 3, n_time))
    data[:, 1, :] = np.sin(2 * np.pi * t / 20)
    data[:, 2, :] = np.sin(2 * np.pi * t / 13 + 0.4)
    _, _, exp_unw, _ = eof(data, n_eof=2)
    _, _, exp_w, _ = eof(data, lat=lat, n_eof=2, weighting="sqrtcos")
    assert (exp_w[0] / exp_w[1]) > (exp_unw[0] / exp_unw[1]) * 1.05


def test_custom_weights_shape():
    data = _rank2_field(n_lon=4, n_lat=3, n_time=15)
    W = np.linspace(0.5, 1.5, 4 * 3).reshape(4, 3)
    maps_a, *_ = eof(data, weights=W, n_eof=2)
    maps_b, *_ = eof(data, weights=W.reshape(-1, order="F"), n_eof=2)
    assert np.allclose(maps_a, maps_b, equal_nan=True)


def test_standardize_scales_pcs():
    data = _rank2_field()
    _, pcs, _, eig = eof(data, n_eof=2, standardize=False)
    _, pcs_s, _, _ = eof(data, n_eof=2, standardize=True)
    ratio = np.std(pcs[0]) / np.std(pcs_s[0])
    assert np.isclose(ratio, np.sqrt(eig[0]), rtol=1e-6)


def test_meof_normalize_balances_variance():
    n_lon, n_lat, n_time = 4, 3, 40
    t = np.arange(n_time)
    a = np.sin(2 * np.pi * t / 10)
    b = np.cos(2 * np.pi * t / 10)
    d1 = np.broadcast_to(a, (n_lon, n_lat, n_time)).copy()
    d2 = 100.0 * np.broadcast_to(b, (n_lon, n_lat, n_time)).copy()
    *_, exp_raw, _ = meof(d1, d2, n_eof=2, normalize=False)
    *_, exp_norm, _ = meof(d1, d2, n_eof=2, normalize=True)
    assert exp_raw[0] / exp_raw[1] > 10
    assert abs(exp_norm[0] - exp_norm[1]) < 5.0


def test_meof_time_mismatch():
    with pytest.raises(ValueError):
        meof(np.zeros((2, 2, 5)), np.zeros((2, 2, 6)))


def test_make_weights_rejects_bad_lat():
    with pytest.raises(ValueError):
        make_weights(2, 2, lat=np.array([-10.0, 100.0]))


def test_eof_class_matches_function():
    data = _rank2_field()
    solver = Eof(data, n_eof=2)
    maps, pcs, expvar, eig = eof(data, n_eof=2)
    assert np.allclose(solver.maps(), maps, equal_nan=True)
    assert np.allclose(solver.pcs(), pcs)
    assert np.allclose(solver.expvar, expvar)
    reco = solver.reconstruct()
    anom = data - data.mean(axis=2, keepdims=True)
    assert np.linalg.norm(reco - anom) / np.linalg.norm(anom) < 1e-8


def test_packed_eof_class():
    rng = np.random.default_rng(3)
    X = rng.standard_normal((8, 20))
    solver = PackedEof(X, n_eof=3)
    L, Y, eig, expvar = eofcore(X, n_eof=3)
    assert np.allclose(solver.modes, L)
    assert np.allclose(solver.pcs, Y)
    assert np.allclose(solver.eigenvalues, eig)
    assert np.allclose(solver.expvar, expvar)


def test_spatial_weights_class():
    lat = np.array([-30.0, 0.0, 30.0])
    W = SpatialWeights.from_latitude(lat, n_lon=5)
    assert W.array.shape == (5, 3)
    assert np.allclose(W.array, latitude_weights(lat, 5))


def test_multivariate_eof_class():
    n_lon, n_lat, n_time = 4, 3, 40
    t = np.arange(n_time)
    d1 = np.broadcast_to(np.sin(2 * np.pi * t / 10), (n_lon, n_lat, n_time)).copy()
    d2 = 100.0 * np.broadcast_to(np.cos(2 * np.pi * t / 10), (n_lon, n_lat, n_time)).copy()
    solver = MultivariateEof(d1, d2, n_eof=2, normalize=True)
    m1, m2, pcs, expvar, eig = meof(d1, d2, n_eof=2, normalize=True)
    cm1, cm2 = solver.maps()
    assert np.allclose(cm1, m1, equal_nan=True)
    assert np.allclose(cm2, m2, equal_nan=True)
    assert np.allclose(solver.pcs(), pcs)
    assert np.allclose(solver.expvar, expvar)
    assert np.allclose(solver.eigenvalues, eig)
