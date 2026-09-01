"""Multivariate (combined) EOF of two gridded fields."""

from __future__ import annotations

import numpy as np

from .core import PackedEof
from .reshape import GridPacker
from .weights import SpatialWeights


class MultivariateEof:
    """Combined EOF of two ``(lon, lat, time)`` fields.

    Fields must share the time length. Horizontal grids may differ.

    ``normalize=True`` (default) divides each packed field by its own
    standard deviation before concatenation so incompatible units do
    not dominate. Spatial maps are scaled back to physical units.
    """

    def __init__(
        self,
        data1: np.ndarray,
        data2: np.ndarray,
        lat1: np.ndarray | None = None,
        lat2: np.ndarray | None = None,
        n_eof: int | None = None,
        *,
        weights1: np.ndarray | None = None,
        weights2: np.ndarray | None = None,
        weighting: str = "sqrtcos",
        normalize: bool = True,
    ) -> None:
        data1 = np.asarray(data1, dtype=float)
        data2 = np.asarray(data2, dtype=float)
        if data1.ndim != 3 or data2.ndim != 3:
            raise ValueError("data1 and data2 must be (n_lon, n_lat, n_time)")
        if data1.shape[2] != data2.shape[2]:
            raise ValueError("size(data1, 3) must equal size(data2, 3)")

        self.n_lon1, self.n_lat1, self.n_time = data1.shape
        self.n_lon2, self.n_lat2, _ = data2.shape
        self.normalize = bool(normalize)
        self.weights1 = SpatialWeights(
            self.n_lon1, self.n_lat1, lat=lat1, weights=weights1, method=weighting
        )
        self.weights2 = SpatialWeights(
            self.n_lon2, self.n_lat2, lat=lat2, weights=weights2, method=weighting
        )
        self._packer1 = GridPacker(self.n_lon1, self.n_lat1)
        self._packer2 = GridPacker(self.n_lon2, self.n_lat2)

        d1, self._loc1 = self._packer1.pack(self.weights1.apply(data1))
        d2, self._loc2 = self._packer2.pack(self.weights2.apply(data2))
        self._n1 = d1.shape[0]
        self._s1 = self._s2 = 1.0
        if self.normalize:
            self._s1 = _field_std(d1)
            self._s2 = _field_std(d2)
            if self._s1 > 0:
                d1 = d1 / self._s1
            if self._s2 > 0:
                d2 = d2 / self._s2
        packed = np.vstack([d1, d2])
        self._solver = PackedEof(packed, n_eof=n_eof)

    @property
    def n_eof(self) -> int:
        return self._solver.n_eof

    @property
    def eigenvalues(self) -> np.ndarray:
        return self._solver.eigenvalues

    @property
    def expvar(self) -> np.ndarray:
        return self._solver.expvar

    def maps(
        self,
        n_eof: int | None = None,
        *,
        unweight: bool = True,
        standardize: bool = False,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return ``(maps1, maps2)`` on each field's grid."""
        k = self.n_eof if n_eof is None else int(n_eof)
        L = self._solver.modes[:, :k].copy()
        L1 = L[: self._n1, :]
        L2 = L[self._n1 :, :]
        if self.normalize:
            L1 = L1 * self._s1
            L2 = L2 * self._s2
        if unweight:
            L1 = _unweight(L1, self.weights1, self._loc1)
            L2 = _unweight(L2, self.weights2, self._loc2)
        if standardize:
            std = np.sqrt(np.abs(self._solver.eigenvalues[:k]))
            std_safe = np.where(std > 0, std, 1.0)
            L1 = L1 * std_safe
            L2 = L2 * std_safe
        return (
            self._packer1.unpack(L1, self._loc1),
            self._packer2.unpack(L2, self._loc2),
        )

    def pcs(self, n_eof: int | None = None, *, standardize: bool = False) -> np.ndarray:
        k = self.n_eof if n_eof is None else int(n_eof)
        Y = self._solver.pcs[:k, :].copy()
        if standardize:
            std = np.sqrt(np.abs(self._solver.eigenvalues[:k]))
            std_safe = np.where(std > 0, std, 1.0)
            Y = Y / std_safe[:, np.newaxis]
        return Y


def meof(
    data1: np.ndarray,
    data2: np.ndarray,
    lat1: np.ndarray | None = None,
    lat2: np.ndarray | None = None,
    n_eof: int | None = None,
    *,
    weights1: np.ndarray | None = None,
    weights2: np.ndarray | None = None,
    weighting: str = "sqrtcos",
    unweight: bool = True,
    normalize: bool = True,
    standardize: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """One-shot wrapper around :class:`MultivariateEof`."""
    solver = MultivariateEof(
        data1,
        data2,
        lat1=lat1,
        lat2=lat2,
        n_eof=n_eof,
        weights1=weights1,
        weights2=weights2,
        weighting=weighting,
        normalize=normalize,
    )
    m1, m2 = solver.maps(unweight=unweight, standardize=standardize)
    return m1, m2, solver.pcs(standardize=standardize), solver.expvar, solver.eigenvalues


def _field_std(d: np.ndarray) -> float:
    d = d - d.mean(axis=1, keepdims=True)
    n = d.size
    if n <= 1:
        return 0.0
    return float(np.sqrt(np.sum(d * d) / (n - 1)))


def _unweight(maps: np.ndarray, weights: SpatialWeights, loc: np.ndarray) -> np.ndarray:
    w = weights.flatten()[loc].astype(float).copy()
    w[w == 0] = np.inf
    out = maps / w[:, np.newaxis]
    out[~np.isfinite(out)] = 0.0
    return out
