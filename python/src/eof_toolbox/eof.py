"""Grid-level EOF solver."""

from __future__ import annotations

from typing import Any

import numpy as np

from .core import PackedEof
from .reshape import GridPacker
from .weights import SpatialWeights, latitude_weights, make_weights


class Eof:
    """EOF analysis of a ``(n_lon, n_lat, n_time)`` field.

    The eigenproblem is solved at construction. Retrieve maps and PCs
    with :meth:`maps` / :meth:`pcs`, or reconstruct anomalies with
    :meth:`reconstruct`.

    Parameters
    ----------
    data :
        Spatiotemporal array. NaNs mark unused cells and must be
        constant in time.
    lat :
        Latitude in degrees, length ``n_lat``. Builds area weights
        ``sqrt(cos(lat))`` unless ``weights`` is given.
    n_eof :
        Number of modes. Default is ``min(n_valid_space, n_time)``.
    weights :
        Multiplied onto ``data`` before the eigenproblem.
    weighting :
        ``'sqrtcos'``, ``'cos'``, or ``'none'`` when ``lat`` is used.
    """

    def __init__(
        self,
        data: np.ndarray,
        lat: np.ndarray | None = None,
        n_eof: int | None = None,
        *,
        weights: np.ndarray | None = None,
        weighting: str = "sqrtcos",
    ) -> None:
        data = np.asarray(data, dtype=float)
        if data.ndim != 3:
            raise ValueError("data must be (n_lon, n_lat, n_time)")
        self.n_lon, self.n_lat, self.n_time = data.shape
        self.spatial_weights = SpatialWeights(
            self.n_lon, self.n_lat, lat=lat, weights=weights, method=weighting
        )
        self._packer = GridPacker(self.n_lon, self.n_lat)
        packed, self._loc = self._packer.pack(self.spatial_weights.apply(data))
        if packed.shape[0] == 0:
            raise ValueError("no finite spatial locations")
        self._solver = PackedEof(packed, n_eof=n_eof)

    @property
    def n_eof(self) -> int:
        return self._solver.n_eof

    @property
    def eigenvalues(self) -> np.ndarray:
        return self._solver.eigenvalues

    @property
    def expvar(self) -> np.ndarray:
        """Percent of total *weighted* variance."""
        return self._solver.expvar

    def maps(self, n_eof: int | None = None, *, unweight: bool = True,
             standardize: bool = False) -> np.ndarray:
        """Spatial patterns, shape ``(n_lon, n_lat, n_eof)``."""
        L = self._scaled_modes(n_eof, unweight=unweight, standardize=standardize)
        k = L.shape[1]
        return self._packer.unpack(L, self._loc)[:, :, :k]

    def pcs(self, n_eof: int | None = None, *, standardize: bool = False) -> np.ndarray:
        """Principal components, shape ``(n_eof, n_time)``."""
        k = self.n_eof if n_eof is None else int(n_eof)
        Y = self._solver.pcs[:k, :].copy()
        if standardize:
            std = np.sqrt(np.abs(self._solver.eigenvalues[:k]))
            std_safe = np.where(std > 0, std, 1.0)
            Y = Y / std_safe[:, np.newaxis]
        return Y

    def reconstruct(self, n_eof: int | None = None) -> np.ndarray:
        """Reconstructed anomalies in physical units, shape of the input."""
        maps = self.maps(n_eof, unweight=True, standardize=False)
        pcs = self.pcs(n_eof, standardize=False)
        weighted_maps = maps * self.spatial_weights.array[:, :, np.newaxis]
        return np.einsum("ijm,mt->ijt", weighted_maps, pcs)

    def _scaled_modes(
        self, n_eof: int | None, *, unweight: bool, standardize: bool
    ) -> np.ndarray:
        k = self.n_eof if n_eof is None else int(n_eof)
        L = self._solver.modes[:, :k].copy()
        if unweight:
            w = self.spatial_weights.flatten()[self._loc].astype(float).copy()
            w[w == 0] = np.inf
            L = L / w[:, np.newaxis]
            L[~np.isfinite(L)] = 0.0
        if standardize:
            std = np.sqrt(np.abs(self._solver.eigenvalues[:k]))
            std_safe = np.where(std > 0, std, 1.0)
            L = L * std_safe
        return L


def eof(
    data: np.ndarray,
    lat: np.ndarray | None = None,
    n_eof: int | None = None,
    *,
    weights: np.ndarray | None = None,
    weighting: str = "sqrtcos",
    unweight: bool = True,
    standardize: bool = False,
    **_: Any,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """One-shot wrapper around :class:`Eof` (MATLAB-compatible returns)."""
    solver = Eof(
        data, lat=lat, n_eof=n_eof, weights=weights, weighting=weighting
    )
    maps = solver.maps(unweight=unweight, standardize=standardize)
    pcs = solver.pcs(standardize=standardize)
    return maps, pcs, solver.expvar, solver.eigenvalues


__all__ = ["Eof", "eof", "latitude_weights", "make_weights", "SpatialWeights"]
