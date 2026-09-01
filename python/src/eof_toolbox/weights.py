"""Latitude and area weights for gridded EOF analysis."""

from __future__ import annotations

import numpy as np


class SpatialWeights:
    """``(n_lon, n_lat)`` multiplicative weights for a lat-lon (or Cartesian) grid.

    Parameters
    ----------
    n_lon, n_lat :
        Horizontal grid size.
    lat :
        Latitude in degrees, length ``n_lat``. Ignored if ``weights`` is set.
    weights :
        User array multiplied as-is. Accepted shapes: ``(n_lat,)``,
        ``(n_lon, n_lat)``, or Fortran-order ``(n_lon * n_lat,)``.
    method :
        How to turn ``lat`` into weights: ``'sqrtcos'`` (default area
        inner product), ``'cos'`` (legacy v1), or ``'none'``.
    """

    def __init__(
        self,
        n_lon: int,
        n_lat: int,
        lat: np.ndarray | None = None,
        weights: np.ndarray | None = None,
        method: str = "sqrtcos",
    ) -> None:
        if n_lon < 1 or n_lat < 1:
            raise ValueError("n_lon and n_lat must be positive")
        self.n_lon = int(n_lon)
        self.n_lat = int(n_lat)
        self.method = method.lower()
        self._array = self._build(lat, weights)

    @classmethod
    def from_latitude(
        cls,
        lat: np.ndarray,
        n_lon: int,
        method: str = "sqrtcos",
    ) -> SpatialWeights:
        lat = np.asarray(lat, dtype=float).reshape(-1)
        return cls(n_lon, lat.size, lat=lat, method=method)

    @property
    def array(self) -> np.ndarray:
        """Weight map with shape ``(n_lon, n_lat)``."""
        return self._array

    def apply(self, data: np.ndarray) -> np.ndarray:
        """Broadcast-multiply onto ``(n_lon, n_lat[, n_time])`` data."""
        data = np.asarray(data, dtype=float)
        if data.ndim == 2:
            return data * self._array
        if data.ndim == 3:
            return data * self._array[:, :, np.newaxis]
        raise ValueError("data must be 2-D or 3-D")

    def flatten(self) -> np.ndarray:
        """Fortran-order vector of length ``n_lon * n_lat``."""
        return np.reshape(self._array, (self.n_lon * self.n_lat,), order="F")

    def _build(self, lat: np.ndarray | None, weights: np.ndarray | None) -> np.ndarray:
        if weights is not None:
            return self._from_explicit(np.asarray(weights, dtype=float))
        if lat is None:
            return np.ones((self.n_lon, self.n_lat))
        lat = np.asarray(lat, dtype=float).reshape(-1)
        if lat.size != self.n_lat:
            raise ValueError("len(lat) must equal n_lat (axis 1 of data)")
        return _latitude_map(lat, self.n_lon, self.method)

    def _from_explicit(self, w: np.ndarray) -> np.ndarray:
        n_loc = self.n_lon * self.n_lat
        if w.ndim == 1 and w.size == self.n_lat:
            return np.broadcast_to(w[np.newaxis, :], (self.n_lon, self.n_lat)).copy()
        if w.shape == (self.n_lon, self.n_lat):
            return w
        if w.size == n_loc:
            return np.reshape(w, (self.n_lon, self.n_lat), order="F")
        raise ValueError(
            "weights must be (n_lat,), (n_lon, n_lat), or (n_lon*n_lat,)"
        )


def _latitude_map(lat: np.ndarray, n_lon: int, method: str) -> np.ndarray:
    if np.any(lat > 90) or np.any(lat < -90):
        raise ValueError("Latitude must lie in [-90, 90].")
    rad = np.deg2rad(lat)
    if method == "sqrtcos":
        wlat = np.sqrt(np.maximum(np.cos(rad), 0.0))
    elif method == "cos":
        wlat = np.maximum(np.cos(rad), 0.0)
    elif method == "none":
        wlat = np.ones_like(lat)
    else:
        raise ValueError("method must be 'sqrtcos', 'cos', or 'none'")
    return np.broadcast_to(wlat[np.newaxis, :], (n_lon, lat.size)).copy()


def latitude_weights(
    lat: np.ndarray,
    n_lon: int,
    method: str = "sqrtcos",
) -> np.ndarray:
    """Return a ``(n_lon, n_lat)`` map from latitudes in degrees."""
    return SpatialWeights.from_latitude(lat, n_lon, method=method).array


def make_weights(
    n_lon: int,
    n_lat: int,
    lat: np.ndarray | None = None,
    weights: np.ndarray | None = None,
    method: str = "sqrtcos",
) -> np.ndarray:
    """Resolve latitudes / explicit weights into an ``(n_lon, n_lat)`` array."""
    return SpatialWeights(n_lon, n_lat, lat=lat, weights=weights, method=method).array
