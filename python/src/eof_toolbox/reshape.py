"""Pack a (lon, lat, time) cube onto valid spatial locations."""

from __future__ import annotations

import numpy as np


class GridPacker:
    """Flatten a lon-lat grid in MATLAB/Fortran order and drop NaN cells.

    A location is kept only if it is finite at every time step.
    Space index = ``lon + lat * n_lon`` (0-based).
    """

    def __init__(self, n_lon: int, n_lat: int) -> None:
        self.n_lon = int(n_lon)
        self.n_lat = int(n_lat)

    def pack(self, data: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Return ``(n_valid, n_time)`` and 0-based kept indices."""
        data = np.asarray(data, dtype=float)
        if data.ndim == 2:
            data = data[:, :, np.newaxis]
        n_lon, n_lat, n_time = data.shape
        if n_lon != self.n_lon or n_lat != self.n_lat:
            raise ValueError("data horizontal size does not match packer")
        flat = np.reshape(data, (n_lon * n_lat, n_time), order="F")
        loc = np.nonzero(np.all(np.isfinite(flat), axis=1))[0]
        return flat[loc, :], loc

    def unpack(self, packed: np.ndarray, loc: np.ndarray) -> np.ndarray:
        """Scatter ``(n_valid, n_mode)`` back to ``(n_lon, n_lat, n_mode)``."""
        packed = np.asarray(packed, dtype=float)
        if packed.ndim == 1:
            packed = packed[:, np.newaxis]
        n_mode = packed.shape[1]
        out = np.full((self.n_lon * self.n_lat, n_mode), np.nan, dtype=float)
        out[loc, :] = packed
        return np.reshape(out, (self.n_lon, self.n_lat, n_mode), order="F")


def reshape_3d_to_2d(data: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    data = np.asarray(data, dtype=float)
    if data.ndim == 2:
        data = data[:, :, np.newaxis]
    return GridPacker(data.shape[0], data.shape[1]).pack(data)


def reshape_2d_to_3d(
    packed: np.ndarray,
    shape3: tuple[int, int, int],
    loc: np.ndarray,
) -> np.ndarray:
    n_lon, n_lat, n_mode = shape3
    maps = GridPacker(n_lon, n_lat).unpack(packed, loc)
    if maps.shape[2] != n_mode:
        raise ValueError("n_mode mismatch")
    return maps
