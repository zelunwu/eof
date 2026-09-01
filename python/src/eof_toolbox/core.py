"""Packed space-by-time EOF (PCA) solver."""

from __future__ import annotations

import numpy as np


class PackedEof:
    """Leading EOFs of a matrix ``X`` with shape ``(n_locations, n_timesteps)``.

    The time mean of each location is removed. The smaller of the
    space-space and time-time covariance matrices is diagonalized
    (Björnsson & Venegas, 1997). Eigenvalues are those of ``X @ X.T / T``.

    Columns of :attr:`modes` are orthonormal. The largest-magnitude
    entry of each column is forced positive so signs are reproducible.
    Reconstruction of anomalies: ``X_anom ≈ modes @ pcs``.
    """

    def __init__(self, X: np.ndarray, n_eof: int | None = None) -> None:
        X = np.array(X, dtype=float, copy=True)
        if X.ndim != 2:
            raise ValueError("X must be 2-D (n_locations, n_timesteps)")
        self.n_locations, self.n_timesteps = X.shape
        X -= X.mean(axis=1, keepdims=True)
        self._anomaly = X
        if n_eof is None:
            n_eof = min(self.n_locations, self.n_timesteps)
        n_eof = int(min(n_eof, self.n_locations, self.n_timesteps))
        if n_eof < 1:
            raise ValueError("n_eof must be >= 1")
        self.n_eof = n_eof
        self._decompose()

    def _decompose(self) -> None:
        X = self._anomaly
        n_loc, n_time = self.n_locations, self.n_timesteps
        total_var = float(np.sum(X * X) / n_time)
        if n_loc <= n_time:
            eig_values, L = _eigensolve(X @ X.T, self.n_eof)
            eig_values = eig_values / n_time
            Y = (X.T @ L).T
        else:
            eig_values, R = _eigensolve(X.T @ X, self.n_eof)
            eig_values = eig_values / n_time
            scale = 1.0 / np.sqrt(
                np.maximum(eig_values * n_time, np.finfo(float).tiny)
            )
            L = np.real(X @ (R * scale))
            Y = L.T @ X
        self._modes, self._pcs = _pin_signs(L, Y)
        self._eig = eig_values
        if total_var > 0:
            self._expvar = eig_values / total_var * 100.0
        else:
            self._expvar = np.zeros_like(eig_values)

    @property
    def modes(self) -> np.ndarray:
        """Spatial modes, shape ``(n_locations, n_eof)``."""
        return self._modes

    @property
    def pcs(self) -> np.ndarray:
        """Principal components, shape ``(n_eof, n_timesteps)``."""
        return self._pcs

    @property
    def eigenvalues(self) -> np.ndarray:
        """Covariance eigenvalues (variance per mode)."""
        return self._eig

    @property
    def expvar(self) -> np.ndarray:
        """Percent of total variance explained by each mode."""
        return self._expvar

    def reconstruct(self, n_eof: int | None = None) -> np.ndarray:
        """Reconstructed anomalies using the first ``n_eof`` modes."""
        k = self.n_eof if n_eof is None else int(n_eof)
        return self._modes[:, :k] @ self._pcs[:k, :]


def eofcore(
    X: np.ndarray, n_eof: int | None = None
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Functional wrapper around :class:`PackedEof`."""
    solver = PackedEof(X, n_eof=n_eof)
    return solver.modes, solver.pcs, solver.eigenvalues, solver.expvar


def _eigensolve(S: np.ndarray, n_eof: int) -> tuple[np.ndarray, np.ndarray]:
    S = 0.5 * (S + S.T)
    eig_values, V = np.linalg.eigh(S)
    order = np.argsort(eig_values)[::-1]
    eig_values = np.maximum(np.real(eig_values[order[:n_eof]]), 0.0)
    V = np.real(V[:, order[:n_eof]])
    return eig_values, V


def _pin_signs(L: np.ndarray, Y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    L = np.array(L, copy=True)
    Y = np.array(Y, copy=True)
    for k in range(L.shape[1]):
        idx = int(np.argmax(np.abs(L[:, k])))
        if L[idx, k] < 0:
            L[:, k] *= -1
            Y[k, :] *= -1
    return L, Y
