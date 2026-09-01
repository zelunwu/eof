"""Empirical Orthogonal Function toolbox (space-by-time PCA)."""

from .core import PackedEof, eofcore
from .eof import Eof, eof
from .meof import MultivariateEof, meof
from .reshape import GridPacker, reshape_2d_to_3d, reshape_3d_to_2d
from .weights import SpatialWeights, latitude_weights, make_weights

__all__ = [
    "Eof",
    "GridPacker",
    "MultivariateEof",
    "PackedEof",
    "SpatialWeights",
    "eof",
    "eofcore",
    "latitude_weights",
    "make_weights",
    "meof",
    "reshape_2d_to_3d",
    "reshape_3d_to_2d",
]
__version__ = "2.0.0"
