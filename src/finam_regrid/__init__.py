"""ESMF/ESMPy regridding components for FINAM"""

from .adapter import Regrid
from .tools import to_esmf

try:
    from ._version import __version__
except ModuleNotFoundError:  # pragma: no cover
    # package is not installed
    __version__ = "0.0.0.dev0"


__all__ = ["to_esmf"]
__all__ = ["Regrid"]
