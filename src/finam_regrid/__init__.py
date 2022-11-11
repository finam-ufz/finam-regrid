"""ESMF/ESMPy regridding components for FINAM"""
from ESMF.api.constants import ExtrapMethod, RegridMethod, UnmappedAction

from .adapter import Regrid

try:
    from ._version import __version__
except ModuleNotFoundError:  # pragma: no cover
    # package is not installed
    __version__ = "0.0.0.dev0"


__all__ = ["Regrid"]
__all__ += ["ExtrapMethod", "RegridMethod", "UnmappedAction"]
