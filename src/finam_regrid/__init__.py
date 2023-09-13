"""
`ESMPy <https://earthsystemmodeling.org/esmpy/>`_
regridding adapter for the `FINAM <https://finam.pages.ufz.de/>`_ model coupling framework.

.. toctree::
   :hidden:

   self

Adapter
=======

.. autosummary::
   :toctree: generated
   :caption: Adapter

    Regrid

Constants
=========

.. autosummary::
   :toctree: generated
   :caption: Constants

    ExtrapMethod
    RegridMethod
    UnmappedAction
"""
from esmpy.api.constants import ExtrapMethod, RegridMethod, UnmappedAction

from .adapter import Regrid

try:
    from ._version import __version__
except ModuleNotFoundError:  # pragma: no cover
    # package is not installed
    __version__ = "0.0.0.dev0"


__all__ = ["Regrid"]
__all__ += ["ExtrapMethod", "RegridMethod", "UnmappedAction"]
