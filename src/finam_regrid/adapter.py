"""ESMF regridding adapters."""

import esmpy
import finam as fm
import numpy as np
from finam.tools.log_helper import ErrorLogger

from .tools import create_transformer, to_esmf


class Regrid(fm.adapters.regrid.ARegridding):
    """
    FINAM adapter for regridding using `ESMPy <https://earthsystemmodeling.org/esmpy/>`_.

    Supports all of ESMPy's  :class:`regridding methods <.RegridMethod>`.
    For parameters passed as ``**regrid_args``, see the ESMPy class
    `Regrid <https://earthsystemmodeling.org/esmpy_doc/release/latest/html/regrid.html>`_


    .. warning::
        Does currently not support masked input data. Raises a ``NotImplementedError`` in that case.


    Examples
    --------

    Simple usage with defaults and grid specifications from connected components:

    .. testcode:: constructor

        import finam_regrid as fmr

        adapter = fmr.Regrid()

    Using a specific regridding method:

    .. testcode:: constructor

        adapter = fmr.Regrid(
            regrid_method=fmr.RegridMethod.CONSERVE_2ND,
        )

    Using a specific regridding method and extrapolation:

    .. testcode:: constructor

        adapter = fmr.Regrid(
            regrid_method=fmr.RegridMethod.CONSERVE_2ND,
            extrap_method=fmr.ExtrapMethod.NEAREST_IDAVG,
        )

    Parameters
    ----------

    in_grid : finam.Grid, optional
        Input grid specification. Will be retrieved from upstream component if not specified.
    out_grid : finam.Grid, optional
        Output grid specification. Will be retrieved from downstream component if not specified.
    zero_region : Region or None, optional
        specify which region of the field indices will be zeroed out before
        adding the values resulting from the interpolation. If None, defaults to Region.TOTAL.
    **regrid_args : Any
        Keyword argument passed to the ESMPy class
        `Regrid <https://earthsystemmodeling.org/esmpy_doc/release/latest/html/regrid.html>`_.

        **Important keyword arguments are:**

    regrid_method : RegridMethod
        Regridding method. See :class:`.RegridMethod`. Defaults to :attr:`.RegridMethod.BILINEAR`.
    extrap_method : ExtrapMethod
        Extrapolation method. See :class:`.ExtrapMethod`. Defaults to ``None``.
    unmapped_action : UnmappedAction
        Action on unmapped cells. See :class:`.UnmappedAction`. Defaults to :attr:`.UnmappedAction.IGNORE`.
    """

    def __init__(self, in_grid=None, out_grid=None, zero_region=None, **regrid_args):
        super().__init__(in_grid, out_grid)
        self.regrid_args = regrid_args
        self.regrid = None
        self.in_grid = None
        self.out_grid = None
        self.in_field = None
        self.out_field = None
        self.zero_region = zero_region

        if "unmapped_action" not in self.regrid_args:
            self.regrid_args["unmapped_action"] = esmpy.UnmappedAction.IGNORE

    def _update_grid_specs(self):
        transformer = create_transformer(self.input_grid.crs, self.output_grid.crs)

        self.in_grid, self.in_field = to_esmf(self.input_grid, transformer)
        self.out_grid, self.out_field = to_esmf(self.output_grid)

        self.regrid = esmpy.Regrid(
            self.in_field,
            self.out_field,
            **self.regrid_args,
        )

    def _get_data(self, time, target):
        in_data = self.pull_data(time, target)

        if fm.data.is_masked_array(in_data):
            with ErrorLogger(self.logger):
                msg = "Regridding is currently not implemented for masked data"
                raise NotImplementedError(msg)

        self.in_field.data[...] = self.input_grid.to_canonical(
            fm.data.strip_time(in_data, self.input_grid).magnitude
        )
        self.out_field.data[...] = np.nan

        self.regrid(self.in_field, self.out_field, zero_region=self.zero_region)

        return self.output_grid.from_canonical(self.out_field.data) * fm.data.get_units(
            in_data
        )

    def _finalize(self):
        self.regrid.destroy()
        self.in_field.destroy()
        self.out_field.destroy()
        self.in_grid.destroy()
        self.out_grid.destroy()

        self.regrid = None
        self.in_field = None
        self.out_field = None
        self.in_grid = None
        self.out_grid = None
