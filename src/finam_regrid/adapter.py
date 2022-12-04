"""ESMF regridding adapters."""
import ESMF
import finam as fm
import numpy as np

from .tools import create_transformer, to_esmf


class Regrid(fm.adapters.regrid.ARegridding):
    """
    FINAM adapter for regridding using `ESMPy <https://earthsystemmodeling.org/esmpy/>`_.

    Supports all of ESMPy's  :class:`regridding methods <.RegridMethod>`.
    For parameters passed as ``**regrid_args``, see the ESMPy class
    `Regrid <https://earthsystemmodeling.org/esmpy_doc/release/latest/html/regrid.html>`_

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

    def __init__(self, in_grid=None, out_grid=None, **regrid_args):
        super().__init__(in_grid, out_grid)
        self.regrid_args = regrid_args
        self.regrid = None
        self.in_grid = None
        self.out_grid = None
        self.in_field = None
        self.out_field = None

        if "unmapped_action" not in self.regrid_args:
            self.regrid_args["unmapped_action"] = ESMF.UnmappedAction.IGNORE

    def _update_grid_specs(self):
        transformer = create_transformer(self.input_grid.crs, self.output_grid.crs)

        self.in_grid, self.in_field = to_esmf(self.input_grid, transformer)
        self.out_grid, self.out_field = to_esmf(self.output_grid)

        self.regrid = ESMF.Regrid(
            self.in_field,
            self.out_field,
            **self.regrid_args,
        )

    def _get_data(self, time, target):
        in_data = self.pull_data(time, target)

        self.in_field.data[:] = fm.data.get_magnitude(fm.data.strip_time(in_data))[:]
        self.out_field.data[:] = np.nan

        self.regrid(self.in_field, self.out_field)

        return self.out_field.data * fm.data.get_units(in_data)

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
