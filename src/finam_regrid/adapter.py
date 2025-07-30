"""ESMF regridding adapters."""

import esmpy
import finam as fm
import numpy as np
from finam.tools.log_helper import ErrorLogger

from .tools import RegridCRS, create_transformer, is_latlon, to_esmf


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
    regrid_crs : RegridCRS or None, optional
        specify which CRS should be used in the regridder.
        Options: (i) RegridCRS.SRC (source grid, default), (ii) RegridCRS.DST (target grid),
        (iii) RegridCRS.SPH (covert both grids to WGS84 and assume it if not present) and
        (iv) a valid CRS specifier for pyproj.
        Using RegridCRS.SPH will use spherical coordinates in the ESMF regridder.
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

    def __init__(
        self,
        in_grid=None,
        out_grid=None,
        zero_region=None,
        regrid_crs=None,
        **regrid_args,
    ):
        super().__init__(in_grid, out_grid)
        self.regrid_args = regrid_args
        self.regrid = None
        self.in_grid = None
        self.out_grid = None
        self.in_field = None
        self.out_field = None
        self.zero_region = zero_region
        self.output_mask = fm.Mask.FLEX
        self.regrid_crs = regrid_crs if regrid_crs is not None else RegridCRS.SRC
        if "unmapped_action" not in self.regrid_args:
            self.regrid_args["unmapped_action"] = esmpy.UnmappedAction.IGNORE

    def _update_grid_specs(self):
        assume_target_crs = False
        if self.regrid_crs == RegridCRS.SRC:
            target_crs = self.input_grid.crs
        elif self.regrid_crs == RegridCRS.DST:
            target_crs = self.output_grid.crs
        elif self.regrid_crs == RegridCRS.SPH:
            target_crs = "WGS84"
            # for spherical we just assume missing crs info as lat-lon
            assume_target_crs = True
        else:
            target_crs = self.regrid_crs
        sph = is_latlon(target_crs) if target_crs is not None else False
        src_transformer = create_transformer(
            self.input_grid.crs, target_crs, assume_target_crs
        )
        dst_transformer = create_transformer(
            self.output_grid.crs, target_crs, assume_target_crs
        )
        self.in_grid, self.in_field = to_esmf(self.input_grid, src_transformer, sph)
        self.out_grid, self.out_field = to_esmf(self.output_grid, dst_transformer, sph)
        self.regrid = esmpy.Regrid(
            self.in_field,
            self.out_field,
            **self.regrid_args,
        )

    def _get_data(self, time, target):
        in_data = self.pull_data(time, target)

        if fm.data.has_masked_values(in_data):
            with ErrorLogger(self.logger):
                msg = "Regridding is currently not implemented for masked data"
                raise NotImplementedError(msg)

        self.in_field.data[...] = self.input_grid.to_canonical(
            fm.data.strip_time(in_data, self.input_grid).magnitude
        )
        self.out_field.data[...] = np.nan

        self.regrid(self.in_field, self.out_field, zero_region=self.zero_region)

        return self.output_grid.from_canonical(self.out_field.data.copy())

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
