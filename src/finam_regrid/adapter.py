"""ESMF regridding adapters."""

import esmpy
import finam as fm
import numpy as np

from .tools import RegridCRS, create_transformer, is_latlon, to_esmf


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

        import finam_regrid as fm_rg

        adapter = fm_rg.Regrid()

    Using a specific regridding method:

    .. testcode:: constructor

        adapter = fm_rg.Regrid(
            regrid_method=fm_rg.RegridMethod.CONSERVE_2ND,
        )

    Using a specific regridding method and extrapolation:

    .. testcode:: constructor

        adapter = fm_rg.Regrid(
            regrid_method=fm_rg.RegridMethod.CONSERVE_2ND,
            extrap_method=fm_rg.ExtrapMethod.NEAREST_IDAVG,
        )

    Parameters
    ----------

    in_grid : finam.Grid, optional
        Input grid specification. Will be retrieved from upstream component if not specified.
    out_grid : finam.Grid, optional
        Output grid specification. Will be retrieved from downstream component if not specified.
    out_mask : :any:`finam.Mask` value or valid boolean mask for :any:`MaskedArray` or None, optional
        masking specification of the regridding output. Options:
            * :any:`finam.Mask.FLEX`: data will be unmasked
            * :any:`finam.Mask.NONE`: data will be unmasked and given as plain numpy array
            * valid boolean mask for MaskedArray
            * None: will be determined by connected target
    regrid_crs : :class:`.RegridCRS`, crs or None, optional
        specify which CRS should be used in the regridder. Options:
            * :attr:`.RegridCRS.SRC` source grid (default),
            * :attr:`.RegridCRS.DST` target grid,
            * :attr:`.RegridCRS.SPH` covert both grids to WGS84
            * a valid CRS specifier for pyproj
    zero_region : :class:`.Region` or None, optional
        specify which region of the field indices will be zeroed out before
        adding the values resulting from the interpolation. If None, defaults to Region.TOTAL.
    **regrid_args : Any
        Keyword argument passed to the ESMPy class
        `Regrid <https://earthsystemmodeling.org/esmpy_doc/release/latest/html/regrid.html>`_.

        Important keyword arguments are documented in the "Other Parameters" section.

    Other Parameters
    ----------------

    regrid_method : :class:`.RegridMethod`
        Regridding method. See :class:`.RegridMethod`. Defaults to :attr:`.RegridMethod.BILINEAR`.
    line_type : :class:`.LineType`
        select the path of the line that connects two points on the surface of a sphere.
        This in turn controls the path along which distances are calculated
        and the shape of the edges that make up a cell.
        If ``None``, defaults to: :attr:`.LineType.GREAT_CIRCLE` for regridmethod == :attr:`.RegridMethod.CONSERVE`,
        or :attr:`.LineType.CART` for regridmethod != :attr:`.RegridMethod.CONSERVE`.
    unmapped_action : :class:`.UnmappedAction`
        Action on unmapped cells. See :class:`.UnmappedAction`. Defaults to :attr:`.UnmappedAction.IGNORE`.
    extrap_method : :class:`.ExtrapMethod`
        Extrapolation method. See :class:`.ExtrapMethod`. Defaults to ``None``.
    extrap_num_src_pnts: int
        The number of source points to use for the extrapolation methods that use more than one source point
        (e.g. :attr:`.ExtrapMethod.NEAREST_IDAVG`). If not specified, defaults to 8.
    extrap_dist_exponent: float
        The exponent to raise the distance to when  calculating weights for the :attr:`.ExtrapMethod.NEAREST_IDAVG`
        extrapolation method. A higher value reduces the influence of more distant points.
        If not specified, defaults to ``2.0``.
    extrap_num_levels: int
        The number of levels to output for the extrapolation methods
        that fill levels (e.g. :attr:`.ExtrapMethod.CREEP`).
        When a method is used that requires this, then an error will be returned if it is not specified.
    pole_method : :class:`.PoleMethod`
        specifies which type of artificial pole to construct on the source Grid for regridding.
        If ``None``, defaults to: :attr:`.PoleMethod.NONE` for regridmethod == :attr:`.RegridMethod.CONSERVE`, or
        :attr:`.PoleMethod.ALLAVG` for regridmethod != :attr:`.RegridMethod.CONSERVE`.
    regrid_pole_npoints: int
        specifies how many points to average over if polemethod == :attr:`.PoleMethod.ALLAVG`.
    ignore_degenerate: bool
        Ignore degenerate cells when checking the input Grids or Meshes for errors.
        If this is set to True, then the regridding proceeds, but degenerate cells will be skipped.
        If set to False, a degenerate cell produces an error.
        This currently only applies to :attr:`.RegridMethod.CONSERVE`,
        other regrid methods currently always skip degenerate cells.
        If ``None``, defaults to ``False``.
    """

    def __init__(
        self,
        in_grid=None,
        out_grid=None,
        out_mask=None,
        regrid_crs=None,
        zero_region=None,
        **regrid_args,
    ):
        super().__init__(in_grid, out_grid, out_mask)
        self.regrid_args = regrid_args
        self.regrid = None
        self.in_grid = None
        self.out_grid = None
        self.in_field = None
        self.out_field = None
        self.zero_region = zero_region
        # self.output_mask = fm.Mask.FLEX
        self.regrid_crs = regrid_crs if regrid_crs is not None else RegridCRS.SRC
        if "unmapped_action" not in self.regrid_args:
            self.regrid_args["unmapped_action"] = esmpy.UnmappedAction.IGNORE

    def _update_grid_specs(self):
        # determine masks for in and output
        self._check_and_set_out_mask()
        src_mask = None
        dst_mask = None
        if self._need_mask(self.input_mask):
            src_mask = self.input_grid.to_canonical(self.input_mask).astype(np.int32)
            self.regrid_args["src_mask_values"] = [1]
        if self._need_mask(self.output_mask):
            dst_mask = self.output_grid.to_canonical(self.output_mask).astype(np.int32)
            self.regrid_args["dst_mask_values"] = [1]
        # determine regrid crs
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
        # create transformer
        src_transformer = create_transformer(
            self.input_grid.crs, target_crs, assume_target_crs
        )
        dst_transformer = create_transformer(
            self.output_grid.crs, target_crs, assume_target_crs
        )
        # create grids and regridder
        sph = is_latlon(target_crs) if target_crs is not None else False
        self.in_grid, self.in_field = to_esmf(
            self.input_grid, src_transformer, sph, src_mask
        )
        self.out_grid, self.out_field = to_esmf(
            self.output_grid, dst_transformer, sph, dst_mask
        )
        self.regrid = esmpy.Regrid(
            self.in_field,
            self.out_field,
            **self.regrid_args,
        )

    def _get_data(self, time, target):
        in_data = self.pull_data(time, target)
        self._check_in_data(in_data)

        self.in_field.data[...] = self.input_grid.to_canonical(
            fm.data.strip_time(in_data, self.input_grid).magnitude
        )
        self.out_field.data[...] = np.nan

        self.regrid(self.in_field, self.out_field, zero_region=self.zero_region)

        data = self.output_grid.from_canonical(self.out_field.data.copy())
        if fm.data.tools.mask_specified(self.output_mask):
            return fm.data.tools.to_masked(data, mask=self.output_mask)
        return data

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
