"""ESMF regridding adapters."""
import ESMF
import finam as fm
import numpy as np

from .tools import create_transformer, to_esmf


class Regrid(fm.adapters.regrid.ARegridding):
    """
    Regrid data between two grid specifications using ESMF.
    """

    def __init__(self, method=None, in_grid=None, out_grid=None):
        super().__init__(in_grid, out_grid)
        self.method = method
        self.regrid = None
        self.in_field = None
        self.out_field = None

    def _update_grid_specs(self):
        transformer = create_transformer(self.input_grid.crs, self.output_grid.crs)

        _g1, self.in_field = to_esmf(self.input_grid, transformer)
        _g2, self.out_field = to_esmf(self.output_grid)

        self.regrid = ESMF.Regrid(
            self.in_field,
            self.out_field,
            regrid_method=self.method,
            unmapped_action=ESMF.UnmappedAction.IGNORE,
        )

    def _get_data(self, time, target):
        in_data = self.pull_data(time, target)

        self.in_field.data[:] = fm.data.get_magnitude(fm.data.strip_time(in_data))[:]
        self.out_field.data[:] = np.nan

        self.regrid(self.in_field, self.out_field)

        return self.out_field.data * fm.data.get_units(in_data)
