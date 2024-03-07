import unittest
from datetime import datetime, timedelta

import finam as fm
import numpy as np

from finam_regrid import Regrid, RegridMethod


class TestAdapter(unittest.TestCase):
    def setup_run(self, in_grid, out_grid, regrid_method, masked=False):
        time = datetime(2000, 1, 1)
        in_info = fm.Info(
            time=time,
            grid=in_grid,
            units="m",
        )

        in_data = np.zeros(shape=in_info.grid.data_shape, order=in_info.grid.order)
        if len(in_data.shape) == 1:
            in_data.data[0] = 1.0
        else:
            in_data.data[0, 0] = 1.0

        if masked:
            in_data = np.ma.masked_where(in_data > 0, in_data)

        self.source = fm.modules.CallbackGenerator(
            callbacks={
                "Output": (
                    lambda t: in_data.copy(),
                    in_info,
                )
            },
            start=datetime(2000, 1, 1),
            step=timedelta(days=1),
        )

        self.sink = fm.modules.DebugConsumer(
            {"Input": fm.Info(None, grid=out_grid, units=None)},
            start=datetime(2000, 1, 1),
            step=timedelta(days=1),
        )

        self.composition = fm.Composition([self.source, self.sink])
        self.composition.initialize()

        (
            self.source.outputs["Output"]
            >> Regrid(regrid_method=regrid_method)
            >> self.sink.inputs["Input"]
        )

    def test_adapter_grid_nearest(self):
        self.setup_run(
            regrid_method=RegridMethod.NEAREST_STOD,
            in_grid=fm.UniformGrid(
                dims=(3, 7),
                spacing=(3.0, 3.0, 3.0),
                data_location=fm.Location.POINTS,
            ),
            out_grid=fm.UniformGrid(dims=(9, 19), data_location=fm.Location.POINTS),
        )
        self.composition.run(end_time=datetime(2000, 1, 5))

        result = self.sink.data["Input"]
        self.assertEqual(result[0, 0, 0], 1.0 * fm.UNITS.meter)
        self.assertEqual(result[0, 0, 1], 1.0 * fm.UNITS.meter)
        self.assertEqual(result[0, 1, 0], 1.0 * fm.UNITS.meter)
        self.assertEqual(result[0, 1, 1], 1.0 * fm.UNITS.meter)
        self.assertEqual(result[0, 0, 2], 0.0 * fm.UNITS.meter)
        self.assertEqual(result[0, 2, 0], 0.0 * fm.UNITS.meter)

    def test_adapter_grid_nearest_masked(self):
        self.setup_run(
            regrid_method=RegridMethod.NEAREST_STOD,
            in_grid=fm.UniformGrid(
                dims=(3, 7),
                spacing=(3.0, 3.0, 3.0),
                data_location=fm.Location.POINTS,
            ),
            out_grid=fm.UniformGrid(dims=(9, 19), data_location=fm.Location.POINTS),
            masked=True,
        )
        with self.assertRaises(NotImplementedError):
            self.composition.run(end_time=datetime(2000, 1, 5))

    def test_adapter_grid_linear(self):
        self.setup_run(
            regrid_method=RegridMethod.BILINEAR,
            in_grid=fm.UniformGrid(
                dims=(5, 10),
                spacing=(2.0, 2.0, 2.0),
                data_location=fm.Location.POINTS,
            ),
            out_grid=fm.UniformGrid(dims=(9, 19), data_location=fm.Location.POINTS),
        )
        self.composition.run(end_time=datetime(2000, 1, 5))

        result = self.sink.data["Input"]
        self.assertEqual(result[0, 0, 0], 1.0 * fm.UNITS.meter)
        self.assertEqual(result[0, 0, 1], 0.5 * fm.UNITS.meter)
        self.assertEqual(result[0, 1, 0], 0.5 * fm.UNITS.meter)
        self.assertEqual(result[0, 1, 1], 0.25 * fm.UNITS.meter)

    def test_adapter_grid_conserve(self):
        self.setup_run(
            regrid_method=RegridMethod.CONSERVE,
            in_grid=fm.UniformGrid(
                dims=(5, 10),
                spacing=(2.0, 2.0, 2.0),
            ),
            out_grid=fm.UniformGrid(dims=(9, 19)),
        )
        self.composition.run(end_time=datetime(2000, 1, 5))

        result = self.sink.data["Input"]
        self.assertEqual(result[0, 0, 0], 1.0 * fm.UNITS.meter)
        self.assertEqual(result[0, 0, 1], 1.0 * fm.UNITS.meter)
        self.assertEqual(result[0, 1, 0], 1.0 * fm.UNITS.meter)
        self.assertEqual(result[0, 1, 1], 1.0 * fm.UNITS.meter)

    def test_adapter_grid_conserve_2nd(self):
        self.setup_run(
            regrid_method=RegridMethod.CONSERVE_2ND,
            in_grid=fm.UniformGrid(
                dims=(5, 10),
                spacing=(2.0, 2.0, 2.0),
            ),
            out_grid=fm.UniformGrid(dims=(9, 19)),
        )
        self.composition.run(end_time=datetime(2000, 1, 5))

        result = self.sink.data["Input"]
        self.assertEqual(result[0, 0, 0], 1.0 * fm.UNITS.meter)
        self.assertEqual(result[0, 0, 1], 1.0 * fm.UNITS.meter)
        self.assertEqual(result[0, 1, 0], 1.0 * fm.UNITS.meter)
        self.assertEqual(result[0, 1, 1], 1.0 * fm.UNITS.meter)

    def test_adapter_grid_crs(self):
        out_grid = fm.UniformGrid(
            dims=(9, 19), data_location=fm.Location.POINTS, crs="EPSG:25832"
        )

        self.setup_run(
            regrid_method=RegridMethod.BILINEAR,
            in_grid=fm.UniformGrid(
                dims=(5, 10),
                spacing=(2.0, 2.0, 2.0),
                data_location=fm.Location.POINTS,
                crs="EPSG:32632",
            ),
            out_grid=out_grid,
        )
        self.composition.run(end_time=datetime(2000, 1, 2))

        self.assertEqual(self.sink.inputs["Input"].info.grid, out_grid)
        self.assertAlmostEqual(
            fm.data.get_magnitude(self.sink.data["Input"])[0, 0, 0], 1.0
        )
        self.assertAlmostEqual(
            fm.data.get_magnitude(self.sink.data["Input"])[0, 0, 1], 0.5
        )
        self.assertAlmostEqual(
            fm.data.get_magnitude(self.sink.data["Input"])[0, 1, 0], 0.5
        )
        self.assertAlmostEqual(
            fm.data.get_magnitude(self.sink.data["Input"])[0, 1, 1], 0.25
        )

    def test_adapter_mesh_nearest(self):
        self.setup_run(
            regrid_method=RegridMethod.NEAREST_STOD,
            in_grid=_create_mesh(),
            out_grid=_create_mesh(),
        )
        self.composition.run(end_time=datetime(2000, 1, 5))
        result = self.sink.data["Input"]

    def test_adapter_mesh_linear(self):
        self.setup_run(
            regrid_method=RegridMethod.BILINEAR,
            in_grid=_create_mesh(),
            out_grid=_create_mesh(),
        )
        self.composition.run(end_time=datetime(2000, 1, 5))
        result = self.sink.data["Input"]

    def test_adapter_mesh_conserve(self):
        self.setup_run(
            regrid_method=RegridMethod.CONSERVE,
            in_grid=_create_mesh(),
            out_grid=_create_mesh(),
        )
        self.composition.run(end_time=datetime(2000, 1, 5))
        result = self.sink.data["Input"]

    def test_adapter_mesh_conserve_2nd(self):
        self.setup_run(
            regrid_method=RegridMethod.CONSERVE_2ND,
            in_grid=_create_mesh(),
            out_grid=_create_mesh(),
        )
        self.composition.run(end_time=datetime(2000, 1, 5))
        result = self.sink.data["Input"]


def _create_mesh():
    grid = fm.UniformGrid((16, 13))
    points = grid.points
    cells = grid.cells
    types = grid.cell_types
    grid = fm.UnstructuredGrid(points, cells, types, data_location=fm.Location.CELLS)

    return grid


if __name__ == "__main__":
    unittest.main()
