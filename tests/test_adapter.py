import unittest
from datetime import datetime, timedelta

import finam as fm
import numpy as np

from finam_regrid import Regrid, RegridMethod


class TestAdapter(unittest.TestCase):
    def test_adapter(self):
        time = datetime(2000, 1, 1)
        in_info = fm.Info(
            time=time,
            grid=fm.UniformGrid(
                dims=(5, 10),
                spacing=(2.0, 2.0, 2.0),
                data_location=fm.Location.POINTS,
            ),
            units="m",
        )
        out_spec = fm.UniformGrid(dims=(9, 19), data_location=fm.Location.POINTS)

        in_data = np.zeros(shape=in_info.grid.data_shape, order=in_info.grid.order)
        in_data.data[0, 0] = 1.0

        source = fm.modules.CallbackGenerator(
            callbacks={
                "Output": (
                    lambda t: in_data,
                    in_info,
                )
            },
            start=datetime(2000, 1, 1),
            step=timedelta(days=1),
        )

        sink = fm.modules.DebugConsumer(
            {"Input": fm.Info(None, grid=out_spec)},
            start=datetime(2000, 1, 1),
            step=timedelta(days=1),
        )

        composition = fm.Composition([source, sink])
        composition.initialize()

        (source.outputs["Output"] >> Regrid() >> sink.inputs["Input"])

        composition.connect()

        composition.run(t_max=datetime(2000, 1, 5))

        result = sink.data["Input"]
        self.assertEqual(result[0, 0, 0], 1.0 * fm.UNITS.meter)
        self.assertEqual(result[0, 0, 1], 0.5 * fm.UNITS.meter)
        self.assertEqual(result[0, 1, 0], 0.5 * fm.UNITS.meter)
        self.assertEqual(result[0, 1, 1], 0.25 * fm.UNITS.meter)

    def test_adapter_conserve(self):
        time = datetime(2000, 1, 1)
        in_info = fm.Info(
            time=time,
            grid=fm.UniformGrid(
                dims=(5, 10),
                spacing=(2.0, 2.0, 2.0),
            ),
            units="m",
        )
        out_spec = fm.UniformGrid(dims=(9, 19))

        in_data = np.zeros(shape=in_info.grid.data_shape, order=in_info.grid.order)
        in_data.data[0, 0] = 1.0

        source = fm.modules.CallbackGenerator(
            callbacks={
                "Output": (
                    lambda t: in_data,
                    in_info,
                )
            },
            start=datetime(2000, 1, 1),
            step=timedelta(days=1),
        )

        sink = fm.modules.DebugConsumer(
            {"Input": fm.Info(None, grid=out_spec)},
            start=datetime(2000, 1, 1),
            step=timedelta(days=1),
        )

        composition = fm.Composition([source, sink])
        composition.initialize()

        (
            source.outputs["Output"]
            >> Regrid(regrid_method=RegridMethod.CONSERVE)
            >> sink.inputs["Input"]
        )

        composition.connect()

        composition.run(t_max=datetime(2000, 1, 5))

        result = sink.data["Input"]
        self.assertEqual(result[0, 0, 0], 1.0 * fm.UNITS.meter)
        self.assertEqual(result[0, 0, 1], 1.0 * fm.UNITS.meter)
        self.assertEqual(result[0, 1, 0], 1.0 * fm.UNITS.meter)
        self.assertEqual(result[0, 1, 1], 1.0 * fm.UNITS.meter)

    def test_adapter_crs(self):
        time = datetime(2000, 1, 1)
        in_info = fm.Info(
            time=time,
            grid=fm.UniformGrid(
                dims=(5, 10),
                spacing=(2.0, 2.0, 2.0),
                data_location=fm.Location.POINTS,
                crs="EPSG:32632",
            ),
            units="m",
        )
        out_spec = fm.UniformGrid(
            dims=(9, 19), data_location=fm.Location.POINTS, crs="EPSG:25832"
        )

        in_data = np.zeros(shape=in_info.grid.data_shape, order=in_info.grid.order)
        in_data.data[0, 0] = 1.0

        source = fm.modules.CallbackGenerator(
            callbacks={"Output": (lambda t: in_data, in_info)},
            start=datetime(2000, 1, 1),
            step=timedelta(days=1),
        )

        sink = fm.modules.DebugConsumer(
            {"Input": fm.Info(None, grid=out_spec)},
            start=datetime(2000, 1, 1),
            step=timedelta(days=1),
        )

        composition = fm.Composition([source, sink])
        composition.initialize()

        (source.outputs["Output"] >> Regrid() >> sink.inputs["Input"])

        composition.run(t_max=datetime(2000, 1, 2))

        self.assertEqual(sink.inputs["Input"].info.grid, out_spec)
        self.assertAlmostEqual(fm.data.get_magnitude(sink.data["Input"])[0, 0, 0], 1.0)
        self.assertAlmostEqual(fm.data.get_magnitude(sink.data["Input"])[0, 0, 1], 0.5)
        self.assertAlmostEqual(fm.data.get_magnitude(sink.data["Input"])[0, 1, 0], 0.5)
        self.assertAlmostEqual(fm.data.get_magnitude(sink.data["Input"])[0, 1, 1], 0.25)