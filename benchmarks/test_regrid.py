import datetime as dt
import gc
import unittest

import finam as fm
import pytest

import finam_regrid as fmr


class TestRegrid(unittest.TestCase):
    @pytest.fixture(autouse=True)
    def setupBenchmark(self, benchmark):
        self.benchmark = benchmark

    def tearDown(self):
        self.adapter.finalize()

    def setup_adapter(self, grid1, grid2, adapter):
        self.adapter = None
        self.source = None
        self.data = None
        gc.collect()

        time = dt.datetime(2000, 1, 1)
        self.data = fm.data.full(1.0, fm.Info(time=time, grid=grid1))

        self.source = fm.components.CallbackGenerator(
            callbacks={"Step": (lambda t: self.data, fm.Info(None, grid=grid1))},
            start=time,
            step=dt.timedelta(1.0),
        )
        self.adapter = adapter
        self.source.initialize()
        self.source.outputs["Step"] >> self.adapter

        self.adapter.get_info(fm.Info(None, grid=grid2))
        self.source.connect(time)
        self.source.connect(time)
        self.source.validate()

    @pytest.mark.benchmark(group="adapters-regrid")
    def test_regrid_nearest_01_32x16(self):
        grid1 = fm.UniformGrid((32, 16))
        grid2 = fm.UniformGrid((32, 16), origin=(0.25, 0.25))

        self.setup_adapter(
            grid1, grid2, fmr.Regrid(regrid_method=fmr.RegridMethod.NEAREST_STOD)
        )
        result = self.benchmark(
            self.adapter.get_data, time=dt.datetime(2000, 1, 1), target=None
        )
        del result
        gc.collect()

    @pytest.mark.benchmark(group="adapters-regrid")
    def test_regrid_nearest_02_512x256(self):
        grid1 = fm.UniformGrid((512, 256))
        grid2 = fm.UniformGrid((512, 256), origin=(0.25, 0.25))

        self.setup_adapter(
            grid1, grid2, fmr.Regrid(regrid_method=fmr.RegridMethod.NEAREST_STOD)
        )
        result = self.benchmark(
            self.adapter.get_data, time=dt.datetime(2000, 1, 1), target=None
        )
        del result
        gc.collect()

    @pytest.mark.benchmark(group="adapters-regrid")
    def test_regrid_nearest_03_1024x512(self):
        grid1 = fm.UniformGrid((1024, 512))
        grid2 = fm.UniformGrid((1024, 512), origin=(0.25, 0.25))

        self.setup_adapter(
            grid1, grid2, fmr.Regrid(regrid_method=fmr.RegridMethod.NEAREST_STOD)
        )
        result = self.benchmark(
            self.adapter.get_data, time=dt.datetime(2000, 1, 1), target=None
        )
        del result
        gc.collect()

    @pytest.mark.benchmark(group="adapters-regrid")
    def test_regrid_linear_01_32x16(self):
        grid1 = fm.UniformGrid((32, 16))
        grid2 = fm.UniformGrid((32, 16), origin=(0.25, 0.25))

        self.setup_adapter(
            grid1, grid2, fmr.Regrid(regrid_method=fmr.RegridMethod.BILINEAR)
        )
        result = self.benchmark(
            self.adapter.get_data, time=dt.datetime(2000, 1, 1), target=None
        )
        del result
        gc.collect()

    @pytest.mark.benchmark(group="adapters-regrid")
    def test_regrid_linear_02_512x256(self):
        grid1 = fm.UniformGrid((512, 256))
        grid2 = fm.UniformGrid((512, 256), origin=(0.25, 0.25))

        self.setup_adapter(
            grid1, grid2, fmr.Regrid(regrid_method=fmr.RegridMethod.BILINEAR)
        )
        result = self.benchmark(
            self.adapter.get_data, time=dt.datetime(2000, 1, 1), target=None
        )
        del result
        gc.collect()

    @pytest.mark.benchmark(group="adapters-regrid")
    def test_regrid_linear_03_1024x512(self):
        grid1 = fm.UniformGrid((1024, 512))
        grid2 = fm.UniformGrid((1024, 512), origin=(0.25, 0.25))

        self.setup_adapter(
            grid1, grid2, fmr.Regrid(regrid_method=fmr.RegridMethod.BILINEAR)
        )
        result = self.benchmark(
            self.adapter.get_data, time=dt.datetime(2000, 1, 1), target=None
        )
        del result
        gc.collect()

    @pytest.mark.benchmark(group="adapters-regrid")
    def test_regrid_conserve_01_32x16(self):
        grid1 = fm.UniformGrid((32, 16))
        grid2 = fm.UniformGrid((32, 16), origin=(0.25, 0.25))

        self.setup_adapter(
            grid1, grid2, fmr.Regrid(regrid_method=fmr.RegridMethod.CONSERVE)
        )
        result = self.benchmark(
            self.adapter.get_data, time=dt.datetime(2000, 1, 1), target=None
        )
        del result
        gc.collect()

    @pytest.mark.benchmark(group="adapters-regrid")
    def test_regrid_conserve_02_512x256(self):
        grid1 = fm.UniformGrid((512, 256))
        grid2 = fm.UniformGrid((512, 256), origin=(0.25, 0.25))

        self.setup_adapter(
            grid1, grid2, fmr.Regrid(regrid_method=fmr.RegridMethod.CONSERVE)
        )
        result = self.benchmark(
            self.adapter.get_data, time=dt.datetime(2000, 1, 1), target=None
        )
        del result
        gc.collect()

    @pytest.mark.benchmark(group="adapters-regrid")
    def test_regrid_conserve_03_1024x512(self):
        grid1 = fm.UniformGrid((1024, 512))
        grid2 = fm.UniformGrid((1024, 512), origin=(0.25, 0.25))

        self.setup_adapter(
            grid1, grid2, fmr.Regrid(regrid_method=fmr.RegridMethod.CONSERVE)
        )
        result = self.benchmark(
            self.adapter.get_data, time=dt.datetime(2000, 1, 1), target=None
        )
        del result
        gc.collect()

    @pytest.mark.benchmark(group="adapters-regrid")
    def test_regrid_conserve_2nd_01_32x16(self):
        grid1 = fm.UniformGrid((32, 16))
        grid2 = fm.UniformGrid((32, 16), origin=(0.25, 0.25))

        self.setup_adapter(
            grid1, grid2, fmr.Regrid(regrid_method=fmr.RegridMethod.CONSERVE_2ND)
        )
        result = self.benchmark(
            self.adapter.get_data, time=dt.datetime(2000, 1, 1), target=None
        )
        del result
        gc.collect()

    @pytest.mark.benchmark(group="adapters-regrid")
    def test_regrid_conserve_2nd_02_512x256(self):
        grid1 = fm.UniformGrid((512, 256))
        grid2 = fm.UniformGrid((512, 256), origin=(0.25, 0.25))

        self.setup_adapter(
            grid1, grid2, fmr.Regrid(regrid_method=fmr.RegridMethod.CONSERVE_2ND)
        )
        result = self.benchmark(
            self.adapter.get_data, time=dt.datetime(2000, 1, 1), target=None
        )
        del result
        gc.collect()

    @pytest.mark.benchmark(group="adapters-regrid")
    def test_regrid_conserve_2nd_03_1024x512(self):
        grid1 = fm.UniformGrid((1024, 512))
        grid2 = fm.UniformGrid((1024, 512), origin=(0.25, 0.25))

        self.setup_adapter(
            grid1, grid2, fmr.Regrid(regrid_method=fmr.RegridMethod.CONSERVE_2ND)
        )
        result = self.benchmark(
            self.adapter.get_data, time=dt.datetime(2000, 1, 1), target=None
        )
        del result
        gc.collect()
