import unittest

import ESMF
import finam as fm
from numpy.testing import assert_allclose

from finam_regrid import to_esmf


class TestTools(unittest.TestCase):
    def test_to_esmf_grid(self):
        grid = fm.UniformGrid((20, 15))

        g, f = to_esmf(grid)

        self.assertIsInstance(g, ESMF.Grid)
        assert_allclose(g.get_coords(0)[:, 0], grid.axes[0])
        assert_allclose(g.get_coords(1)[0, :], grid.axes[1])

    def test_to_esmf_mesh(self):
        points = [
            [0, 0],
            [1, 0],
            [0, 1],
            [1, 1],
        ]
        cells = [
            [0, 1, 3],
            [0, 3, 2],
        ]
        types = [
            fm.CellType.TRI,
            fm.CellType.TRI,
        ]
        grid = fm.UnstructuredGrid(
            points, cells, types, data_location=fm.Location.CELLS
        )

        g, f = to_esmf(grid)

        self.assertIsInstance(g, ESMF.Mesh)

    def test_to_esmf_points(self):
        points = [
            [0, 0],
            [1, 0],
            [0, 1],
            [1, 1],
        ]
        grid = fm.UnstructuredPoints(points)

        g, f = to_esmf(grid)

        self.assertIsInstance(g, ESMF.LocStream)
