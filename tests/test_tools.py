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
        assert_allclose(
            g.get_coords(0, staggerloc=ESMF.StaggerLoc.CORNER)[:, 0], grid.axes[0]
        )
        assert_allclose(
            g.get_coords(1, staggerloc=ESMF.StaggerLoc.CORNER)[0, :], grid.axes[1]
        )

    def test_to_esmf_grid_point(self):
        grid = fm.UniformGrid((20, 15), data_location=fm.Location.POINTS)

        g, f = to_esmf(grid)

        self.assertIsInstance(g, ESMF.Grid)
        assert_allclose(
            g.get_coords(0, staggerloc=ESMF.StaggerLoc.CORNER)[:, 0], grid.axes[0]
        )
        assert_allclose(
            g.get_coords(1, staggerloc=ESMF.StaggerLoc.CORNER)[0, :], grid.axes[1]
        )

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

    def test_regrid(self):
        grid1 = fm.UniformGrid((21, 17))
        grid2 = fm.UniformGrid((11, 9), spacing=(2.0, 2.0))

        g1, f1 = to_esmf(grid1)
        g2, f2 = to_esmf(grid2)

        regrid = ESMF.Regrid(
            f1,
            f2,
            regrid_method=ESMF.RegridMethod.BILINEAR,
            unmapped_action=ESMF.UnmappedAction.IGNORE,
        )

        f1.data[:] = 2.0

        regrid(f1, f2)

        self.assertTrue(all(f == 2 for f in f2))
