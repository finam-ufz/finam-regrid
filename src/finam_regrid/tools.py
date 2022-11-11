"""Tools for ESMF regridding."""
from __future__ import annotations

import ESMF
import finam as fm
import numpy as np

ESMF_CELLTYPES = {
    fm.CellType.TRI: ESMF.api.constants.MeshElemType.TRI,
    fm.CellType.QUAD: ESMF.api.constants.MeshElemType.QUAD,
}

ESMF_DIM_NAMES = ["ESMF:X", "ESMF:Y", "ESMF:Z"]

ESMF_STAGGER_LOC = {
    fm.Location.CELLS: ESMF.StaggerLoc.CENTER,
    fm.Location.POINTS: ESMF.StaggerLoc.CORNER,
}

ESMF_MESH_LOC = {
    fm.Location.CELLS: ESMF.MeshLoc.ELEMENT,
    fm.Location.POINTS: ESMF.MeshLoc.NODE,
}


def to_esmf(grid):
    """Converts a FINAM grid specification to the corresponding ESMF type."""
    if isinstance(grid, fm.data.grid_tools.StructuredGrid):
        return _to_esmf_grid(grid)
    if isinstance(grid, fm.UnstructuredPoints):
        return _to_esmf_points(grid)
    if isinstance(grid, fm.UnstructuredGrid):
        return _to_esmf_mesh(grid)

    raise ValueError(f"Grid type '{grid.__class__.__name__}' not supported")


def _to_esmf_grid(grid: fm.data.grid_tools.StructuredGrid):
    dims = np.array(grid.dims, dtype=np.int32)
    loc = ESMF_STAGGER_LOC[grid.data_location]

    g = ESMF.Grid(
        dims,
        num_peri_dims=1,
        staggerloc=[loc],
        coord_sys=ESMF.CoordSys.CART,
    )

    axes = grid.axes
    for i, ax in enumerate(axes):
        coords = np.asarray(ax[g.lower_bounds[loc][i] : g.upper_bounds[loc][i]])
        size = [1] * len(axes)
        size[i] = coords.size

        gridX = g.get_coords(i)
        gridX[...] = coords.reshape(size)

    field = ESMF.Field(g, name=grid.name)
    field.data[:] = np.nan

    return g, field


def _to_esmf_mesh(grid: fm.UnstructuredGrid):
    loc = ESMF_MESH_LOC[grid.data_location]

    mesh = ESMF.Mesh(
        parametric_dim=grid.mesh_dim, spatial_dim=grid.dim, coord_sys=ESMF.CoordSys.CART
    )

    node_ids = np.arange(grid.point_count)
    # Does for some reason create weird coordinates with `parametric_dim=2, spatial_dim=3`
    # Therefore removing the z coordinate
    node_coords = np.array([p[:2] for p in grid.points]).flatten()
    node_owner = np.zeros(grid.point_count)

    mesh.add_nodes(grid.point_count, node_ids, node_coords, node_owner)

    elems = grid.cells
    num_elems = elems.shape[0]
    elem_ids = np.arange(num_elems)
    elem_types = np.asarray([ESMF_CELLTYPES[e] for e in grid.cell_types])

    mesh.add_elements(
        num_elems,
        elem_ids,
        elem_types,
        elems.flatten(),
        element_coords=grid.cell_centers,
    )

    field = ESMF.Field(mesh, name=grid.name, meshloc=loc)
    field.data[:] = np.nan

    return mesh, field


def _to_esmf_points(grid: fm.UnstructuredPoints):
    locstream = ESMF.LocStream(grid.point_count, coord_sys=ESMF.CoordSys.CART)

    for i in range(grid.dim):
        locstream[ESMF_DIM_NAMES[i]] = grid.points[:, i]

    field = ESMF.Field(locstream, name=grid.name)
    field.data[:] = np.nan

    return locstream, field
