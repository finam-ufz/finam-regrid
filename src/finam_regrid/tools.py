"""Tools for ESMF regridding."""
from __future__ import annotations

import ESMF
import finam as fm
import numpy as np
from pyproj import Transformer, crs

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


def create_transformer(in_crs, out_crs):
    """Creates a transformer for conversion between different CRS.

    Returns
    -------
    Transformer or None
        Return None if no transform is required.
    """
    in_crs = None if in_crs is None else crs.CRS(in_crs)
    out_crs = None if out_crs is None else crs.CRS(out_crs)
    transformer = (
        None
        if (in_crs is None and out_crs is None) or in_crs == out_crs
        else Transformer.from_crs(in_crs, out_crs)
    )
    return transformer


def _transform_points(transformer, points):
    if transformer is None:
        return points
    return np.asarray(list(transformer.itransform(points)))


def to_esmf(grid, transformer=None):
    """Converts a FINAM grid specification to the corresponding ESMF type."""
    if isinstance(grid, fm.data.grid_tools.StructuredGrid):
        return _to_esmf_grid(grid, transformer)
    if isinstance(grid, fm.UnstructuredPoints):
        return _to_esmf_points(grid, transformer)
    if isinstance(grid, fm.UnstructuredGrid):
        return _to_esmf_mesh(grid, transformer)

    raise ValueError(f"Grid type '{grid.__class__.__name__}' not supported")


def _to_esmf_grid(grid: fm.data.grid_tools.StructuredGrid, transformer):
    dims = np.array([d - 1 for d in grid.dims], dtype=np.int32)
    loc = ESMF_STAGGER_LOC[grid.data_location]

    g = ESMF.Grid(
        dims,
        staggerloc=[],
        coord_sys=ESMF.CoordSys.CART,
    )
    g.add_coords(staggerloc=ESMF.StaggerLoc.CENTER)
    g.add_coords(staggerloc=ESMF.StaggerLoc.CORNER)

    axes = zip(grid.axes, grid.cell_axes)
    for i, (ax_corner, ax) in enumerate(axes):
        size = [1] * grid.dim
        size_corner = [1] * grid.dim
        size[i] = ax.size
        size_corner[i] = ax_corner.size

        grid_corner = g.get_coords(i, staggerloc=ESMF.StaggerLoc.CORNER)
        grid_center = g.get_coords(i, staggerloc=ESMF.StaggerLoc.CENTER)

        grid_corner[...] = ax_corner.reshape(size_corner)
        grid_center[...] = ax.reshape(size)

    if transformer is not None:
        pass

    field = ESMF.Field(g, name=grid.name, staggerloc=loc)
    field.data[:] = np.nan

    return g, field


def _to_esmf_mesh(grid: fm.UnstructuredGrid, transformer):
    loc = ESMF_MESH_LOC[grid.data_location]

    mesh = ESMF.Mesh(
        parametric_dim=grid.mesh_dim, spatial_dim=grid.dim, coord_sys=ESMF.CoordSys.CART
    )

    node_ids = np.arange(grid.point_count)
    points = _transform_points(transformer, grid.points)

    # Does for some reason create weird coordinates with `parametric_dim=2, spatial_dim=3`
    # Therefore removing the z coordinate
    node_coords = np.array([p[:2] for p in points]).flatten()
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


def _to_esmf_points(grid: fm.UnstructuredPoints, transformer):
    locstream = ESMF.LocStream(grid.point_count, coord_sys=ESMF.CoordSys.CART)

    points = _transform_points(transformer, grid.points)

    for i in range(grid.dim):
        locstream[ESMF_DIM_NAMES[i]] = points[:, i]

    field = ESMF.Field(locstream, name=grid.name)
    field.data[:] = np.nan

    return locstream, field
