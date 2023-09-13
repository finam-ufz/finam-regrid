"""Tools for ESMF regridding."""
from __future__ import annotations

import esmpy
import finam as fm
import numpy as np
from pyproj import Transformer, crs

ESMF_CELLTYPES = {
    fm.CellType.TRI: esmpy.api.constants.MeshElemType.TRI,
    fm.CellType.QUAD: esmpy.api.constants.MeshElemType.QUAD,
}

ESMF_DIM_NAMES = ["ESMF:X", "ESMF:Y", "ESMF:Z"]

ESMF_STAGGER_LOC = {
    fm.Location.CELLS: esmpy.StaggerLoc.CENTER,
    fm.Location.POINTS: esmpy.StaggerLoc.CORNER,
}

ESMF_MESH_LOC = {
    fm.Location.CELLS: esmpy.MeshLoc.ELEMENT,
    fm.Location.POINTS: esmpy.MeshLoc.NODE,
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
    if isinstance(grid, fm.data.StructuredGrid):
        return _to_esmf_grid(grid, transformer)
    if isinstance(grid, fm.UnstructuredPoints):
        return _to_esmf_points(grid, transformer)
    if isinstance(grid, fm.UnstructuredGrid):
        return _to_esmf_mesh(grid, transformer)

    raise ValueError(f"Grid type '{grid.__class__.__name__}' not supported")


def _to_esmf_grid(grid: fm.data.StructuredGrid, transformer):
    dims = np.array([d - 1 for d in grid.dims], dtype=np.int32)
    loc = ESMF_STAGGER_LOC[grid.data_location]

    g = esmpy.Grid(
        dims,
        staggerloc=[],
        coord_sys=esmpy.CoordSys.CART,
    )
    g.add_coords(staggerloc=esmpy.StaggerLoc.CENTER)
    g.add_coords(staggerloc=esmpy.StaggerLoc.CORNER)

    points = _transform_points(transformer, grid.points)
    cells = _transform_points(transformer, grid.cell_centers)

    for i in range(len(grid.axes)):
        grid_corner = g.get_coords(i, staggerloc=esmpy.StaggerLoc.CORNER)
        grid_center = g.get_coords(i, staggerloc=esmpy.StaggerLoc.CENTER)

        grid_corner[...] = points[:, i].reshape(grid.dims, order=grid.order)
        grid_center[...] = cells[:, i].reshape(dims, order=grid.order)

    field = esmpy.Field(g, name=grid.name, staggerloc=loc)
    field.data[:] = np.nan

    return g, field


def _to_esmf_mesh(grid: fm.UnstructuredGrid, transformer):
    loc = ESMF_MESH_LOC[grid.data_location]

    mesh = esmpy.Mesh(
        parametric_dim=grid.mesh_dim,
        spatial_dim=grid.dim,
        coord_sys=esmpy.CoordSys.CART,
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

    field = esmpy.Field(mesh, name=grid.name, meshloc=loc)
    field.data[:] = np.nan

    return mesh, field


def _to_esmf_points(grid: fm.UnstructuredPoints, transformer):
    locstream = esmpy.LocStream(grid.point_count, coord_sys=esmpy.CoordSys.CART)

    points = _transform_points(transformer, grid.points)

    for i in range(grid.dim):
        locstream[ESMF_DIM_NAMES[i]] = points[:, i]

    field = esmpy.Field(locstream, name=grid.name)
    field.data[:] = np.nan

    return locstream, field
