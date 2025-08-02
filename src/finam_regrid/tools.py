"""Tools for ESMF regridding."""

from __future__ import annotations

from enum import Enum

import esmpy
import finam as fm
import numpy as np
from finam.data.grid_tools import ESMF_TYPE_MAP
from pyproj import Transformer, crs

ESMF_DIM_NAMES = ["ESMF:X", "ESMF:Y", "ESMF:Z"]

ESMF_STAGGER_LOC_2D = {
    fm.Location.CELLS: esmpy.StaggerLoc.CENTER,
    fm.Location.POINTS: esmpy.StaggerLoc.CORNER,
}

ESMF_STAGGER_LOC_3D = {
    fm.Location.CELLS: esmpy.StaggerLoc.CENTER_VCENTER,
    fm.Location.POINTS: esmpy.StaggerLoc.CORNER_VFACE,
}

ESMF_STAGGER_LOC = {
    2: ESMF_STAGGER_LOC_2D,
    3: ESMF_STAGGER_LOC_3D,
}

ESMF_MESH_LOC = {
    fm.Location.CELLS: esmpy.MeshLoc.ELEMENT,
    fm.Location.POINTS: esmpy.MeshLoc.NODE,
}


class RegridCRS(Enum):
    """Selector for CRS used to regrid on."""

    SRC = 0
    """CRS of the source grid"""
    DST = 1
    """CRS of the target grid"""
    SPH = 2
    """Use lat-lon representation for both grids."""


def is_latlon(coordsys):
    """
    Check if given CRS is using lat-lon coordinates.

    Parameters
    ----------
    coordsys : any
        CRS specifier

    Returns
    -------
    bool
        Whether the CRS is using lat-lon coordinates.
    """
    coordsys = crs.CRS(coordsys)
    return coordsys.is_geographic and all(
        ax.unit_name == "degree" for ax in coordsys.axis_info
    )


def _shp(i, dim=3):
    res = dim * [1]
    res[i] = -1
    return res


def create_transformer(in_crs, out_crs, assume_target_crs=False):
    """Creates a transformer for conversion between different CRS.

    Returns
    -------
    Transformer or None
        Return None if no transform is required.
    """
    in_crs = None if in_crs is None else crs.CRS(in_crs)
    out_crs = None if out_crs is None else crs.CRS(out_crs)
    # if in_crs is not given, we will assume it to be the same as out_crs
    in_crs = out_crs if (assume_target_crs and in_crs is None) else in_crs
    transformer = (
        None
        if (in_crs is None and out_crs is None) or in_crs == out_crs
        else Transformer.from_crs(in_crs, out_crs, always_xy=True)
    )
    return transformer


def _transform_points(transformer, points):
    if transformer is None:
        return points
    return np.asarray(transformer.transform(*points.T)).T


def to_esmf(grid, transformer=None, spherical=False, mask=None):
    """Converts a FINAM grid specification to the corresponding ESMF type."""
    if isinstance(grid, fm.data.StructuredGrid):
        return _to_esmf_grid(grid, transformer, spherical, mask)
    if isinstance(grid, fm.UnstructuredPoints) or (
        isinstance(grid, fm.UnstructuredGrid) and np.all(grid.cell_types == 0)
    ):
        return _to_esmf_points(grid, transformer, spherical, mask)
    if isinstance(grid, fm.UnstructuredGrid):
        return _to_esmf_mesh(grid, transformer, spherical, mask)

    raise ValueError(f"Grid type '{grid.__class__.__name__}' not supported")


def _to_esmf_grid(grid: fm.data.StructuredGrid, transformer, spherical, mask):
    dims = np.array([d - 1 for d in grid.dims], dtype=np.int32)
    grid_dim = grid.mesh_dim
    loc = ESMF_STAGGER_LOC[grid_dim][grid.data_location]
    p_loc = ESMF_STAGGER_LOC[grid_dim][fm.Location.POINTS]
    c_loc = ESMF_STAGGER_LOC[grid_dim][fm.Location.CELLS]
    g = esmpy.Grid(
        dims,
        staggerloc=[p_loc, c_loc],
        coord_sys=esmpy.CoordSys.SPH_DEG if spherical else esmpy.CoordSys.CART,
    )
    if transformer is None:
        for i in range(grid.dim):
            grid_corner = g.get_coords(i, staggerloc=p_loc)
            grid_center = g.get_coords(i, staggerloc=c_loc)
            grid_corner[...] = grid.axes[i].reshape(*_shp(i, grid.dim))
            grid_center[...] = grid.cell_axes[i].reshape(*_shp(i, grid.dim))
    else:
        points = fm.data.grid_tools.gen_points(grid.axes, order="F")
        points = _transform_points(transformer, points)
        cell_centers = fm.data.grid_tools.gen_points(grid.cell_axes, order="F")
        cell_centers = _transform_points(transformer, cell_centers)
        for i in range(grid.dim):
            grid_corner = g.get_coords(i, staggerloc=p_loc)
            grid_center = g.get_coords(i, staggerloc=c_loc)
            grid_corner[...] = points[:, i].reshape(grid.dims, order="F")
            grid_center[...] = cell_centers[:, i].reshape(dims, order="F")

    if mask is not None:
        mask_item = g.add_item(esmpy.GridItem.MASK, staggerloc=loc)
        mask_item[...] = mask

    field = esmpy.Field(g, name=grid.name, staggerloc=loc)
    field.data[:] = np.nan
    return g, field


def _to_esmf_mesh(grid: fm.UnstructuredGrid, transformer, spherical, mask):
    loc = ESMF_MESH_LOC[grid.data_location]
    mesh = esmpy.Mesh(
        parametric_dim=grid.mesh_dim,
        spatial_dim=grid.dim,
        coord_sys=esmpy.CoordSys.SPH_DEG if spherical else esmpy.CoordSys.CART,
    )
    num_node = grid.point_count
    points = _transform_points(transformer, grid.points)
    # deal with mask
    pnt_mask = None
    ele_mask = None
    if grid.data_location == fm.Location.POINTS:
        pnt_mask = mask
    else:
        ele_mask = mask
    # Does for some reason create weird coordinates with `parametric_dim=2, spatial_dim=3`
    mesh.add_nodes(
        node_count=num_node,
        node_ids=np.arange(num_node, dtype=int) + 1,
        node_coords=points.ravel().astype(float),
        node_owners=np.zeros(num_node, dtype=int),
        node_mask=pnt_mask,
    )

    elem_types = ESMF_TYPE_MAP[grid.cell_types]
    if np.any(elem_types == -1):
        # this should only occure for line elements in 1D
        # vertices are covered by the UnstructuredPoints class
        raise ValueError("ESMF can't be used to regrid 1D data.")

    num_elem = grid.cell_count
    cell_centers = _transform_points(transformer, grid.cell_centers)
    mesh.add_elements(
        element_count=num_elem,
        element_ids=np.arange(num_elem, dtype=int) + 1,
        element_types=elem_types,
        element_conn=grid.cells_connectivity.astype(float),
        element_coords=cell_centers.ravel().astype(float),
        element_mask=ele_mask,
    )
    field = esmpy.Field(mesh, name=grid.name, meshloc=loc)
    field.data[:] = np.nan
    return mesh, field


def _to_esmf_points(grid: fm.UnstructuredPoints, transformer, spherical, mask):
    locstream = esmpy.LocStream(
        grid.point_count,
        coord_sys=esmpy.CoordSys.SPH_DEG if spherical else esmpy.CoordSys.CART,
    )

    points = _transform_points(transformer, grid.points)

    for i in range(grid.dim):
        locstream[ESMF_DIM_NAMES[i]] = points[:, i]

    if mask is not None:
        locstream["ESMF:Mask"] = mask

    field = esmpy.Field(locstream, name=grid.name)
    field.data[:] = np.nan

    return locstream, field
