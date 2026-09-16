"""Conversions between legacy PMAP indexing and NumPy grid indexing."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np


def validate_mesh_dimensions(mesh_x: int, mesh_y: int) -> tuple[int, int]:
    """Validate and return positive mesh dimensions as integers."""
    if isinstance(mesh_x, bool) or not isinstance(mesh_x, int):
        raise TypeError("mesh_x must be an integer")
    if isinstance(mesh_y, bool) or not isinstance(mesh_y, int):
        raise TypeError("mesh_y must be an integer")
    if mesh_x <= 0 or mesh_y <= 0:
        raise ValueError("mesh_x and mesh_y must be positive")
    return mesh_x, mesh_y


def pmap_index_to_row_col(
    pmap_index: int, mesh_x: int, mesh_y: int
) -> tuple[int, int]:
    """Convert a one-based legacy PMAP index to zero-based row and column."""
    mesh_x, mesh_y = validate_mesh_dimensions(mesh_x, mesh_y)
    if isinstance(pmap_index, bool) or not isinstance(pmap_index, int):
        raise TypeError("pmap_index must be an integer")
    if not 1 <= pmap_index <= mesh_x * mesh_y:
        raise IndexError("pmap_index is outside the mesh")
    zero_based = pmap_index - 1
    return divmod(zero_based, mesh_x)


def row_col_to_pmap_index(
    row: int, column: int, mesh_x: int, mesh_y: int
) -> int:
    """Convert zero-based NumPy row and column to a one-based PMAP index."""
    mesh_x, mesh_y = validate_mesh_dimensions(mesh_x, mesh_y)
    if isinstance(row, bool) or not isinstance(row, int):
        raise TypeError("row must be an integer")
    if isinstance(column, bool) or not isinstance(column, int):
        raise TypeError("column must be an integer")
    if not 0 <= row < mesh_y or not 0 <= column < mesh_x:
        raise IndexError("row and column are outside the mesh")
    return row * mesh_x + column + 1


def pmap_to_array(
    pmap_linear: Sequence[object] | np.ndarray, mesh_x: int, mesh_y: int
) -> np.ndarray:
    """Convert flat PMAP data to shape ``(mesh_y, mesh_x)`` without narrowing dtype."""
    mesh_x, mesh_y = validate_mesh_dimensions(mesh_x, mesh_y)
    flat = np.asarray(pmap_linear)
    expected_size = mesh_x * mesh_y
    if flat.ndim != 1 or flat.size != expected_size:
        raise ValueError(f"pmap_linear must have shape ({expected_size},)")
    return flat.reshape((mesh_y, mesh_x), order="C").copy()


def array_to_pmap(array: np.ndarray, mesh_x: int, mesh_y: int) -> np.ndarray:
    """Convert a ``(mesh_y, mesh_x)`` array to flat PMAP order."""
    mesh_x, mesh_y = validate_mesh_dimensions(mesh_x, mesh_y)
    values = np.asarray(array)
    expected_shape = (mesh_y, mesh_x)
    if values.shape != expected_shape:
        raise ValueError(f"array must have shape {expected_shape}")
    return values.reshape(-1, order="C").copy()


def round_trip_preserves_data(
    pmap_linear: Sequence[object] | np.ndarray, mesh_x: int, mesh_y: int
) -> bool:
    """Return whether flat-to-grid-to-flat conversion preserves values and dtype."""
    original = np.asarray(pmap_linear)
    restored = array_to_pmap(pmap_to_array(original, mesh_x, mesh_y), mesh_x, mesh_y)
    return np.array_equal(original, restored) and original.dtype == restored.dtype
