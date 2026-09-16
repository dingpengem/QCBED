import numpy as np
import pytest

from qcbed.grids import (
    array_to_pmap,
    pmap_index_to_row_col,
    pmap_to_array,
    oshift,
    reciprocal_index_grid,
    round_trip_preserves_data,
    row_col_to_pmap_index,
    validate_mesh_dimensions,
)


@pytest.mark.parametrize(
    "pmap_index, expected",
    [(1, (0, 0)), (2, (0, 1)), (5, (1, 0)), (12, (2, 3))],
)
def test_pmap_indices_on_rectangular_mesh(pmap_index: int, expected: tuple[int, int]) -> None:
    assert pmap_index_to_row_col(pmap_index, mesh_x=4, mesh_y=3) == expected
    row, column = expected
    assert row_col_to_pmap_index(row, column, mesh_x=4, mesh_y=3) == pmap_index


def test_mesh_validation_rejects_invalid_dimensions() -> None:
    with pytest.raises(ValueError):
        validate_mesh_dimensions(4, 0)
    with pytest.raises(TypeError):
        validate_mesh_dimensions(4.0, 3)


@pytest.mark.parametrize("dtype", [np.float64, np.complex128])
def test_flat_and_array_round_trip_preserves_dtype(dtype: np.dtype) -> None:
    values = np.arange(12, dtype=dtype)
    if np.issubdtype(dtype, np.complexfloating):
        values = values + 1j * values

    array = pmap_to_array(values, mesh_x=4, mesh_y=3)
    restored = array_to_pmap(array, mesh_x=4, mesh_y=3)

    assert array.shape == (3, 4)
    assert array.dtype == dtype
    assert restored.dtype == dtype
    np.testing.assert_array_equal(restored, values)
    assert round_trip_preserves_data(values, mesh_x=4, mesh_y=3)


def test_flat_to_array_uses_row_major_numpy_mapping() -> None:
    values = np.arange(12, dtype=np.float64)
    array = pmap_to_array(values, mesh_x=4, mesh_y=3)

    assert array[0, 0] == values[0]
    assert array[1, 0] == values[4]
    assert array[2, 3] == values[11]


def test_index_validation_rejects_out_of_range_values() -> None:
    with pytest.raises(IndexError):
        pmap_index_to_row_col(0, mesh_x=4, mesh_y=3)
    with pytest.raises(IndexError):
        pmap_index_to_row_col(13, mesh_x=4, mesh_y=3)
    with pytest.raises(IndexError):
        row_col_to_pmap_index(3, 0, mesh_x=4, mesh_y=3)


def test_reciprocal_grid_shape_dtype_and_centered_origin() -> None:
    h_grid, k_grid = reciprocal_index_grid(mesh_x=6, mesh_y=4, centered=True)

    assert h_grid.shape == (4, 6)
    assert k_grid.shape == (4, 6)
    assert h_grid.dtype == np.float64
    assert k_grid.dtype == np.float64
    assert h_grid[2, 3] == 0.0
    assert k_grid[2, 3] == 0.0
    np.testing.assert_array_equal(h_grid[2], np.array([-3, -2, -1, 0, 1, 2], dtype=np.float64))
    np.testing.assert_array_equal(k_grid[:, 3], np.array([-2, -1, 0, 1], dtype=np.float64))


def test_reciprocal_grid_unshifted_order_has_origin_at_first_element() -> None:
    h_grid, k_grid = reciprocal_index_grid(mesh_x=6, mesh_y=4, centered=False)

    np.testing.assert_array_equal(h_grid[0], np.array([0, 1, 2, -3, -2, -1], dtype=np.float64))
    np.testing.assert_array_equal(k_grid[:, 0], np.array([0, 1, -2, -1], dtype=np.float64))
    assert h_grid[0, 0] == 0.0
    assert k_grid[0, 0] == 0.0


@pytest.mark.parametrize("dtype", [np.float64, np.complex128])
def test_oshift_matches_legacy_quadrant_swaps_on_rectangular_mesh(dtype) -> None:
    values = np.arange(24, dtype=np.float64).reshape((4, 6))
    expected = np.array(
        [
            [15, 16, 17, 12, 13, 14],
            [21, 22, 23, 18, 19, 20],
            [3, 4, 5, 0, 1, 2],
            [9, 10, 11, 6, 7, 8],
        ],
        dtype=np.float64,
    )

    values = values.astype(dtype)
    expected = expected.astype(dtype)
    if dtype == np.complex128:
        values += 1j * (100 + 2 * values.real)
        expected += 1j * (100 + 2 * expected.real)
    original = values.copy()
    shifted = oshift(values)

    assert shifted.shape == (4, 6)
    assert shifted.dtype == dtype
    np.testing.assert_array_equal(values, original)
    np.testing.assert_array_equal(shifted, expected)
    np.testing.assert_array_equal(oshift(shifted), values)


def test_oshift_moves_centered_origin_to_first_flat_element() -> None:
    values = np.zeros((4, 6), dtype=np.float64)
    values[2, 3] = 1.0

    shifted = oshift(values)

    assert shifted[0, 0] == 1.0
    assert np.count_nonzero(shifted) == 1


def test_reciprocal_grid_coordinates_follow_flat_pmap_indices() -> None:
    h_grid, k_grid = reciprocal_index_grid(mesh_x=6, mesh_y=4, centered=True)
    pmap_index = 1 + 3 * 6 + 4
    row, column = pmap_index_to_row_col(pmap_index, mesh_x=6, mesh_y=4)

    assert (row, column) == (3, 4)
    assert h_grid[row, column] == 1.0
    assert k_grid[row, column] == 1.0


@pytest.mark.parametrize(
    "centered, h_expected, k_expected",
    [(True, [-2, -1, 0, 1, 2], [-1, 0, 1]),
     (False, [0, 1, 2, -2, -1], [0, 1, -1])],
)
def test_odd_reciprocal_indices(centered, h_expected, k_expected) -> None:
    h, k = reciprocal_index_grid(5, 3, centered=centered)
    assert h.shape == k.shape == (3, 5)
    assert h.dtype == k.dtype == np.float64
    for row in range(3):
        np.testing.assert_array_equal(h[row], h_expected)
    for column in range(5):
        np.testing.assert_array_equal(k[:, column], k_expected)


@pytest.mark.parametrize("shape", [(3, 5), (4, 5), (3, 6), (1, 2)])
def test_legacy_oshift_rejects_odd_dimensions(shape) -> None:
    with pytest.raises(ValueError, match="OSHIFT requires even"):
        oshift(np.zeros(shape, dtype=np.complex128))


@pytest.mark.parametrize("mesh_x, mesh_y", [(5, 3), (6, 4)])
def test_pmap_complex_rectangular_orientation(mesh_x, mesh_y) -> None:
    flat = np.arange(mesh_x * mesh_y, dtype=np.float64).astype(np.complex128)
    flat += 1j * (100 - flat.real)
    array = pmap_to_array(flat, mesh_x, mesh_y)
    assert array.shape == (mesh_y, mesh_x)
    assert array.dtype == np.complex128
    for row in range(mesh_y):
        for column in range(mesh_x):
            index = row * mesh_x + column
            assert array[row, column] == complex(index, 100 - index)
            assert pmap_index_to_row_col(index + 1, mesh_x, mesh_y) == (row, column)
            assert row_col_to_pmap_index(row, column, mesh_x, mesh_y) == index + 1
    np.testing.assert_array_equal(array_to_pmap(array, mesh_x, mesh_y), flat)
