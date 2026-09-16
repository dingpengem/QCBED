import numpy as np
import pytest

from qcbed.grids import (
    array_to_pmap,
    pmap_index_to_row_col,
    pmap_to_array,
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
