from dataclasses import FrozenInstanceError
import math
from pathlib import Path

import numpy as np
import pytest

from qcbed import elastic_atomic_scattering_factor, parse_for007
from qcbed.grids import array_to_pmap, pmap_index_to_row_col, reciprocal_index_grid
from qcbed.reciprocal import (
    DirectLattice2D, ReciprocalLattice2D, ReciprocalIndexGrid, ReciprocalVectorGrid,
)


# Analytic geometry uses exact mathematical angles. Allow the documented
# ~1e-11 angular difference from the literal 3.14159265358 in current source.
ATOL = 2e-11


@pytest.mark.parametrize('a,b', [(2.0, 2.0), (2.0, 4.0)])
def test_orthogonal_cells(a, b):
    lattice = ReciprocalLattice2D(DirectLattice2D(a, b, 90))
    assert lattice.direct.area_angstrom_squared == pytest.approx(a*b, abs=ATOL)
    assert lattice.a_star_inverse_angstrom == pytest.approx(1/a, abs=ATOL)
    assert lattice.b_star_inverse_angstrom == pytest.approx(1/b, abs=ATOL)
    assert lattice.cos_gamma_star == pytest.approx(0, abs=ATOL)
    np.testing.assert_allclose(lattice.a_star_vector_inverse_angstrom, [1/a, 0], atol=ATOL, rtol=0)
    np.testing.assert_allclose(lattice.b_star_vector_inverse_angstrom, [0, 1/b], atol=ATOL, rtol=0)


@pytest.mark.parametrize('gamma, cosine', [(60, .5), (120, -.5)])
def test_oblique_cell_and_dual_basis(gamma, cosine):
    direct = DirectLattice2D(2, 3, gamma)
    reciprocal = ReciprocalLattice2D(direct)
    # Independent analytic triangular lattice with sin(gamma)=sqrt(3)/2.
    assert direct.area_angstrom_squared == pytest.approx(3*math.sqrt(3), rel=0, abs=ATOL)
    assert reciprocal.a_star_inverse_angstrom == pytest.approx(1/math.sqrt(3), rel=0, abs=ATOL)
    assert reciprocal.b_star_inverse_angstrom == pytest.approx(2/(3*math.sqrt(3)), rel=0, abs=ATOL)
    assert reciprocal.cos_gamma_star == pytest.approx(-cosine, rel=0, abs=ATOL)
    expected_direct = np.array([[2, 0], [3*cosine, 1.5*math.sqrt(3)]])
    stars = np.array([reciprocal.a_star_vector_inverse_angstrom,
                      reciprocal.b_star_vector_inverse_angstrom])
    np.testing.assert_allclose(expected_direct @ stars.T, np.eye(2), atol=ATOL, rtol=0)
    actual_direct = np.array([direct.a_vector_angstrom, direct.b_vector_angstrom])
    np.testing.assert_allclose(actual_direct @ stars.T, np.eye(2), atol=5e-16, rtol=0)


@pytest.mark.parametrize('centered', [True, False])
@pytest.mark.parametrize('mesh_x,mesh_y', [(6, 4), (5, 3)])
def test_grids_preserve_order_units_and_flattening(centered, mesh_x, mesh_y):
    indices = ReciprocalIndexGrid(mesh_x, mesh_y, centered)
    grid = ReciprocalVectorGrid(ReciprocalLattice2D(DirectLattice2D(2, 4, 90)), indices)
    h, k = reciprocal_index_grid(mesh_x, mesh_y, centered=centered)
    np.testing.assert_array_equal(indices.h, h)
    np.testing.assert_array_equal(indices.k, k)
    for values in (indices.h, indices.k, grid.gx_inverse_angstrom,
                   grid.gy_inverse_angstrom, grid.magnitude_inverse_angstrom, grid.stl_inverse_angstrom):
        assert values.shape == (mesh_y, mesh_x)
        assert values.dtype == np.float64
    origin = (mesh_y//2, mesh_x//2) if centered else (0, 0)
    assert grid.gx_inverse_angstrom[origin] == grid.gy_inverse_angstrom[origin] == 0
    assert grid.magnitude_inverse_angstrom[origin] == grid.stl_inverse_angstrom[origin] == 0
    flat_x = array_to_pmap(grid.gx_inverse_angstrom, mesh_x, mesh_y)
    flat_y = array_to_pmap(grid.gy_inverse_angstrom, mesh_x, mesh_y)
    for pmap in range(1, mesh_x*mesh_y+1):
        row, col = pmap_index_to_row_col(pmap, mesh_x, mesh_y)
        # Explicit ideal rectangular reciprocal metric, not production formula.
        gx, gy = h[row, col]/2, k[row, col]/4
        assert flat_x[pmap-1] == pytest.approx(gx, rel=0, abs=ATOL)
        assert flat_y[pmap-1] == pytest.approx(gy, rel=0, abs=ATOL)
        assert grid.magnitude_inverse_angstrom[row, col] == pytest.approx(math.hypot(gx, gy), rel=0, abs=ATOL)
        assert grid.stl_inverse_angstrom[row, col] == pytest.approx(math.hypot(gx, gy)/2, rel=0, abs=ATOL)


def test_oblique_vectors_metric_and_signs():
    grid = ReciprocalVectorGrid(ReciprocalLattice2D(DirectLattice2D(2, 3, 60)), ReciprocalIndexGrid(5, 3))
    # H=K=1: g=(1/2, sqrt(3)/18), |g|=sqrt(7/27).
    assert grid.gx_inverse_angstrom[2, 3] == pytest.approx(.5, rel=0, abs=ATOL)
    assert grid.gy_inverse_angstrom[2, 3] == pytest.approx(math.sqrt(3)/18, rel=0, abs=ATOL)
    assert grid.magnitude_inverse_angstrom[2, 3] == pytest.approx(math.sqrt(7/27), rel=0, abs=ATOL)
    assert grid.stl_inverse_angstrom[2, 3] == pytest.approx(math.sqrt(7/27)/2, rel=0, abs=ATOL)
    assert grid.gx_inverse_angstrom[0, 1] == -grid.gx_inverse_angstrom[2, 3]
    assert grid.gy_inverse_angstrom[0, 1] == -grid.gy_inverse_angstrom[2, 3]
    assert grid.magnitude_inverse_angstrom[0, 1] == grid.magnitude_inverse_angstrom[2, 3]
    # |g| from Cartesian components must agree with the source metric.
    np.testing.assert_allclose(np.hypot(grid.gx_inverse_angstrom, grid.gy_inverse_angstrom),
                               grid.magnitude_inverse_angstrom, rtol=5e-16, atol=5e-16)


def test_dilation_is_explicit_reciprocal_scaling():
    lattice = ReciprocalLattice2D(DirectLattice2D(2, 3, 60))
    indices = ReciprocalIndexGrid(6, 4)
    original = ReciprocalVectorGrid(lattice, indices)
    dilated = ReciprocalVectorGrid(lattice, indices, dilation=2)
    for name in ('gx_inverse_angstrom', 'gy_inverse_angstrom', 'magnitude_inverse_angstrom', 'stl_inverse_angstrom'):
        np.testing.assert_array_equal(getattr(dilated, name), getattr(original, name)/2)


def test_stl_scattering_path_has_no_two_pi_or_full_g_factor():
    grid = ReciprocalVectorGrid(ReciprocalLattice2D(DirectLattice2D(20, 10, 90)), ReciprocalIndexGrid(4, 2, False))
    assert grid.magnitude_inverse_angstrom[0, 1] == .05
    assert grid.stl_inverse_angstrom[0, 1] == .025
    assert elastic_atomic_scattering_factor(13, grid.stl_inverse_angstrom[0, 1]) == pytest.approx(
        275.340445007875, rel=2e-12, abs=1e-12
    )


@pytest.mark.parametrize('args', [(0, 2, 90), (2, -1, 90), (2, 3, 0), (2, 3, 180), (2, 3, 181),
                                 (math.nan, 2, 90), (2, math.inf, 90), (2, 3, math.nan),
                                 (2, 3, 1e-20), (1e300, 1e300, 90)])
def test_invalid_geometry(args):
    with pytest.raises(ValueError):
        DirectLattice2D(*args)


@pytest.mark.parametrize('args', [(True, 2, 90), (2, '3', 90), (2, 3, False)])
def test_invalid_geometry_types(args):
    with pytest.raises(TypeError):
        DirectLattice2D(*args)


@pytest.mark.parametrize('args', [(0, 2), (3, -1), (True, 2), (3.0, 2)])
def test_invalid_mesh(args):
    with pytest.raises((TypeError, ValueError)):
        ReciprocalIndexGrid(*args)


@pytest.mark.parametrize('dilation', [0, -1, math.nan, math.inf, 1e300, 1e-300])
def test_invalid_dilation(dilation):
    with pytest.raises(ValueError):
        ReciprocalVectorGrid(ReciprocalLattice2D(DirectLattice2D(2, 3, 60)), ReciprocalIndexGrid(2, 2), dilation)


def test_representations_and_buffers_are_immutable():
    direct = DirectLattice2D(2, 3, 60)
    reciprocal = ReciprocalLattice2D(direct)
    indices = ReciprocalIndexGrid(4, 2)
    grid = ReciprocalVectorGrid(reciprocal, indices)
    for obj, name in [(direct, 'a_angstrom'), (reciprocal, 'cos_gamma_star'), (indices, 'centered'), (grid, 'dilation')]:
        with pytest.raises(FrozenInstanceError):
            setattr(obj, name, 0)
    for values in (indices.h, indices.k, grid.gx_inverse_angstrom, grid.gy_inverse_angstrom,
                   grid.magnitude_inverse_angstrom, grid.stl_inverse_angstrom):
        with pytest.raises(ValueError):
            values[0, 0] = 99
        with pytest.raises(ValueError):
            values.setflags(write=True)


def test_reference_input_geometry():
    parsed = parse_for007(Path(__file__).parents[1] / 'tests/fixtures/legacy_reference' / 'FOR007.DAT')
    geometry = parsed.slices[0]
    reciprocal = ReciprocalLattice2D(DirectLattice2D(
        geometry.cell_a_angstrom, geometry.cell_b_angstrom, geometry.lattice_gamma_degrees))
    # Independent 50-digit Decimal reciprocal of input a=b=2.8630754 A.
    assert reciprocal.a_star_inverse_angstrom == pytest.approx(0.34927476936164517, rel=1e-13, abs=0)
    assert reciprocal.b_star_inverse_angstrom == reciprocal.a_star_inverse_angstrom
    # Actual FOR0013.DAT reports 0.34927 at five decimal places. This is a
    # rounded-output cross-check only, not executable numerical equivalence.
    assert reciprocal.a_star_inverse_angstrom == pytest.approx(.34927, rel=0, abs=5e-6)
