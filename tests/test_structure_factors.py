from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from qcbed import (
    AtomSite, Species, SliceGeometry, DirectLattice2D, ReciprocalLattice2D,
    ReciprocalIndexGrid, ReciprocalVectorGrid, elastic_structure_factors,
    parse_for007,
)
from qcbed.grids import array_to_pmap


AL = (Species("AL", 13),)
ORIGIN = AtomSite(1, 0.0, 0.0, 1.0, "S")
ATOL = 3e-12


def grid(a=20, b=10, gamma=90, centered=True, mesh=(6, 4), dilation=1):
    return ReciprocalVectorGrid(ReciprocalLattice2D(DirectLattice2D(a, b, gamma)),
                                ReciprocalIndexGrid(*mesh, centered=centered), dilation)


def evaluate(atoms=(ORIGIN,), *, species=AL, reciprocal=None, thermal=None, vacancy=0):
    return elastic_structure_factors(
        species, atoms, grid() if reciprocal is None else reciprocal,
        applied_debye_waller_angstrom_squared=(0,) * len(species) if thermal is None else thermal,
        vacancy_fraction=vacancy,
    )


# Decimal (60 digits) independent evaluation: source table Lagrange polynomial,
# Taylor sine/cosine at the printed TWOPI, Decimal.exp for thermal factors.
# At H=1,K=0: STL=.025; Al fe=275.340445007875, area=200.
# These literal expectations never call the implementation to generate values.
HAND_CASES = [
    ((ORIGIN,), (0,), 0, 1.376702225039375 + 0j),
    ((replace(ORIGIN, x_fractional=.25),), (0,), 0,
     6.180927559250617e-11 + 1.376702225039375j),
    ((ORIGIN, ORIGIN), (0,), 0, 2.75340445007875 + 0j),
    ((ORIGIN, replace(ORIGIN, x_fractional=.5)), (0,), 0,
     5.55005502248115e-21 + 1.2361855118501233e-10j),
    ((replace(ORIGIN, starting_occupancy=.25),), (0,), 0, .34417555625984375 + 0j),
    ((ORIGIN,), (0,), .25, 1.03252666877953125 + 0j),
    # B=1600 A^2, STL^2=.000625 => exp(-1), no extra 2*pi or 1/4.
    ((ORIGIN,), (1600,), 0, .5064604452069665 + 0j),
]


@pytest.mark.parametrize("atoms,thermal,vacancy,expected", HAND_CASES)
def test_hand_calculated_factors(atoms, thermal, vacancy, expected):
    actual = evaluate(atoms, thermal=thermal, vacancy=vacancy)
    assert actual[2, 4] == pytest.approx(expected, abs=ATOL, rel=0)


def test_origin_and_two_species_lookup():
    atoms = (ORIGIN, replace(ORIGIN, species_index_one_based=2))
    result = evaluate(atoms, species=(Species("AL", 13), Species("H ", 1)))
    # (281.953595001 + 25.327466761)/200 at g=0.
    assert result[2, 3] == pytest.approx(1.53640530881, abs=ATOL, rel=0)
    assert result[2, 4] == pytest.approx(1.500825339960299475, abs=ATOL, rel=0)
    # Species-specific B: H is fully occupied but independently damped.
    first_only = evaluate((ORIGIN,), thermal=(1600,))
    mixed = evaluate(atoms, species=(Species("AL", 13), Species("H ", 1)), thermal=(1600, 0))
    assert mixed[2, 4] - first_only[2, 4] == pytest.approx(.124123114920924475, abs=ATOL, rel=0)


def test_zero_coordinate_has_unit_phase_and_thermal_factor():
    result = evaluate((replace(ORIGIN, x_fractional=.37, y_fractional=-.81),), thermal=(1600,))
    assert result[2, 3] == pytest.approx(1.409767975005, abs=ATOL, rel=0)


def test_rectangular_layout_flattening_phase_sign_and_friedel():
    reciprocal = grid()
    result = evaluate((replace(ORIGIN, x_fractional=.25),), reciprocal=reciprocal)
    assert result.shape == (4, 6)
    assert result.dtype == np.complex128
    assert result.flags.c_contiguous
    flat = array_to_pmap(result, 6, 4)
    # Zero-based row=2,col=4 means Fortran PMAP=17, H=1,K=0.
    assert flat[16] == pytest.approx(6.180927559250617e-11 + 1.376702225039375j, abs=ATOL, rel=0)
    assert result[2, 2] == pytest.approx(6.180927559250617e-11 - 1.376702225039375j, abs=ATOL, rel=0)
    # H=0,K=1 has zero phase: x must not be interpreted as a row coordinate.
    assert result[3, 3].imag == 0
    for row, col in [(1, 2), (1, 3), (1, 4), (2, 2)]:
        assert result[row, col] == pytest.approx(result[4-row, 6-col].conjugate(), abs=ATOL, rel=0)
    y_result = evaluate((replace(ORIGIN, y_fractional=.25),), reciprocal=grid(a=10, b=20))
    assert y_result[3, 3] == pytest.approx(6.180927559250617e-11 + 1.376702225039375j, abs=ATOL, rel=0)


@pytest.mark.parametrize("mesh", [(6, 4), (5, 3)])
def test_unshifted_and_generic_odd_grid(mesh):
    centered = evaluate(reciprocal=grid(mesh=mesh))
    unshifted = evaluate(reciprocal=grid(mesh=mesh, centered=False))
    np.testing.assert_array_equal(unshifted, np.roll(centered, (-(mesh[1]//2), -(mesh[0]//2)), (0, 1)))
    assert unshifted[0, 0] == pytest.approx(1.409767975005, abs=ATOL, rel=0)
    # Raw formula evaluates boundary points too; no hidden legacy map mask.
    assert np.all(np.abs(centered) > 0)


def test_oblique_cell():
    result = evaluate(reciprocal=grid(a=20, b=20, gamma=60))
    # Independent Decimal metric with printed source pi; Lagrange stencil
    # (.01,.02,.03,.04) -> (5.867,5.800,5.692,5.547), multiplied by 47.878009.
    assert result[2, 3] == pytest.approx(.8139299198655844, abs=ATOL, rel=0)
    assert result[2, 4] == pytest.approx(.7886596996082075, abs=ATOL, rel=0)


def test_dilation_scales_area_and_reciprocal_coordinate():
    result = evaluate(reciprocal=grid(a=10, b=5, dilation=2))
    assert result[2, 3] == pytest.approx(1.409767975005, abs=ATOL, rel=0)
    assert result[2, 4] == pytest.approx(1.376702225039375, abs=ATOL, rel=0)


def test_slice_shifts_and_thickness_do_not_enter_atomic_sum():
    geometry = SliceGeometry(2, cell_a_angstrom=20, cell_b_angstrom=10, lattice_gamma_degrees=90)
    shifted = replace(geometry, thickness_angstrom=7, x_shift_fractional=.25, y_shift_fractional=-.375)
    outputs = [evaluate(reciprocal=grid(a=s.cell_a_angstrom, b=s.cell_b_angstrom,
                                      gamma=s.lattice_gamma_degrees)) for s in (geometry, shifted)]
    # SFMAP has no XSHIFT/YSHIFT term. The API deliberately takes reciprocal
    # geometry alone; the nonzero translations belong to a future PROPG call.
    np.testing.assert_array_equal(*outputs)
    assert outputs[1][2, 4] == pytest.approx(1.376702225039375, abs=ATOL, rel=0)


def test_input_immutability_and_empty_or_vacant_slice():
    reciprocal = grid()
    before = [a.copy() for a in (reciprocal.indices.h, reciprocal.indices.k, reciprocal.stl_inverse_angstrom)]
    species, atoms, thermal = list(AL), [ORIGIN], [0.0]
    evaluate(atoms, species=species, thermal=thermal, reciprocal=reciprocal)
    assert species == list(AL) and atoms == [ORIGIN] and thermal == [0.0]
    for original, current in zip(before, (reciprocal.indices.h, reciprocal.indices.k, reciprocal.stl_inverse_angstrom)):
        np.testing.assert_array_equal(original, current)
        assert not current.flags.writeable
    np.testing.assert_array_equal(evaluate(()), np.zeros((4, 6), dtype=np.complex128))
    np.testing.assert_array_equal(evaluate(vacancy=1), np.zeros((4, 6), dtype=np.complex128))


@pytest.mark.parametrize("index", [0, -1, 2, .5, True])
def test_invalid_species_index(index):
    with pytest.raises((TypeError, ValueError), match="species_index_one_based"):
        evaluate((replace(ORIGIN, species_index_one_based=index),))


@pytest.mark.parametrize("field", ["x_fractional", "y_fractional", "starting_occupancy"])
@pytest.mark.parametrize("value", [float("nan"), float("inf"), -float("inf"), True, "0.5"])
def test_invalid_atom_values(field, value):
    with pytest.raises((TypeError, ValueError), match=field):
        evaluate((replace(ORIGIN, **{field: value}),))


def test_negative_occupancy_and_unresolved_g():
    with pytest.raises(ValueError, match="starting_occupancy.*nonnegative"):
        evaluate((replace(ORIGIN, starting_occupancy=-.1),))
    for kind in ("G", "s", "?"):
        with pytest.raises(ValueError, match="record_type.*unresolved G"):
            evaluate((replace(ORIGIN, record_type=kind),))


@pytest.mark.parametrize("value", [-1, float("nan"), float("inf"), True])
def test_invalid_thermal(value):
    with pytest.raises((TypeError, ValueError), match="applied_debye_waller"):
        evaluate(thermal=(value,))


@pytest.mark.parametrize("value", [-.01, 1.01, float("nan"), float("inf"), True])
def test_invalid_vacancy(value):
    with pytest.raises((TypeError, ValueError), match="vacancy_fraction"):
        evaluate(vacancy=value)


def test_invalid_species_thermal_count_grid_and_scattering_domain():
    with pytest.raises(ValueError, match="one value per species"):
        evaluate(thermal=())
    with pytest.raises(ValueError, match="species must not be empty"):
        evaluate(species=())
    for z in (0, 99, True):
        with pytest.raises((TypeError, ValueError), match="atomic_number"):
            evaluate(species=(Species("?", z),))
    with pytest.raises(TypeError, match="reciprocal_grid"):
        evaluate(reciprocal="wrong")
    with pytest.raises(ValueError, match="stl_inverse_angstrom"):
        evaluate(reciprocal=grid(a=.01, b=.01))


def test_unrepresentable_phase_rejected():
    with pytest.raises(ValueError, match="atomic phase.*finite"):
        evaluate((replace(ORIGIN, x_fractional=1e308),))


def test_supplied_fortran_elastic_listing():
    root = Path(__file__).resolve().parents[1] / "tests/fixtures/legacy_reference"
    parsed = parse_for007(root / "FOR007.DAT")
    geometry = parsed.slices[0]
    reciprocal = grid(a=geometry.cell_a_angstrom, b=geometry.cell_b_angstrom,
                      gamma=geometry.lattice_gamma_degrees,
                      mesh=(parsed.mesh.mesh_x, parsed.mesh.mesh_y))
    # FOR007: B=.8600000, constant vacancy=0, no surface correction;
    # FOR0013 confirms the two effective occupancies remain .5 each.
    result = evaluate(parsed.atoms_by_slice[0], species=parsed.species,
                      reciprocal=reciprocal, thermal=(.86,))
    lines = (root / "FOR0013.DAT").read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if "STRUCTURE FACTOR LISTING FOR SLICE" in line)
    rows = lines[start+2:start+27]
    assert len(rows) == 25
    for line in rows:
        h, k, amplitude, phase, *_ = line.split()
        # Compare complex coefficients against the recorded amplitude/phase,
        # rather than intensity alone. The listing thresholds tiny amplitudes.
        expected = float(amplitude) * np.exp(1j * np.deg2rad(float(phase)))
        actual = result[int(k)+parsed.mesh.mesh_y//2, int(h)+parsed.mesh.mesh_x//2]
        assert actual == pytest.approx(expected, rel=2e-7, abs=5e-7)
