from dataclasses import FrozenInstanceError, replace
from pathlib import Path

import numpy as np
import pytest

from qcbed.absorption import (
    AbsorptionCoefficients, absorptive_atomic_scattering_factor,
    absorptive_structure_factors, combine_elastic_absorptive_structure_factors,
)
from qcbed import (AtomSite, Species, DirectLattice2D, ReciprocalLattice2D,
                   ReciprocalIndexGrid, ReciprocalVectorGrid, parse_for007)
from qcbed.grids import array_to_pmap


COEFFICIENTS = AbsorptionCoefficients(0, 1, 1)
ATOM = AtomSite(1, 0, 0, 1, "S")
SPECIES = (Species("AL", 13),)


def atomic(z=13, s=.25, **kwargs):
    options = dict(voltage_kv=200, applied_debye_waller_angstrom_squared=.86,
                   coefficients=COEFFICIENTS)
    options.update(kwargs)
    return absorptive_atomic_scattering_factor(z, s, **options)


def grid(a=20, b=10, gamma=90, dilation=1, centered=True, mesh=(6, 4)):
    return ReciprocalVectorGrid(ReciprocalLattice2D(DirectLattice2D(a, b, gamma)),
                                ReciprocalIndexGrid(*mesh, centered=centered), dilation)


def summed(atoms=(ATOM,), species=SPECIES, reciprocal=None, **kwargs):
    options = dict(voltage_kv=200, applied_debye_waller_angstrom_squared=(.86,)*len(species),
                   coefficients=(COEFFICIENTS,)*len(species), vacancy_fraction=0,
                   absorption_raw="YES_ABSORPTION")
    options.update(kwargs)
    return absorptive_structure_factors(species, atoms, grid() if reciprocal is None else reciprocal, **options)


# Independent 65-digit Decimal calculations, with coefficients (2,.5,1.5).
# Used energy=(1+q)/(2*(1-q)) and distributed B**p so that the C term cancels
# that power; this does not repeat the production expression grouping.
ATOMIC_CASES = [
    (1, 0, .0011644736797559489),
    (1, 1e-6, .0005756034771997381),
    (1, .25, .00021455287057936757),
    (1, 1.5, .00006491914736578576),
    (13, 0, .11615494071626605),
    (13, 1e-6, .057415799886858502),
    (13, .25, .021531909702254684),
    (13, 6, .0006275195646752409),
    (98, 0, 3.4351027768917547),
    (98, 1e-6, 1.6979835073148182),
    (98, .25, .6567876701842111),
    (98, 6, .03531580129877369),
]


@pytest.mark.parametrize("z,s,expected", ATOMIC_CASES)
def test_independent_atomic_factors(z, s, expected):
    value = atomic(z, s, coefficients=AbsorptionCoefficients(2, .5, 1.5))
    assert type(value) is float
    assert value == pytest.approx(expected, rel=3e-14, abs=0)


@pytest.mark.parametrize("kv,b,expected", [
    (100, .86, .05252864259107706),
    (300, .86, .03712155961696001),
    (200, 1.72, .059337303277793465),
])
def test_energy_and_applied_thermal_dependence(kv, b, expected):
    assert atomic(voltage_kv=kv, applied_debye_waller_angstrom_squared=b) == pytest.approx(expected, rel=3e-14)


def test_finite_width_origin_term_and_signed_weights():
    coeff = AbsorptionCoefficients(1, 0, 0)
    assert atomic(s=0, coefficients=coeff) == pytest.approx(.04646197628650642, rel=3e-14)
    assert atomic(s=1e-6, coefficients=coeff) == pytest.approx(.017092405872000787, rel=3e-14)
    assert atomic(s=.001, coefficients=coeff) == 0.0
    assert atomic(s=0, coefficients=AbsorptionCoefficients(-1, 0, 0)) == pytest.approx(-.04646197628650642, rel=3e-14)
    assert atomic(coefficients=AbsorptionCoefficients(0, 0, 0)) == 0.0
    with pytest.raises(FrozenInstanceError):
        coeff.inela = 2


SUM_CASES = [
    (ATOM, 0, -9.270486844482137 + 0j),
    (replace(ATOM, starting_occupancy=.25), 0, -2.317621711120534 + 0j),
    (ATOM, .25, -6.952865133361602 + 0j),
    (replace(ATOM, x_fractional=.25), 0, -4.1621351794569144e-10 - 9.270486844482137j),
]


@pytest.mark.parametrize("atom,vacancy,expected", SUM_CASES)
def test_independent_structure_factors(atom, vacancy, expected):
    # area=200 multiplies FAATM; external damping exp(-.86*.025^2).
    result = summed((atom,), vacancy_fraction=vacancy)
    assert result[2, 4] == pytest.approx(expected, abs=5e-13, rel=0)


def test_origin_area_dilation_and_oblique_cell():
    assert summed()[2, 3] == pytest.approx(-9.292395257301284, rel=3e-14)
    np.testing.assert_allclose(summed(reciprocal=grid(a=10, b=5, dilation=2)), summed(), rtol=3e-14)
    # Area 20*10*sin(30deg) is 100 to the source angle-literal precision.
    assert summed(reciprocal=grid(gamma=30))[2, 3] == pytest.approx(-4.646197628650642, abs=3e-11, rel=0)


def test_rectangular_orientation_friedel_and_flattening():
    result = summed((replace(ATOM, x_fractional=.25),))
    assert result.dtype == np.complex128 and result.shape == (4, 6)
    assert result.flags.c_contiguous
    expected = -4.1621351794569144e-10 - 9.270486844482137j
    assert array_to_pmap(result, 6, 4)[16] == pytest.approx(expected, abs=5e-13, rel=0)
    assert result[2, 2] == pytest.approx(expected.conjugate(), abs=5e-13, rel=0)
    assert result[3, 3].imag == 0
    for row, col in [(1, 2), (1, 3), (1, 4)]:
        assert result[row, col] == pytest.approx(result[4-row, 6-col].conjugate(), abs=5e-13, rel=0)
    y_result = summed((replace(ATOM, y_fractional=.25),), reciprocal=grid(a=10, b=20))
    assert y_result[3, 3] == pytest.approx(expected, abs=5e-13, rel=0)


@pytest.mark.parametrize("mesh", [(6, 4), (5, 3)])
def test_unshifted_order_and_no_masks(mesh):
    centered = summed(reciprocal=grid(mesh=mesh))
    unshifted = summed(reciprocal=grid(mesh=mesh, centered=False))
    np.testing.assert_array_equal(unshifted, np.roll(centered, (-(mesh[1]//2), -(mesh[0]//2)), (0, 1)))
    assert np.all(centered.real < 0)


def test_multiple_species_independent_parameters_and_atom_order():
    atoms = (ATOM, replace(ATOM, species_index_one_based=2))
    species = (Species("AL", 13), Species("H ", 1))
    result = summed(atoms, species, coefficients=(AbsorptionCoefficients(0, 0, 0),
                                                  AbsorptionCoefficients(2, .5, 1.5)))
    # Isolate H with nonzero coefficients; index 2 must not select Al.
    assert result[2, 3] == pytest.approx(-.23289473595118978, rel=3e-14)
    doubled = summed((ATOM, ATOM))
    assert doubled[2, 3] == pytest.approx(-18.58479051460257, rel=3e-14)


@pytest.mark.parametrize("raw,enabled", [("YES_ABSORPTION", True), ("YEanything", True),
    ("YE", True), ("yes", False), (" YE", False), ("NO", False), ("Y", False), ("", False)])
def test_literal_absorption_flag(raw, enabled):
    result = summed(absorption_raw=raw)
    assert bool(np.any(result)) is enabled


def test_immutability_empty_slice_and_vacancy():
    reciprocal = grid()
    before = reciprocal.stl_inverse_angstrom.copy()
    species, atoms, coeff, thermal = list(SPECIES), [ATOM], [COEFFICIENTS], [.86]
    first = summed(atoms, species, reciprocal, coefficients=coeff,
                   applied_debye_waller_angstrom_squared=thermal)
    np.testing.assert_array_equal(first, summed())
    assert species == list(SPECIES) and atoms == [ATOM] and coeff == [COEFFICIENTS] and thermal == [.86]
    np.testing.assert_array_equal(reciprocal.stl_inverse_angstrom, before)
    assert not reciprocal.stl_inverse_angstrom.flags.writeable
    assert not np.any(summed(()))
    assert not np.any(summed(vacancy_fraction=1))


def test_combination_sign_and_input_immutability():
    elastic = np.array([[1+2j, 3-4j, 5j], [6, -7j, 8+9j]], dtype=np.complex128)
    absorptive = np.array([[-2-3j, 4+5j, -6j], [7j, 8, -9-10j]], dtype=np.complex128)
    e, a = elastic.copy(), absorptive.copy()
    expected = [[4+0j, -2+0j, 6+5j], [ -1+0j, 1j, 18+0j]]
    actual = combine_elastic_absorptive_structure_factors(elastic, absorptive)
    np.testing.assert_array_equal(actual, expected)
    assert actual.dtype == np.complex128 and actual.shape == (2, 3)
    np.testing.assert_array_equal(e, elastic)
    np.testing.assert_array_equal(a, absorptive)


@pytest.mark.parametrize("z", [0, -1, 99, 1.5, True, float("nan"), float("inf")])
def test_invalid_atomic_numbers(z):
    with pytest.raises((TypeError, ValueError), match="atomic_number"):
        atomic(z=z)


@pytest.mark.parametrize("s", [-1, 6.01, float("nan"), float("inf"), -float("inf"), True, "1"])
def test_invalid_coordinates(s):
    with pytest.raises((TypeError, ValueError), match="stl_inverse_angstrom"):
        atomic(s=s)


def test_hydrogen_coordinate_bound():
    with pytest.raises(ValueError, match="0..1.5"):
        atomic(z=1, s=1.50001)


@pytest.mark.parametrize("field", ["voltage_kv", "applied_debye_waller_angstrom_squared"])
@pytest.mark.parametrize("value", [0, -1, float("nan"), float("inf"), True, "1"])
def test_invalid_positive_parameters(field, value):
    with pytest.raises((TypeError, ValueError), match=field):
        atomic(**{field: value})


@pytest.mark.parametrize("field", ["inela", "inelb", "inelc"])
@pytest.mark.parametrize("value", [float("nan"), float("inf"), True, "1"])
def test_invalid_coefficients(field, value):
    with pytest.raises((TypeError, ValueError), match=field):
        AbsorptionCoefficients(**({"inela": 0, "inelb": 1, "inelc": 1} | {field: value}))


@pytest.mark.parametrize("field", ["x_fractional", "y_fractional", "starting_occupancy"])
@pytest.mark.parametrize("value", [float("nan"), float("inf"), True])
def test_invalid_atoms(field, value):
    with pytest.raises((TypeError, ValueError), match=field):
        summed((replace(ATOM, **{field: value}),))


@pytest.mark.parametrize("index", [0, 2, True, .5])
def test_invalid_species_index(index):
    with pytest.raises((TypeError, ValueError), match="species_index_one_based"):
        summed((replace(ATOM, species_index_one_based=index),))


@pytest.mark.parametrize("value", [-.1, 1.1, float("nan"), float("inf"), True])
def test_invalid_vacancy(value):
    with pytest.raises((TypeError, ValueError), match="vacancy_fraction"):
        summed(vacancy_fraction=value)


def test_invalid_structure_parameters():
    with pytest.raises(ValueError, match="unresolved G"):
        summed((replace(ATOM, record_type="G"),))
    with pytest.raises(ValueError, match="starting_occupancy"):
        summed((replace(ATOM, starting_occupancy=-.1),))
    with pytest.raises(ValueError, match="one value per species"):
        summed(coefficients=())
    with pytest.raises(ValueError, match="species must not be empty"):
        summed(species=())
    with pytest.raises(TypeError, match="absorption_raw"):
        summed(absorption_raw=True)
    with pytest.raises(TypeError, match="reciprocal_grid"):
        summed(reciprocal="bad")
    with pytest.raises(TypeError, match="coefficients"):
        atomic(coefficients=(0, 1, 1))
    with pytest.raises(ValueError, match="positive"):
        summed(applied_debye_waller_angstrom_squared=(0,))
    with pytest.raises(ValueError, match="stl_inverse_angstrom"):
        summed(reciprocal=grid(a=.01))


@pytest.mark.parametrize("bad", [np.ones((2, 3)), np.ones((2, 3), dtype=np.complex64), [[0j]]])
def test_combination_rejects_implicit_dtype_conversion(bad):
    with pytest.raises(TypeError, match="complex128"):
        combine_elastic_absorptive_structure_factors(bad, np.ones((2, 3), dtype=np.complex128))


def test_combination_rejects_shapes_and_nonfinite():
    good = np.zeros((2, 3), dtype=np.complex128)
    for bad in (np.zeros((3, 2), dtype=np.complex128), np.zeros(6, dtype=np.complex128),
                np.zeros((0, 3), dtype=np.complex128), np.full((2, 3), np.nan, dtype=np.complex128)):
        with pytest.raises(ValueError, match="shape|finite"):
            combine_elastic_absorptive_structure_factors(good, bad)


def test_unrepresentable_inputs_fail_explicitly():
    with pytest.raises(ValueError, match="voltage in volts"):
        atomic(voltage_kv=1e308)
    with pytest.raises(ValueError, match="atomic phase"):
        summed((replace(ATOM, x_fractional=1e308),))
    huge = np.full((2, 3), 1e308+1e308j, dtype=np.complex128)
    with pytest.raises(ValueError, match="finite"):
        combine_elastic_absorptive_structure_factors(huge, huge)


def test_saved_fortran_absorption_listing():
    root = Path(__file__).resolve().parents[1] / "reference_data"
    parsed = parse_for007(root / "FOR007.DAT")
    s = parsed.slices[0]
    reciprocal = grid(a=s.cell_a_angstrom, b=s.cell_b_angstrom, gamma=s.lattice_gamma_degrees,
                      mesh=(parsed.mesh.mesh_x, parsed.mesh.mesh_y))
    result = summed(parsed.atoms_by_slice[0], parsed.species, reciprocal,
                    voltage_kv=parsed.energy_kv, absorption_raw=parsed.absorption_raw)
    # Explicit source-input values: INELA/B/C=(0,1,1), B=.86, VAC=0,d=1.
    lines = (root / "FOR0013.DAT").read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if "STRUCTURE FACTOR LISTING FOR SLICE" in line)
    rows = lines[start+2:start+27]
    assert len(rows) == 25
    for line in rows:
        h, k, _, _, amplitude, phase = line.split()
        expected = float(amplitude) * np.exp(1j*np.deg2rad(float(phase)))
        actual = result[int(k)+parsed.mesh.mesh_y//2, int(h)+parsed.mesh.mesh_x//2]
        assert actual == pytest.approx(expected, rel=3e-7, abs=5e-8)
