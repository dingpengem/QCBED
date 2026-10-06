import cmath
from dataclasses import FrozenInstanceError, replace
import math
from pathlib import Path

import numpy as np
import pytest

from qcbed import AtomSite, Species, MeshConfig, SliceGeometry
from qcbed.legacy_io import LegacyInput, parse_for007
from qcbed.simulation import (
    SliceAppliedParameters, SamplePreparationParameters, PreparedSliceType,
    PreparedMultisliceSample, prepare_multislice_sample, build_tilt_propagators, simulate_one_tilt,
)
from qcbed.reciprocal import DirectLattice2D, ReciprocalLattice2D, ReciprocalIndexGrid, ReciprocalVectorGrid
from qcbed.structure_factors import elastic_structure_factors
from qcbed.absorption import AbsorptionCoefficients, absorptive_structure_factors
from qcbed.phase_grating import phase_grating_from_structure_factors
from qcbed.propagator import PropagationParameters, propagation_kernel
from qcbed.multislice import simulate_slice_sequence


def model(n=1, mesh=(4, 2), gamma=90, absorption="NO", two_species=False):
    species = (Species("Al", 13), Species("H ", 1)) if two_species else (Species("Al", 13),)
    atoms = tuple((AtomSite(1, .25+i*.125, .125, .5, "S"), AtomSite(2, i*.5, .375, .25, "S"))
                  if two_species else (AtomSite(1, 0, 0, 1, "S"),) for i in range(n))
    geometries = tuple(SliceGeometry(2+i, .125*i, -.25*i, 20, 10, gamma) for i in range(n))
    return LegacyInput("Synthetic", 200, False, (0, 0, 1), MeshConfig(*mesh), geometries,
                       species, atoms, "out.dat", "in.dat", title_raw="Synthetic", absorption_raw=absorption)


def config(m, *, sequence=None, method="direct", ordering="unshifted", dsap=.5):
    return SamplePreparationParameters(
        applied_voltage_kv=200,
        slices=tuple(SliceAppliedParameters((.86, 1.2) if len(m.species) == 2 else (.86,), .2, 1)
                     for _ in m.slices),
        dsap_inverse_angstrom=dsap,
        absorption_coefficients=tuple(AbsorptionCoefficients(0, 1, 1) for _ in m.species)
            if m.absorption_raw[:2] == "YE" else None,
        slice_sequence_zero_based=tuple(range(len(m.slices))) if sequence is None else sequence,
        beam_indices_hk=((0, 0), (1, 0), (-1, 0), (1, 0)),
        reciprocal_ordering=ordering, phase_grating_method=method,
        boundary_policy="legacy", reflection_multipliers="unity",
    )


def simulate(prepared, h=0, k=0):
    return simulate_one_tilt(prepared, tilt_h_undilated_fractional=h, tilt_k_undilated_fractional=k)


@pytest.mark.parametrize("method", ["direct", "fft"])
@pytest.mark.parametrize("sequence", [(0,), (0, 0, 0)])
def test_one_type_uniform_origin_and_no_rmt(method, sequence):
    m = model()
    p = replace(config(m, sequence=sequence, method=method, dsap=0),
                slices=(SliceAppliedParameters((0,), 0, 1),), beam_indices_hk=((0, 0),))
    prepared = prepare_multislice_sample(m, p)
    # Hard-checked source Al origin coefficient 281.953595001 / area 200.
    # Independent printed-source sigma at 200 kV from the phase-grating audit.
    phase = -.0007288360746343612*1.409767975005
    expected = cmath.exp(1j*phase*len(sequence))
    result = simulate(prepared)
    np.testing.assert_allclose(result.amplitudes, [expected], atol=8e-16, rtol=0)
    np.testing.assert_allclose(result.intensities, [1.], atol=1e-15, rtol=0)
    assert prepared.wavelength_angstrom == pytest.approx(.0250794525129636, rel=1e-14)


@pytest.mark.parametrize("raw, enabled", [("NO", False), ("yes", False), (" YES", False),
                                         ("'YES'", False), ("YEanything  ", True)])
def test_raw_absorption_semantics_and_independent_dc(raw, enabled):
    m = model(absorption=raw)
    p = replace(config(m, dsap=0), slices=(SliceAppliedParameters((.86,), 0, 1),))
    prepared = prepare_multislice_sample(m, p)
    assert prepared.absorption_raw == raw and prepared.absorption_enabled == enabled
    # Independent Al DC grating for 200 kV, area=200, B=.86, weights=(0,1,1).
    expected = (.9932497254026199-.0010205542781868584j if enabled
                else cmath.exp(-1j*.0007288360746343612*1.409767975005))
    np.testing.assert_allclose(prepared.slice_types[0].phase_grating, expected, atol=4e-16, rtol=0)
    np.testing.assert_allclose(simulate(prepared).amplitudes, [expected, 0, 0, 0], atol=4e-16, rtol=0)


def nested_dft(values, to_reciprocal):
    ny, nx = values.shape
    result = np.zeros_like(values)
    for r in range(ny):
        for c in range(nx):
            result[r, c] = sum(values[j, i]*cmath.exp((-1 if to_reciprocal else 1)*2j*math.pi*(i*c/nx+j*r/ny))
                               for j in range(ny) for i in range(nx))/(nx*ny if to_reciprocal else 1)
    return result


def tiny_reference(m, p, h0, k0):
    """Independent retained Fourier modes, printed atomic values, scalar DFTs.

    On the 2x4 grid, legacy masks leave only K=0 and H=-1,0,1. Their
    STL is 0 or .025: no production scattering or grid routine is called.
    """
    fe_zero = (281.953595001, 25.327466761)
    fe_one = (275.340445007875, 24.824622984184895)
    sigma, wavelength = .0007288360746343612, .0250794525129636
    cosine = math.cos(3.14159265358/2)
    gratings, kernels = [], []
    for geometry, atoms, state in zip(m.slices, m.atoms_by_slice, p.slices):
        assert state.dilation == 1 and geometry.lattice_gamma_degrees == 90
        coefficients = {}
        for h in (-1, 0, 1):
            total = 0j
            for atom in atoms:
                z = atom.species_index_one_based-1
                f = fe_zero[z] if h == 0 else fe_one[z]
                total += f/200*cmath.exp(1j*6.283185307*h*atom.x_fractional)*math.exp(
                    -state.debye_waller_angstrom_squared[z]*(h*.025)**2)*atom.starting_occupancy*(1-state.vacancy_fraction)
            coefficients[h] = total
        t = np.empty((2, 4), dtype=np.complex128)
        for r in range(2):
            for c in range(4):
                potential = sum(f*cmath.exp(-2j*math.pi*h*c/4) for h, f in coefficients.items())
                t[r, c] = cmath.exp(-1j*sigma*potential)
        kernel = np.zeros((2, 4), dtype=np.complex128)
        for h in (-1, 0, 1):
            # Cartesian dual vectors, independent of the production metric grouping.
            gx, gy = h/20, -h*cosine/20
            tx, ty = h0/20, -h0*cosine/20+k0/10
            phase = 6.2831853*(wavelength/2*(gx*gx+gy*gy+2*(tx*gx+ty*gy))*geometry.thickness_angstrom
                              -h*geometry.x_shift_fractional)
            kernel[0, h % 4] = cmath.exp(1j*phase)
        gratings.append(t)
        kernels.append(kernel)
    wave = np.ones((2, 4), dtype=np.complex128)
    for occurrence, index in enumerate(p.slice_sequence_zero_based):
        reciprocal = nested_dft(wave*gratings[index], True)*kernels[index]
        if occurrence != len(p.slice_sequence_zero_based)-1:
            wave = nested_dft(reciprocal, False)
    return gratings, reciprocal


@pytest.mark.parametrize("sequence", [(0,), (0, 1), (0, 1, 0, 1)])
@pytest.mark.parametrize("tilt", [(0, 0), (.2, 0), (0, -.3), (.2, -.3)])
@pytest.mark.parametrize("method", ["direct", "fft"])
def test_complete_independent_small_model(sequence, tilt, method):
    m = model(n=2, two_species=True)
    p = config(m, sequence=sequence, method=method)
    prepared = prepare_multislice_sample(m, p)
    expected_gratings, expected_wave = tiny_reference(m, p, *tilt)
    result = simulate(prepared, *tilt)
    for stored, expected in zip(prepared.slice_types, expected_gratings):
        np.testing.assert_allclose(stored.phase_grating, expected, atol=5e-16, rtol=0)
    np.testing.assert_allclose(result.final_reciprocal_wave, expected_wave, atol=3e-15, rtol=0)
    expected_beams = expected_wave[0, [0, 1, 3, 1]]
    np.testing.assert_allclose(result.amplitudes, expected_beams, atol=3e-15, rtol=0)
    np.testing.assert_allclose(result.intensities, abs(expected_beams)**2, atol=4e-15, rtol=0)
    assert result.amplitudes.dtype == np.complex128 and result.intensities.dtype == np.float64


@pytest.mark.parametrize("mesh,gamma", [((4, 4), 90), ((6, 4), 90), ((6, 4), 60)])
@pytest.mark.parametrize("ordering", ["centered", "unshifted"])
def test_direct_fft_and_manual_lower_level_composition(mesh, gamma, ordering):
    m = model(n=2, mesh=mesh, gamma=gamma, two_species=True, absorption="YES")
    p = config(m, sequence=(0, 1, 0, 1), ordering=ordering)
    p = replace(p, slices=(p.slices[0], replace(p.slices[1], dilation=.9, vacancy_fraction=.1)),
                beam_indices_hk=((0, 0), (1, 1), (-1, 0), (0, -1), (1, 1)))
    direct = prepare_multislice_sample(m, p)
    fft = prepare_multislice_sample(m, replace(p, phase_grating_method="fft"))
    manual_gratings, manual_kernels = [], []
    for geometry, atoms, state in zip(m.slices, m.atoms_by_slice, p.slices):
        g = ReciprocalVectorGrid(ReciprocalLattice2D(DirectLattice2D(20, 10, gamma)),
                                  ReciprocalIndexGrid(*mesh, centered=ordering == "centered"), state.dilation)
        e = elastic_structure_factors(m.species, atoms, g,
            applied_debye_waller_angstrom_squared=state.debye_waller_angstrom_squared, vacancy_fraction=state.vacancy_fraction)
        a = absorptive_structure_factors(m.species, atoms, g, voltage_kv=200,
            applied_debye_waller_angstrom_squared=state.debye_waller_angstrom_squared,
            vacancy_fraction=state.vacancy_fraction, coefficients=p.absorption_coefficients, absorption_raw="YES")
        manual_gratings.append(phase_grating_from_structure_factors(e, a, g, voltage_kv=200,
                                dsap_inverse_angstrom=p.dsap_inverse_angstrom, method="direct"))
        pp = PropagationParameters(.0250794525129636, geometry.thickness_angstrom,
                                   .2/state.dilation, -.3/state.dilation,
                                   geometry.x_shift_fractional, geometry.y_shift_fractional)
        manual_kernels.append(propagation_kernel(g, pp, ordering=ordering))
    expected = simulate_slice_sequence(manual_gratings, manual_kernels, p.slice_sequence_zero_based,
        beam_indices_hk=p.beam_indices_hk, reciprocal_ordering=ordering, propagator_ordering=ordering)
    for ds, fs, manual in zip(direct.slice_types, fft.slice_types, manual_gratings):
        np.testing.assert_array_equal(ds.phase_grating, manual)
        np.testing.assert_allclose(fs.phase_grating, manual, atol=4e-16, rtol=0)
    a, b = simulate(direct, .2, -.3), simulate(fft, .2, -.3)
    np.testing.assert_allclose(a.final_reciprocal_wave, expected.final_reciprocal_wave, atol=2e-15, rtol=0)
    np.testing.assert_allclose(b.final_reciprocal_wave, a.final_reciprocal_wave, atol=2e-15, rtol=0)
    assert a.beam_indices_hk == p.beam_indices_hk and a.amplitudes[1] == a.amplitudes[-1]
    other = prepare_multislice_sample(m, replace(p, reciprocal_ordering="unshifted" if ordering == "centered" else "centered"))
    np.testing.assert_allclose(simulate(other, .2, -.3).amplitudes, a.amplitudes, atol=2e-15, rtol=0)


def test_preparation_reuse_and_per_type_propagator_construction(monkeypatch):
    import qcbed.simulation as module
    counts = {"grating": 0, "propagator": 0}
    old_grating, old_prop = module.phase_grating_from_structure_factors, module.propagation_kernel
    tilts = []
    def grating(*args, **kwargs):
        counts["grating"] += 1
        return old_grating(*args, **kwargs)
    def prop(grid, parameters, **kwargs):
        counts["propagator"] += 1
        tilts.append((parameters.tilt_h_fractional, parameters.tilt_k_fractional))
        return old_prop(grid, parameters, **kwargs)
    monkeypatch.setattr(module, "phase_grating_from_structure_factors", grating)
    monkeypatch.setattr(module, "propagation_kernel", prop)
    m = model(n=2)
    p = config(m, sequence=(0, 1)*6)
    p = replace(p, slices=(p.slices[0], replace(p.slices[1], dilation=.5)))
    prepared = prepare_multislice_sample(m, p)
    assert counts == {"grating": 2, "propagator": 0}
    simulate(prepared, .2, -.3)
    assert counts == {"grating": 2, "propagator": 2}
    assert tilts == [(.2, -.3), (.4, -.6)]  # Source caller's per-slice division, once.
    simulate(prepared, -.1, .1)
    assert counts == {"grating": 2, "propagator": 4}


def test_shifts_only_change_propagators_and_immutable_prepared_state():
    m = model()
    p = config(m)
    original = prepare_multislice_sample(m, p)
    changed = prepare_multislice_sample(replace(m, slices=(replace(m.slices[0], x_shift_fractional=.25),)), p)
    np.testing.assert_array_equal(original.slice_types[0].phase_grating, changed.slice_types[0].phase_grating)
    kwargs = dict(tilt_h_undilated_fractional=0, tilt_k_undilated_fractional=0)
    a, b = build_tilt_propagators(original, **kwargs), build_tilt_propagators(changed, **kwargs)
    assert b[0][0, 1]/a[0][0, 1] == pytest.approx(1.7948965149208059e-9-1j, abs=3e-16, rel=0)
    before = original.slice_types[0].phase_grating.copy()
    simulate(original, .3, -.2)
    np.testing.assert_array_equal(before, original.slice_types[0].phase_grating)
    for values in (original.slice_types[0].phase_grating, original.slice_types[0].reciprocal_grid.stl_inverse_angstrom, a[0]):
        with pytest.raises(ValueError):
            values.setflags(write=True)
    with pytest.raises(FrozenInstanceError):
        original.wavelength_angstrom = 1
    with pytest.raises(FrozenInstanceError):
        p.slices[0].dilation = 2
    assert m.slices[0].x_shift_fractional == 0


def test_reference_input_preparation_and_zero_tilt():
    path = Path(__file__).parents[1]/"tests/fixtures/legacy_reference"/"FOR007.DAT"
    m = parse_for007(path)
    lines = path.read_text().splitlines()
    # Explicit initial applied values are supported by these input records and
    # SFMAP/UPDPARS: VAC=0 -> d=1 and B=DW, and all reflection factors are unity.
    assert list(map(float, lines[10].split())) == [840., 1.]
    assert list(map(float, lines[13].split())) == [0., 0.]
    assert list(map(float, lines[14].split())) == [0., 0.]
    assert list(map(float, lines[16].split())) == [0., 0., 1., 1., 1., 1.]
    assert float(lines[17]) == .86 and float(lines[7]) == 2.5
    assert all(list(map(float, line.split()[2:])) == [1., 1.] for line in lines[46:54])
    assert 414*m.slices[0].thickness_angstrom < 840 <= 415*m.slices[0].thickness_angstrom
    beam_start = next(i for i, line in enumerate(lines) if line.strip() == "out.dat")+1
    beams = tuple(tuple(map(int, line.split())) for line in lines[beam_start+1:beam_start+10])
    p = replace(config(m, method="fft", sequence=(0,)*415, dsap=2.5),
                slices=(SliceAppliedParameters((.86,), 0, 1),), beam_indices_hk=beams)
    prepared = prepare_multislice_sample(m, p)
    assert len(prepared.slice_types) == 1
    assert prepared.parameters.slice_sequence_zero_based == (0,)*415
    assert prepared.slice_types[0].phase_grating.shape == (16, 16)
    assert prepared.beam_positions_row_col == ((0, 0), (1, 1), (15, 1), (15, 15), (1, 15), (0, 2), (14, 0), (0, 14), (2, 0))
    assert np.all(np.isfinite(prepared.slice_types[0].phase_grating))
    kernels = build_tilt_propagators(prepared, tilt_h_undilated_fractional=0, tilt_k_undilated_fractional=0)
    assert len(kernels) == 1 and np.all(np.isfinite(kernels[0]))
    result = simulate(prepared)
    assert result.intensities.shape == (9,) and np.all(np.isfinite(result.intensities))
    assert np.all(result.intensities >= 0)


@pytest.mark.parametrize("field,bad", [("dilation", 0), ("dilation", -1), ("dilation", np.inf),
                                      ("vacancy_fraction", -.1), ("vacancy_fraction", 1.1),
                                      ("vacancy_fraction", True), ("vacancy_fraction", np.nan)])
def test_invalid_applied_scalars(field, bad):
    with pytest.raises((TypeError, ValueError), match=field):
        replace(SliceAppliedParameters((.86,), 0, 1), **{field: bad})


@pytest.mark.parametrize("bad", [-1, np.nan, np.inf, True, ".86"])
def test_invalid_thermal_scalar(bad):
    with pytest.raises((TypeError, ValueError), match="debye_waller"):
        SliceAppliedParameters((bad,), 0, 1)


@pytest.mark.parametrize("field,bad", [("applied_voltage_kv", 0), ("applied_voltage_kv", np.inf),
                                      ("applied_voltage_kv", True), ("dsap_inverse_angstrom", -1),
                                      ("dsap_inverse_angstrom", np.nan), ("phase_grating_method", "auto"),
                                      ("reciprocal_ordering", "FFT"), ("boundary_policy", "none"),
                                      ("reflection_multipliers", "refined"), ("slice_sequence_zero_based", (1,)),
                                      ("slice_sequence_zero_based", ()), ("slice_sequence_zero_based", (True,)),
                                      ("absorption_coefficients", (1,))])
def test_invalid_configuration(field, bad):
    with pytest.raises((TypeError, ValueError), match="voltage|dsap|method|ordering|boundary|reflection|sequence|coefficients"):
        replace(config(model()), **{field: bad})


@pytest.mark.parametrize("case", ["cell", "mesh", "count", "atoms_count", "species", "G", "thermal", "beams", "domain"])
def test_preflight_before_expensive_sums(monkeypatch, case):
    import qcbed.simulation as module
    def forbidden(*args, **kwargs):
        raise AssertionError("structure-factor calculation started before validation")
    monkeypatch.setattr(module, "elastic_structure_factors", forbidden)
    m = model(n=2)
    p = config(m)
    if case == "cell":
        m = replace(m, slices=(m.slices[0], replace(m.slices[1], cell_b_angstrom=11)))
    elif case == "mesh":
        m = replace(m, mesh=MeshConfig(5, 3))
    elif case == "count":
        m = replace(m, slices=m.slices[:1])
    elif case == "atoms_count":
        m = replace(m, atoms_by_slice=m.atoms_by_slice[:1])
    elif case in ("species", "G"):
        bad = replace(m.atoms_by_slice[1][0], **({"species_index_one_based": 2} if case == "species" else {"record_type": "G"}))
        m = replace(m, atoms_by_slice=(m.atoms_by_slice[0], (bad,)))
    elif case == "thermal":
        p = replace(p, slices=(p.slices[0], replace(p.slices[1], debye_waller_angstrom_squared=())))
    elif case == "beams":
        p = replace(p, beam_indices_hk=((99, 0),))
    elif case == "domain":
        m = replace(m, slices=tuple(replace(s, cell_a_angstrom=.01, cell_b_angstrom=.01) for s in m.slices))
    with pytest.raises(ValueError) as error:
        prepare_multislice_sample(m, p)
    if case in ("cell", "species", "G", "thermal"):
        assert "slice type 1" in str(error.value)


@pytest.mark.parametrize("enabled", [True, False])
def test_absorption_configuration_not_guessed(enabled):
    m = model(absorption="YE" if enabled else "NO")
    p = config(m)
    with pytest.raises(ValueError, match="absorption_coefficients"):
        prepare_multislice_sample(m, replace(p, absorption_coefficients=None if enabled else (AbsorptionCoefficients(0, 1, 1),)))
    if enabled:
        with pytest.raises(ValueError, match="one entry per global species"):
            prepare_multislice_sample(m, replace(p, absorption_coefficients=()))
        with pytest.raises(ValueError, match="slice type 0.*positive for absorption"):
            prepare_multislice_sample(m, replace(p, slices=(replace(p.slices[0], debye_waller_angstrom_squared=(0,)),)))


@pytest.mark.parametrize("field,bad", [("x_fractional", np.nan), ("y_fractional", np.inf),
                                      ("starting_occupancy", -.1), ("species_index_one_based", 0)])
def test_invalid_atoms_identify_slice(field, bad):
    m = model()
    m = replace(m, atoms_by_slice=((replace(m.atoms_by_slice[0][0], **{field: bad}),),))
    with pytest.raises(ValueError, match="slice type 0"):
        prepare_multislice_sample(m, config(m))


def test_invalid_wavelength_and_prepared_cross_slice_mesh():
    m = model()
    with pytest.raises(ValueError, match="wavelength"):
        prepare_multislice_sample(m, replace(config(m), applied_voltage_kv=1e308))
    prepared = prepare_multislice_sample(m, config(m))
    with pytest.raises(ValueError, match="mesh shape"):
        replace(prepared, mesh=MeshConfig(6, 4))
    with pytest.raises(ValueError, match="ordering"):
        replace(prepared, parameters=replace(prepared.parameters, reciprocal_ordering="centered"))
    with pytest.raises(ValueError, match="phase_grating shape"):
        replace(prepared.slice_types[0], phase_grating=np.ones((4, 6), dtype=np.complex128))
    with pytest.raises(ValueError, match="geometry cell"):
        replace(prepared.slice_types[0], geometry=replace(m.slices[0], cell_a_angstrom=30))


def test_cross_slice_prepared_mesh_rejected_and_input_lists_snapshotted():
    m = model(n=2)
    p = config(m)
    prepared = prepare_multislice_sample(m, p)
    grid = ReciprocalVectorGrid(ReciprocalLattice2D(DirectLattice2D(20, 10, 90)),
                                ReciprocalIndexGrid(6, 4, centered=False))
    different_mesh = PreparedSliceType(m.slices[1], grid, np.ones((4, 6), dtype=np.complex128))
    with pytest.raises(ValueError, match="slice type 1.*mesh shape"):
        replace(prepared, slice_types=(prepared.slice_types[0], different_mesh))
    atoms = [list(group) for group in m.atoms_by_slice]
    sequence, beams = [0, 1], [[0, 0]]
    p = replace(p, slice_sequence_zero_based=sequence, beam_indices_hk=beams)
    snapshot = prepare_multislice_sample(replace(m, atoms_by_slice=atoms), p)
    before = simulate(snapshot).amplitudes.copy()
    atoms[0].clear()
    sequence[0], beams[0][0] = 1, 1
    np.testing.assert_array_equal(simulate(snapshot).amplitudes, before)
    assert snapshot.parameters.slice_sequence_zero_based == (0, 1)
    assert snapshot.parameters.beam_indices_hk == ((0, 0),)


@pytest.mark.parametrize("bad", [np.nan, np.inf, True, "0"])
def test_invalid_tilt_inputs(bad):
    prepared = prepare_multislice_sample(model(), config(model()))
    with pytest.raises((TypeError, ValueError), match="tilt_h_undilated_fractional"):
        simulate(prepared, bad, 0)
    with pytest.raises((TypeError, ValueError), match="tilt_k_undilated_fractional"):
        simulate(prepared, 0, bad)
