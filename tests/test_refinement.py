from dataclasses import FrozenInstanceError, replace
from pathlib import Path
import cmath

import numpy as np
import pytest

from qcbed import AtomSite, Species, MeshConfig, SliceGeometry, LegacyInput, ExperimentalIntensityData
from qcbed.refinement import (RefinementRegion, RefinementReflection, RefinementState,
    build_refinement_parameter_layout as layout, pack_refinement_parameters as pack,
    apply_refinement_parameters as apply)
from qcbed.refinement_io import parse_refinement_state
from qcbed.refinement_trial import evaluate_refinement_trial as evaluate, _gradient, _reflection_scales
from qcbed.scan import LegacyScanDefinition, generate_scan_tilts
from qcbed.simulation import SamplePreparationParameters, SliceAppliedParameters, prepare_multislice_sample
from qcbed.tilt_series import simulate_tilt_series_serial
from qcbed.objective import compute_legacy_objective
from qcbed.reciprocal import ReciprocalIndexGrid


def state():
    model = LegacyInput("synthetic", 200., False, (0, 0, 1), MeshConfig(4, 2),
        (SliceGeometry(2., 0., 0., 20., 10., 90.),), (Species("Al", 13),),
        ((AtomSite(1, 0., 0., 1., "S"),),), "out", "in", absorption_raw="NO")
    region = RefinementRegion(2., 0, 1, 0, (0.,), (0,), (0.,), (0,), (0, 0., 0., 1.))
    scan = LegacyScanDefinition((1, 0, 0), (0, 1, 0), (0, 0), (0, 0), (0, 0), 2, 2)
    return RefinementState(model, (region,), 0., 0, (), (), scan, ("NN",)*3,
                           ((),), ("NONE",), 0., 1.)


def data(gradient=0):
    return ExperimentalIntensityData(2, 2, gradient, ((0, 0),),
        np.array([[1.], [2.], [3.], [4.]]), np.ones((4, 1)), np.ones((4, 1)))


def test_empty_and_single_parameter():
    s = state()
    assert layout(s).parameters == ()
    assert pack(s).shape == (0,) and pack(s).dtype == np.float64
    assert apply(s, pack(s)) == s
    s = replace(s, regions=(replace(s.regions[0], thickness_mask=1),))
    assert [(p.category, p.target) for p in layout(s).parameters] == [("depth", (0,))]
    mapped = apply(s, np.array([-5.]))
    assert mapped.regions[0].depth_angstrom == 5
    assert s.regions[0].depth_angstrom == 2


def all_categories():
    s = state()
    model = replace(s.model, energy_refined=True, absorption_raw="YES",
        atoms_by_slice=((replace(s.model.atoms_by_slice[0][0], x_mask=1, y_mask=1,
                                x_fractional=.25, y_fractional=-.5),),))
    r = replace(s.regions[0], thickness_mask=1, vacancy_flag=1,
        vacancy_values=(-.2, -.3), vacancy_masks=(1, 1), debye_waller=(-.8,), debye_waller_masks=(1,))
    return replace(s, model=model, regions=(r,), dilation=-2., dilation_mask=1,
        absorption=((-1., -2., -3.),), absorption_masks=((1, 1, 1),),
        scan_masks=("YY", "NY", "YN"), reflections=((RefinementReflection(1, 0, -2., 3., "YY"),),))


def test_exact_order_and_literal_pack_unpack_oracle():
    s = all_categories()
    expected = [("voltage", ()), ("depth", (0,)), ("vacancy_coefficient", (0, 0)),
        ("vacancy_coefficient", (0, 1)), ("dilation", ()), ("absorption", (0, 0)),
        ("absorption", (0, 1)), ("absorption", (0, 2)), ("debye_waller", (0, 0)),
        ("atom_coordinate", (0, 0, 0)), ("atom_coordinate", (0, 0, 1)),
        ("scan_coordinate", (0, 0)), ("scan_coordinate", (0, 1)),
        ("scan_coordinate", (1, 1)), ("scan_coordinate", (2, 0)),
        ("reflection_multiplier", (0, 0, 0)), ("reflection_multiplier", (0, 0, 1))]
    assert [(p.category, p.target) for p in layout(s).parameters] == expected
    np.testing.assert_array_equal(pack(s), [200, 2, .2, .3, .5, 1, 2, 3, -.8, .25, -.5,
                                           1e-10, 1e-10, 1e-10, 1e-10, -2, 3])
    # Explicit source-ordered cursor, independent of production dispatch.
    v = -np.arange(1., 18.)
    out = apply(s, v)
    assert out.model.energy_kv == 1
    assert out.regions[0].depth_angstrom == 2
    assert out.regions[0].vacancy_values == (3, 4)
    assert out.dilation == .2
    assert out.absorption == ((6, 7, 8),)
    assert out.regions[0].debye_waller == (9,)
    assert (out.model.atoms_by_slice[0][0].x_fractional, out.model.atoms_by_slice[0][0].y_fractional) == (-10, -11)
    assert out.scan.corner_a == (-12, -13)
    assert out.scan.corner_b == (0, -14) and out.scan.corner_c == (-15, 0)
    assert (out.reflections[0][0].elastic, out.reflections[0][0].absorptive) == (-16, -17)
    np.testing.assert_array_equal(v, -np.arange(1., 18.))


@pytest.mark.parametrize("v,expected", [(-2., 2.), (2., 1.), (-.3, .3), (0., 0.)])
def test_constant_vacancy_operation_order(v, expected):
    s = state()
    s = replace(s, regions=(replace(s.regions[0], vacancy_masks=(1,)),))
    assert apply(s, np.array([v])).regions[0].vacancy_values == (expected,)


def test_linked_absorption_and_inactive_values():
    s = state()
    s = replace(s, model=replace(s.model, absorption_raw="YES"),
                absorption=((-4., 2., 99.),), absorption_masks=((0, 1, -1),))
    out = apply(s, np.array([-3.]))
    assert out.absorption == ((-4., 3., 3.),)
    assert s.absorption == ((-4., 2., 99.),)


def test_simplex_affine_conversion_no_clipping_and_roundtrip():
    s = state()
    s = replace(s, regions=(replace(s.regions[0], thickness_mask=1, depth_angstrom=10.),), ranges=((1.2, .8),))
    # a=1/(.8-1.2)/10=-.25; b=12. q=1.5 maps to physical 6, outside [8,12].
    assert apply(s, np.array([1.5]), coordinate_mode="legacy_simplex").regions[0].depth_angstrom == pytest.approx(6)
    assert apply(s, np.array([1.5])).regions[0].depth_angstrom == 1.5
    assert apply(s, pack(s)) == s
    assert layout(s).parameters[0].range_factors == (1.2, .8)
    with pytest.raises(ValueError, match="singular"):
        apply(replace(s, ranges=((1., 1.),)), np.array([0.]), coordinate_mode="legacy_simplex")
    with pytest.raises(ValueError, match="range"):
        apply(replace(s, ranges=()), np.array([0.]), coordinate_mode="legacy_simplex")


@pytest.mark.parametrize("bad", [np.ones(1), np.ones((0, 1)), np.array([], dtype=np.float32), [], np.array([np.nan])])
def test_invalid_vector(bad):
    with pytest.raises((TypeError, ValueError), match="trial_vector"):
        apply(state(), bad)


def test_deep_immutability_and_dynamic_count():
    s = state()
    refs = [RefinementReflection(i, 0, 1, 1, "YY") for i in range(501)]
    s = replace(s, reflections=(refs,))
    refs.clear()
    assert len(layout(s).parameters) == 1002
    packed = pack(s)
    with pytest.raises(ValueError):
        packed.setflags(write=True)
    with pytest.raises(FrozenInstanceError):
        layout(s).parameters[0].category = "depth"


@pytest.mark.parametrize("filename", ["FOR007.DAT", "STARTINGFOR007.DAT"])
def test_reference_exact_layout(filename):
    s = parse_refinement_state(Path(__file__).parents[1]/"tests/fixtures/legacy_reference"/filename)
    expected = [("depth", (0,)), ("absorption", (0, 1)), ("absorption", (0, 2)), ("debye_waller", (0, 0))]
    expected += [("scan_coordinate", (i, j)) for i in range(3) for j in range(2)]
    expected += [("reflection_multiplier", (0, i, 0)) for i in range(8)]
    assert [(p.category, p.target) for p in layout(s).parameters] == expected
    assert len(pack(s)) == 18
    assert all(p.range_factors == (1.05, .95) for p in layout(s).parameters)
    assert [(r.h, r.k) for r in s.reflections[0]] == [(1, 1), (1, -1), (-1, 1), (-1, -1), (2, 0), (0, -2), (-2, 0), (2, 0)]
    if filename == "FOR007.DAT":
        np.testing.assert_array_equal(pack(s), [840, 1, 1, .86, -.508156, .407641, .407641,
                                              .508156, -.407641, -.508156, *([1.]*8)])


@pytest.mark.parametrize("mode,vector", [("physical", [4.]), ("legacy_simplex", [-1.])])
def test_trial_manual_forward_oracle_and_determinism(mode, vector):
    s = state()
    s = replace(s, regions=(replace(s.regions[0], thickness_mask=1),), ranges=((1.5, .5),))
    # a=-.5, b=3: q=-.5 maps to 4. Use q=-1 ->5, three slices.
    depth = 4 if mode == "physical" else 5
    count = 2 if depth == 4 else 3
    pars = SamplePreparationParameters(200., (SliceAppliedParameters((0.,), 0., 1.),), 0., None,
        (0,)*count, ((0, 0),), "centered", "fft", "legacy", "unity")
    manual = simulate_tilt_series_serial(prepare_multislice_sample(s.model, pars), generate_scan_tilts(s.scan))
    expected = compute_legacy_objective(manual.intensities, data(), iref=1.)
    trial = np.array(vector)
    first = evaluate(s, trial, data(), coordinate_mode=mode)
    second = evaluate(s, trial, data(), coordinate_mode=mode)
    np.testing.assert_array_equal(first.amplitudes, manual.amplitudes)
    np.testing.assert_array_equal(first.amplitudes, second.amplitudes)
    assert first.objective == expected.objective and first.scaling == expected.scaling
    assert s.regions[0].depth_angstrom == 2
    np.testing.assert_array_equal(trial, vector)
    with pytest.raises(ValueError):
        first.amplitudes.setflags(write=True)


def test_nonunit_origin_reflection_complex_reference():
    s = state()
    s = replace(s, reflections=((RefinementReflection(0, 0, 1., 1., "YN"),),))
    result = evaluate(s, np.array([2.]), data())
    # Source-checked Al origin 281.953595001 / area 200; interaction phase
    # per dz=2 from independently checked phase-grating test values.
    expected = cmath.exp(-1j*.0007288360746343612*1.409767975005*2)
    np.testing.assert_allclose(result.amplitudes, expected, rtol=0, atol=1e-15)


@pytest.mark.parametrize("sym,point", [("INVE", (-1, 0)), ("RECT", (-1, 0)),
    ("SQUA", (0, 1)), ("HEXA", (0, -1))])
def test_symmetry_first_exact_match_and_duplicate_precedence(sym, point):
    grid = ReciprocalIndexGrid(4, 4, centered=True)
    refs = (RefinementReflection(1, 0, 2, 3, "YY"), RefinementReflection(1, 0, 8, 9, "YY"),
            RefinementReflection(*point, 4, 5, "YY"))
    es, ab = _reflection_scales(grid, refs, sym)
    assert es[2, 3] == 2 and ab[2, 3] == 3
    assert es[point[1]+2, point[0]+2] == 4  # exact beats earlier symmetry match
    es, ab = _reflection_scales(grid, refs[:2], sym)
    assert es[point[1]+2, point[0]+2] == 2


def test_gradient_literal_stencil():
    v = np.array([[1.], [2.], [3.], [4.]])
    np.testing.assert_array_equal(_gradient(v, 2, 2), [[.125], [-1.], [-2.125], [-3.25]])
    result = evaluate(state(), np.empty(0), data(1))
    np.testing.assert_array_equal(result.objective_intensities, _gradient(result.calculated_intensities, 2, 2))


def test_rejection_does_not_corrupt_baseline():
    s = state()
    bad = replace(s, regions=(replace(s.regions[0], vacancy_values=(.1,)),))
    with pytest.raises(ValueError, match="hazards"):
        evaluate(bad, np.empty(0), data())
    assert np.isfinite(evaluate(s, np.empty(0), data()).objective)
    assert bad.regions[0].vacancy_values == (.1,)


def test_multiple_regions_species_slices_source_order():
    s = state()
    model = replace(s.model, slices=s.model.slices*2, species=(Species("Al", 13), Species("Si", 14)),
        atoms_by_slice=((replace(s.model.atoms_by_slice[0][0], y_mask=1),),
                        (replace(s.model.atoms_by_slice[0][0], x_mask=1),)))
    r0 = replace(s.regions[0], thickness_mask=1, vacancy_masks=(1,), debye_waller=(.1, .2), debye_waller_masks=(0, 1))
    r1 = replace(r0, vacancy_flag=2, vacancy_values=(.2, .3), vacancy_masks=(1, 0), debye_waller_masks=(1, 1))
    s = replace(s, model=model, regions=(r0, r1), reflections=((), ()), symmetry=("NONE", "NONE"))
    assert [(p.category, p.target) for p in layout(s).parameters] == [
        ("depth", (0,)), ("depth", (1,)), ("constant_vacancy", (0, 0)),
        ("vacancy_coefficient", (1, 0)), ("debye_waller", (0, 1)),
        ("debye_waller", (1, 0)), ("debye_waller", (1, 1)),
        ("atom_coordinate", (0, 0, 1)), ("atom_coordinate", (1, 0, 0))]
    out = apply(s, np.arange(1., 10.))
    assert out.regions[0].debye_waller == (.1, 5.)
    assert out.regions[1].debye_waller == (6., 7.)
    assert out.regions[1].vacancy_values == (4., .3)


@pytest.mark.parametrize("value,expected", [(0., 1e-10), (-1e-11, 1e-10), (-1e-10, -1e-10)])
def test_pack_near_zero_boundary_without_baseline_mutation(value, expected):
    s = state()
    s = replace(s, scan=replace(s.scan, corner_a=(value, 0)), scan_masks=("YN", "NN", "NN"))
    assert pack(s)[0] == expected
    assert s.scan.corner_a[0] == value
    assert apply(s, np.array([value])).scan.corner_a[0] == value  # no repack on apply


def test_trial_all_scaling_branches_and_direct_fft():
    s = state()
    for mode in (0., -1., 1.25):
        sample = replace(s, iref=mode)
        fft = evaluate(sample, np.empty(0), data())
        direct = evaluate(sample, np.empty(0), data(), phase_grating_method="direct")
        expected = compute_legacy_objective(fft.objective_intensities, data(), iref=mode)
        assert fft.objective == expected.objective and fft.scaling == expected.scaling
        np.testing.assert_allclose(fft.amplitudes, direct.amplitudes, rtol=0, atol=1e-15)


def test_trial_preflight_rejections():
    s = state()
    with pytest.raises(ValueError, match="coordinate_mode"):
        evaluate(s, np.empty(0), data(), coordinate_mode="hidden")
    with pytest.raises(ValueError, match="dimensions"):
        evaluate(replace(s, scan=replace(s.scan, scan_rows=3)), np.empty(0), data())
    with pytest.raises(ValueError, match="G-atom"):
        m = replace(s.model, atoms_by_slice=((replace(s.model.atoms_by_slice[0][0], record_type="G"),),))
        evaluate(replace(s, model=m), np.empty(0), data())
    with pytest.raises(ValueError, match="surface"):
        evaluate(replace(s, regions=(replace(s.regions[0], free_surface=(1, 0, 0, 1)),)), np.empty(0), data())
    with pytest.raises(ValueError, match="VOLTS"):
        m = replace(s.model, absorption_raw="YES", energy_refined=True)
        evaluate(replace(s, model=m, absorption=((0., 1., 1.),), absorption_masks=((0, 0, 0),)), np.array([300.]), data())


def test_parser_range_error_has_record_context(tmp_path):
    source = Path(__file__).parents[1]/"tests/fixtures/legacy_reference"/"FOR007.DAT"
    lines = source.read_text().splitlines()
    lines[-1] = "NaN 0.95"
    target = tmp_path/"bad.dat"
    target.write_text("\n".join(lines))
    with pytest.raises(ValueError, match="bad.dat: section 'parameter 18 range factors', record"):
        parse_refinement_state(target, experimental_header_path=source.parent/"in.dat")


def test_region_sequence_cycles_until_depth_is_reached():
    s = state()
    model = replace(s.model, slices=tuple(replace(s.model.slices[0], thickness_angstrom=v) for v in (2., 3., 5.)),
                    atoms_by_slice=s.model.atoms_by_slice*3)
    regions = (replace(s.regions[0], depth_angstrom=4., slice_count=2),
               replace(s.regions[0], depth_angstrom=11.))
    s = replace(s, model=model, regions=regions, reflections=((), ()), symmetry=("NONE", "NONE"))
    parameters = SamplePreparationParameters(200., (SliceAppliedParameters((0.,), 0., 1.),)*3,
        0., None, (0, 1, 2, 2, 2), ((0, 0),), "centered", "fft", "legacy", "unity")
    manual = simulate_tilt_series_serial(prepare_multislice_sample(model, parameters), generate_scan_tilts(s.scan))
    np.testing.assert_array_equal(evaluate(s, np.empty(0), data()).amplitudes, manual.amplitudes)


def test_reference_problem_short_trial_supported_path():
    # Keep all 18 source-selected categories but shorten geometry/scan for a
    # focused integration test. This is not a reference-executable comparison.
    original = parse_refinement_state(Path(__file__).parents[1]/"tests/fixtures/legacy_reference"/"FOR007.DAT")
    s = replace(original, dsap=.1,
        regions=(replace(original.regions[0], depth_angstrom=2.0245),),
        scan=replace(original.scan, scan_rows=2, scan_columns=2))
    trial = pack(s).copy()
    trial[10] = 1.2
    result = evaluate(s, trial, data(1))
    assert result.amplitudes.shape == (4, 1)
    assert result.amplitudes.dtype == np.complex128
    assert np.isfinite(result.objective)
    assert original.regions[0].depth_angstrom == 840.
    assert original.reflections[0][0].elastic == 1.
