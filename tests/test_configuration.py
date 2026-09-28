"""Canonical compiler equivalence against the unchanged legacy input path."""
from dataclasses import FrozenInstanceError, replace
import builtins
import importlib.util
import io
from pathlib import Path

import numpy as np
import pytest

from qcbed import compile_config
from qcbed.configuration import (
    ParameterCategory as C, ParameterTarget, ParameterBinding, BoundsConfig,
    AbsorptionConfig, VacancyConfig, VacancyKind, SpeciesConfig, ThermalFactorConfig,
    LegacyExportConfig, MicroscopeConfig, TargetBound, RelativeBound,
)
from qcbed.refinement_driver import parse_legacy_refinement
from qcbed.refinement import (pack_refinement_parameters, build_refinement_parameter_layout,
                              apply_refinement_parameters)
from qcbed.refinement_subset import select_refinement_categories
from qcbed.refinement_bounds import RefinementBoundRules, resolve_refinement_bounds
from qcbed.refinement_trial import _prepare_trial, evaluate_refinement_trial
from qcbed.experimental import parse_experimental_intensities
from qcbed.scan import generate_scan_tilts

ROOT = Path(__file__).resolve().parents[1]


def example(name):
    spec = importlib.util.spec_from_file_location(name, ROOT/'examples'/f'{name}.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def raw():
    return example('reference_configuration').reference_configuration(legacy_header=True)


def legacy_full():
    runner = example('run_full_refinement')
    state, controls = parse_legacy_refinement(ROOT/'reference_data/FOR007.DAT', seed=2026)
    state, controls = select_refinement_categories(state, controls, runner.REFINEMENT_CATEGORIES)
    rules = RefinementBoundRules(runner.SPECIFIC_BOUNDS, runner.CATEGORY_BOUNDS,
        runner.CATEGORY_RELATIVE_BOUNDS, runner.CATEGORY_RELATIVE_FALLBACK)
    settings = runner.make_config(controls, refinement_bounds=resolve_refinement_bounds(state, rules),
                                  bound_rules=rules)
    return state, settings


def test_complete_reference_equivalence(raw):
    state, runner = legacy_full()
    result = compile_config(raw)
    assert result.state == state  # Includes raw fields, atoms, geometry and masks.
    assert result.runner == runner
    assert result.optimizer == runner.refinement
    assert result.execution == runner.refinement.execution
    assert result.bounds == runner.refinement_bounds
    assert result.bound_rules == runner.bound_rules
    assert result.layout == build_refinement_parameter_layout(state)
    assert result.range_factors == state.ranges == ((1.05, .95),)*7
    assert result.plamda == (.5,) + (.0005,)*6
    np.testing.assert_array_equal(result.packed_starting_vector, pack_refinement_parameters(state))
    np.testing.assert_array_equal(generate_scan_tilts(result.state.scan), generate_scan_tilts(state.scan))
    data = parse_experimental_intensities(ROOT/'reference_data/in.dat',
        beam_indices_hk=result.experiment.beam_indices_hk, scan=result.state.scan)
    assert data.gradient_option == result.experiment.gradient_option == 1
    assert result.experiment.beam_indices_hk == ((0,0),(1,1),(1,-1),(-1,-1),(-1,1),(2,0),(0,-2),(-2,0),(0,2))
    assert result.configuration.legacy_export.dsmax_inverse_angstrom == .5


def test_full_18_parameter_reference_and_binding_order(raw):
    state, controls = parse_legacy_refinement(ROOT/'reference_data/FOR007.DAT', seed=2026)
    raw = replace(raw, refinement=replace(raw.refinement, active_categories=None))
    forward = compile_config(raw)
    backward = compile_config(replace(raw, refinement=replace(raw.refinement,
        parameters=tuple(reversed(raw.refinement.parameters)))))
    assert forward.state == backward.state
    assert forward.state == state
    assert forward.layout == backward.layout == build_refinement_parameter_layout(state)
    assert forward.plamda == backward.plamda == controls.plamda
    assert forward.range_factors == backward.range_factors == state.ranges
    np.testing.assert_array_equal(forward.packed_starting_vector,
        [840,1,1,.86,-.508156,.407641,.407641,.508156,-.407641,-.508156,*([1.]*8)])
    assert [(r.h,r.k) for r in forward.state.reflections[0]][-1] == (2,0)
    assert len(forward.state.reflections[0]) == 8


def test_pure_compile_and_deep_immutability(raw, monkeypatch):
    raw = replace(raw, refinement=replace(raw.refinement, parameters=list(raw.refinement.parameters)))
    before = repr(raw)
    def forbidden(*args, **kwargs):
        raise AssertionError('compiler attempted filesystem I/O')
    monkeypatch.setattr(builtins, 'open', forbidden)
    monkeypatch.setattr(io, 'open', forbidden)
    resolved = compile_config(raw)
    assert repr(raw) == before
    assert isinstance(raw.refinement.parameters, tuple)
    with pytest.raises(FrozenInstanceError):
        raw.title = 'changed'
    with pytest.raises(FrozenInstanceError):
        resolved.state = None
    with pytest.raises(ValueError):
        resolved.packed_starting_vector.setflags(write=True)


def test_duplicate_beams_and_distinct_duplicate_reflection_targets(raw):
    experiment = replace(raw.experiment, beam_indices_hk=((1,1),(0,0),(1,1)))
    raw = replace(raw, experiment=experiment, refinement=replace(raw.refinement, active_categories=None))
    r = compile_config(raw)
    assert r.experiment.beam_indices_hk == ((1,1),(0,0),(1,1))
    entries = [p for p in r.layout.parameters if p.category == 'reflection_multiplier']
    assert entries[4].target == (0,4,0) and entries[7].target == (0,7,0)


def test_starting_values_not_projected_and_linked_absorption(raw):
    region = replace(raw.sample.regions[0], depth_angstrom=-840.,
        debye_waller=(ThermalFactorConfig('al', -.86),))
    coeff = replace(raw.sample.absorption.coefficients[0], b=-1., c=-7., c_mode='linked_to_b')
    sample = replace(raw.sample, regions=(region,), dilation=-2.,
        absorption=replace(raw.sample.absorption, coefficients=(coeff,)))
    params = tuple(p for p in raw.refinement.parameters if p.target != ParameterTarget(C.ABSORPTION, ('al','c')))
    raw = replace(raw, sample=sample, refinement=replace(raw.refinement,
        parameters=params, active_categories=None))
    resolved = compile_config(raw)
    s = resolved.state
    assert s.regions[0].depth_angstrom == -840.
    assert s.regions[0].debye_waller == (-.86,)
    assert s.dilation == -2.
    assert s.absorption == ((0., -1., -7.),)
    assert s.absorption_masks == ((0,1,-1),)
    mapped = apply_refinement_parameters(s, resolved.packed_starting_vector)
    assert mapped.regions[0].depth_angstrom == 840.
    assert mapped.absorption == ((0.,1.,1.),)
    assert resolved.bounds.entries[0].lower == 600.
    assert resolved.configuration.sample.regions[0].depth_angstrom == -840.


@pytest.mark.parametrize('kind', ['vacancy', 'gas'])
def test_unsupported_physics_still_rejected_at_existing_boundary(raw, kind):
    region = raw.sample.regions[0]
    if kind == 'vacancy':
        region = replace(region, vacancy=VacancyConfig(VacancyKind.LEGACY_PROFILE_1, (.1,.2)))
    else:
        region = replace(region, slices=(replace(region.slices[0], record_type='G'),))
    config = replace(raw, sample=replace(raw.sample, regions=(region,)))
    state = compile_config(config).state
    with pytest.raises(ValueError, match='unsupported'):
        _prepare_trial(state, original_voltage=200.)


def test_multiple_species_regions_and_stored_slice_order(raw):
    first = raw.sample.regions[0]
    stored = first.slices[0]
    # Map by stable species ID, regardless of thermal declaration order.
    thermals = (ThermalFactorConfig('al', .86), ThermalFactorConfig('si', .3))
    first = replace(first, debye_waller=thermals)
    second_slice = replace(stored, id='second_slice', atoms=(replace(stored.atoms[1], species='si'),))
    second = replace(first, id='second', slices=(second_slice,), reflections=(), depth_angstrom=10.)
    sample = replace(raw.sample, species=(SpeciesConfig('si','SI',14), *raw.sample.species),
        regions=(first,second), absorption=AbsorptionConfig())
    params = (ParameterBinding(ParameterTarget(C.DEBYE_WALLER,('second','al')), (.9,1.2), .07),
              ParameterBinding(ParameterTarget(C.DEPTH,('second',)), (1.3,.8), .09),
              ParameterBinding(ParameterTarget(C.ATOM_COORDINATE,('second_slice','a1','y')), (1.1,.9), .02))
    config = replace(raw, sample=sample, refinement=replace(raw.refinement,
        parameters=params, active_categories=None, bounds=BoundsConfig()))
    r = compile_config(config)
    assert [s.atomic_number for s in r.state.model.species] == [14,13]
    assert r.state.regions[0].debye_waller == (.3,.86)
    assert r.state.regions[1].debye_waller_masks == (0,1)
    assert r.state.model.atoms_by_slice[0][0].species_index_one_based == 2
    assert r.state.model.atoms_by_slice[1][0].species_index_one_based == 1
    assert [(p.category,p.target) for p in r.layout.parameters] == [
        ('depth',(1,)), ('debye_waller',(1,1)), ('atom_coordinate',(1,0,1))]
    assert r.plamda == (.09,.07,.02)
    assert r.range_factors == ((1.3,.8),(.9,1.2),(1.1,.9))


@pytest.mark.parametrize('change', ['duplicate_binding','unknown_target','duplicate_id','linked_refined','unknown_species','raw_flag'])
def test_bad_references_and_conflicts(raw, change):
    if change == 'duplicate_binding':
        raw = replace(raw, refinement=replace(raw.refinement, parameters=raw.refinement.parameters*2))
    elif change == 'unknown_target':
        binding = replace(raw.refinement.parameters[0], target=ParameterTarget(C.DEPTH,('absent',)))
        raw = replace(raw, refinement=replace(raw.refinement,parameters=(binding,)))
    elif change == 'duplicate_id':
        raw = replace(raw, sample=replace(raw.sample, regions=raw.sample.regions*2))
    elif change == 'unknown_species':
        region = raw.sample.regions[0]
        raw = replace(raw, sample=replace(raw.sample, regions=(replace(region,
            debye_waller=(ThermalFactorConfig('missing', .1),)),)))
    else:
        absorption = raw.sample.absorption
        if change == 'linked_refined':
            absorption = replace(absorption, coefficients=(replace(absorption.coefficients[0], c_mode='linked_to_b'),))
        else:
            absorption = replace(absorption, legacy_raw='yes')
        raw = replace(raw, sample=replace(raw.sample,absorption=absorption))
    with pytest.raises(ValueError):
        compile_config(raw)


@pytest.mark.parametrize('value', [True, float('nan'), float('inf')])
def test_strict_numbers(value):
    with pytest.raises((ValueError,TypeError)):
        MicroscopeConfig(value, (0,0,1))
    with pytest.raises((ValueError,TypeError)):
        SpeciesConfig('s','S',value)


def test_small_forward_equivalence_without_new_numerical_path(raw):
    state, _ = legacy_full()
    compiled = compile_config(raw)
    # Restrict only this test's trial: same existing evaluator for both paths.
    state = replace(state, scan=replace(state.scan,scan_rows=2,scan_columns=2))
    other = replace(compiled.state, scan=state.scan)
    from qcbed.experimental import ExperimentalIntensityData
    data = ExperimentalIntensityData(2,2,0,compiled.experiment.beam_indices_hk,
        np.ones((4,9)), np.ones((4,9)), np.ones((4,9)))
    vector = compiled.packed_starting_vector.copy()
    vector[0] = 4.
    a = evaluate_refinement_trial(state,vector,data,execution=compiled.execution)
    b = evaluate_refinement_trial(other,vector,data,execution=compiled.execution)
    for field in ('amplitudes','calculated_intensities','objective_intensities'):
        np.testing.assert_array_equal(getattr(a,field),getattr(b,field))
    assert a.objective == b.objective and a.scaling == b.scaling


def test_remaining_categories_and_bounds_precedence(raw):
    # Every category has a source-index mapping, including unsupported physics
    # that must remain representable without being executed during compilation.
    region = raw.sample.regions[0]
    region = replace(region, vacancy=VacancyConfig(VacancyKind.LEGACY_PROFILE_2, (-.2,.3)))
    definitions = ((C.VOLTAGE,()), (C.DEPTH,('specimen',)),
        (C.VACANCY_COEFFICIENT,('specimen','a')), (C.VACANCY_COEFFICIENT,('specimen','b')),
        (C.DILATION,()), (C.ABSORPTION,('al','a')), (C.ABSORPTION,('al','b')),
        (C.ABSORPTION,('al','c')), (C.DEBYE_WALLER,('specimen','al')),
        (C.ATOM_COORDINATE,('al_slice','a0','x')), (C.ATOM_COORDINATE,('al_slice','a1','y')),
        (C.SCAN_COORDINATE,('B','second')), (C.REFLECTION_MULTIPLIER,('specimen','r0','absorptive')))
    bindings = tuple(ParameterBinding(ParameterTarget(c,t),(1.05,.95),float(i))
                     for i,(c,t) in enumerate(definitions))
    config = replace(raw, sample=replace(raw.sample,regions=(region,)),
        refinement=replace(raw.refinement,parameters=tuple(reversed(bindings)),active_categories=None,
            bounds=BoundsConfig(specific=(TargetBound(ParameterTarget(C.DEPTH,('specimen',)),700.,900.),),
                relative=(RelativeBound(C.DEPTH,.5,.1),))))
    resolved = compile_config(config)
    assert resolved.plamda == tuple(float(i) for i in range(len(bindings)))
    assert resolved.state.model.energy_refined
    assert resolved.state.regions[0].vacancy_values == (-.2,.3)
    assert resolved.state.regions[0].vacancy_masks == (1,1)
    assert resolved.state.dilation_mask == 1
    assert resolved.state.absorption_masks == ((1,1,1),)
    assert resolved.state.model.atoms_by_slice[0][0].x_mask == 1
    assert resolved.state.model.atoms_by_slice[0][1].y_mask == 1
    assert resolved.state.scan_masks == ('NN','NY','NN')
    assert resolved.state.reflections[0][0].mask == 'NY'
    assert (resolved.bounds.entries[0].lower,resolved.bounds.entries[0].upper) == (700.,900.)
    constant = replace(region,vacancy=VacancyConfig(values=(-.3,)))
    binding = ParameterBinding(ParameterTarget(C.CONSTANT_VACANCY,('specimen',)),(1.05,.95),.2)
    config = replace(config,sample=replace(config.sample,regions=(constant,)),
        refinement=replace(config.refinement,parameters=(binding,),bounds=BoundsConfig()))
    resolved = compile_config(config)
    assert resolved.state.regions[0].vacancy_values == (-.3,)
    assert resolved.state.regions[0].vacancy_masks == (1,)


def test_empty_layout_and_optional_export_metadata(raw):
    config = replace(raw,legacy_export=LegacyExportConfig(),refinement=replace(raw.refinement,
        parameters=(),active_categories=None,optimizer=replace(raw.refinement.optimizer,loops=0),bounds=BoundsConfig()))
    resolved = compile_config(config)
    assert resolved.plamda == resolved.range_factors == resolved.layout.parameters == ()
    assert resolved.packed_starting_vector.shape == (0,)
