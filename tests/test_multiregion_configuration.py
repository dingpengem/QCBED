"""Independent region depths across compilation, MSEQ and runner persistence."""
from dataclasses import replace
import json

import numpy as np
import pytest

from qcbed import compile_config
from qcbed.configuration import ParameterBinding, ParameterTarget, ParameterCategory as C, BoundsConfig, TargetBound
from qcbed.refinement import apply_refinement_parameters, pack_refinement_parameters
from qcbed.refinement_trial import _prepare_trial
from qcbed.refinement_driver import run_legacy_refinement
from qcbed.full_runner import run_full_refinement, load_refinement_checkpoint
from test_dynamic_configuration import small_config, experimental
from test_full_runner import cheap_render


def configuration(patterns, depths, active=()):
    raw = small_config()
    template = raw.sample.regions[0]
    regions = tuple(replace(template, id=f'region-{i}', depth_angstrom=depth,
        slices=tuple(replace(template.slices[0], id=f'slice-{i}-{j}',
            geometry=replace(template.slices[0].geometry, thickness_angstrom=step))
            for j, step in enumerate(pattern)))
        for i, (pattern, depth) in enumerate(zip(patterns, depths)))
    bindings = tuple(ParameterBinding(ParameterTarget(C.DEPTH, (regions[i].id,)),
        (1.1+i*.01, .9), .1+i*.01) for i in reversed(active))
    return replace(raw, sample=replace(raw.sample, regions=regions),
        refinement=replace(raw.refinement, parameters=bindings))


def sequence(state):
    return _prepare_trial(state, original_voltage=state.model.energy_kv)[1]


@pytest.mark.parametrize('patterns,depths,active,expected', [
    (((2.,),), (4.,), (0,), (0,0)),
    (((2.,),(3.,)), (4.,6.), (), (0,0,1,1)),
    (((2.,),(3.,)), (4.,6.), (1,), (0,0,1,1)),
    (((2.,),(3.,)), (4.,6.), (0,1), (0,0,1,1)),
    (((2.,),(1.,3.),(1.,2.,4.)), (4.,5.,8.), (0,2), (0,0,1,2,1,3,4,5,3)),
])
def test_region_patterns_and_independent_mapping(patterns, depths, active, expected):
    compiled = compile_config(configuration(patterns, depths, active))
    assert sequence(compiled.state) == expected
    assert tuple(r.depth_angstrom for r in compiled.state.regions) == depths
    assert tuple(r.slice_count for r in compiled.state.regions) == tuple(map(len, patterns))
    assert [p.target for p in compiled.layout.parameters] == [(i,) for i in active]
    assert compiled.parameter_count == len(active)
    np.testing.assert_array_equal(compiled.packed_starting_vector, [depths[i] for i in active])
    assert compiled.range_factors == tuple((1.1+i*.01,.9) for i in active)
    assert compiled.plamda == tuple(.1+i*.01 for i in active)
    for position, region in enumerate(active):
        vector = compiled.packed_starting_vector.copy()
        vector[position] += .25
        mapped = apply_refinement_parameters(compiled.state, vector)
        assert mapped.regions[region].depth_angstrom == depths[region]+.25
        for other in range(len(depths)):
            if other != region:
                assert mapped.regions[other] == compiled.state.regions[other]
        assert mapped.model == compiled.state.model
        np.testing.assert_array_equal(pack_refinement_parameters(mapped), vector)


def test_dynamic_region_count():
    compiled = compile_config(configuration(((2.,),)*23, tuple(2.*(i+1) for i in range(23)), tuple(range(23))))
    assert compiled.parameter_count == 23
    assert sequence(compiled.state) == tuple(i for i in range(23) for _ in range(i+1))


@pytest.mark.parametrize('region', [0,1,2])
def test_exact_sequence_boundary_changes_only_selected_region(region):
    compiled = compile_config(configuration(((2.,),(1.,3.),(1.,2.,4.)), (4.,4.,7.), (0,1,2)))
    baseline = ((0,0),(1,2),(3,4,5))
    offsets = (0,1,3,6)
    for direction in (-np.inf,np.inf):
        vector = compiled.packed_starting_vector.copy()
        vector[region] = np.nextafter(vector[region],direction)
        mapped = apply_refinement_parameters(compiled.state,vector)
        parts = tuple(tuple(k for k in sequence(mapped) if offsets[i] <= k < offsets[i+1]) for i in range(3))
        expected = list(baseline)
        if direction == np.inf:
            expected[region] += (offsets[region],)
        assert parts == tuple(expected)
        assert sequence(mapped) == sum(parts, ())
        for other in range(3):
            if other != region:
                assert mapped.regions[other] == compiled.state.regions[other]


def test_independent_stable_id_bounds():
    raw = configuration(((2.,),(3.,)), (4.,6.), (0,1))
    raw = replace(raw, refinement=replace(raw.refinement, bounds=BoundsConfig('reject', specific=(
        TargetBound(ParameterTarget(C.DEPTH,('region-1',)),5.,8.),
        TargetBound(ParameterTarget(C.DEPTH,('region-0',)),3.,5.)))))
    compiled = compile_config(raw)
    assert {(b.target,b.lower,b.upper) for b in compiled.bounds.entries} == {((0,),3.,5.),((1,),5.,8.)}
    for values, target in (([2.,6.],(0,)),([4.,9.],(1,))):
        violations = compiled.bounds.violations(apply_refinement_parameters(compiled.state,np.array(values)))
        assert [v.target for v in violations] == [target]


@pytest.mark.parametrize('active', [(1,), (0,1)])
def test_callbacks_reports_checkpoints_and_resume(active, cheap_render, tmp_path):
    compiled = compile_config(configuration(((2.,),(1.,3.)), (4.,7.), active))
    controls = replace(compiled.optimizer, loops=1)
    observed = []
    def observe(kind,stage,name,loop,payload):
        if kind == 'evaluation':
            observed.append(payload)
    direct = run_legacy_refinement(compiled.state,experimental(),controls,progress_callback=observe)
    assert len(observed) == direct.evaluations and observed
    for row in observed:
        assert len(row.physical_state.regions) == 2
        if active == (1,):
            assert row.physical_state.regions[0] == compiled.state.regions[0]
    folder = tmp_path/'run'
    settings = replace(compiled.runner,refinement=controls)
    result = run_full_refinement(compiled.state,experimental(),settings,output_directory=folder)
    cp = load_refinement_checkpoint(folder)
    assert [p['target'] for p in cp['layout']['parameters']] == [[i] for i in active]
    assert len(set(cp['parameter_labels'])) == len(active)
    mapped = apply_refinement_parameters(compiled.state,np.array(cp['retained_vector']))
    np.testing.assert_array_equal(pack_refinement_parameters(mapped),cp['retained_vector'])
    if active == (1,):
        assert mapped.regions[0] == compiled.state.regions[0]
    summary = json.loads((folder/'summary.json').read_text())
    assert len(summary['parameter_changes']) == len(active)
    saved = json.loads((folder/'final_parameters.json').read_text())
    assert [r['depth_angstrom'] for r in saved['regions']] == [r.depth_angstrom for r in mapped.regions]
    assert cp['parameter_labels'] == [f'{j+1}: depth {(i,)}' for j,i in enumerate(active)]
    resumed = run_full_refinement(compiled.state,experimental(),settings,output_directory=folder,resume=True)
    assert resumed.status == result.status == 'completed'
    assert resumed.summary['final_objective'] == result.summary['final_objective']
    assert load_refinement_checkpoint(folder)['retained_vector'] == cp['retained_vector']


def test_interrupted_resume_retains_all_region_depths(cheap_render, tmp_path):
    compiled = compile_config(configuration(((2.,),(1.,3.),(1.,2.,4.)), (4.,7.,8.), (0,2)))
    settings = replace(compiled.runner,refinement=replace(compiled.optimizer,loops=2))
    def stop(row):
        if row['classification'] == 'stage_retained':
            raise KeyboardInterrupt
    folder = tmp_path/'interrupted'
    first = run_full_refinement(compiled.state,experimental(),settings,
        output_directory=folder,progress_hook=stop)
    assert first.status == 'interrupted'
    resumed = run_full_refinement(compiled.state,experimental(),settings,
        output_directory=folder,resume=True)
    direct_folder = tmp_path/'direct'
    direct = run_full_refinement(compiled.state,experimental(),settings,output_directory=direct_folder)
    assert resumed.status == direct.status == 'completed'
    np.testing.assert_array_equal(np.load(folder/'final_active_vector.npy'),
                                  np.load(direct_folder/'final_active_vector.npy'))
    saved = json.loads((folder/'final_parameters.json').read_text())
    assert saved == json.loads((direct_folder/'final_parameters.json').read_text())
    assert saved['regions'][1]['depth_angstrom'] == 7.
    changed = replace(compiled.state,regions=(compiled.state.regions[0],
        replace(compiled.state.regions[1],depth_angstrom=9.),compiled.state.regions[2]))
    with pytest.raises(ValueError,match='fingerprint mismatch'):
        run_full_refinement(changed,experimental(),settings,output_directory=folder,resume=True)
