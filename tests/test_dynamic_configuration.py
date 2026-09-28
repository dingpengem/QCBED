"""Variable physical instances, dynamic vector plumbing and source ordering."""
from dataclasses import replace
import json

import numpy as np
import pytest

from qcbed import compile_config
from qcbed.configuration import (ParameterCategory as C, ParameterTarget, ParameterBinding,
    ReflectionConfig, ThermalFactorConfig, SpeciesConfig, AbsorptionConfig, BoundsConfig,
    CategoryBound, LegacyExportConfig)
from qcbed.models import MeshConfig, SliceGeometry
from qcbed.refinement import apply_refinement_parameters
from qcbed.refinement_categories import PARAMETER_SEMANTICS
from qcbed.refinement_driver import run_legacy_refinement
from qcbed.full_runner import run_full_refinement, load_refinement_checkpoint
from qcbed.experimental import ExperimentalIntensityData
from test_configuration import example
from test_full_runner import cheap_render


def small_config():
    raw = example('reference_configuration').reference_configuration()
    region = raw.sample.regions[0]
    stored = replace(region.slices[0], geometry=SliceGeometry(2.,0.,0.,20.,10.,90.))
    region = replace(region,depth_angstrom=4.,slices=(stored,),reflections=())
    return replace(raw, sample=replace(raw.sample,mesh=MeshConfig(4,2),
        regions=(region,),absorption=AbsorptionConfig()),
        scan=replace(raw.scan,scan_rows=2,scan_columns=2),
        experiment=replace(raw.experiment,beam_indices_hk=((0,0),),gradient_option=0),
        refinement=replace(raw.refinement,parameters=(),active_categories=None,
            bounds=BoundsConfig(),optimizer=replace(raw.refinement.optimizer,
                method='ANNE',loops=0,maxfsimp=0,maxfpara=0)),
        legacy_export=LegacyExportConfig())


def reflection_config(count, both=True):
    raw = small_config()
    region = raw.sample.regions[0]
    records = tuple(ReflectionConfig(f'f{i}',0,0,1.,1.) for i in range(count))
    params = tuple(ParameterBinding(ParameterTarget(C.REFLECTION_MULTIPLIER,
        (region.id,r.id,component)), (1.05+i*.00001,.95), .0001*(i+1))
        for i,r in enumerate(records) for component in (('elastic','absorptive') if both else ('elastic',)))
    return replace(raw,sample=replace(raw.sample,regions=(replace(region,reflections=records),)),
        refinement=replace(raw.refinement,parameters=tuple(reversed(params)),
            bounds=BoundsConfig('reject',absolute=(CategoryBound(C.REFLECTION_MULTIPLIER,.5,1.5),))))


@pytest.mark.parametrize('count', [0,1,3,11])
def test_reflection_instances_drive_all_vector_lengths(count):
    compiled = compile_config(reflection_config(count))
    n = 2*count
    assert compiled.parameter_count == n
    assert len(compiled.layout.parameters) == len(compiled.plamda) == len(compiled.range_factors) == n
    assert len(compiled.bounds.entries) == n
    assert compiled.packed_starting_vector.shape == (n,)
    assert [p.target for p in compiled.layout.parameters] == [(0,i,j) for i in range(count) for j in range(2)]
    assert len(compiled.state.reflections[0]) == count  # duplicate H/K retained
    values = 1.+np.arange(n)*.001
    mapped = apply_refinement_parameters(compiled.state,values)
    assert [(r.elastic,r.absorptive) for r in mapped.reflections[0]] == list(zip(values[::2],values[1::2]))
    assert compiled.plamda == tuple(.0001*(i+1) for i in range(count) for _ in range(2))


def test_compiler_does_not_inherit_optimizer_guard():
    compiled = compile_config(reflection_config(1003,both=False))
    assert len(compiled.layout.parameters) == 1003
    assert compiled.packed_starting_vector.shape == (1003,)
    assert len(compiled.bounds.entries) == 1003
    with pytest.raises(ValueError,match='1000'):
        run_legacy_refinement(compiled.state,experimental(),replace(compiled.optimizer,loops=1))


@pytest.mark.parametrize('region_count,species_count,subset', [(1,2,False),(3,4,False),(3,4,True)])
def test_thermal_targets_follow_region_species_not_slots(region_count,species_count,subset):
    raw = small_config()
    species = tuple(SpeciesConfig(f's{i}',f'{i:02}',13+i) for i in range(species_count))
    regions = []
    bindings = []
    expected = []
    for r in range(region_count):
        original = raw.sample.regions[0]
        stored = replace(original.slices[0],id=f'slice{r}',atoms=tuple(
            replace(a,species='s0') for a in original.slices[0].atoms))
        thermal = tuple(ThermalFactorConfig(sp.id,.2+i*.1) for i,sp in enumerate(species))
        regions.append(replace(original,id=f'r{r}',slices=(stored,),debye_waller=tuple(reversed(thermal))))
        for i,sp in enumerate(species):
            if subset and (r+i)%2:
                continue
            expected.append((r,i))
            bindings.append(ParameterBinding(ParameterTarget(C.DEBYE_WALLER,(f'r{r}',sp.id)),(1.1,.9),.003))
    raw = replace(raw,sample=replace(raw.sample,species=species,regions=tuple(regions)),
        refinement=replace(raw.refinement,parameters=tuple(reversed(bindings)),
            bounds=BoundsConfig(absolute=(CategoryBound(C.DEBYE_WALLER,.01,2.),))))
    compiled = compile_config(raw)
    assert [p.target for p in compiled.layout.parameters] == expected
    assert len(compiled.bounds.entries) == len(expected)
    values = np.arange(len(expected),dtype=np.float64)*.01+.4
    mapped = apply_refinement_parameters(compiled.state,values)
    for (r,i),value in zip(expected,values):
        assert mapped.regions[r].debye_waller[i] == value
    assert sum(sum(r.debye_waller_masks) for r in compiled.state.regions) == len(expected)
    assert len(compiled.plamda) == len(compiled.range_factors) == len(expected)


def experimental():
    return ExperimentalIntensityData(2,2,0,((0,0),),np.ones((4,1)),np.ones((4,1)),np.ones((4,1)))


@pytest.mark.parametrize('method', ['ANNE','SIMP','PARA','SPAR'])
@pytest.mark.parametrize('count', [1,5])
def test_optimizer_callbacks_use_dynamic_layout(method,count):
    compiled = compile_config(reflection_config(count,both=False))
    controls = replace(compiled.optimizer,method=method,loops=1)
    observed = []
    def observe(kind,stage,name,loop,payload):
        if kind == 'evaluation':
            observed.append(payload)
    result = run_legacy_refinement(compiled.state,experimental(),controls,progress_callback=observe,
        bounds_policy='reject',refinement_bounds=compiled.bounds)
    assert result.initial_vector.shape == result.final_vector.shape == (count,)
    assert len(observed) == result.evaluations
    assert all(len(row.physical_state.reflections[0]) == count for row in observed)


@pytest.mark.parametrize('count', [0,5])
def test_checkpoint_and_report_lengths_follow_layout(count,cheap_render,tmp_path):
    compiled = compile_config(reflection_config(count,both=False))
    # loops=0 is the existing supported empty-layout/evaluation-only path.
    folder = tmp_path/f'count{count}'
    result = run_full_refinement(compiled.state,experimental(),compiled.runner,output_directory=folder)
    checkpoint = load_refinement_checkpoint(folder)
    assert len(checkpoint['layout']['parameters']) == count
    assert len(checkpoint['parameter_labels']) == count
    assert len(checkpoint['retained_vector']) == count
    summary = json.loads((folder/'summary.json').read_text())
    assert len(summary['parameter_changes']) == count
    resumed = run_full_refinement(compiled.state,experimental(),compiled.runner,
        output_directory=folder,resume=True)
    assert resumed.status == result.status


def test_taxonomy_is_not_a_vector_length():
    assert {c.value for c in C} == set(PARAMETER_SEMANTICS)
    raw = example('reference_configuration').reference_configuration()
    for selected,n in ((raw,7),(replace(raw,refinement=replace(raw.refinement,active_categories=None)),18)):
        compiled = compile_config(selected)
        assert len(compiled.layout.parameters) == n
        ranks = [PARAMETER_SEMANTICS[p.category].source_order_group for p in compiled.layout.parameters]
        assert ranks == sorted(ranks)
