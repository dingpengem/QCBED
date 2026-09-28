"""Frontend APIs operate on canonical snapshots with stable physical identities."""
from dataclasses import replace, FrozenInstanceError
import math

import numpy as np
import pytest

from qcbed import (compile_config, validate_config, validate_fields, ValidationPhase as Phase,
    available_refinement_targets, Id, get_object, add_object, delete_object,
    deletion_dependencies, DependencyError, duplicate_object, reorder_objects,
    update_object, replace_object, set_refinement)
from qcbed.configuration import (MicroscopeConfig, ParameterTarget, ParameterCategory as C,
    BoundsConfig, TargetBound, ThermalFactorConfig, SpeciesConfig)
from test_dynamic_configuration import small_config, reflection_config
from test_multiregion_configuration import configuration
from test_configuration import example


def test_field_diagnostics_before_construction():
    issues = validate_fields(MicroscopeConfig,{'energy_kv':math.nan,'zone_axis':(True,0,1),'other':0},
                             section='microscope',object_path=('microscope',))
    assert {v.field for v in issues} == {'energy_kv','zone_axis','other'}
    assert all(v.phase == Phase.FIELD and v.section == 'microscope' for v in issues)
    assert validate_fields(MicroscopeConfig,{'zone_axis':(0,0,1)})[0].code == 'missing_field'
    assert validate_config({})[0].code == 'invalid_config_type'
    assert validate_config(small_config(),phases=(Phase.FIELD,)) == ()


def test_duplicate_ids_have_both_scoped_locations():
    raw = configuration(((2.,),(3.,)),(4.,6.))
    r0,r1 = raw.sample.regions
    # Same atom ID across different slices is valid; same slice ID globally is not.
    assert validate_config(raw) == ()
    r1 = replace(r1,slices=(replace(r1.slices[0],id=r0.slices[0].id),))
    raw = replace(raw,sample=replace(raw.sample,regions=(r0,r1),species=raw.sample.species*2))
    issues = validate_config(raw)
    duplicates = [v for v in issues if v.code == 'duplicate_id']
    assert len(duplicates) == 2
    assert all(v.field == 'id' and v.related_locations for v in duplicates)
    slice_issue = next(v for v in duplicates if 'slices' in v.object_path)
    assert Id(r1.id) in slice_issue.object_path
    assert Id(r0.id) in slice_issue.related_locations[0].object_path
    species_issue = next(v for v in duplicates if 'species' in v.object_path)
    assert species_issue.location.occurrence == 1
    assert species_issue.related_locations[0].occurrence == 0


def test_reference_and_target_errors_locate_objects():
    raw = small_config()
    region = raw.sample.regions[0]; stored = region.slices[0]
    atom = replace(stored.atoms[0],species='unknown')
    bad = replace(raw,sample=replace(raw.sample,regions=(replace(region,slices=(replace(stored,atoms=(atom,)),)),)))
    issue = next(v for v in validate_config(bad) if v.code == 'unknown_species')
    assert issue.field == 'species' and Id(atom.id) in issue.object_path and Id(stored.id) in issue.object_path
    from qcbed.configuration import ParameterBinding
    target = ParameterTarget(C.DEPTH,('missing',))
    bad = replace(raw,refinement=replace(raw.refinement,parameters=(ParameterBinding(target,(1.1,.9),.1),)))
    issue = validate_config(bad)[0]
    assert issue.code == 'unknown_target' and target in issue.object_path


def test_execution_and_export_phases_are_distinct():
    raw = small_config()
    raw = replace(raw,scan=replace(raw.scan,scan_rows=1))
    assert validate_config(raw,phases=(Phase.COMPILE,)) == ()
    issues = validate_config(raw)
    assert [(i.code,i.phase,i.field) for i in issues] == [('singleton_scan',Phase.EXECUTION,'scan_rows')]
    export = validate_config(raw,phases=(Phase.EXPORT,))
    assert any(i.code == 'missing_dsmax' and i.phase == Phase.EXPORT for i in export)
    assert any(i.code == 'legacy_omission' and i.severity == 'warning' for i in export)


@pytest.mark.parametrize('count',[0,1,8,23])
def test_catalogue_dynamic_and_matches_all_compiler_targets(count):
    raw = reflection_config(count)
    infos = available_refinement_targets(raw)
    assert len(infos) == 15+2*count
    assert len({v.target for v in infos}) == len(infos)
    assert sum(v.selected for v in infos) == sum(v.active for v in infos) == 2*count
    assert all(v.label and v.units and math.isfinite(v.starting_value) for v in infos)
    assert available_refinement_targets(compile_config(raw)) == infos
    # Activate each supported target through the canonical API; compiler resolves it.
    base = small_config()
    for info in available_refinement_targets(base):
        if info.available:
            changed = set_refinement(base,info.target,enabled=True,range_factors=(1.1,.9),plamda=.01)
            assert compile_config(changed).parameter_count == 1
        else:
            assert info.reason


def test_catalogue_multiregion_and_filtered_binding_state():
    raw = configuration(((2.,),(1.,3.),(1.,2.,4.)),(4.,5.,8.),(0,2))
    infos = available_refinement_targets(raw)
    assert len(infos) == 41
    depth = [v for v in infos if v.category == C.DEPTH]
    assert [v.starting_value for v in depth] == [4.,5.,8.]
    assert [v.selected for v in depth] == [True,False,True]
    reference = example('reference_configuration').reference_configuration()
    infos = available_refinement_targets(reference)
    assert sum(v.selected for v in infos) == 18
    assert sum(v.active for v in infos) == 7


def test_linked_and_gas_targets_remain_visible_with_reasons():
    raw = example('reference_configuration').reference_configuration()
    coeff = replace(raw.sample.absorption.coefficients[0],c=-7.,c_mode='linked_to_b')
    target = ParameterTarget(C.ABSORPTION,('al','c'))
    raw = replace(raw,sample=replace(raw.sample,absorption=replace(raw.sample.absorption,coefficients=(coeff,))),
        refinement=replace(raw.refinement,parameters=tuple(b for b in raw.refinement.parameters if b.target != target)))
    info = next(v for v in available_refinement_targets(raw) if v.target == target)
    assert not info.available and 'linked' in info.reason and info.starting_value == -7.
    with pytest.raises(ValueError,match='linked'):
        set_refinement(raw,target,enabled=True,range_factors=(1.1,.9),plamda=.1)
    raw = small_config()
    r = raw.sample.regions[0]
    raw = replace(raw,sample=replace(raw.sample,regions=(replace(r,slices=(replace(r.slices[0],record_type='G'),)),)))
    infos = [v for v in available_refinement_targets(raw) if v.category == C.ATOM_COORDINATE]
    assert infos and all(not v.available and v.reason for v in infos)
    assert any(i.code == 'unsupported_gas_atoms' for i in validate_config(raw))


def test_pure_add_update_replace_and_reorder():
    raw = configuration(((2.,),(3.,)),(4.,6.),(0,1))
    original = repr(raw)
    reordered = reorder_objects(raw,('sample','regions'),('region-1','region-0'))
    assert [r.id for r in reordered.sample.regions] == ['region-1','region-0']
    assert reordered.refinement == raw.refinement
    np.testing.assert_array_equal(compile_config(reordered).packed_starting_vector,[6.,4.])
    path = ('sample','regions',Id('region-0'))
    updated = update_object(raw,path,depth_angstrom=9.)
    assert updated.sample.regions[1] == raw.sample.regions[1]
    assert get_object(updated,path).depth_angstrom == 9.
    with pytest.raises(ValueError,match='IDs'):
        replace_object(raw,path,replace(raw.sample.regions[0],id='renamed'))
    with pytest.raises(ValueError,match='permutation'):
        reorder_objects(raw,('sample','regions'),('region-0','region-0'))
    assert repr(raw) == original
    with pytest.raises(FrozenInstanceError): updated.sample.regions[0].depth_angstrom = 3.


def test_deletion_reports_bounds_and_bindings_then_explicit_resolution():
    raw = configuration(((2.,),(3.,)),(4.,6.),(0,))
    target = ParameterTarget(C.DEPTH,('region-0',))
    raw = replace(raw,refinement=replace(raw.refinement,bounds=BoundsConfig(specific=(TargetBound(target,1.,10.),))))
    path = ('sample','regions',Id('region-0'))
    deps = deletion_dependencies(raw,path)
    assert len(deps) == 2
    with pytest.raises(DependencyError) as error: delete_object(raw,path)
    assert error.value.dependencies == deps
    changed = set_refinement(raw,target,enabled=False)
    changed = update_object(changed,('refinement','bounds'),specific=())
    changed = delete_object(changed,path)
    assert [r.id for r in changed.sample.regions] == ['region-1']
    assert compile_config(changed).parameter_count == 0


def test_duplication_remaps_only_owned_ids_and_optional_bindings():
    raw = configuration(((2.,),),(4.,),(0,))
    path = ('sample','regions',Id('region-0'))
    copied = duplicate_object(raw,path,copy_refinement=True)
    first,second = copied.sample.regions
    assert first == raw.sample.regions[0]
    assert second.id != first.id and second.slices[0].id != first.slices[0].id
    assert second.slices[0].atoms[0].id != first.slices[0].atoms[0].id
    assert second.slices[0].atoms[0].species == first.slices[0].atoms[0].species
    assert compile_config(copied).parameter_count == 2
    assert compile_config(duplicate_object(raw,path)).parameter_count == 1
    assert duplicate_object(raw,path) == duplicate_object(raw,path)  # deterministic pure ID allocation


def test_atom_reflection_and_species_operations():
    raw = reflection_config(2)
    r = raw.sample.regions[0]; s = r.slices[0]
    rp = ('sample','regions',Id(r.id))
    ap = rp+('slices',Id(s.id),'atoms')
    copied = duplicate_object(raw,ap+(Id(s.atoms[0].id),))
    assert len(get_object(copied,ap)) == 3
    fresh = get_object(copied,ap)[-1]
    assert delete_object(copied,ap+(Id(fresh.id),)) == raw
    copied = duplicate_object(raw,rp+('reflections',Id(r.reflections[0].id)),copy_refinement=True)
    assert len(copied.sample.regions[0].reflections) == 3
    assert compile_config(copied).parameter_count == 6
    assert copied.sample.regions[0].reflections[-1].h == r.reflections[0].h
    species = SpeciesConfig('new','SI',14)
    added = add_object(raw,('sample','species'),species)
    assert any(i.code == 'thermal_species_mismatch' for i in validate_config(added))
    added = update_object(added,rp,debye_waller=r.debye_waller+(ThermalFactorConfig('new',.2),))
    added = update_object(added,ap+(Id(s.atoms[0].id),),species='new')
    assert not validate_config(added)
    deps = deletion_dependencies(added,('sample','species',Id('new')))
    assert {d.field for d in deps} == {'species'} and len(deps) == 2


def test_binding_changes_require_explicit_values_and_preserve_filters():
    raw = small_config(); target = ParameterTarget(C.DEPTH,(raw.sample.regions[0].id,))
    with pytest.raises(ValueError,match='explicit'):
        set_refinement(raw,target,enabled=True)
    changed = set_refinement(raw,target,enabled=True,range_factors=(1.05,.95),plamda=-.1)
    same = set_refinement(changed,target,enabled=True)
    assert same == changed
    assert set_refinement(same,target,enabled=False) == raw
    with pytest.raises(ValueError,match='unknown'):
        set_refinement(raw,ParameterTarget(C.DEPTH,('absent',)),enabled=True,range_factors=(1.,2.),plamda=.1)


def test_target_limit_is_execution_readiness_not_editor_capacity():
    raw = reflection_config(1003,both=False)
    assert compile_config(raw).parameter_count == 1003
    assert sum(i.selected for i in available_refinement_targets(raw)) == 1003
    assert any(i.code == 'optimizer_parameter_count' and i.phase == Phase.EXECUTION for i in validate_config(raw))


def test_duplicate_region_copies_atom_and_bound_targets_by_scope():
    from qcbed.configuration import ParameterBinding
    raw = configuration(((2.,),(3.,)),(4.,6.))
    r = raw.sample.regions[0]; s = r.slices[0]; a = s.atoms[0]
    target = ParameterTarget(C.ATOM_COORDINATE,(s.id,a.id,'x'))
    raw = replace(raw,refinement=replace(raw.refinement,
        parameters=(ParameterBinding(target,(1.1,.9),.1),),
        bounds=BoundsConfig(specific=(TargetBound(target,-1.,1.),))))
    result = duplicate_object(raw,('sample','regions',Id(r.id)),copy_refinement=True)
    copy = result.sample.regions[-1]
    copied_target = ParameterTarget(C.ATOM_COORDINATE,(copy.slices[0].id,copy.slices[0].atoms[0].id,'x'))
    assert result.refinement.parameters[-1].target == copied_target
    assert result.refinement.bounds.specific[-1].target == copied_target
    assert compile_config(result).parameter_count == 2
    assert result.sample.regions[1] == raw.sample.regions[1]  # Same atom IDs in other scope unchanged.
    with pytest.raises(ValueError,match='globally unique'):
        add_object(raw,('sample','regions',Id(raw.sample.regions[1].id),'slices'),s)
    info = next(v for v in available_refinement_targets(result) if v.target == copied_target)
    assert info.stable_ids == (copy.slices[0].id,copy.slices[0].atoms[0].id)
    assert 'x' not in info.stable_ids


def test_duplicates_in_atom_scope_and_metadata_mismatch():
    raw = small_config(); r = raw.sample.regions[0]; s = r.slices[0]
    bad = replace(raw,sample=replace(raw.sample,regions=(replace(r,slices=(replace(s,atoms=s.atoms+s.atoms[:1]),)),)))
    issue = next(i for i in validate_config(bad) if i.code == 'duplicate_id')
    assert issue.field == 'id' and 'atoms' in issue.object_path
    assert issue.location.occurrence == 2 and issue.related_locations[0].occurrence == 0
    raw = replace(raw,legacy_export=replace(raw.legacy_export,title_raw='other'))
    issue = next(i for i in validate_config(raw) if i.code == 'legacy_metadata_mismatch')
    assert issue.section == 'legacy_export' and issue.field == 'title_raw'


def test_local_constructor_domain_errors_are_structured():
    from qcbed.models import MeshConfig
    issues = validate_fields(MeshConfig,{'mesh_x':0,'mesh_y':2},section='sample',object_path=('sample','mesh'))
    assert len(issues) == 1
    assert issues[0].code == 'invalid_field_value' and issues[0].phase == Phase.FIELD
    assert issues[0].object_path == ('sample','mesh')
    assert 'positive' in issues[0].message


def test_nonpositive_fixed_voltage_execution_diagnostic():
    raw = small_config()
    raw = replace(raw,microscope=replace(raw.microscope,energy_kv=-200.))
    assert validate_config(raw,phases=(Phase.COMPILE,)) == ()
    issue = next(v for v in validate_config(raw) if v.code == 'nonpositive_voltage')
    assert issue.phase == Phase.EXECUTION and issue.field == 'energy_kv'
    assert raw.microscope.energy_kv == -200.
