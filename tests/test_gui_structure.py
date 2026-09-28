"""Canonical structural transactions/history without Qt."""
from dataclasses import replace
from pathlib import Path
import pytest
from qcbed import (Id, get_object, available_refinement_targets, load_yaml_config,
                   compile_config, prepare_run_snapshot)
from qcbed.config_editing import DependencyError
from qcbed.configuration import (SpeciesConfig, ThermalFactorConfig, ReflectionConfig,
    AbsorptionSpeciesConfig, ParameterTarget, ParameterBinding, ParameterCategory as C,
    TargetBound, CategoryBound, RelativeBound)
from qcbed.gui.document import EditorDocument
from qcbed.gui.structure import StructuralDraft, removal_plan, resolved_bound_rows
from test_gui_run import run_config_fixture


@pytest.fixture
def document(tmp_path):
    return EditorDocument(run_config_fixture(tmp_path),tmp_path/'config.yml')


def region_path(config): return ('sample','regions',Id(config.sample.regions[0].id))


def test_species_transaction_dependencies_history_save(document,tmp_path):
    original = document.config
    draft = StructuralDraft(original)
    path = region_path(original)
    draft.add(('sample','species'),SpeciesConfig('cr','Cr',24))
    assert not document.commit_structure(draft.config)  # Thermal completion required.
    assert not document.undo_stack
    draft.add(path+('debye_waller',),ThermalFactorConfig('cr',.003))
    atom = original.sample.regions[0].slices[0].atoms[0]
    atom_path = path+('slices',Id(original.sample.regions[0].slices[0].id),'atoms',Id(atom.id))
    draft.update(atom_path,species='cr',x_fractional=-1.123456789012345,starting_occupancy=.8123456789012345)
    draft.move(('sample','species',Id('cr')),-1)
    assert document.commit_structure(draft.config)
    assert document.config.sample.species[0].id == 'cr'
    assert get_object(document.config,atom_path).x_fractional == -1.123456789012345
    target = ParameterTarget(C.ATOM_COORDINATE,(original.sample.regions[0].slices[0].id,atom.id,'x'))
    # Use the catalogue's authoritative coordinate tokens.
    target = next(t.target for t in available_refinement_targets(document.config)
                  if t.category == C.ATOM_COORDINATE and atom.id in t.target.identifiers)
    draft = StructuralDraft(document.config)
    binding = ParameterBinding(target,(1.1,.9),.1)
    draft.add(('refinement','parameters'),binding)
    draft.add(('refinement','bounds','specific'),TargetBound(target,-2.,2.))
    assert document.commit_structure(draft.config)
    before_delete = document.config
    plan = removal_plan(before_delete,('sample','species',Id('cr')))
    assert any(d.field == 'species' for d in plan)
    assert any('parameters' in d.object_path for d in plan)
    assert any('specific' in d.object_path for d in plan)
    with pytest.raises(DependencyError): draft.delete(('sample','species',Id('cr')))
    assert draft.config == before_delete
    draft.delete(('sample','species',Id('cr')),remove_dependencies=True)
    assert document.commit_structure(draft.config)
    assert not document.config.refinement.parameters
    assert not document.config.refinement.bounds.specific
    assert document.undo()
    assert document.config == before_delete
    assert document.redo()
    after_delete = document.config
    assert document.undo()
    document.save(tmp_path/'structural.yml')
    assert load_yaml_config(tmp_path/'structural.yml') == before_delete
    assert document.redo()
    assert document.config == after_delete


def test_regions_slices_reflections_targets_and_snapshot(document):
    original = document.config
    snapshot = prepare_run_snapshot(original,base_directory=document.base_directory)
    draft = StructuralDraft(original)
    path = region_path(original)
    draft.duplicate(path)
    second = draft.config.sample.regions[1]
    second_path = ('sample','regions',Id(second.id))
    draft.update(second_path,depth_angstrom=7.25,symmetry='INVE')
    draft.duplicate(second_path+('slices',Id(second.slices[0].id)))
    second = get_object(draft.config,second_path)
    last_slice = second.slices[-1]
    draft.update(second_path+('slices',Id(last_slice.id),'geometry'),thickness_angstrom=3.5)
    draft.move(second_path+('slices',Id(last_slice.id)),-1)
    draft.add(second_path+('reflections',),ReflectionConfig('r1',1,1,.9,1.))
    draft.duplicate(second_path+('reflections',Id('r1')))
    draft.move(second_path,-1)
    assert document.commit_structure(draft.config)
    assert document.config.sample.regions[1] == original.sample.regions[0]
    assert document.config.sample.regions[0].slices[0].id == last_slice.id
    assert [r.h for r in document.config.sample.regions[0].reflections] == [1,1]
    catalogue = available_refinement_targets(document.config)
    assert sum(t.category == C.DEPTH for t in catalogue) == 2
    assert sum(t.category == C.DEBYE_WALLER for t in catalogue) == 2
    assert sum(t.category == C.REFLECTION_MULTIPLIER for t in catalogue) == 4
    assert snapshot.compiled.configuration == original
    draft.delete(second_path)
    assert document.commit_structure(draft.config)
    assert document.config == original


def test_bounds_provenance_and_linked_absorption(document):
    draft = StructuralDraft(document.config)
    path = region_path(document.config)
    species = document.config.sample.species[0].id
    draft.update(('sample','absorption'),enabled=True)
    draft.add(('sample','absorption','coefficients'),AbsorptionSpeciesConfig(species,.1,.2,.3))
    target = ParameterTarget(C.ABSORPTION,(species,'c'))
    draft.add(('refinement','parameters'),ParameterBinding(target,(1.1,.9),.01))
    assert document.commit_structure(draft.config)
    draft.update(('sample','absorption','coefficients',Id(species,'species')),c_mode='linked_to_b')
    assert not document.commit_structure(draft.config)
    draft.delete(('refinement','parameters',Id(target,'target')))
    depth = ParameterTarget(C.DEPTH,(path[-1].value,))
    draft.add(('refinement','parameters'),ParameterBinding(depth,(1.1,.9),.1))
    draft.add(('refinement','bounds','relative'),RelativeBound(C.DEPTH,.1,.01))
    draft.add(('refinement','bounds','absolute'),CategoryBound(C.DEPTH,1.,10.))
    draft.add(('refinement','bounds','specific'),TargetBound(depth,2.,8.))
    draft.update(('refinement','bounds'),policy='reject')
    assert document.commit_structure(draft.config)
    rows = resolved_bound_rows(document.config)
    assert rows[0][1:3] == (2.,8.)
    assert 'specific target' in rows[0][3]
    catalogue = available_refinement_targets(document.config)
    assert not next(t for t in catalogue if t.target == target).available
    assert any('Absorption entry' in d.message for d in removal_plan(document.config,('sample','species',Id(species))))
    draft.delete(('refinement','bounds','specific',Id(depth,'target')))
    assert document.commit_structure(draft.config)
    assert resolved_bound_rows(document.config)[0][1:3] == (1.,10.)
    draft.delete(('refinement','bounds','absolute',Id(C.DEPTH,'category')))
    assert document.commit_structure(draft.config)
    assert 'category relative' in resolved_bound_rows(document.config)[0][3]


def test_common_edits_history_and_invalid_text(document):
    before = document.config
    assert not document.apply({('scan','scan_rows'):'bad'})
    assert not document.undo_stack
    assert document.apply({('title',):'changed'})
    assert document.undo()
    assert document.config == before
    assert document.redo()
    assert document.config.title == 'changed'
    assert document.undo()
    assert document.apply({('title',):'branch'})
    assert not document.redo_stack


def test_species_duplicate_and_dynamic_thermal_binding(document):
    from qcbed import set_refinement
    draft = StructuralDraft(document.config)
    species = document.config.sample.species[0]
    draft.duplicate(('sample','species',Id(species.id)))
    copied = draft.config.sample.species[-1]
    assert copied.id != species.id
    assert copied.atomic_number == species.atomic_number
    path = region_path(document.config)
    draft.add(path+('debye_waller',),ThermalFactorConfig(copied.id,.009))
    assert document.commit_structure(draft.config)
    config = document.config
    for info in available_refinement_targets(config):
        if info.category == C.DEBYE_WALLER:
            config = set_refinement(config,info.target,enabled=True,range_factors=(1.1,.9),plamda=.001)
    assert document.commit_structure(config)
    assert compile_config(document.config).parameter_count == 2
    original_targets = tuple(b.target for b in config.refinement.parameters)
    draft = StructuralDraft(config)
    draft.move(('sample','species',Id(copied.id)),-1)
    assert document.commit_structure(draft.config)
    assert tuple(b.target for b in document.config.refinement.parameters) == original_targets
    assert compile_config(document.config).parameter_count == 2


def test_unbounded_parameter_does_not_shift_bound_labels(document):
    from qcbed import set_refinement
    config = document.config
    targets = available_refinement_targets(config)
    voltage = next(t.target for t in targets if t.category == C.VOLTAGE)
    depth = next(t.target for t in targets if t.category == C.DEPTH)
    for target in (voltage,depth):
        config = set_refinement(config,target,enabled=True,range_factors=(1.1,.9),plamda=.1)
    draft = StructuralDraft(config)
    draft.add(('refinement','bounds','specific'),TargetBound(depth,2.,8.))
    rows = resolved_bound_rows(draft.config)
    assert rows[0][1:3] == (None,None)
    assert rows[1][1:3] == (2.,8.)
    assert 'depth' in rows[1][0]
