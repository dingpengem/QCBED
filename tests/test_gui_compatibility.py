"""GUI-independent Stage 5 adapters retain canonical state and numerical fidelity."""
from dataclasses import replace
from pathlib import Path
import numpy as np
import pytest
from qcbed import compile_config, load_yaml_config, load_experimental_data
from qcbed.gui.compatibility import import_document, export_document, compatibility_report
from qcbed.gui.document import EditorDocument
from qcbed.configuration import LegacyExportConfig,AbsorptionConfig,AbsorptionSpeciesConfig
from qcbed.refinement_io import _parse_refinement_input
from qcbed.refinement_trial import evaluate_refinement_trial
from test_legacy_config import policy, equal_state, ROOT
from test_dynamic_configuration import reflection_config
from test_multiregion_configuration import configuration


def test_import_unsaved_export_roundtrip_and_readonly_reports(tmp_path):
    source = ROOT/'tests/fixtures/legacy_reference/FOR007.DAT'
    experimental = ROOT/'tests/fixtures/legacy_reference/in.dat'
    original_data = experimental.read_bytes()
    document = import_document(source,experimental,policy(tmp_path))
    assert document.dirty and document.needs_save_as
    before = document.config
    history = (list(document.undo_stack),list(document.redo_stack))
    report = compatibility_report(before,base_directory=document.base_directory)
    assert not report['errors']
    assert {'execution','seed','bounds','run','objective_mode','invalid_trial_policy'} <= report['omissions'].keys()
    export_document(document,tmp_path/'standalone.dat',filename_limits=(4096,4096))
    export_document(document,tmp_path/'bundle',bundle=True)
    assert (tmp_path/'bundle/compatibility.json').exists()
    assert experimental.read_bytes() == original_data
    assert document.config is before
    assert (document.undo_stack,document.redo_stack) == history
    assert document.dirty
    migrated = import_document(tmp_path/'bundle/FOR007.dat',tmp_path/'bundle/in.dat',policy(tmp_path))
    first,second = compile_config(before),compile_config(migrated.config)
    equal_state(second.state,first.state)
    assert second.plamda == first.plamda
    assert second.range_factors == first.range_factors
    vector = first.packed_starting_vector.copy(); vector[0] = 4.
    results = []
    for compiled in (first,second):
        data = load_experimental_data(compiled,base_directory=tmp_path)
        results.append(evaluate_refinement_trial(compiled.state,vector,data,execution=compiled.execution,objective_mode='masked'))
    for field in ('amplitudes','calculated_intensities','objective_intensities'):
        np.testing.assert_array_equal(getattr(results[0],field),getattr(results[1],field))
    assert results[0].scaling == results[1].scaling
    assert results[0].objective == results[1].objective
    document.save(tmp_path/'migrated.yml')
    assert not document.needs_save_as and not document.dirty
    assert load_yaml_config(tmp_path/'migrated.yml') == document.config


@pytest.mark.parametrize('count',[0,3,11])
def test_dynamic_multiregion_linked_export(tmp_path,count):
    configs = (reflection_config(count),configuration(((2.,),(1.,3.),(1.,2.,4.)),(4.,7.,8.),(0,2)))
    for i,raw in enumerate(configs):
        body = tmp_path/f'body{i}.dat'; body.write_text('1 1\n1 1\n'*3)
        raw = replace(raw,experiment=replace(raw.experiment,path=str(body),beam_indices_hk=((0,0),)),
            legacy_export=LegacyExportConfig(.12345678901234567),
            sample=replace(raw.sample,absorption=AbsorptionConfig(True,
                (AbsorptionSpeciesConfig('al',.1,-2.,-7.,'linked_to_b'),))))
        document = EditorDocument(raw,tmp_path/'current.yml')
        export_document(document,tmp_path/f'bundle{i}',bundle=True)
        imported = import_document(tmp_path/f'bundle{i}/FOR007.dat',tmp_path/f'bundle{i}/in.dat',policy(tmp_path))
        equal_state(compile_config(imported.config).state,compile_config(raw).state)
        assert imported.config.sample.absorption.coefficients[0].c_mode == 'linked_to_b'
        assert imported.config.legacy_export.dsmax_inverse_angstrom == raw.legacy_export.dsmax_inverse_angstrom
        assert document.config is raw


def test_export_protection_and_unexportable_report(tmp_path):
    document = import_document(ROOT/'tests/fixtures/legacy_reference/FOR007.DAT',ROOT/'tests/fixtures/legacy_reference/in.dat',policy(tmp_path))
    existing = tmp_path/'existing.dat'; existing.write_text('keep me')
    with pytest.raises(FileExistsError): export_document(document,existing,filename_limits=(4096,4096))
    assert existing.read_text() == 'keep me'
    with pytest.raises(FileExistsError): export_document(document,tmp_path,bundle=True)
    no_metadata = replace(document.config,legacy_export=LegacyExportConfig())
    report = compatibility_report(no_metadata,base_directory=tmp_path)
    assert any(i.code == 'missing_dsmax' for i in report['errors'])
    report = compatibility_report(document.config,base_directory=tmp_path,bundle=False)
    assert report['errors']  # Never guess filename capacities.


def test_ambiguous_input_rejected(tmp_path):
    source = tmp_path/'FOR007.dat'
    source.write_text((ROOT/'tests/fixtures/legacy_reference/FOR007.DAT').read_text()+'\nUNKNOWN LEGACY SETTING\n')
    with pytest.raises(ValueError,match='trailing'):
        import_document(source,ROOT/'tests/fixtures/legacy_reference/in.dat',policy(tmp_path))


def test_missing_canonical_data_path_is_not_an_export_destination(tmp_path):
    document = import_document(ROOT/'tests/fixtures/legacy_reference/FOR007.DAT',ROOT/'tests/fixtures/legacy_reference/in.dat',policy(tmp_path))
    missing = tmp_path/'missing-input.dat'
    assert document.apply({('experiment','path'):str(missing)})
    with pytest.raises(ValueError,match='canonical experimental'):
        export_document(document,missing,filename_limits=(4096,4096))
    assert not missing.exists()
