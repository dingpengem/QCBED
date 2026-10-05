"""Editor behaviour without importing Qt or requiring the GUI extra."""
from dataclasses import replace
from pathlib import Path
import subprocess
import sys

import pytest

from qcbed import Id, available_refinement_targets, load_yaml_config, update_object
from qcbed.gui.document import EditorDocument, editable_fields, display_value, target_label

REFERENCE = Path(__file__).resolve().parent / 'data/reference_config.yml'


@pytest.fixture
def document():
    return EditorDocument.open(REFERENCE)


def test_reference_unchanged_roundtrip(document, tmp_path):
    original = document.config
    from qcbed import get_object
    drafts = {s.key: display_value(getattr(get_object(original, s.path), s.name)) for s in editable_fields(original)}
    assert document.apply(drafts)
    assert document.config == original
    assert document.summary()['active_parameters'] == 7
    destination = tmp_path / 'saved.yml'
    document.save(destination)
    assert load_yaml_config(destination) == original
    assert not document.dirty
    document.save(overwrite=True)
    assert load_yaml_config(destination) == original


@pytest.mark.parametrize('key,text', [(('scan','scan_rows'),'2.0'),
    (('microscope','energy_kv'),'nan'), (('scan','corner_a'),'1, 2, 3'),
    (('microscope','zone_axis'),'1, 0.5, 0')])
def test_invalid_drafts_preserve_snapshot(document, key, text):
    original = document.config
    assert not document.apply({key:text})
    assert document.config is original
    assert document.issues[0].severity == 'error'
    assert document.issues[0].field == key[-1]


def test_region_scan_path_and_unedited_fields(document, tmp_path):
    old = document.config
    region = old.sample.regions[0]
    assert document.apply({('sample','regions',Id(region.id),'depth_angstrom'):'901.1234567890123',
        ('scan','scan_rows'):'3', ('scan','scan_columns'):'4',
        ('scan','basis_1'):'1, 0, 0', ('scan','corner_a'):'-0.5123456789012345, 0.4',
        ('experiment','path'):str(tmp_path/'data.dat')})
    new = document.config
    assert new.sample.regions[0].depth_angstrom == 901.1234567890123
    assert new.sample.regions[0].slices == region.slices
    assert new.sample.species == old.sample.species
    assert new.refinement == old.refinement
    assert new.execution == old.execution
    assert new.scan.corner_a == (-0.5123456789012345, 0.4)
    assert document.summary()['tilts'] == 12
    assert new.experiment.path == str(tmp_path/'data.dat')
    document.validate(check_files=True)
    assert any(i.phase.value == 'files' for i in document.issues)


def test_binding_selection_and_explicit_settings(document):
    original = document.config
    binding = original.refinement.parameters[0]
    target = binding.target
    assert document.apply({}, [(target, False, '', '', '')])
    assert all(b.target != target for b in document.config.refinement.parameters)
    assert not document.apply({}, [(target, True, '', '', '')])
    assert document.apply({}, [(target, True, '-0.25', '0.125', '0.012345678901234567')])
    result = next(b for b in document.config.refinement.parameters if b.target == target)
    assert result.range_factors == (-0.25, 0.125)
    assert result.plamda == 0.012345678901234567
    assert document.config.refinement.active_categories == original.refinement.active_categories


def test_readiness_distinct_from_draft_validity(document):
    assert document.apply({('scan','scan_rows'):'1'})
    assert document.summary()['execution'] == 'Blocked'
    assert not document.apply({('scan','scan_columns'):'invalid'})
    assert document.summary()['execution'] == 'Blocked'


def test_dynamic_regions_catalogue(document):
    from qcbed import duplicate_object
    config = duplicate_object(document.config, ('sample','regions',Id(document.config.sample.regions[0].id)))
    # The public helper returns a configuration snapshot with fresh nested IDs.
    other = EditorDocument(config, REFERENCE)
    depths = [s for s in editable_fields(config) if s.name == 'depth_angstrom']
    assert len(depths) == 2
    assert depths[0].key != depths[1].key
    assert other.apply({depths[1].key:'123.456'})
    assert other.config.sample.regions[0] == config.sample.regions[0]
    catalogue = available_refinement_targets(other.config)
    assert sum(i.category.value == 'depth' for i in catalogue) == 2
    assert all(target_label(other.config, i) for i in catalogue)


def test_gui_package_and_help_do_not_import_qt():
    result = subprocess.run([sys.executable, '-c',
        'import sys; import qcbed.gui; import qcbed.gui.document; assert "PySide6" not in sys.modules'],
        capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    result = subprocess.run([sys.executable, '-m', 'qcbed.gui', '--help'], capture_output=True, text=True)
    assert result.returncode == 0


def test_other_common_fields_and_raw_metadata(document, tmp_path):
    old = document.config
    assert document.apply({('title',):'New project', ('microscope','energy_kv'):'201.12345678901234',
        ('microscope','zone_axis'):'1, 1, 0', ('scan','basis_2'):'0, 1, 1',
        ('scan','corner_c'):'-0.2, -0.6', ('execution','execution_backend'):'serial',
        ('execution','batch_size'):'', ('run','output_directory'):str(tmp_path/'out'),
        ('run','output_filename'):str(tmp_path/'result.dat')})
    new = document.config
    assert new.title == 'New project'
    assert new.microscope.energy_kv == 201.12345678901234
    assert new.microscope.zone_axis == (1,1,0)
    assert new.execution.execution_backend == 'serial'
    assert new.execution.batch_size is None
    assert new.sample == old.sample
    assert new.refinement == old.refinement
    assert new.legacy_export == replace(old.legacy_export, title_raw=None, output_filename_raw=None)
    assert not document.apply({('execution','execution_backend'):'batched', ('execution','batch_size'):''})
    assert document.apply({('execution','execution_backend'):'batched', ('execution','batch_size'):'13'})
    assert document.config.execution.batch_size == 13
    document.save(tmp_path/'common.yml')
    assert load_yaml_config(tmp_path/'common.yml') == document.config


def test_relative_paths_use_document_base(document, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert document.apply({('experiment','path'):'relative-data.dat'})
    assert document.config.experiment.path == str(REFERENCE.parent/'relative-data.dat')
    document.save(tmp_path/'relocated.yml')
    assert load_yaml_config(tmp_path/'relocated.yml') == document.config


def test_invalid_yaml_is_a_presentable_error(tmp_path):
    path = tmp_path/'broken.yml'
    path.write_text('sample: [unterminated')
    with pytest.raises(ValueError, match='Invalid YAML'):
        EditorDocument.open(path)


def test_invalid_path_draft_does_not_escape_or_commit(document):
    original = document.config
    assert not document.apply({('experiment','path'):'bad\x00path'})
    assert document.config is original
    assert document.issues


def test_existing_output_is_readiness_only(document, tmp_path):
    from qcbed.gui.document import diagnostic_presentation
    from qcbed import validate_config, ValidationPhase
    assert document.apply({('run','output_directory'):str(tmp_path)})
    document.validate(check_files=True)
    assert not any(i.code=='output_exists' for i in document.issues)
    assert not validate_config(document.config, phases=(ValidationPhase.FIELD, ValidationPhase.COMPILE))
    summary = document.summary()
    assert summary['configuration'].startswith('Valid')
    assert summary['compile'].startswith('Valid')
    assert summary['execution'] != 'Blocked'
    document.save(tmp_path/'editable.yml')
    assert load_yaml_config(tmp_path/'editable.yml') == document.config
