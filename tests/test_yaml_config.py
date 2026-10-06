from dataclasses import asdict, replace
from enum import Enum
from pathlib import Path

import numpy as np
import pytest
import yaml

from qcbed import compile_config, load_yaml_config, load_experimental_data
from qcbed.refinement_trial import evaluate_refinement_trial
from test_configuration import example
from test_multiregion_configuration import configuration, sequence
from test_dynamic_configuration import reflection_config

ROOT = Path(__file__).resolve().parents[1]
REFERENCE_YAML = Path(__file__).resolve().parent/'data/reference_config.yml'


def plain(value):
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, dict):
        return {k:plain(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):
        return [plain(v) for v in value]
    return value


def write(tmp_path, raw):
    path = tmp_path/'config.yml'
    path.write_text(yaml.safe_dump(raw,sort_keys=False))
    return path


def expected_paths(raw, directory):
    return replace(raw,
        experiment=replace(raw.experiment,path=str((directory/raw.experiment.path).resolve())),
        run=replace(raw.run,output_directory=str((directory/raw.run.output_directory).resolve()),
                    output_filename=str((directory/raw.run.output_filename).resolve())),
        legacy_export=replace(raw.legacy_export,
            output_filename_raw=None if raw.legacy_export.output_filename_raw is None else
                str((directory/raw.run.output_filename).resolve())+'   '))


def compare(a,b):
    assert a == b
    x,y = compile_config(a),compile_config(b)
    assert x == y
    np.testing.assert_array_equal(x.packed_starting_vector,y.packed_starting_vector)
    assert x.range_factors == y.range_factors and x.plamda == y.plamda
    assert x.parameter_count == y.parameter_count
    return x,y


def test_reference_complete_equivalence_and_forward(tmp_path,monkeypatch):
    direct = example('reference_configuration').reference_configuration()
    direct = replace(direct,experiment=replace(direct.experiment,path='../tmp/configuration_reference/in_body.dat'),
                     run=replace(direct.run,output_directory='../tmp/full_refinement'))
    monkeypatch.chdir(tmp_path)
    loaded = load_yaml_config(REFERENCE_YAML)
    x,y = compare(loaded,expected_paths(direct,REFERENCE_YAML.parent))
    # Use a temporary body, never alter the reference input or require tmp artifacts.
    body = tmp_path/'body.dat'
    body.write_text(''.join((ROOT/'tests/fixtures/legacy_reference/in.dat').read_text().splitlines(keepends=True)[1:]))
    states = []
    for config in (loaded,expected_paths(direct,REFERENCE_YAML.parent)):
        compiled = compile_config(replace(config,experiment=replace(config.experiment,path=str(body))))
        data = load_experimental_data(compiled,base_directory=tmp_path)
        states.append((compiled,data))
    for name in ('intensities','variances','weights'):
        np.testing.assert_array_equal(getattr(states[0][1],name),getattr(states[1][1],name))
    vector = x.packed_starting_vector.copy()
    vector[0] = 4.
    results = [evaluate_refinement_trial(c.state,vector,d,execution=c.execution,objective_mode='masked')
               for c,d in states]
    for name in ('amplitudes','calculated_intensities','objective_intensities'):
        np.testing.assert_array_equal(getattr(results[0],name),getattr(results[1],name))
    assert results[0].scaling == results[1].scaling
    assert results[0].objective == results[1].objective


@pytest.mark.parametrize('count',[0,1,5,19])
def test_variable_regions_and_reflection_counts(tmp_path,count):
    for raw in (reflection_config(count), configuration(((2.,),(1.,3.),(1.,2.,4.)),(4.,7.,8.),(0,2))):
        loaded = load_yaml_config(write(tmp_path,plain(asdict(raw))))
        x,y = compare(loaded,expected_paths(raw,tmp_path))
        assert sequence(x.state) == sequence(y.state)


@pytest.fixture
def document():
    return plain(asdict(example('reference_configuration').reference_configuration()))


@pytest.mark.parametrize('path,value',[
    (('extra',),1), (('sample','bogus'),0),
    (('scan','scan_columns'),True), (('scan','scan_rows'),2.0),
    (('microscope','energy_kv'),'200'), (('microscope','energy_kv'),float('nan')),
    (('sample','dilation'),float('inf')), (('execution','execution_backend'),'gpu'),
    (('scan','corner_a'),[1,2,3]), (('sample','mesh','mesh_x'),0),
    (('schema_version',),2), (('refinement','active_categories'),['bogus']),
])
def test_invalid_fields(tmp_path,document,path,value):
    node = document
    for key in path[:-1]:
        node = node[key]
    node[path[-1]] = value
    with pytest.raises((TypeError,ValueError)):
        load_yaml_config(write(tmp_path,document))


@pytest.mark.parametrize('kind',['missing','duplicate_id','unresolved','target','thermal','slice_string','bool_enum'])
def test_invalid_semantics(tmp_path,document,kind):
    region = document['sample']['regions'][0]
    if kind == 'missing': del document['microscope']['energy_kv']
    if kind == 'duplicate_id': document['sample']['species'] *= 2
    if kind == 'unresolved': region['slices'][0]['atoms'][0]['species'] = 'missing'
    if kind == 'target': document['refinement']['parameters'][0]['target']['identifiers'] = ['missing']
    if kind == 'thermal': region['debye_waller'][0]['species'] = 'missing'
    if kind == 'slice_string': region['slices'][0]['geometry']['thickness_angstrom'] = '2.0'
    if kind == 'bool_enum': region['vacancy']['kind'] = False
    with pytest.raises((TypeError,ValueError)):
        load_yaml_config(write(tmp_path,document))


@pytest.mark.parametrize('text',[
    'title: first\ntitle: second\n',
    'sample:\n  mesh: {mesh_x: 2, mesh_x: 4}\n',
    '!!python/object/apply:os.system ["echo forbidden"]',
    'a: &a [*a]', 'a: &a {}\nb: {<<: *a}',
    '---\ntitle: one\n---\ntitle: two', '[]', '',
])
def test_unsafe_or_ambiguous_yaml(tmp_path,text):
    path = tmp_path/'bad.yml'
    path.write_text(text)
    with pytest.raises((ValueError,TypeError,yaml.YAMLError)):
        load_yaml_config(path)


def test_no_substitution_or_physical_normalization(tmp_path,document):
    document['experiment']['path'] = '${HOME}/in.dat'
    document['sample']['regions'][0]['depth_angstrom'] = -840.
    document['sample']['regions'][0]['slices'][0]['atoms'][0]['x_fractional'] = -2.5
    document['refinement']['bounds'] = {'policy':'none'}
    config = load_yaml_config(write(tmp_path,document))
    assert config.experiment.path == str(tmp_path/'${HOME}/in.dat')
    assert config.sample.regions[0].depth_angstrom == -840.
    assert config.sample.regions[0].slices[0].atoms[0].x_fractional == -2.5


def test_optional_export_and_defaults(tmp_path):
    raw = configuration(((2.,),),(4.,))
    value = plain(asdict(raw))
    value['legacy_export'] = {'filename':'export/FOR007.dat'}
    value['execution'] = {}
    loaded = load_yaml_config(write(tmp_path,value))
    assert loaded.legacy_export.filename == str(tmp_path/'export/FOR007.dat')
    assert loaded.execution.execution_backend == 'serial'
    assert not (tmp_path/'export').exists()


@pytest.mark.parametrize('count',[1,3,11])
def test_dynamic_thermal_targets_and_reversed_bindings(tmp_path,count):
    from qcbed.configuration import ParameterBinding, ParameterTarget, ParameterCategory as C
    raw = configuration(((2.,),)*count,(4.,)*count)
    bindings = tuple(ParameterBinding(ParameterTarget(C.DEBYE_WALLER,(r.id,'al')),(1.1,.9),.01)
                     for r in reversed(raw.sample.regions))
    raw = replace(raw,refinement=replace(raw.refinement,parameters=bindings))
    loaded = load_yaml_config(write(tmp_path,plain(asdict(raw))))
    compiled,_ = compare(loaded,expected_paths(raw,tmp_path))
    assert compiled.parameter_count == count
    assert [p.target for p in compiled.layout.parameters] == [(i,0) for i in range(count)]


def test_python_compiler_needs_no_yaml(monkeypatch):
    import sys
    monkeypatch.setitem(sys.modules,'yaml',None)
    assert compile_config(configuration(((2.,),),(4.,),(0,))).parameter_count == 1


def test_raw_filename_resolution_and_literal_paths(tmp_path,document):
    document['experiment']['path'] = 'data/body.dat'
    document['legacy_export']['experimental_filename_raw'] = 'data/body.dat   '
    document['run']['output_filename'] = 'reports/out.dat'
    document['legacy_export']['output_filename_raw'] = 'reports/out.dat   '
    config = load_yaml_config(write(tmp_path,document))
    assert config.legacy_export.experimental_filename_raw == str(tmp_path/'data/body.dat')+'   '
    assert config.legacy_export.output_filename_raw == str(tmp_path/'reports/out.dat')+'   '
    assert config.run.output_filename == str(tmp_path/'reports/out.dat')
    assert not (tmp_path/'data').exists()


def test_title_only_change_preserves_compilation_and_numerical_outputs(tmp_path):
    original = load_yaml_config(REFERENCE_YAML)
    # A canonical title edit clears its associated raw legacy spelling, as in the GUI.
    renamed = replace(original,title='QCBED_test',
                      legacy_export=replace(original.legacy_export,title_raw=None))
    first,second = compile_config(original),compile_config(renamed)
    assert second.state.model.title == 'QCBED_test'
    # Normalize only the two presentation fields; every other compiled field,
    # including model ordering, masks, layout, bounds and optimizer, remains exact.
    normalized_model = replace(second.state.model,title=first.state.model.title,
                               title_raw=first.state.model.title_raw)
    assert first == replace(second,configuration=original,
                           state=replace(second.state,model=normalized_model))
    np.testing.assert_array_equal(first.packed_starting_vector,second.packed_starting_vector)
    assert first.range_factors == second.range_factors
    assert first.plamda == second.plamda
    body = tmp_path/'body.dat'
    body.write_text(''.join((ROOT/'tests/fixtures/legacy_reference/in.dat').read_text().splitlines(keepends=True)[1:]))
    results = []
    for config in (original,renamed):
        compiled = compile_config(replace(config,experiment=replace(config.experiment,path=str(body))))
        data = load_experimental_data(compiled,base_directory=tmp_path)
        vector = compiled.packed_starting_vector.copy()
        vector[0] = 4.
        results.append(evaluate_refinement_trial(compiled.state,vector,data,
                       execution=compiled.execution,objective_mode='masked'))
    for field in ('amplitudes','calculated_intensities','objective_intensities'):
        np.testing.assert_array_equal(getattr(results[0],field),getattr(results[1],field))
    assert results[0].scaling == results[1].scaling  # YNORM
    assert results[0].objective == results[1].objective
