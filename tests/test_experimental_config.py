"""Header authority and exact body/forward equivalence across all reader modes."""
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from qcbed import compile_config, load_experimental_data
from qcbed.experimental import parse_experimental_intensities
from qcbed.legacy_io import LegacyInputParseError
from qcbed.refinement_trial import evaluate_refinement_trial
from qcbed.scan import generate_scan_tilts
from test_configuration import example, legacy_full
from test_dynamic_configuration import small_config

ROOT = Path(__file__).resolve().parents[1]


def configured(path, *, columns=3, rows=2, option=1, header='absent', beams=((0,0),(0,0))):
    raw = small_config()
    return compile_config(replace(raw,
        scan=replace(raw.scan,scan_columns=columns,scan_rows=rows),
        experiment=replace(raw.experiment,path=str(path),beam_indices_hk=beams,
                           gradient_option=option,header=header)))


def assert_data(a,b):
    assert (a.scan_columns,a.scan_rows,a.gradient_option,a.beam_indices_hk) == (
        b.scan_columns,b.scan_rows,b.gradient_option,b.beam_indices_hk)
    for field in ('intensities','variances','weights'):
        np.testing.assert_array_equal(getattr(a,field),getattr(b,field))
        with pytest.raises(ValueError):
            getattr(b,field).setflags(write=True)


def body_text():
    # Each logical row uses continuation and a discarded surplus value.
    return ''.join(f'{a},\n{b}, {c}, 999\n' for a,b,c in (
        (1,2,3),(4,5,6),(10,20,30),(40,50,60),
        (0,2,3),(4,5,6),(1,2,3),(4,5,6),
        (0,1,1),(1,1,1),(1,1,0),(1,1,1)))


@pytest.mark.parametrize('option', [-5,0,1,2])
def test_continued_records_rectangular_scan_duplicates_and_gradient(tmp_path,option):
    plain, headed = tmp_path/'body.dat', tmp_path/'legacy.dat'
    plain.write_text(body_text())
    headed.write_text(f'3,\n2,\n{option},888\n'+body_text())
    absent = configured(plain,option=option)
    validate = configured(headed,option=option,header='validate')
    old = parse_experimental_intensities(headed,beam_indices_hk=absent.experiment.beam_indices_hk)
    new = load_experimental_data(absent,base_directory=tmp_path)
    transitional = load_experimental_data(validate,base_directory=tmp_path)
    assert_data(old,new)
    assert_data(old,transitional)
    np.testing.assert_array_equal(new.intensities,[[1,10],[2,20],[3,30],[4,40],[5,50],[6,60]])
    assert new.variances[0,0] == new.weights[0,0] == 0
    tilts = generate_scan_tilts(absent.state.scan)
    assert tilts.shape == (6,2)
    np.testing.assert_array_equal(tilts,generate_scan_tilts(validate.state.scan))
    for mode in ('legacy','masked'):
        results = [evaluate_refinement_trial(absent.state,absent.packed_starting_vector,data,
                    execution=absent.execution,objective_mode=mode) for data in (old,new,transitional)]
        for result in results[1:]:
            for name in ('amplitudes','calculated_intensities','objective_intensities'):
                np.testing.assert_array_equal(getattr(results[0],name),getattr(result,name))
            np.testing.assert_equal(result.scaling,results[0].scaling)
            np.testing.assert_equal(result.objective,results[0].objective)
        if option != 1:
            np.testing.assert_array_equal(results[0].calculated_intensities,results[0].objective_intensities)


@pytest.mark.parametrize('header,field,record', [('4\n2\n1\n','NUMP',1),
    ('3\n3\n1\n','NUML',2),('3\n2\n0\n','GRADIENTOPT',3)])
def test_header_never_overrides_configuration(tmp_path,header,field,record):
    path = tmp_path/'mismatch.dat'
    path.write_text(header+body_text())
    config = configured(path,header='validate')
    with pytest.raises(LegacyInputParseError) as caught:
        load_experimental_data(config,base_directory=tmp_path)
    message = str(caught.value)
    assert field in message and f'record {record}' in message and 'canonical configuration' in message
    assert config.state.scan.scan_columns == 3 and config.experiment.gradient_option == 1


def test_no_header_autodetection_and_relative_paths(tmp_path,monkeypatch):
    folder = tmp_path/'inputs'
    folder.mkdir()
    # First data row deliberately looks exactly like a valid metadata header.
    text = body_text().replace('1,\n2, 3, 999\n','3 2 1\n',1)
    (folder/'body.dat').write_text(text)
    config = configured('body.dat')
    monkeypatch.chdir(tmp_path)
    data = load_experimental_data(config,base_directory=folder)
    np.testing.assert_array_equal(data.intensities[:3,0],[3,2,1])
    with pytest.raises(LegacyInputParseError):
        load_experimental_data(configured('body.dat',header='validate'),base_directory=folder)
    (folder/'headed.dat').write_text('3 2 1\n'+text)
    with pytest.raises(LegacyInputParseError):
        load_experimental_data(configured('headed.dat'),base_directory=folder)


@pytest.mark.parametrize('change', ['short','extra','invalid'])
@pytest.mark.parametrize('header', ['absent','validate'])
def test_body_failures_retain_record_diagnostics(tmp_path,change,header):
    body = body_text()
    if change == 'short':
        body = '\n'.join(body.splitlines()[:-2])
    elif change == 'extra':
        body += '123\n'
    else:
        body = body.replace('1,','NaN,',1)
    path = tmp_path/'invalid.dat'
    path.write_text(('3 2 1\n' if header == 'validate' else '')+body)
    with pytest.raises(LegacyInputParseError,match='record'):
        load_experimental_data(configured(path,header=header),base_directory=tmp_path)


def test_empty_beams_and_nonboolean_gradient(tmp_path):
    path = tmp_path/'empty.dat'
    path.write_text('')
    config = configured(path,columns=1,rows=1,option=-5,beams=())
    data = load_experimental_data(config,base_directory=tmp_path)
    assert data.intensities.shape == (1,0) and data.gradient_option == -5
    with pytest.raises(TypeError,match='CompiledConfig'):
        load_experimental_data(config.configuration,base_directory=tmp_path)


def test_reference_three_modes_setup_and_objective_equivalence(tmp_path):
    # Reference header is one line; production reader tests above separately
    # exercise multi-line logical headers. Reference files are never rewritten.
    source = ROOT/'tests/fixtures/legacy_reference/in.dat'
    body = tmp_path/'body.dat'
    body.write_text(''.join(source.read_text().splitlines(keepends=True)[1:]))
    raw = example('reference_configuration').reference_configuration()
    assert raw.experiment.header == 'absent'
    absent = compile_config(replace(raw,experiment=replace(raw.experiment,path=str(body))))
    validate = compile_config(example('reference_configuration').reference_configuration(legacy_header=True))
    state,runner = legacy_full()
    old = parse_experimental_intensities(source,beam_indices_hk=absent.experiment.beam_indices_hk,scan=state.scan)
    new = load_experimental_data(absent,base_directory=tmp_path)
    transitional = load_experimental_data(validate,base_directory=ROOT/'tests/fixtures/legacy_reference')
    assert_data(old,new)
    assert_data(old,transitional)
    assert new.intensities.shape == (1521,9)
    assert np.count_nonzero(new.variances == 0) == 117
    assert np.count_nonzero(new.weights == 0) == 6140
    for compiled in (absent,validate):
        assert replace(compiled.state,model=replace(compiled.state.model,
            experimental_intensity_filename=state.model.experimental_intensity_filename,
            experimental_intensity_filename_raw=state.model.experimental_intensity_filename_raw)) == state
        assert compiled.runner == runner
        np.testing.assert_array_equal(generate_scan_tilts(compiled.state.scan),generate_scan_tilts(state.scan))
    vector = absent.packed_starting_vector.copy()
    vector[0] = 4.  # Two occurrences, all original tilts/beams and GRADIENT=1.
    results = [evaluate_refinement_trial(s,vector,d,execution=absent.execution,objective_mode='masked')
               for s,d in ((state,old),(absent.state,new),(validate.state,transitional))]
    for value in results[1:]:
        for name in ('amplitudes','calculated_intensities','objective_intensities'):
            np.testing.assert_array_equal(getattr(results[0],name),getattr(value,name))
        assert results[0].scaling == value.scaling and results[0].objective == value.objective
