from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from qcbed import compile_config, load_yaml_config, load_experimental_data
from qcbed.legacy_config import (MigrationPolicy,import_legacy_config,export_for007,write_yaml_config)
from qcbed.configuration import RunConfig,BoundsConfig,LegacyExportConfig
from qcbed.execution import ForwardExecutionConfig
from qcbed.refinement_io import _parse_refinement_input
from qcbed.refinement import build_refinement_parameter_layout
from qcbed.refinement_trial import evaluate_refinement_trial
from test_multiregion_configuration import configuration
from test_dynamic_configuration import reflection_config

ROOT = Path(__file__).resolve().parents[1]


def policy(tmp_path):
    return MigrationPolicy(123,ForwardExecutionConfig('batched',8),
        RunConfig(str(tmp_path/'run'),str(tmp_path/'out.dat'),None,1,50,60.),'legacy','raise',BoundsConfig())


def equal_state(a,b):
    # Export relocation changes filename metadata, never physical state/masks.
    assert replace(a,model=replace(a.model,output_filename=b.model.output_filename,
        output_filename_raw=b.model.output_filename_raw,
        experimental_intensity_filename=b.model.experimental_intensity_filename,
        experimental_intensity_filename_raw=b.model.experimental_intensity_filename_raw)) == b
    assert build_refinement_parameter_layout(a) == build_refinement_parameter_layout(b)


def test_reference_roundtrip(tmp_path):
    source = ROOT/'reference_data/FOR007.DAT'
    exp = ROOT/'reference_data/in.dat'
    old,controls,_,_ = _parse_refinement_input(source,experimental_header_path=exp)
    raw = import_legacy_config(source,experimental_path=exp,policy=policy(tmp_path))
    compiled = compile_config(raw)
    equal_state(compiled.state,old)
    assert raw.legacy_export.dsmax_inverse_angstrom == .5
    bundle = tmp_path/'bundle'
    report = export_for007(compiled,bundle,bundle=True)
    assert {'execution','seed','bounds','objective_mode','invalid_trial_policy','run'} <= report.keys()
    parsed,actual,_,_ = _parse_refinement_input(bundle/'FOR007.dat')
    equal_state(parsed,old)
    assert actual == controls
    imported = import_legacy_config(bundle/'FOR007.dat',experimental_path=bundle/'in.dat',policy=policy(tmp_path))
    equal_state(compile_config(imported).state,old)
    yaml = tmp_path/'migrated.yml'
    write_yaml_config(imported,yaml)
    assert load_yaml_config(yaml) == imported
    data = load_experimental_data(compiled,base_directory=tmp_path)
    other = load_experimental_data(compile_config(imported),base_directory=tmp_path)
    for name in ('intensities','variances','weights'):
        np.testing.assert_array_equal(getattr(data,name),getattr(other,name))
    vector = compiled.packed_starting_vector.copy(); vector[0] = 4.
    results = [evaluate_refinement_trial(s,vector,d,execution=compiled.execution,objective_mode='masked')
               for s,d in ((old,data),(parsed,other))]
    for name in ('amplitudes','calculated_intensities','objective_intensities'):
        np.testing.assert_array_equal(getattr(results[0],name),getattr(results[1],name))
    assert results[0].scaling == results[1].scaling
    assert results[0].objective == results[1].objective


@pytest.mark.parametrize('count',[0,1,7,21])
def test_dynamic_roundtrip(tmp_path,count):
    configs = [reflection_config(count),configuration(((2.,),(1.,3.),(1.,2.,4.)),(4.,7.,8.),(0,2))]
    for i,raw in enumerate(configs):
        body = tmp_path/f'body{i}.dat'; body.write_text('1 1\n1 1\n'*3)
        raw = replace(raw,experiment=replace(raw.experiment,path=str(body)),legacy_export=LegacyExportConfig(.12345678901234567))
        compiled = compile_config(raw)
        bundle = tmp_path/f'bundle{i}'
        export_for007(compiled,bundle,bundle=True)
        parsed,controls,_,_ = _parse_refinement_input(bundle/'FOR007.dat')
        equal_state(parsed,compiled.state)
        assert controls['plamda'] == compiled.plamda
        imported = import_legacy_config(bundle/'FOR007.dat',experimental_path=bundle/'in.dat',policy=policy(tmp_path))
        equal_state(compile_config(imported).state,compiled.state)
        assert imported.legacy_export.dsmax_inverse_angstrom == raw.legacy_export.dsmax_inverse_angstrom
        migrated = compile_config(imported)
        values = compiled.packed_starting_vector
        data = load_experimental_data(compiled,base_directory=tmp_path)
        other = load_experimental_data(migrated,base_directory=tmp_path)
        before = evaluate_refinement_trial(compiled.state,values,data,execution=compiled.execution)
        after = evaluate_refinement_trial(migrated.state,values,other,execution=compiled.execution)
        np.testing.assert_array_equal(before.amplitudes,after.amplitudes)
        np.testing.assert_array_equal(before.calculated_intensities,after.calculated_intensities)
        assert before.scaling == after.scaling and before.objective == after.objective


def test_refuse_missing_metadata_and_overwrite(tmp_path):
    raw = reflection_config(1)
    with pytest.raises(ValueError,match='DSMAX'):
        export_for007(raw,tmp_path/'bad',bundle=True)
    raw = replace(raw,legacy_export=LegacyExportConfig(.5))
    with pytest.raises(ValueError,match='headerless'):
        export_for007(raw,tmp_path/'bad')
    with pytest.raises(ValueError,match='protected'):
        export_for007(raw,ROOT/'reference_data/no-write',bundle=True)
    target = tmp_path/'existing'; target.write_text('keep')
    with pytest.raises(FileExistsError): export_for007(raw,target,bundle=True)
    assert target.read_text() == 'keep'


def test_linked_empty_g_and_duplicate_beams(tmp_path):
    from qcbed.configuration import AbsorptionConfig,AbsorptionSpeciesConfig
    raw = configuration(((2.,),(3.,)),(4.,6.),())
    regions = list(raw.sample.regions)
    regions[1] = replace(regions[1],slices=(replace(regions[1].slices[0],atoms=(),record_type='G',record_type_raw='G'),))
    body = tmp_path/'body.dat'; body.write_text('1 2\n3 4\n'*6)
    raw = replace(raw,sample=replace(raw.sample,regions=tuple(regions),
        absorption=AbsorptionConfig(True,(AbsorptionSpeciesConfig('al',.12345678901234567,-2.,-7.,'linked_to_b'),))),
        experiment=replace(raw.experiment,path=str(body),beam_indices_hk=((0,0),(0,0))),legacy_export=LegacyExportConfig(.5))
    export_for007(raw,tmp_path/'bundle',bundle=True)
    imported = import_legacy_config(tmp_path/'bundle/FOR007.dat',experimental_path=tmp_path/'bundle/in.dat',policy=policy(tmp_path))
    assert imported.sample.regions[1].slices[0].record_type == 'G'
    assert imported.sample.regions[1].slices[0].atoms == ()
    assert imported.sample.absorption.coefficients[0].c_mode == 'linked_to_b'
    assert imported.sample.absorption.coefficients[0].c == -7.
    assert imported.experiment.beam_indices_hk == ((0,0),(0,0))
    equal_state(compile_config(imported).state,compile_config(raw).state)


def test_standalone_limits_and_no_truncation(tmp_path):
    raw = import_legacy_config(ROOT/'reference_data/FOR007.DAT',experimental_path=ROOT/'reference_data/in.dat',policy=policy(tmp_path))
    with pytest.raises(ValueError,match='filename_limits'):
        export_for007(raw,tmp_path/'out')
    with pytest.raises(ValueError,match='exceed'):
        export_for007(raw,tmp_path/'out',filename_limits=(9,10))
    export_for007(raw,tmp_path/'out',experimental_filename='in.dat',output_filename='out.dat',filename_limits=(9,10))
    parsed,_,_,_ = _parse_refinement_input(tmp_path/'out',experimental_header_path=ROOT/'reference_data/in.dat')
    equal_state(parsed,compile_config(raw).state)
    for invalid in (replace(raw,scan=replace(raw.scan,basis_1=(.5,1.,0.))),
                    replace(raw,sample=replace(raw.sample,regions=(replace(raw.sample.regions[0],depth_angstrom=-1.),)))):
        with pytest.raises(ValueError): export_for007(invalid,tmp_path/'bad',bundle=True)
        assert not (tmp_path/'bad').exists()


def test_full_parameter_categories(tmp_path):
    from qcbed.configuration import ParameterBinding,ParameterTarget,ParameterCategory as C,VacancyConfig,VacancyKind,SurfaceConfig,SurfaceKind
    raw = import_legacy_config(ROOT/'reference_data/FOR007.DAT',experimental_path=ROOT/'reference_data/in.dat',policy=policy(tmp_path))
    for flag in (0,2):
        region = replace(raw.sample.regions[0],vacancy=VacancyConfig(VacancyKind(flag),(-.2,) if flag == 0 else (-.2,.3)),
                         free_surface=SurfaceConfig(SurfaceKind.LEGACY_1,(.1,.2,.3)))
        targets = [(C.VOLTAGE,()),(C.DILATION,()),(C.ATOM_COORDINATE,('slice-0','atom-0','x')),
                   (C.ABSORPTION,('species-0','a')),(C.REFLECTION_MULTIPLIER,('region-0','reflection-0','absorptive'))]
        targets += [(C.CONSTANT_VACANCY,('region-0',))] if flag == 0 else [(C.VACANCY_COEFFICIENT,('region-0',x)) for x in ('a','b')]
        bindings = raw.refinement.parameters+tuple(ParameterBinding(ParameterTarget(category,ids),(1.03,.97),.12345678901234567) for category,ids in targets)
        config = replace(raw,sample=replace(raw.sample,regions=(region,)),refinement=replace(raw.refinement,parameters=bindings))
        export_for007(config,tmp_path/f'b{flag}',bundle=True)
        imported = import_legacy_config(tmp_path/f'b{flag}/FOR007.dat',experimental_path=tmp_path/f'b{flag}/in.dat',policy=policy(tmp_path))
        equal_state(compile_config(imported).state,compile_config(config).state)
        assert compile_config(imported).plamda == compile_config(config).plamda


def test_cli_migration_and_export(tmp_path,capsys):
    from qcbed.legacy import main
    yaml = tmp_path/'migration.yml'
    assert main(['import','--for007',str(ROOT/'reference_data/FOR007.DAT'),
        '--experimental',str(ROOT/'reference_data/in.dat'),'--yaml',str(yaml),
        '--seed','123','--execution','serial','--objective','masked','--invalid-trial-policy','reject',
        '--bounds','none','--output-directory',str(tmp_path/'run'),'--output-filename',str(tmp_path/'out.dat'),
        '--snapshot-every','0','--checkpoint-every','1','--max-evaluations','50','--max-seconds','60']) == 0
    assert main(['export','--config',str(yaml),'--destination',str(tmp_path/'bundle'),'--bundle']) == 0
    assert 'masked' in capsys.readouterr().out
    assert (tmp_path/'bundle/compatibility.json').exists()
    with pytest.raises(FileExistsError): write_yaml_config(load_yaml_config(yaml),yaml)


def test_variable_species_thermal_targets(tmp_path):
    from qcbed.configuration import SpeciesConfig,ThermalFactorConfig,ParameterBinding,ParameterTarget,ParameterCategory as C
    raw = configuration(((2.,),(1.,3.),(2.,)),(4.,5.,6.))
    species = tuple(SpeciesConfig(f's{i}',f'{i:02}',13+i,f'{i:02}') for i in range(4))
    regions = tuple(replace(r,slices=tuple(replace(s,atoms=tuple(replace(a,species='s0') for a in s.atoms)) for s in r.slices),
        debye_waller=tuple(ThermalFactorConfig(sp.id,.1+j*.02) for j,sp in enumerate(species))) for r in raw.sample.regions)
    bindings = tuple(ParameterBinding(ParameterTarget(C.DEBYE_WALLER,(r.id,sp.id)),(1.1,.9),.01)
                     for r in regions for sp in species)
    body = tmp_path/'body.dat'; body.write_text('1 1\n1 1\n'*3)
    raw = replace(raw,sample=replace(raw.sample,species=species,regions=regions),
        experiment=replace(raw.experiment,path=str(body)),refinement=replace(raw.refinement,parameters=bindings),
        legacy_export=LegacyExportConfig(.5))
    export_for007(raw,tmp_path/'bundle',bundle=True)
    imported = import_legacy_config(tmp_path/'bundle/FOR007.dat',experimental_path=tmp_path/'bundle/in.dat',policy=policy(tmp_path))
    equal_state(compile_config(imported).state,compile_config(raw).state)
    assert compile_config(imported).parameter_count == 12


def test_import_rejects_unknown_trailing_records(tmp_path):
    path = tmp_path/'legacy.dat'
    path.write_text((ROOT/'reference_data/FOR007.DAT').read_text()+'unknown extension\n')
    with pytest.raises(ValueError,match='trailing'):
        import_legacy_config(path,experimental_path=ROOT/'reference_data/in.dat',policy=policy(tmp_path))
