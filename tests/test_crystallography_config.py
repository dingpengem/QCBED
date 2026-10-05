from dataclasses import replace
import builtins
from pathlib import Path
import numpy as np
import pytest
import yaml
from qcbed import load_yaml_config, save_yaml_config, compile_config, compile_config_from_sources
from qcbed.configuration import CbedReflectionBasisConfig, ExperimentalBeamConfig, CrystalScanConfig
from qcbed.crystallography_config import load_config_structure
from qcbed.crystallography import build_zone_axis_simulation_sample
from qcbed.config_validation import validate_config, ValidationPhase

ROOT=Path(__file__).resolve().parents[1]
EXAMPLE=ROOT/'examples/crystallography_config.yml'


def config():
    c=load_yaml_config(EXAMPLE)
    return replace(c,run=replace(c.run,output_filename=str(Path(c.run.output_directory)/'final_out.dat')))


def test_compile_direct_equivalence_and_nine_beams():
    raw=config(); structure=load_config_structure(raw)
    compiled=compile_config(raw,structure=structure)
    bundle=compiled.crystallography
    assert compiled.source_configuration is raw
    assert compiled.configuration.crystallography is None
    assert len(bundle.indexed_beams)==9
    assert [(b.stable_id,b.experimental_index) for b in bundle.indexed_beams]==[(b.id,b.lattice_coordinates) for b in raw.crystallography.beams]
    for beam in bundle.indexed_beams:
        assert beam.original_hkl==(*beam.experimental_index,0)
        assert beam.internal_hk.as_HKL()==beam.original_hkl
    applied=bundle.simulation_sample.depth_sample.original_structure
    c=raw.crystallography
    direct=build_zone_axis_simulation_sample(applied,bundle.zone_cell,
        slice_boundaries_fractional=c.slice_boundaries_fractional,depth_angstrom=c.depth_angstrom,
        mesh=c.mesh,indexed_beams=bundle.indexed_beams,energy_kv=raw.microscope.energy_kv,
        dsap_inverse_angstrom=c.dsap_inverse_angstrom,absorption_coefficients=None,
        reciprocal_ordering='unshifted',phase_grating_method='direct')
    np.testing.assert_array_equal(bundle.simulation_cell.direct_lattice.vectors,direct.depth_sample.simulation_cell.direct_lattice.vectors)
    assert direct.parameters==bundle.simulation_sample.parameters
    assert compiled.state.model.atoms_by_slice==tuple(tuple(replace(a,record_type_raw="S") for a in row) for row in direct.model.atoms_by_slice)
    from qcbed.refinement_trial import _prepare_trial
    applied,sequence,_,absorption=_prepare_trial(compiled.state,original_voltage=raw.microscope.energy_kv)
    assert applied==direct.parameters.slices
    assert sequence==direct.parameters.slice_sequence_zero_based
    assert absorption==direct.parameters.absorption_coefficients
    for a,b in zip(direct.prepared.slice_types,bundle.simulation_sample.prepared.slice_types):
        np.testing.assert_array_equal(a.phase_grating,b.phase_grating)
    assert compiled.state.scan.corner_a==raw.scan.corner_a
    assert compiled.state.scan.basis_2==(0.,-1.,0.)


def test_pure_compiler_and_round_trip(tmp_path,monkeypatch):
    raw=config(); structure=load_config_structure(raw)
    def fail(*args,**kwargs): raise AssertionError('filesystem I/O in pure compiler')
    with monkeypatch.context() as m:
        m.setattr(builtins,'open',fail)
        m.setattr(Path,'open',fail)
        compile_config(raw,structure=structure)
    destination=tmp_path/'saved.yml'
    saved=save_yaml_config(raw,destination,base_directory=tmp_path)
    assert load_yaml_config(destination)==saved==raw
    save_yaml_config(raw,destination,base_directory=tmp_path,overwrite=True)
    assert load_yaml_config(destination)==raw


@pytest.mark.parametrize('hkl',[(0,0,0),(0,0,1),(2,0,0)])
def test_invalid_basis(hkl):
    raw=config();structure=load_config_structure(raw)
    bad=replace(raw,crystallography=replace(raw.crystallography,cbed_basis=CbedReflectionBasisConfig((1,0,0),hkl)))
    with pytest.raises(ValueError,match='cbed_basis'):
        compile_config(bad,structure=structure)
    issues=validate_config(bad,phases=(ValidationPhase.COMPILE,),structure=structure)
    assert issues[0].section=='crystallography'
    assert issues[0].location.field=='cbed_basis'


def test_strict_rejection_and_no_mixed_modes():
    raw=config();structure=load_config_structure(raw)
    bad=replace(raw,microscope=replace(raw.microscope,zone_axis=(1,1,0)),
        crystallography=replace(raw.crystallography,zone_axis=(1,1,0),
            cbed_basis=CbedReflectionBasisConfig((1,-1,0),(0,0,1))))
    with pytest.raises(ValueError,match='Unsupported depth-resolved'):
        compile_config(bad,structure=structure)
    derived=compile_config(raw,structure=structure).configuration
    with pytest.raises(ValueError,match='mutually exclusive'):
        compile_config(replace(raw,sample=derived.sample),structure=structure)


@pytest.mark.parametrize('text', ['bogus', '1\nno lattice\nSi 0 0 0\n'])
def test_unreadable_or_invalid_lattice(tmp_path,text):
    raw=config();p=tmp_path/'bad.xyz';p.write_text(text)
    bad=replace(raw,crystallography=replace(raw.crystallography,structure=replace(raw.crystallography.structure,path=str(p),format='xyz')))
    with pytest.raises(ValueError,match='cannot read structure'):
        compile_config_from_sources(bad)


def test_missing_file_and_missing_ase(tmp_path,monkeypatch):
    raw=config()
    bad=replace(raw,crystallography=replace(raw.crystallography,structure=replace(raw.crystallography.structure,path=str(tmp_path/'absent.cif'))))
    with pytest.raises(ValueError,match='not found'):
        compile_config_from_sources(bad)
    original=builtins.__import__
    def blocked(name,*a,**kw):
        if name=='ase.io': raise ImportError('simulated absent ASE')
        return original(name,*a,**kw)
    monkeypatch.setattr(builtins,'__import__',blocked)
    with pytest.raises(ValueError,match=r'qcbed\[crystallography\]'):
        compile_config_from_sources(raw)


def test_strict_yaml_indices(tmp_path):
    data=yaml.safe_load(EXAMPLE.read_text())
    for value,match in [([1,0,0,0],'4-index'),([1.5,0,0],'integer')]:
        data['crystallography']['cbed_basis']['g1_hkl']=value
        p=tmp_path/'bad.yml';p.write_text(yaml.safe_dump(data))
        with pytest.raises((ValueError,TypeError),match=match):load_yaml_config(p)


def test_prepared_dataset_bridge_no_array_or_corner_changes():
    from test_prepared_dataset_integration import assembled
    from qcbed.data_preparation import SignalMode
    _,conversion,data=assembled(SignalMode.CONVENTIONAL,True)
    raw=config()
    raw=replace(raw,crystallography=replace(raw.crystallography,
        beams=tuple(ExperimentalBeamConfig(b.id,n) for b,n in zip(data.beams,data.lattice_coordinates))),
        scan=CrystalScanConfig(data.scan_geometry.scan_columns,data.scan_geometry.scan_rows,
            data.scan_geometry.corner_a,data.scan_geometry.corner_b,data.scan_geometry.corner_c))
    compiled=compile_config_from_sources(raw)
    bundle=compiled.crystallography
    snapshots=[a.copy() for a in (data.intensities,data.variances,data.weights)]
    indexed=bundle.index_dataset(data)
    assert indexed.dataset is data
    adapted=bundle.adapt_dataset(data,objective_mode='masked')
    assert adapted.scan.corner_a==compiled.state.scan.corner_a==conversion.scan_geometry.corner_a
    assert adapted.scan.corner_b==compiled.state.scan.corner_b
    assert adapted.scan.corner_c==compiled.state.scan.corner_c
    assert adapted.experimental.beam_indices_hk==compiled.experiment.beam_indices_hk
    for actual,old in zip((data.intensities,data.variances,data.weights),snapshots):
        np.testing.assert_array_equal(actual,old)


@pytest.mark.parametrize('rows', [((3,0,0),(0,4,0),(0,0,5)),
    ((3,0,0),(-1.5,2.598076211353316,0),(0,0,5)),
    ((3,0,5),(1,4,-5),(0,0,5))])
def test_noncubic_compile_uses_same_backend(rows):
    from qcbed.crystallography import DirectLattice
    raw=config();structure=replace(load_config_structure(raw),direct_lattice=DirectLattice(rows))
    result=compile_config(raw,structure=structure)
    assert result.crystallography.diagnostics['strict_b5_compatible']
    assert result.crystallography.diagnostics['simulation_atom_count']==2
    for beam in result.crystallography.indexed_beams:
        g=np.array(beam.original_hkl)@structure.direct_lattice.reciprocal().vectors
        gi=np.array(beam.internal_hk.as_HKL())@result.crystallography.simulation_cell.direct_lattice.reciprocal().vectors
        np.testing.assert_allclose(g,gi,rtol=1e-14,atol=1e-14)


def test_basis_sign_order_duplicates_and_invalid_transform():
    raw=config();structure=load_config_structure(raw)
    c=replace(raw.crystallography,cbed_basis=CbedReflectionBasisConfig((0,-2,0),(1,0,0)),
        beams=raw.crystallography.beams+(ExperimentalBeamConfig('duplicate',(1,0)),))
    compiled=compile_config(replace(raw,crystallography=c),structure=structure)
    assert compiled.crystallography.indexed_beams[1].original_hkl==(0,-2,0)
    assert compiled.experiment.beam_indices_hk[-1]==compiled.experiment.beam_indices_hk[1]
    bad=replace(raw,crystallography=replace(c,scan_basis_transform=((1,0),(0,2))))
    with pytest.raises(ValueError,match='unimodular'):compile_config(bad,structure=structure)


def test_compile_snapshot_without_objective(tmp_path,monkeypatch):
    from qcbed.canonical_run import prepare_run_snapshot,canonical_checkpoint_identity
    raw=config()
    snapshot=prepare_run_snapshot(raw,base_directory=tmp_path)
    before=canonical_checkpoint_identity(snapshot.compiled)
    assert snapshot.experimental.beam_indices_hk==snapshot.compiled.experiment.beam_indices_hk
    assert snapshot.experimental.intensities.shape==(1,9)
    assert canonical_checkpoint_identity(snapshot.compiled)==before
    from qcbed.gui.document import EditorDocument
    document=EditorDocument(raw,EXAMPLE)
    assert document.config==raw
    assert document.summary()['compile']=='Valid (committed snapshot)'


def test_no_ase_for_legacy_compile(monkeypatch):
    raw=load_yaml_config(ROOT/'examples/qcbed_config.yml')
    original=builtins.__import__
    def blocked(name,*a,**kw):
        if name.startswith('ase'): raise AssertionError('legacy path imported ASE')
        return original(name,*a,**kw)
    monkeypatch.setattr(builtins,'__import__',blocked)
    assert compile_config_from_sources(raw).crystallography is None


@pytest.mark.parametrize('fmt,suffix',[('cif','.cif'),('vasp','.vasp')])
def test_ase_common_structure_files(tmp_path,fmt,suffix):
    from ase import Atoms
    from ase.io import write
    raw=config()
    path=tmp_path/('structure'+suffix)
    write(path,Atoms('Si',scaled_positions=[(.1,.2,.3)],cell=[4,4,4],pbc=True),format=fmt)
    raw=replace(raw,crystallography=replace(raw.crystallography,
        structure=replace(raw.crystallography.structure,path=str(path),format=fmt)))
    result=compile_config_from_sources(raw)
    assert len(result.crystallography.original_structure.sites)==1
    assert result.crystallography.simulation_cell is not None


def test_source_relative_paths_not_cwd_and_save_protection(tmp_path,monkeypatch):
    raw=load_yaml_config(EXAMPLE)  # Explicitly exercise the old output-filename reader.
    monkeypatch.chdir(tmp_path)
    assert load_yaml_config(EXAMPLE)==raw
    with pytest.raises(ValueError,match='structure source'):
        save_yaml_config(raw,raw.crystallography.structure.path,base_directory=tmp_path,overwrite=True)
    with pytest.raises(ValueError,match='unresolved structure'):
        compile_config(raw)
