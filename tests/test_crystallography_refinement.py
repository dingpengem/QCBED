"""B7a groups/aliases/application/bounds, not optimizer integration."""
from dataclasses import replace, FrozenInstanceError, asdict
import hashlib
import json
import numpy as np
import pytest
from examples.b7a_refinement_layout import real_config, run_diagnostics, ROOT
from qcbed import compile_config, compile_config_from_sources, load_yaml_config
from qcbed.configuration import AtomicSiteGroup, ParameterBinding, TargetBound, BoundsConfig
from qcbed.crystallography import AtomicSite, AtomicStructureModel, DirectLattice
from qcbed.crystallography.site_groups import resolve_site_groups, group_target
from qcbed.crystallography_refinement import apply_crystallographic_parameters, physical_target_layout, site_occurrence_provenance
from qcbed.config_targets import available_refinement_targets
from qcbed.canonical_run import canonical_checkpoint_identity, run_snapshot, RunSnapshot


@pytest.fixture
def compiled():return compile_config_from_sources(real_config())


@pytest.mark.parametrize('kw', [dict(id=''),dict(element='XX'),dict(site_indices=(True,)),
    dict(site_indices=(1.,)),dict(site_indices=(-1,)),dict(site_indices=(0,0)),dict(site_indices=()),
    dict(occupancy=0),dict(occupancy=1.1),dict(debye_waller_b=-1),dict(debye_waller_b=float('nan')),
    dict(absorption=(0.,float('inf'),1.))])
def test_group_invalid_fields(kw):
    args=dict(id='g',element='Al',site_indices=(0,),debye_waller_b=.86,occupancy=1.)
    args.update(kw)
    with pytest.raises((ValueError,TypeError)):AtomicSiteGroup(**args)


def test_group_membership_validation_and_override():
    source=AtomicStructureModel(DirectLattice(np.eye(3)*4),
        (AtomicSite('Al',(0,0,0),.5,.1),AtomicSite('Si',(.5,.5,.5),1.,.2)))
    al=AtomicSiteGroup('al','Al',(0,),.86,1.,(0.,1.,1.))
    si=AtomicSiteGroup('si','Si',(1,),.3,.8,(0.,1.,1.))
    result=resolve_site_groups(source,(si,al))
    assert [s.symbol for s in result.sites]==['Al','Si'] # Group declaration never sorts atoms.
    assert result.sites[0].occupancy==1 and result.sites[0].debye_waller_b==.86
    assert source.sites[0].occupancy==.5 and source.sites[0].debye_waller_b==.1
    for groups in ((al,), (al,replace(si,id='al')), (al,replace(si,site_indices=(0,))),
                   (al,replace(si,site_indices=(2,))), (al,replace(si,element='Al'))):
        with pytest.raises(ValueError):resolve_site_groups(source,groups)
    with pytest.raises(FrozenInstanceError):al.id='changed'


def test_nonempty_layout_catalogue_and_provenance(compiled):
    assert compiled.parameter_count==1
    assert compiled.packed_starting_vector.tolist()==[.86]
    assert [(p.category,p.target) for p in compiled.layout.parameters]==[('debye_waller',(0,0))]
    active=[t for t in available_refinement_targets(compiled) if t.active]
    assert len(active)==1 and active[0].target==group_target('Al_fcc') and active[0].available
    assert active[0].stable_ids==('Al_fcc',)
    assert physical_target_layout(compiled)[0]['lower']==.5
    p=site_occurrence_provenance(compiled)
    assert {i for row in p['stored_slices'] for i,g in row}=={0,1,2,3}
    assert sum(len(p['stored_slices'][i]) for i in p['sequence'])==830
    assert all(g=='Al_fcc' for row in p['stored_slices'] for _,g in row)


@pytest.mark.parametrize('value',[.5,1.5,.8602])
def test_bounds_accepted_and_scattering_updated(compiled,value):
    changed=apply_crystallographic_parameters(compiled,np.array([value]))
    assert compiled.packed_starting_vector.tolist()==[.86]
    assert changed.packed_starting_vector.tolist()==[value]
    sample=changed.crystallography.simulation_sample
    for sites in (sample.depth_sample.original_structure.sites,sample.depth_sample.recast_structure.sites,
                  tuple(s for slab in sample.depth_sample.slices for s in slab.sites)):
        assert all(s.debye_waller_b==value and s.physical_group_id=='Al_fcc' for s in sites)
        assert all(s.occupancy==1 for s in sites)
    assert all(p.debye_waller_angstrom_squared==(value,) for p in sample.parameters.slices)
    assert any(not np.array_equal(a.phase_grating,b.phase_grating) for a,b in zip(
        sample.prepared.slice_types,compiled.crystallography.simulation_sample.prepared.slice_types))
    assert sample.parameters.slice_sequence_zero_based==compiled.crystallography.simulation_sample.parameters.slice_sequence_zero_based
    assert changed.state.scan==compiled.state.scan
    assert changed.state.absorption==compiled.state.absorption


@pytest.mark.parametrize('value',[np.nextafter(.5,-np.inf),np.nextafter(1.5,np.inf),-.86,np.nan,np.inf,-np.inf])
def test_bounds_rejected_before_rebuild(compiled,value,monkeypatch):
    import qcbed.crystallography_refinement as api
    def forbidden(*a,**kw):raise AssertionError('rebuild must follow validation')
    monkeypatch.setattr(api,'compile_config',forbidden)
    with pytest.raises(ValueError):apply_crystallographic_parameters(compiled,np.array([value]))


def test_two_groups_order_and_unrelated_state(compiled):
    cfg=compiled.source_configuration
    groups=(AtomicSiteGroup('Si_g','Si',(1,),.3,1.),AtomicSiteGroup('Al_g','Al',(0,),.86,1.))
    source=AtomicStructureModel(DirectLattice(np.eye(3)*4.049),
        (AtomicSite('Al',(0,0,0)),AtomicSite('Si',(.5,.5,.5))))
    targets=(group_target('Si_g'),group_target('Al_g'))
    config=replace(cfg,crystallography=replace(cfg.crystallography,site_groups=groups),
        refinement=replace(cfg.refinement,parameters=tuple(ParameterBinding(t,(1.05,.95),.0002) for t in targets),
            bounds=BoundsConfig('reject',tuple(TargetBound(t,.1,1.5) for t in targets))))
    c=compile_config(config,structure=source)
    assert c.packed_starting_vector.tolist()==[.86,.3] # Original species order, not binding/group order.
    v=c.packed_starting_vector.copy();v[0]=.87
    altered=apply_crystallographic_parameters(c,v)
    assert altered.crystallography.original_structure.sites[1]==c.crystallography.original_structure.sites[1]
    assert altered.crystallography.site_groups[0]==c.crystallography.site_groups[0]
    reordered=replace(config,refinement=replace(config.refinement,parameters=config.refinement.parameters[::-1]))
    other=compile_config(reordered,structure=source)
    assert other.layout==c.layout
    assert canonical_checkpoint_identity(other)==canonical_checkpoint_identity(c)


def test_yaml_snapshot_and_checkpoint_identity(compiled,tmp_path):
    from qcbed.config_serialization import save_yaml_config
    config=compiled.source_configuration
    # Modern Save owns a stable filename inside the output directory.
    from pathlib import Path
    config=replace(config,run=replace(config.run,output_filename=str(Path(config.run.output_directory)/'final_out.dat')))
    compiled=compile_config_from_sources(config)
    saved=save_yaml_config(config,tmp_path/'snapshot.yml',base_directory=ROOT)
    loaded=load_yaml_config(tmp_path/'snapshot.yml')
    assert loaded==saved
    restored=compile_config_from_sources(loaded)
    identical=apply_crystallographic_parameters(compiled,compiled.packed_starting_vector)
    assert canonical_checkpoint_identity(restored)==canonical_checkpoint_identity(compiled)==canonical_checkpoint_identity(identical)
    assert asdict(restored.state)==asdict(compiled.state)==asdict(identical.state)
    assert canonical_checkpoint_identity(compiled)['version']==2
    renamed=replace(config,crystallography=replace(config.crystallography,
        site_groups=(replace(config.crystallography.site_groups[0],id='renamed'),)),
        refinement=replace(config.refinement,
            parameters=(replace(config.refinement.parameters[0],target=group_target('renamed')),),
            bounds=BoundsConfig('reject',(TargetBound(group_target('renamed'),.5,1.5),))))
    assert canonical_checkpoint_identity(compile_config_from_sources(renamed))!=canonical_checkpoint_identity(compiled)
    legacy=compile_config_from_sources(load_yaml_config(ROOT/'examples/fcc_al_real_preflight.yml'))
    payload=dict(state=asdict(legacy.state),runner=asdict(legacy.runner),layout=asdict(legacy.layout),
        experiment=asdict(legacy.experiment),run=asdict(legacy.configuration.run))
    assert canonical_checkpoint_identity(legacy)==dict(format='qcbed-canonical-run',version=1,
        digest=hashlib.sha256(json.dumps(payload,sort_keys=True,allow_nan=False).encode()).hexdigest())
    from examples.a9_real_variance_closure import prepare_reference_variance
    data,_=prepare_reference_variance()
    exp=compiled.crystallography.adapt_dataset(data['conventional'],objective_mode='masked').experimental
    with pytest.raises(ValueError,match='checkpoint'):run_snapshot(RunSnapshot(compiled,exp),resume=True)


def test_real_baseline_roundtrip_and_perturbation_without_optimizer(monkeypatch):
    import qcbed.refinement_driver as driver
    def forbidden(*a,**kw):raise AssertionError('No optimizer allowed')
    monkeypatch.setattr(driver,'run_legacy_refinement',forbidden)
    report,baseline,reapplied,perturbed=run_diagnostics()
    assert report['round_trip_exact'] and report['b6b_baseline_match']
    assert report['baseline']['ynorm']==14870390.0181119
    assert report['baseline']['objective']==.939675432747781
    assert report['baseline']['determinism_repeat_evaluations']==0
    assert report['perturbation']['result']['determinism_repeat_evaluations']==0
    assert report['perturbation']['max_phase_grating_change']>0
    assert report['perturbation']['max_amplitude_change']>0
    assert report['perturbation']['result']['objective']!=report['baseline']['objective']
    assert baseline.source_configuration==reapplied.source_configuration


def test_unsupported_same_element_groups_and_duplicate_authority(compiled):
    source=compiled.crystallography.original_structure
    groups=(AtomicSiteGroup('a','Al',(0,1),.86,1.),AtomicSiteGroup('b','Al',(2,3),.86,1.))
    with pytest.raises(ValueError,match='same-element'):resolve_site_groups(source,groups)
    from qcbed.configuration import ThermalFactorConfig
    config=compiled.source_configuration
    bad=replace(config,crystallography=replace(config.crystallography,
        debye_waller=(ThermalFactorConfig('Al',.86),)))
    with pytest.raises(ValueError,match='conflicting authorities'):compile_config_from_sources(bad)


def test_existing_edit_helper_selects_group_target(compiled):
    from qcbed.config_editing import set_refinement
    raw=compiled.source_configuration
    disabled=set_refinement(raw,group_target('Al_fcc'),enabled=False)
    assert compile_config_from_sources(disabled).parameter_count==0
    enabled=set_refinement(disabled,group_target('Al_fcc'),enabled=True,
        range_factors=(1.05,.95),plamda=.0002)
    assert enabled==raw
