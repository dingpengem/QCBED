"""B7c fixed depth candidates, with unchanged one-dimensional PARAMIN for B."""
from dataclasses import replace, asdict
from pathlib import Path
import json
import numpy as np
from .configuration import ParameterBinding, TargetBound
from .crystallography.site_groups import requested_depth_target
from .crystallography_refinement import apply_crystallographic_parameters, physical_target_layout
from .depth_geometry import requested_depth_to_sample_geometry
from .canonical_run import canonical_checkpoint_identity
from .objective_diagnostics import evaluate_preflight_objectives
from .single_target_refinement import audit_binary_masks, _digest, _run_journal, load_single_target_checkpoint
from .data_preparation import prepared_dataset_manifest
from .full_runner import _atomic_json

STAGE_SETTINGS=dict(optimizer='PARAMIN',max_cycles=0,tolerance=0.,random_order=False,
                    maximum_optimizer_proposals=5,bounds=[.5,1.5])


def enable_requested_depth(config):
    """Add physical requested depth with existing reference limits and D1/D2."""
    target=requested_depth_target()
    return replace(config,refinement=replace(config.refinement,
        parameters=config.refinement.parameters+(ParameterBinding(target,(1.05,.95),.5),),
        bounds=replace(config.refinement.bounds,specific=config.refinement.bounds.specific+(TargetBound(target,600.,1200.),))))


def sample_extent(compiled):
    sample=compiled.crystallography.simulation_sample
    geometry=requested_depth_to_sample_geometry(sample.requested_depth_angstrom,
        (s.thickness_angstrom for s in sample.model.slices))
    if geometry.sequence!=sample.parameters.slice_sequence_zero_based or geometry.effective_depth!=sample.effective_depth_angstrom:
        raise ValueError('compiled geometry differs from canonical depth expansion')
    return dict(requested_depth=geometry.requested_depth,effective_depth=geometry.effective_depth,
        occurrence_count=geometry.occurrence_count,stored_slice_count=geometry.stored_slice_count,
        sequence_digest=_digest(geometry.sequence),plateau_lower_exclusive=geometry.plateau_lower_exclusive)


def evaluate_depth_point(compiled,datasets,depth,b,*,signal_mode="conventional"):
    """Apply the source-ordered vector [requested depth,B]; validate invariants."""
    expected=['sample.requested_depth','site_group.Al_fcc.debye_waller_b']
    if [v['target_id'] for v in physical_target_layout(compiled)]!=expected:
        raise ValueError('B7c requires exactly source-ordered requested depth then group B')
    current=apply_crystallographic_parameters(compiled,np.array([depth,b],dtype=np.float64))
    base=compiled.source_configuration;cfg=current.source_configuration
    restored=replace(cfg,crystallography=replace(cfg.crystallography,
        depth_angstrom=base.crystallography.depth_angstrom,
        site_groups=tuple(replace(g,debye_waller_b=old.debye_waller_b) for g,old in
                         zip(cfg.crystallography.site_groups,base.crystallography.site_groups))))
    if restored!=base:raise ValueError('unrelated configuration changed')
    a,bundle=compiled.crystallography,current.crystallography
    np.testing.assert_array_equal(a.original_structure.direct_lattice.vectors,bundle.original_structure.direct_lattice.vectors)
    for old,site in zip(a.original_structure.sites,bundle.original_structure.sites,strict=True):
        if replace(site,debye_waller_b=old.debye_waller_b)!=old:raise ValueError('unrelated site state changed')
    if a.simulation_sample.model.slices!=bundle.simulation_sample.model.slices:
        raise ValueError('stored slice geometry changed')
    if signal_mode not in ('conventional','gradient'):raise ValueError('unsupported signal mode')
    count=7974 if signal_mode=='conventional' else 7281
    audit=audit_binary_masks(datasets)
    if audit[signal_mode]['active_count']!=count:raise ValueError('mask changed')
    manifest=prepared_dataset_manifest(datasets[signal_mode])
    if manifest['signal_mode']!=signal_mode:raise ValueError('dataset signal mode mismatch')
    r,_=evaluate_preflight_objectives(current,{signal_mode:datasets[signal_mode]},
        baseline_manifests={signal_mode:manifest},input_kind='confirmed_real',
        physical_input_confirmation='Confirmed B7c FCC Al depth/B',verify_repeat=False)
    if prepared_dataset_manifest(datasets[signal_mode])!=manifest:raise ValueError('experimental arrays changed')
    mode=r['modes'][signal_mode]
    return dict(B=float(b),**sample_extent(current),ynorm=mode['ynorm'],objective=mode['objective'],
                active_count=count,active_unique_weights=[1.],invariants='unchanged',
                zero_weight_excluded=mode['decomposition']['zero_weight_contribution_exact'])


def run_depth_profile(compiled,datasets,audit,output_directory,*,resume=False,stop_after=None,cancel_requested=None):
    """15-proposal candidate-depth profile with per-stage journal replay.

    Outer checkpoint never substitutes for child journals; on resume the completed
    prefix is replay-verified by the same PARAMIN engine without physical reruns.
    """
    if audit.get('status')!='PASS' or audit.get('depth_gate')!='READY':raise ValueError('depth audit gate required')
    depths=audit['candidate_depths']
    thicknesses=tuple(s.thickness_angstrom for s in compiled.crystallography.simulation_sample.model.slices)
    base_geometry=requested_depth_to_sample_geometry(840.,thicknesses)
    expected=[base_geometry.plateau_lower_exclusive,840.,float(np.nextafter(base_geometry.effective_depth,np.inf))]
    if depths!=expected:raise ValueError('candidate depths must be audited adjacent geometry states')
    layout=physical_target_layout(compiled)
    if [(p['target_id'],p['lower'],p['upper']) for p in layout]!=[
        ('sample.requested_depth',600.,1200.),('site_group.Al_fcc.debye_waller_b',.5,1.5)]:
        raise ValueError('unexpected two-target layout/bounds')
    manifests={m:prepared_dataset_manifest(d) for m,d in datasets.items()}
    identity=json.loads(json.dumps(dict(canonical=canonical_checkpoint_identity(compiled),
        data_digest=_digest(manifests),target_layout=layout,depths=depths,stage_settings=STAGE_SETTINGS)))
    out=Path(output_directory)
    if resume:
        cp=load_depth_checkpoint(out)
        if cp['identity']!=identity:raise ValueError('B7c checkpoint identity mismatch')
        for row in cp['history']:
            child=out/f"depth_{row['depth_candidate_index']}"/'checkpoint_latest.json'
            if not child.is_file():raise ValueError('missing B7c child journal')
            load_single_target_checkpoint(child)
    else:
        if out.exists():raise FileExistsError(out)
        # Baseline and sensitivity grid precede candidate search.
        initial=evaluate_depth_point(compiled,datasets,840.,.86)
        if (initial['ynorm'],initial['objective'])!=(14870390.0181119,.939675432747781):
            raise ValueError('baseline differs from B6b')
        grid=[evaluate_depth_point(compiled,datasets,d,b) for d in depths for b in (.84,.86,.88)]
        out.mkdir(parents=True)
        cp=dict(format='qcbed-b7c-depth-profile',version=1,identity=identity,status='running',
            baseline=initial,audit=audit,sensitivity_grid=grid,history=[],best=None,last_evaluated=None,
            source_configuration=asdict(compiled.source_configuration),stage=0)
        _save(out,cp)
    prior=len(cp['history']); remaining=stop_after
    for index,d in enumerate(depths):
        stage=out/f'depth_{index}'
        stage_identity=dict(parent=identity,requested_depth=d,
            geometry=asdict(requested_depth_to_sample_geometry(d,thicknesses)))
        stage_identity=json.loads(json.dumps(stage_identity))
        old_count=len(load_single_target_checkpoint(stage)['history']) if stage.exists() else 0
        result=_run_journal(.86,.95,1.05,lambda b:evaluate_depth_point(compiled,datasets,d,b),
            stage,stage_identity,resume=stage.exists(),stop_after=remaining,
            cancel_requested=cancel_requested,settings=STAGE_SETTINGS)
        histories=[]
        for j in range(len(depths)):
            if not (out/f'depth_{j}'/'checkpoint_latest.json').exists():continue
            records=load_single_target_checkpoint(out/f'depth_{j}')['history']
            for row in records:
                # Rejected proposals still retain requested/effective geometry.
                extent=requested_depth_to_sample_geometry(depths[j],thicknesses)
                histories.append(dict(row,evaluation_index=len(histories),depth_candidate_index=j,
                    requested_depth=depths[j],effective_depth=extent.effective_depth,
                    occurrence_count=extent.occurrence_count,stored_slice_count=extent.stored_slice_count))
        valid=[r for r in histories if r['bound_status']=='accepted']
        cp.update(history=histories,best=min(valid,key=lambda r:r['objective']) if valid else None,
            last_evaluated=valid[-1] if valid else None,stage=index,
            status=result['status'] if index==len(depths)-1 or result['status']!='completed' else 'running')
        cp['stage_retained_final']=result.get('retained_final')
        _save(out,cp)
        if remaining is not None:
            remaining-=len(result['history'])-old_count
        if result['status']!='completed':break
        if remaining is not None and remaining<=0 and index<len(depths)-1:
            cp['status']='interrupted';_save(out,cp);break
    if cp['status']=='completed' and cp['stage']==len(depths)-1:
        cp['retained_final']=cp['best'].copy() # Explicit finite-profile selection policy.
    _save(out,cp)
    if cp['best'] is not None:
        from .config_serialization import save_yaml_config
        best=cp['best'];snapshot=apply_crystallographic_parameters(compiled,np.array([best['requested_depth'],best['B']]))
        save_yaml_config(replace(snapshot.source_configuration,title='B7c diagnostic depth/B profile output'),
            out.resolve()/'refined_snapshot.yml',base_directory=out.resolve(),overwrite=True)
        _atomic_json(out/'refined_geometry.json',sample_extent(snapshot))
    return cp


def _save(out,cp):
    body={k:v for k,v in cp.items() if k!='journal_digest'}
    cp['journal_digest']=_digest(body);_atomic_json(out/'checkpoint_latest.json',cp)


def load_depth_checkpoint(path):
    p=Path(path);p=p/'checkpoint_latest.json' if p.is_dir() else p
    try:cp=json.loads(p.read_text())
    except (OSError,json.JSONDecodeError) as exc:raise ValueError(f'Cannot read B7c checkpoint: {exc}') from exc
    if cp.get('format')!='qcbed-b7c-depth-profile' or cp.get('version')!=1:raise ValueError('unsupported B7c checkpoint')
    if cp.get('journal_digest')!=_digest({k:v for k,v in cp.items() if k!='journal_digest'}):raise ValueError('B7c checksum mismatch')
    return cp
