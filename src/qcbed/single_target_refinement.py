"""B7b controlled PARAMIN integration, separate from general production runs.

Unchanged PARAMIN algorithm; physical group application; binary inclusion; durable
proposal replay checkpoints. No multi-target, GUI, or general runner enablement.
"""
from dataclasses import replace, asdict
import hashlib
import json
from pathlib import Path
import numpy as np
from .canonical_run import canonical_checkpoint_identity
from .crystallography_refinement import apply_crystallographic_parameters, physical_target_layout
from .objective_diagnostics import evaluate_preflight_objectives
from .data_preparation import prepared_dataset_manifest
from .full_runner import _atomic_json
from .paramin import _run as _paramin

TARGET = 'site_group.Al_fcc.debye_waller_b'
SETTINGS = dict(optimizer='PARAMIN', max_cycles=1, tolerance=0., random_order=False,
                maximum_optimizer_proposals=9, bounds=[.5,1.5])


def _digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,allow_nan=False).encode()).hexdigest()


def audit_binary_masks(datasets):
    """B7b input contract only; do not change generic legacy weighted objectives."""
    result={}
    for mode,data in datasets.items():
        w=data.weights
        if data.variances is None:
            raise ValueError(f"{mode}: validated variance is required")
        if not np.isfinite(w).all() or not np.all((w==0)|(w==1)):
            raise ValueError(f'{mode}: support must be binary inclusion, not fractional confidence')
        active=w>0
        if not np.any(active) or not np.all(np.isfinite(data.variances[active])) or np.any(data.variances[active]<=0):
            raise ValueError(f'{mode}: active variance must be finite and positive')
        result[mode]=dict(active_count=int(active.sum()),active_unique_weights=np.unique(w[active]).tolist(),
                         zero_weight_count=int((w==0).sum()),meaning='binary inclusion; variance is denominator')
    return result


class _Pause(Exception):
    pass


def load_single_target_checkpoint(path):
    """Validate this explicit B7b journal; never reinterpret old runner checkpoints."""
    path=Path(path)
    if path.is_dir():path=path/'checkpoint_latest.json'
    try:
        cp=json.loads(path.read_text())
        if cp['format']!='qcbed-b7b-paramin-replay' or cp['version']!=2:
            raise ValueError('unsupported B7b checkpoint format/version')
        body={k:v for k,v in cp.items() if k!='journal_digest'}
        if cp['journal_digest']!=_digest(body):raise ValueError('B7b checkpoint checksum mismatch')
        for i,row in enumerate(cp['history']):
            if row['evaluation_index']!=i or not np.isfinite(row['B']):raise ValueError('invalid history order/value')
            accepted=.5<=row['B']<=1.5
            if (row['bound_status']=='accepted')!=accepted:raise ValueError('invalid bound status')
            if accepted and (not np.isfinite(row['objective']) or row['objective']<0 or not np.isfinite(row['ynorm'])):
                raise ValueError('invalid history objective')
            if not accepted and (row['objective'] is not None or row['ynorm'] is not None):
                raise ValueError('rejected proposal must not contain fabricated physics')
        return cp
    except (OSError,KeyError,TypeError,json.JSONDecodeError) as exc:
        raise ValueError(f'Cannot read valid B7b checkpoint: {exc}') from exc


def _write_checkpoint(directory,cp):
    body={k:v for k,v in cp.items() if k!='journal_digest'}
    cp['journal_digest']=_digest(body)
    _atomic_json(directory/'checkpoint_latest.json',cp)


def _run_journal(initial,bottom,top,evaluate,directory,identity,*,resume=False,stop_after=None,cancel_requested=None,
                 initial_diagnostics=None, settings=None):
    """Restart deterministic optimizer and replay validated completed proposals.

    The journal stores all answers needed to reconstruct PARAMIN's local state.
    No physics is recomputed for replay entries. Cancellation checks occur before
    a new evaluation, after the last complete result has been durably written.
    """
    settings = SETTINGS if settings is None else settings
    directory=Path(directory)
    if stop_after is not None and (type(stop_after) is not int or stop_after<1):
        raise ValueError('stop_after must be a positive number of new proposals')
    if resume:
        cp=load_single_target_checkpoint(directory)
        if cp['identity']!=identity or cp['settings']!=settings:
            raise ValueError('B7b checkpoint configuration/data/ownership identity mismatch')
    else:
        directory.mkdir(parents=True,exist_ok=False)
        cp=dict(format='qcbed-b7b-paramin-replay',version=2,identity=identity,settings=settings,
            history=[],best=None,last_evaluated=None,status='running',initial_diagnostics=initial_diagnostics,
            resume_count=0)
        _write_checkpoint(directory,cp)
    history=cp['history']; replay_count=len(history); cursor=0; new_count=0
    if resume:cp['resume_count']+=1
    def visit(vector):
        nonlocal cursor,new_count
        b=float(vector[0])
        if cursor<replay_count:
            row=history[cursor]
            if b!=row['B']:raise ValueError(f'PARAMIN replay proposal mismatch at {cursor}')
            cursor+=1
            return row['objective'] if row['bound_status']=='accepted' else float('inf')
        if cancel_requested is not None and cancel_requested():
            cp['status']='cancelled';raise _Pause()
        if stop_after is not None and new_count>=stop_after:
            cp['status']='interrupted';raise _Pause()
        if cursor>=settings['maximum_optimizer_proposals']:
            raise RuntimeError('PARAMIN exceeded declared controlled budget')
        if not np.isfinite(b):raise ValueError('nonfinite PARAMIN proposal')
        accepted=.5<=b<=1.5
        row=dict(evaluation_index=cursor,B=b,bound_status='accepted' if accepted else 'rejected',
                 objective=None,ynorm=None)
        if accepted:
            result=evaluate(b)
            if not np.isfinite(result['objective']) or result['objective']<0 or not np.isfinite(result['ynorm']):
                raise ValueError('invalid objective; no finite substitute allowed')
            row.update(result)
        history.append(row);cursor+=1;new_count+=1
        if accepted:
            cp['last_evaluated']=row.copy()
            if cp['best'] is None or row['objective']<=cp['best']['objective']:cp['best']=row.copy()
        cp['status']='running'
        _write_checkpoint(directory,cp)
        return row['objective'] if accepted else float('inf')
    try:
        # Same source-derived loop used by run_legacy_paramin. Physical proposals
        # are not abs-folded: B7a bounds reject before any sample rebuild.
        with np.errstate(divide='raise',invalid='raise',over='raise'):
            result=_paramin(np.array([initial]),np.array([bottom]),np.array([top]),
                lambda v:v.copy(),visit,max_cycles=settings['max_cycles'],
                tolerance=settings['tolerance'],random_order=False,allow_invalid=True)
        if cursor!=len(history):raise ValueError('checkpoint contains unreplayed extra proposals')
        cp.update(status='completed',retained_final=dict(B=float(result.best_vector[0]),objective=float(result.best_objective)),
                  optimizer_cycles=result.cycles,optimizer_converged=result.converged)
    except _Pause:
        pass
    except Exception as exc:
        cp['status']='failed';cp['error']=f'{type(exc).__name__}: {exc}'
        _write_checkpoint(directory,cp)
        raise
    _write_checkpoint(directory,cp)
    return cp


def _problem(compiled,datasets):
    layout=physical_target_layout(compiled)
    if len(layout)!=1 or layout[0]['target_id']!=TARGET or (layout[0]['lower'],layout[0]['upper'])!=(.5,1.5):
        raise ValueError('B7b requires exactly Al_fcc B with bounds[0.5,1.5]')
    if compiled.packed_starting_vector.tolist()!=[.86]:raise ValueError('B7b baseline must remain B0.86')
    audit=audit_binary_masks(datasets)
    if set(audit)!={'conventional','gradient'} or audit['conventional']['active_count']!=7974 or audit['gradient']['active_count']!=7281:
        raise ValueError('real conventional/gradient support counts changed')
    manifests={m:prepared_dataset_manifest(d) for m,d in datasets.items()}
    identity=dict(canonical=canonical_checkpoint_identity(compiled),data_digest=_digest(manifests),
                  target_layout=list(layout),signal_mode='conventional')
    # Normalize tuples so written/read checkpoint identities compare exactly.
    return audit,manifests,json.loads(json.dumps(identity))


def run_single_target_refinement(compiled,datasets,output_directory,*,resume=False,stop_after=None,cancel_requested=None):
    """Controlled backend Run/Resume/Cancel for B7b only, never general run_config."""
    audit,manifests,identity=_problem(compiled,datasets)
    destination=Path(output_directory)
    if not resume and destination.exists():raise FileExistsError(destination)
    if resume:
        existing=load_single_target_checkpoint(destination)
        if existing['identity']!=identity:raise ValueError('B7b checkpoint configuration/data/ownership identity mismatch')
    baseline_config=compiled.source_configuration
    baseline_sample=compiled.crystallography.simulation_sample
    def physical_evaluation(b,*,arrays=False):
        current=apply_crystallographic_parameters(compiled,np.array([b],dtype=np.float64))
        # Undo only B in the source snapshot to prove all other config fields fixed.
        cfg=current.source_configuration
        restored=replace(cfg,crystallography=replace(cfg.crystallography,
            site_groups=tuple(replace(g,debye_waller_b=original.debye_waller_b)
                for g,original in zip(cfg.crystallography.site_groups,baseline_config.crystallography.site_groups))))
        if restored!=baseline_config:raise ValueError('unrelated canonical state changed')
        sample=current.crystallography.simulation_sample
        for changed,original in zip(current.crystallography.original_structure.sites,
                                    compiled.crystallography.original_structure.sites,strict=True):
            if replace(changed,debye_waller_b=original.debye_waller_b)!=original:
                raise ValueError('unrelated original site state changed')
        np.testing.assert_array_equal(current.crystallography.original_structure.direct_lattice.vectors,
                                      compiled.crystallography.original_structure.direct_lattice.vectors)
        sites=sample.depth_sample.original_structure.sites+sample.depth_sample.recast_structure.sites+tuple(
            s for slab in sample.depth_sample.slices for s in slab.sites)
        if any(s.debye_waller_b!=b or s.physical_group_id!='Al_fcc' for s in sites):
            raise ValueError('group B propagation failure')
        if sample.parameters.slice_sequence_zero_based!=baseline_sample.parameters.slice_sequence_zero_based:
            raise ValueError('depth sequence changed')
        if sum(len(sample.depth_sample.slices[i].sites) for i in sample.parameters.slice_sequence_zero_based)!=830:
            raise ValueError('simulation occurrence count changed')
        if audit_binary_masks(datasets)!=audit or _digest({m:prepared_dataset_manifest(d) for m,d in datasets.items()})!=identity['data_digest']:
            raise ValueError('experimental arrays/support changed')
        report,values=evaluate_preflight_objectives(current,{'conventional':datasets['conventional']},
            baseline_manifests={'conventional':manifests['conventional']},input_kind='confirmed_real',
            physical_input_confirmation='B7b confirmed FCC Al [001] (200)/(020)',verify_repeat=False)
        r=report['modes']['conventional']
        result=dict(ynorm=r['ynorm'],objective=r['objective'],active_count=r['active_count'],
            active_unique_weights=[1.],zero_weight_excluded=r['decomposition']['zero_weight_contribution_exact'],
            decomposition=r['decomposition'],unrelated_state_unchanged=True,
            original_sites_affected=4,simulation_atom_occurrences_affected=830)
        return (result,values['conventional']) if arrays else result
    diagnostics=None
    if not resume:
        baseline,a=physical_evaluation(.86,arrays=True)
        if not np.allclose([baseline['ynorm'],baseline['objective']],
                          [14870390.0181119,.939675432747781],rtol=1e-14,atol=0):
            raise ValueError('baseline differs from B6b; optimizer not started')
        repeat,z=physical_evaluation(.86,arrays=True)
        errors={name:float(np.max(np.abs(a[name]-z[name]))) for name in
                ('complex_amplitudes','calculated_intensities','objective_intensities')}
        if any(errors.values()) or baseline!=repeat:raise ValueError('deterministic repeat failed')
        delta=.0002
        local=[dict(B=.86-delta,**physical_evaluation(.86-delta)),dict(B=.86,**baseline),
               dict(B=.86+delta,**physical_evaluation(.86+delta))]
        diagnostics=dict(source_configuration=asdict(baseline_config),
            original_structure=dict(lattice=compiled.crystallography.original_structure.direct_lattice.vectors.tolist(),
                sites=[asdict(s) for s in compiled.crystallography.original_structure.sites]),
            site_groups=[asdict(g) for g in compiled.crystallography.site_groups],
            slice_sequence=list(baseline_sample.parameters.slice_sequence_zero_based),
            mask_audit=audit,baseline=baseline,local_sensitivity=local,repeat_errors=errors,
                         pre_optimizer_physical_evaluations=4)
    bottom=compiled.range_factors[0][1];top=compiled.range_factors[0][0]
    cp=_run_journal(.86,bottom,top,physical_evaluation,destination,identity,resume=resume,
        stop_after=stop_after,cancel_requested=cancel_requested,initial_diagnostics=diagnostics)
    cp['physical_evaluation_count']=sum(r['bound_status']=='accepted' for r in cp['history'])
    cp['proposal_count']=len(cp['history'])
    _write_checkpoint(destination,cp)
    if cp['best'] is not None:
        from .config_serialization import save_yaml_config
        best=apply_crystallographic_parameters(compiled,np.array([cp['best']['B']]))
        output=replace(best.source_configuration,title='B7b diagnostic refinement output — best recorded B')
        save_yaml_config(output,destination.resolve()/'refined_snapshot.yml',
                         base_directory=destination.resolve(),overwrite=True)
    return cp
