"""B7d workflow adapters; unchanged QCBED ANNEAL and AMOEBA engines.

Checkpoints replay the exact completed prefix from the original seed/initial
simplex, not from best-so-far. No numerical objective is rerun during replay.
"""
from dataclasses import replace, asdict
from pathlib import Path
import json
import numpy as np
from .anneal import _run as anneal
from .amoeba import _run as amoeba, build_legacy_simplex
from .refinement import _physical_vector
from .canonical_run import canonical_checkpoint_identity
from .crystallography_refinement import physical_target_layout
from .depth_refinement import evaluate_depth_point
from .depth_geometry import requested_depth_to_sample_geometry
from .data_preparation import prepared_dataset_manifest
from .single_target_refinement import _digest, audit_binary_masks
from .full_runner import _atomic_json

WINDOW=(1000.,1100.)
SEEDS=(2026,2027,2028)


def case_settings(config):
    """Explicit diagnostic settings; never modify the baseline config on disk.

    Relative initial widths .05 mean sigma(depth)=52.5 Å at midpoint1050.
    Minimum .0001 permits sub-percent adaptation; original cooling .9 retained.
    MAXITS59 permits60 proposals plus source initial/reset evaluations.
    """
    return replace(config.refinement.optimizer,maxtemp=.05,mintemp=.0001,
                   maxfsimp=59,seed=2026)


def boundary_hit(row,window=WINDOW):
    """Flag an endpoint plateau or its immediately adjacent slice interval."""
    d=row['requested_depth'];lo=row['plateau_lower_exclusive'];hi=row['effective_depth']
    # One local slice thickness defines adjacency; no arbitrary Å tolerance.
    width=hi-lo
    return bool(d<=window[0] or d>=window[1] or
                lo<=window[0]+width or hi>=window[1]-width)


def select_global_best(runs):
    if not runs or any(r['status']!='completed' for r in runs):
        raise ValueError('global selection requires completed runs')
    return min((r['best'] for r in runs),key=lambda r:r['objective']).copy()


def _save(path,cp):
    cp['journal_digest']=_digest({k:v for k,v in cp.items() if k!='journal_digest'})
    _atomic_json(path/'checkpoint_latest.json',cp)


def load_stage_checkpoint(path):
    p=Path(path);p=p/'checkpoint_latest.json' if p.is_dir() else p
    try:cp=json.loads(p.read_text())
    except (OSError,json.JSONDecodeError) as e:raise ValueError(f'Cannot read B7d checkpoint: {e}') from e
    if cp.get('format')!='qcbed-b7d-replay' or cp.get('version') not in (1,2):
        raise ValueError('unsupported B7d checkpoint')
    if cp.get('journal_digest')!=_digest({k:v for k,v in cp.items() if k!='journal_digest'}):
        raise ValueError('B7d checksum mismatch')
    return cp


class _Pause(Exception):pass


def _run_stage(compiled,evaluate,path,*,optimizer,initial,seed,settings,identity,
               resume=False,stop_after=None,cancel_requested=None,local_budget=40,window=WINDOW,parameter_bounds=None,progress_hook=None,dynamic=False):
    if optimizer not in (('ANNEAL','AMOEBA','PARAMIN','SIMPARA') if dynamic else ('ANNEAL','AMOEBA')):raise ValueError('unsupported stage')
    # Reuse the production control validator before invoking private engine seams.
    replace(compiled.optimizer, **asdict(settings))
    if type(local_budget) is not int or local_budget < 0:
        raise ValueError('local_budget must be a nonnegative integer')
    custom_bounds=parameter_bounds is not None
    if dynamic:
        limits = np.asarray(parameter_bounds, dtype=np.float64)
        initial = np.asarray(initial, dtype=np.float64)
        if limits.shape != (compiled.parameter_count, 2) or not np.isfinite(limits).all() or np.any(limits[:,0] >= limits[:,1]):
            raise ValueError('Explicit finite ordered bounds required for every active parameter')
        if initial.shape != (compiled.parameter_count,) or not np.isfinite(initial).all() or np.any(initial < limits[:,0]) or np.any(initial > limits[:,1]):
            raise ValueError('Initial vector must match layout and lie within bounds')
    else:
        depth_bounds,b_bounds=((600.,1200.),(.5,1.5)) if parameter_bounds is None else parameter_bounds
        if any(len(b)!=2 or not all(np.isfinite(x) for x in b) or not 0<b[0]<b[1] for b in (depth_bounds,b_bounds)):
            raise ValueError('explicit positive finite ordered physical bounds required')
        window=tuple(window)
        if len(window)!=2 or not all(np.isfinite(x) for x in window) or not depth_bounds[0]<=window[0]<window[1]<=depth_bounds[1]:
            raise ValueError('case window must lie within legal depth bounds')
        initial=np.asarray(initial,dtype=np.float64)
        if initial.shape!=(2,) or not np.isfinite(initial).all() or not (
            window[0]<=initial[0]<=window[1] and b_bounds[0]<=initial[1]<=b_bounds[1]):
            raise ValueError('initial state must be within case window/B bounds')
    if stop_after is not None and (type(stop_after) is not int or stop_after<1):
        raise ValueError('stop_after must be a positive integer')
    if type(seed) is not int or seed<0:raise ValueError('seed must be a nonnegative integer')
    if dynamic:
        identity = json.loads(json.dumps(dict(problem=identity,optimizer=optimizer,initial=initial.tolist(),
            seed=seed,settings=asdict(settings),parameter_bounds=limits.tolist(),local_budget=local_budget,
            layout=physical_target_layout(compiled))))
    else:
        identity=json.loads(json.dumps(dict(problem=identity,optimizer=optimizer,initial=initial.tolist(),
            seed=seed,settings=asdict(settings),window=window,b_bounds=list(b_bounds),local_budget=local_budget)))
        if custom_bounds:identity['parameter_bounds']=[list(depth_bounds),list(b_bounds)]
    path=Path(path)
    if resume:
        cp=load_stage_checkpoint(path)
        if cp['identity']!=identity:raise ValueError('B7d checkpoint identity mismatch')
    else:
        path.mkdir(parents=True,exist_ok=False)
        cp=dict(format='qcbed-b7d-replay',version=2 if dynamic else 1,identity=identity,status='running',
            history=[],best=None,last_evaluated=None,optimizer_acceptance='not exposed by engine; best-so-far recorded')
        _save(path,cp)
    history=cp['history'];prefix=len(history);cursor=0;new=0
    rng=np.random.Generator(np.random.PCG64(seed))
    thicknesses=tuple(s.thickness_angstrom for s in compiled.crystallography.simulation_sample.model.slices)
    def visit(vector):
        nonlocal cursor,new
        v=np.asarray(vector,dtype=np.float64)
        if not np.isfinite(v).all():raise ValueError('nonfinite optimizer proposal')
        rng_state=json.loads(json.dumps(rng.bit_generator.state))
        if cursor<prefix:
            row=history[cursor]
            if v.tolist()!=(row['physical_active_vector'] if dynamic else [row['requested_depth'],row['B']]) or row['rng_state']!=rng_state:
                raise ValueError(f'optimizer replay mismatch at {cursor}')
            cursor+=1
            return row['objective'] if row['bound_status']=='accepted' else float('inf')
        if cancel_requested is not None and cancel_requested():
            cp['status']='cancelled';raise _Pause()
        if stop_after is not None and new>=stop_after:
            cp['status']='interrupted';raise _Pause()
        if dynamic:
            accepted = bool(np.all(v >= limits[:,0]) and np.all(v <= limits[:,1]))
            row = dict(evaluation_index=cursor,stage=optimizer,physical_active_vector=v.tolist(),
                objective=None,ynorm=None,bound_status='accepted' if accepted else 'rejected',
                rng_state=rng_state,optimizer_acceptance=None)
        else:
            accepted=window[0]<=v[0]<=window[1] and b_bounds[0]<=v[1]<=b_bounds[1]
            geometry=requested_depth_to_sample_geometry(v[0],thicknesses) if v[0]>0 else None
            row=dict(evaluation_index=cursor,stage=optimizer,requested_depth=float(v[0]),B=float(v[1]),
                effective_depth=None if geometry is None else geometry.effective_depth,
                occurrence_count=None if geometry is None else geometry.occurrence_count,
                stored_slice_count=len(thicknesses),sequence_digest=None if geometry is None else _digest(geometry.sequence),
                plateau_lower_exclusive=None if geometry is None else geometry.plateau_lower_exclusive,
                objective=None,ynorm=None,bound_status='accepted' if accepted else 'rejected',
                rng_state=rng_state,optimizer_acceptance=None)
        if accepted:
            row.update(evaluate(v.copy()) if dynamic else evaluate(float(v[0]),float(v[1])))
            if not np.isfinite(row['objective']) or row['objective']<0 or not np.isfinite(row['ynorm']):
                raise ValueError('invalid physical objective')
        row['best_so_far']=bool(accepted and (cp['best'] is None or row['objective']<=cp['best']['objective']))
        history.append(row);cursor+=1;new+=1
        if row['best_so_far']:cp['best']=row.copy()
        if accepted:cp['last_evaluated']=row.copy()
        cp['last_proposal']=row.copy();cp['rng_state']=rng_state;cp['status']='running'
        _save(path,cp)
        if progress_hook is not None:progress_hook(cp)
        return row['objective'] if accepted else float('inf')
    try:
        if optimizer=='ANNEAL':
            result=anneal(initial,lambda v:v.copy(),visit,max_iterations=settings.maxfsimp,
                maximum_temperature=settings.maxtemp,minimum_temperature=settings.mintemp,
                cooling_rate=settings.coolrate,rng=rng,allow_invalid=True)
            retained=result.best_vector
            detail=dict(iterations=result.iterations,evaluations=result.evaluations,
                temperature=float(result.temperature),amplitudes=result.amplitudes.tolist(),
                convergence='ANNEAL has no convergence flag')
        elif optimizer=='AMOEBA':
            simplex=build_legacy_simplex(compiled.state,initial,np.asarray(compiled.plamda,dtype=np.float64))
            result=amoeba(simplex,lambda p:visit(_physical_vector(compiled.layout,p,'legacy_simplex')),
                ftol=settings.ftol,max_evaluations=local_budget,fmin=settings.chimin,allow_invalid=True)
            retained=_physical_vector(compiled.layout,result.best_vector,'legacy_simplex')
            detail=dict(iterations=result.iterations,evaluations=result.evaluations,
                stop_reason=result.stop_reason,converged=result.converged,
                final_simplex=result.final_simplex.tolist(),final_objectives=[float(x) if np.isfinite(x) else None for x in result.final_objectives])
        elif optimizer == 'PARAMIN':
            from .paramin import _run as paramin
            from .refinement import apply_refinement_parameters, pack_refinement_parameters
            def project(v):
                return v.copy()  # Modern explicit bounds reject proposals before any physical folding.
            result=paramin(initial,np.array([p.range_factors[1] for p in compiled.layout.parameters]),
                np.array([p.range_factors[0] for p in compiled.layout.parameters]),project,visit,
                max_cycles=settings.maxfpara,tolerance=settings.ftol,random_order=bool(settings.paraopt),rng=rng,allow_invalid=True)
            retained=result.best_vector
            detail=dict(evaluations=result.evaluations,converged=result.converged)
        else:
            from .simpara import _run as simpara
            result=simpara(initial,np.asarray(compiled.plamda,dtype=np.float64),visit,
                max_evaluations=settings.maxfsimp,high_temperature=settings.maxtemp,
                low_temperature=settings.mintemp,cooling=settings.coolrate,
                random_order=bool(settings.paraopt),rng=rng,allow_invalid=True)
            retained=result.best_vector
            detail=dict(evaluations=result.evaluations,stop_reason=result.stop_reason)
        if cursor!=len(history):raise ValueError('unreplayed history')
        matching=[r for r in history if (r['physical_active_vector'] if dynamic else [r['requested_depth'],r['B']])==retained.tolist() and r['objective'] is not None]
        if not matching:raise ValueError('retained state lacks recorded evaluation')
        cp.update(status='completed',retained_final=matching[-1].copy(),optimizer_result=detail)
    except _Pause:pass
    except Exception as exc:
        cp.update(status='failed',error=f'{type(exc).__name__}: {exc}');_save(path,cp);raise
    _save(path,cp)
    return cp


def run_optimizer_stage(compiled,datasets,path,*,optimizer,initial,seed,settings=None,signal_mode="conventional",**kwargs):
    """Canonical physical two-target adapter, preserving B7b/B7c checkpoints."""
    layout=physical_target_layout(compiled)
    if [(p['target_id'],p['lower'],p['upper']) for p in layout]!=[
        ('sample.requested_depth',600.,1200.),('site_group.Al_fcc.debye_waller_b',.5,1.5)]:
        raise ValueError('requires B7c physical layout and legal bounds')
    if signal_mode not in ('conventional','gradient'):raise ValueError('unsupported signal mode')
    audit=audit_binary_masks(datasets)
    if audit[signal_mode]['active_count']!=(7974 if signal_mode=='conventional' else 7281):raise ValueError('mask changed')
    manifests={k:prepared_dataset_manifest(v) for k,v in datasets.items()}
    identity=dict(canonical=canonical_checkpoint_identity(compiled),data=_digest(manifests),layout=layout)
    if signal_mode=='gradient':
        # Opt-in identity extension; historical conventional checkpoint identities stay unchanged.
        from .crystallography_refinement import site_occurrence_provenance
        identity.update(signal_mode=signal_mode,plamda=compiled.plamda,
            range_factors=compiled.range_factors,site_groups=site_occurrence_provenance(compiled))
    def evaluate(d,b):
        if {k:prepared_dataset_manifest(v) for k,v in datasets.items()}!=manifests:raise ValueError('experimental data changed')
        result=evaluate_depth_point(compiled,datasets,d,b,signal_mode=signal_mode)
        if signal_mode=='gradient':result['signal_mode']=signal_mode
        return result
    return _run_stage(compiled,evaluate,path,optimizer=optimizer,initial=initial,seed=seed,
        settings=case_settings(compiled.source_configuration) if settings is None else settings,
        identity=identity,**kwargs)
