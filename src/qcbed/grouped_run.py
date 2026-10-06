"""Canonical dynamic physical-vector execution using durable stage journals.

No optimizer, propagation or objective mathematics live here. This orchestrator
accepts an immutable RunSnapshot and uses its explicit experimental body, never
an example script, temporary YAML, or a reconstructed experimental dataset.
"""
from dataclasses import asdict
from pathlib import Path
from types import SimpleNamespace
import json
import os
import time
import numpy as np
from .anneal_local_refinement import _run_stage, load_stage_checkpoint
from .canonical_run import canonical_checkpoint_identity
from .crystallography_refinement import physical_target_layout, apply_crystallographic_parameters, site_occurrence_provenance
from .depth_refinement import sample_extent
from .single_target_refinement import _digest
from .full_runner import _atomic_json, _digest as problem_digest

FORMAT='qcbed-grouped-run'


def check_grouped(compiled, experimental=None):
    """Explicit execution capability gate, separate from configuration validity."""
    crystal=compiled.crystallography
    if crystal is None or not crystal.site_groups:
        raise ValueError('Modern refinement requires explicit physical site groups; membership is never inferred.')
    from .config_targets import available_refinement_targets
    for info in available_refinement_targets(compiled):
        if info.active and not info.available: raise ValueError(info.reason)
    from .projected_elastic import check_ready
    check_ready(compiled)
    layout=physical_target_layout(compiled)
    elastic_adapter=any(r['category']=='projected_elastic' for r in layout)
    if not layout: raise ValueError('Select at least one complete supported refinement binding.')
    if compiled.runner.objective_mode!='masked' or compiled.runner.bounds_policy!='reject':
        raise ValueError('Grouped execution requires masked objective and explicit reject bounds')
    if compiled.execution.execution_backend!='batched':
        raise ValueError('Grouped execution currently supports the validated batched backend only')
    for row in layout:
        if row['lower'] is None or row['upper'] is None or not np.isfinite([row['lower'],row['upper']]).all() or row['lower']>=row['upper']:
            raise ValueError(f"{row['target_id']}: explicit finite ordered lower/upper bounds required")
        if not row['lower']<=row['initial_value']<=row['upper']:
            raise ValueError(f"{row['target_id']}: starting value outside bounds; never clipped")
        if row['category'] not in ('scan_coordinate','projected_elastic') and row['lower'] < 0:
            raise ValueError(f"{row['target_id']}: lower bound must be nonnegative")
        if (row['category'] in ('depth','voltage') or (row['category']=='debye_waller' and compiled.state.model.absorption_enabled)) and row['lower'] <= 0:
            raise ValueError(f"{row['target_id']}: lower bound must be positive")
    if compiled.optimizer.loops < 1: raise ValueError('At least one optimizer loop is required')
    if 'AMOEBA' in compiled.optimizer.stages and not elastic_adapter:
        if any(p['initial_value']==0 for p in layout): raise ValueError('AMOEBA scaling requires nonzero physical starting values')
        from .amoeba import build_legacy_simplex
        build_legacy_simplex(compiled.state,compiled.packed_starting_vector,np.asarray(compiled.plamda,dtype=np.float64))
    if 'PARAMIN' in compiled.optimizer.stages and not elastic_adapter:
        if any(row['initial_value']==0 or 0 in p.range_factors for row,p in zip(layout,compiled.layout.parameters)):
            raise ValueError('PARAMIN requires nonzero starting values and D1/D2 divisors')
    if 'SIMPARA' in compiled.optimizer.stages and any(v < 0 for v in compiled.plamda):
        raise ValueError('SIMPARA requires nonnegative absolute uncertainty widths (PLAMDA)')
    if experimental is not None:
        w=experimental.weights;active=w!=0
        if experimental.gradient_option not in (0,1):raise ValueError('Select conventional or gradient explicitly')
        if not np.any(active) or not np.all((w==0)|(w==1)):
            raise ValueError('Grouped experiments require binary inclusion weights; support fractions are not continuous objective weights')
        if not np.isfinite(experimental.variances[active]).all() or np.any(experimental.variances[active]<=0):
            raise ValueError('Positive finite active variances required as statistical denominator')
        if tuple(experimental.beam_indices_hk)!=tuple((b.internal_hk.H,b.internal_hk.K) for b in crystal.indexed_beams):
            raise ValueError('Experimental beam ordering differs from compiled original-hkl mapping')
    return layout


def grouped_identity(snapshot):
    c=snapshot.compiled
    return json.loads(json.dumps(dict(canonical=canonical_checkpoint_identity(c),
        problem=problem_digest(c.state,snapshot.experimental), optimizer=asdict(c.optimizer),
        bounds=asdict(c.bounds), execution=asdict(c.execution),
        original_provenance=site_occurrence_provenance(c))))


def load_grouped_checkpoint(directory, snapshot=None):
    path=Path(directory)/'checkpoint_latest.json'
    try:
        cp=json.loads(path.read_text())
        if not isinstance(cp,dict):raise ValueError('Grouped checkpoint must contain an object')
        if cp.get('format')!=FORMAT or cp.get('version') not in (1,2):raise ValueError('Not a grouped canonical checkpoint; old checkpoints are not reinterpreted')
        if cp.get('digest')!=_digest({k:v for k,v in cp.items() if k!='digest'}):raise ValueError('Grouped checkpoint checksum mismatch')
        if snapshot is not None and cp['identity']!=grouped_identity(snapshot):raise ValueError('Grouped checkpoint configuration/data/ownership/target-layout identity mismatch')
        return cp
    except (OSError,KeyError,TypeError,json.JSONDecodeError) as exc:
        raise ValueError(f'Cannot read grouped checkpoint: {exc}') from exc


def _save(directory,cp):
    cp['digest']=_digest({k:v for k,v in cp.items() if k!='digest'})
    _atomic_json(Path(directory)/'checkpoint_latest.json',cp)


def grouped_history(directory,cp):
    """Read committed child journals, including a child saved before parent update."""
    history=[]
    for index,stage in enumerate(cp['plan']):
        path=Path(directory)/'checkpoints'/f'{index:04d}_{stage}'
        if not (path/'checkpoint_latest.json').exists():break
        child=load_stage_checkpoint(path)
        if child['identity']['problem']!=cp['identity']:raise ValueError('Child journal ownership differs from parent')
        for row in child['history']:
            history.append(dict(row,global_evaluation_index=len(history)+1,
                classification='trial' if row['bound_status']=='accepted' else 'bound_rejection',
                physical_active_vector=row.get('physical_active_vector', [row.get('requested_depth'),row.get('B')]),loop_index=index//len(cp['stages']),
                stage=stage,stage_index=index))
    return history


def evaluate_grouped(snapshot,vector):
    """Existing original-cell rebuild, attached B5 propagation and objective."""
    from .scan import generate_scan_tilts
    from .tilt_series import simulate_tilt_series_batched
    from .refinement_trial import _gradient
    from .masked_objective import compute_selected_objective
    from .objective_diagnostics import residual_diagnostics
    c=apply_crystallographic_parameters(snapshot.compiled,np.asarray(vector,dtype=np.float64))
    exp=snapshot.experimental
    series=simulate_tilt_series_batched(c.crystallography.simulation_sample.prepared,
        generate_scan_tilts(c.state.scan),batch_size=c.execution.batch_size)
    values=_gradient(series.intensities,exp.scan_rows,exp.scan_columns) if exp.gradient_option else series.intensities
    result=compute_selected_objective(values,exp,iref=c.state.iref,objective_mode=c.runner.objective_mode)
    decomposition,arrays=residual_diagnostics(values,exp,result.scaling,result.objective)
    arrays.update(amplitudes=series.amplitudes,calculated=series.intensities,
        scaled=arrays['scaled_calculated'],beams=np.asarray(exp.beam_indices_hk),
        mask=exp.weights!=0,scan_shape=np.array([exp.scan_rows,exp.scan_columns]),
        physical_vector=np.asarray(vector),objective=result.objective,ynorm=result.scaling)
    row=dict(**sample_extent(c),physical_active_vector=np.asarray(vector).tolist(),objective=float(result.objective),ynorm=float(result.scaling),
        active_count=int((exp.weights!=0).sum()),active_unique_weights=[1.],
        signal_mode='gradient' if exp.gradient_option else 'conventional',per_beam=decomposition['beams'])
    return row,arrays


def run_grouped(snapshot,*,resume=False,progress_hook=None,cancel_requested=None):
    c=snapshot.compiled;layout=check_grouped(c,snapshot.experimental)
    from .config_targets import available_refinement_targets
    from .modern_yaml import target_name
    labels={target_name(i.target):i.label for i in available_refinement_targets(c)}
    stored_layout=tuple(dict(p,label=p.get('label',labels.get(p['target_id'],p['target_id']))) for p in layout)
    out=Path(c.configuration.run.output_directory)
    export=Path(c.configuration.run.output_filename)
    root=Path(__file__).resolve().parents[2]
    if not out.is_absolute() or not export.is_absolute():raise ValueError('Resolve output paths against an explicit project base first')
    if any(p.is_relative_to(root/n) for p in (out,export) for n in ('legacy_fortran','reference_data','tests/fixtures')):
        raise ValueError('Output is inside protected source/reference directories')
    if export in (out/'checkpoint_latest.json',out/'summary.json',out/'final_parameters.json',out/'report.html',out/'.run.lock') or export.is_relative_to(out/'checkpoints') or export.is_relative_to(out/'patterns'):
        raise ValueError('Output filename conflicts with internal run artifacts')
    if not resume and export.exists():raise FileExistsError(export)
    if resume:cp=load_grouped_checkpoint(out,snapshot)
    else:
        out.mkdir(parents=True,exist_ok=False)
        cp=dict(format=FORMAT,version=2,identity=grouped_identity(snapshot),status='running',
            stages=list(c.optimizer.stages),plan=list(c.optimizer.stages)*c.optimizer.loops,
            initial=[p['initial_value'] for p in layout],layout=stored_layout,elapsed_seconds=0.,restart_count=0,
            signal_mode='gradient' if snapshot.experimental.gradient_option else 'conventional',
            beam_metadata=[dict(stable_id=b.stable_id,lattice_coordinates=list(b.experimental_index),
                original_hkl=list(b.original_hkl),internal_hk=[b.internal_hk.H,b.internal_hk.K]) for b in c.crystallography.indexed_beams])
        _save(out,cp)
        (out/'patterns').mkdir()
    # A stale lock requires explicit investigation; never bypass concurrent writer protection.
    lock=out/'.run.lock'
    with lock.open('x') as f:f.write(str(os.getpid()))
    started=time.monotonic();prior_elapsed=cp['elapsed_seconds'];new_count=0
    if resume:cp['restart_count']+=1
    initial=cp['initial'];paused_reason=None;paused_stop=None
    def stop():
        nonlocal paused_reason,paused_stop
        if cancel_requested is not None and cancel_requested():paused_reason='cancelled';paused_stop='cancelled';return True
        budget=c.configuration.run.max_total_evaluations
        seconds=c.configuration.run.max_wall_time_seconds
        if (budget is not None and new_count>=budget) or (seconds is not None and time.monotonic()-started>=seconds):
            paused_reason='interrupted';paused_stop='execution_evaluations' if budget is not None and new_count>=budget else 'execution_time';return True
        return False
    def publish(child=None):
        nonlocal new_count
        if child is not None:new_count+=1
        history=grouped_history(out,cp);valid=[r for r in history if r['objective'] is not None]
        cp.update(elapsed_seconds=prior_elapsed+time.monotonic()-started,evaluations=len(history),
            best=min(valid,key=lambda r:r['objective']) if valid else None,
            last_evaluated=valid[-1] if valid else None)
        _save(out,cp)
        if child is not None and progress_hook is not None:
            event=history[-1].copy();event['elapsed_seconds']=cp['elapsed_seconds'];event['restart_count']=cp['restart_count']
            progress_hook(event)
    try:
        cp['status']='running';cp.pop('stop_reason',None)
        cp['execution_limits']=dict(max_total_evaluations=c.configuration.run.max_total_evaluations,max_wall_time_seconds=c.configuration.run.max_wall_time_seconds,max_stage_iterations=c.configuration.run.max_stage_iterations)
        _save(out,cp)
        for index,stage in enumerate(cp['plan']):
            path=out/'checkpoints'/f'{index:04d}_{stage}'
            def evaluate(*values):
                vector = values[0] if cp['version']==2 else values
                row,arrays=evaluate_grouped(snapshot,vector)
                if cp['version']==1: row['B']=float(vector[1])
                cursor=len(load_stage_checkpoint(path)['history'])
                filename=f'patterns/{index:04d}_{cursor:06d}.npz'
                temporary=out/(filename+'.tmp')
                with temporary.open('wb') as f:
                    np.savez_compressed(f,**arrays);f.flush();os.fsync(f.fileno())
                os.replace(temporary,out/filename)
                row['pattern_file']=filename
                return row
            child=_run_stage(c,evaluate,path,optimizer=stage,initial=initial,seed=c.optimizer.seed,
                settings=c.source_configuration.refinement.optimizer,identity=cp['identity'],
                resume=path.exists(),cancel_requested=stop,local_budget=c.optimizer.maxfsimp,
                window=(layout[0]['lower'],layout[0]['upper']),
                parameter_bounds=tuple((p['lower'],p['upper']) for p in layout),progress_hook=publish,dynamic=cp['version']==2)
            cp['stage_result']=child.get('optimizer_result')
            if child['status']!='completed':
                cp['status']=paused_reason or child['status'];cp['stop_reason']=paused_stop or child['status'];publish();break
            cp['retained_final']=child['retained_final']
            selected=child['best'] if stage=='ANNEAL' else child['retained_final']
            initial=selected['physical_active_vector'] if cp['version']==2 else [selected['requested_depth'],selected['B']]
            publish()
        else:
            cp['status']='completed';publish()
            _atomic_json(out/'final_parameters.json',cp['retained_final'])
            from .output_io import LegacySimulatedOutput,write_legacy_out_dat
            with np.load(out/cp['retained_final']['pattern_file'],allow_pickle=False) as a:
                from .refinement import apply_refinement_parameters
                retained_scan=(apply_refinement_parameters(c.state,np.asarray([v for p,v in zip(c.layout.parameters,cp['retained_final']['physical_active_vector']) if p.category!='projected_elastic'],dtype=np.float64)).scan
                    if cp['version']==2 else c.state.scan)
                output=LegacySimulatedOutput(retained_scan,snapshot.experimental.beam_indices_hk,a['scaled'])
            write_legacy_out_dat(out/'final_out.dat',output)
            if export!=out/'final_out.dat':
                export.parent.mkdir(parents=True,exist_ok=True);write_legacy_out_dat(export,output)
        summary=dict(status=cp['status'],final_objective=(cp.get('retained_final') or cp.get('last_evaluated') or {}).get('objective'),
            output_directory=str(out),evaluations=cp.get('evaluations',0),runtime_seconds=cp['elapsed_seconds'])
        from .termination import termination_summary
        summary.update(termination_summary(out,cp))
        _atomic_json(out/'summary.json',summary)
        from html import escape
        (out/'report.html').write_text('<!doctype html><meta charset="utf-8"><title>QCBED grouped results</title><h1>QCBED grouped results</h1><pre>'+escape(json.dumps(dict(summary=summary,signal_mode=cp['signal_mode'],best=cp.get('best'),last_evaluated=cp.get('last_evaluated'),retained_final=cp.get('retained_final')),indent=2))+'</pre>')
        return SimpleNamespace(status=cp['status'],summary=summary)
    except BaseException as exc:
        cp['status']='failed';cp['error']=f'{type(exc).__name__}: {exc}';publish();raise
    finally:lock.unlink()
