"""Read-only modern journal/NPZ adapter. Never calls a forward evaluator."""
from pathlib import Path
from zipfile import BadZipFile
import numpy as np
from qcbed.grouped_run import load_grouped_checkpoint,grouped_history
from .results import RunResults


def load_grouped_results(directory,*,evaluation=None):
    directory=Path(directory).resolve();cp=load_grouped_checkpoint(directory)
    history=grouped_history(directory,cp)
    valid=[r for r in history if r['objective'] is not None]
    selected=next((r for r in valid if r['global_evaluation_index']==evaluation),None) if evaluation is not None else (valid[-1] if valid else None)
    if evaluation is not None and selected is None:raise ValueError('No completed evaluation with that identity')
    best=min(valid,key=lambda r:r['objective']) if valid else None
    final=cp.get('retained_final')
    current=selected or {};vector=current.get('physical_active_vector',cp['initial'])
    summary=dict(status=cp['status'],partial=cp['status']!='completed',output_directory=str(directory),
        checkpoint_status='Ownership-aware grouped journal; compatible resume required',
        final_state_status='Present' if (directory/'final_parameters.json').exists() else 'Unavailable / partial',
        retained_objective=current.get('objective'),retained_ynorm=current.get('ynorm'),
        displayed_evaluation=current.get('global_evaluation_index'),latest_valid_evaluation=valid[-1]['global_evaluation_index'] if valid else None,
        historical=bool(valid and selected is not valid[-1]),recorded_evaluations=len(history),
        recorded_elapsed_seconds=cp['elapsed_seconds'],recorded_stage=current.get('stage'),restart_count=cp['restart_count'],
        signal_mode=cp['signal_mode'],initial_parameters=cp['initial'],best=best,
        last_evaluated=valid[-1] if valid else None,retained_final=final,
        requested_depth=current.get('requested_depth'),effective_depth=current.get('effective_depth'),
        optimizer_result=cp.get('stage_result'),convergence=('Converged' if (cp.get('stage_result') or {}).get('converged') else 'Not converged')+' (optimizer criterion only; no physical convergence claim)')
    from qcbed.termination import termination_summary
    summary.update(termination_summary(directory,cp))
    parameters=tuple(dict(label=p.get('label',p['target_id']).replace('Linked coefficient group ','Shared refinement parameter group '),target_id=p['target_id'],category=p['category'],units=p['units'],
        initial=p['initial_value'],final=v if selected else None,delta=v-p['initial_value'] if selected else None,
        latest=valid[-1]['physical_active_vector'][i] if valid else None,
        best=best['physical_active_vector'][i] if best else None,
        retained=(final.get('physical_active_vector', [final.get('requested_depth'),final.get('B')])[i] if final else None),
        bounds=dict(lower=p['lower'],upper=p['upper'])) for i,(p,v) in enumerate(zip(cp['layout'],vector)))
    from qcbed.projected_elastic import results_coefficients
    def vector_of(row):return row.get('physical_active_vector') if row else None
    summary['projected_elastic']=results_coefficients(cp['layout'],cp['initial'],displayed=vector_of(selected),
        latest=vector_of(valid[-1] if valid else None),best=vector_of(best),retained=vector_of(final))
    comparison=None;warnings=[]
    if current.get('pattern_file'):
        try:
            path=directory/current['pattern_file']
            if not path.resolve().is_relative_to(directory):raise ValueError('Pattern path escapes run directory')
            with np.load(path,allow_pickle=False) as a:
                if not np.array_equal(a['physical_vector'],vector) or float(a['objective'])!=current['objective'] or float(a['ynorm'])!=current['ynorm']:
                    raise ValueError('Pattern and journal state do not match')
                mask=a['mask']
                comparison=dict(experimental=np.where(mask,a['experimental'],np.nan),
                    scaled=np.where(mask,a['scaled'],np.nan),mask=mask.copy(),beams=a['beams'].copy(),scan_shape=a['scan_shape'].copy(),
                    beam_metadata=cp['beam_metadata'],per_beam=current['per_beam'])
            for a in comparison.values():
                if isinstance(a,np.ndarray):a.setflags(write=False)
        except (OSError,ValueError,KeyError,BadZipFile,EOFError) as exc:warnings.append(str(exc))
    if comparison is None:warnings.append('No stored pattern for this state; no simulation was recomputed.')
    return RunResults(directory,summary,parameters,tuple(history),comparison,tuple(warnings))
