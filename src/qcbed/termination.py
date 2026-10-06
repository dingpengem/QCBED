"""Read-only termination reporting; optimizer cycles are not objective calls."""
from pathlib import Path

LABELS = dict(max_evaluations='Evaluation budget reached',max_iterations='Maximum iterations reached',
    max_cycles='Stage cycle budget reached',objective_spread='Converged',objective_threshold='Converged',converged='Converged',
    temperature='Temperature schedule completed',temperature_schedule='Temperature schedule completed',
    cancelled='Cancelled by user',execution_evaluations='Execution evaluation limit reached',execution_time='Execution wall-time limit reached',
    interrupted='Interrupted — historical record does not specify the limit',failed='Failed',running='Running')


def termination_summary(directory, checkpoint):
    from .anneal_local_refinement import load_stage_checkpoint
    stages=[]
    for i,name in enumerate(checkpoint['plan']):
        path=Path(directory)/'checkpoints'/f'{i:04d}_{name}'
        if not (path/'checkpoint_latest.json').exists():continue
        child=load_stage_checkpoint(path);result=child.get('optimizer_result') or {}
        reason=result.get('stop_reason') or child['status']
        if child['status']!='completed':reason=checkpoint.get('stop_reason') or checkpoint['status']
        count=result.get('iterations',child.get('optimizer_iterations'))
        history=child['history']
        stages.append(dict(stage=i+1,optimizer=name,iterations=count,objective_evaluations=sum(r['objective'] is not None for r in history),
            attempted_evaluations=len(history),stop_reason=LABELS.get(reason,reason),converged=result.get('converged',False)))
    reason=checkpoint.get('stop_reason') or (stages[-1]['stop_reason'] if stages else checkpoint['status'])
    return dict(iterations=sum(s['iterations'] for s in stages) if all(s['iterations'] is not None for s in stages) else None,
        objective_evaluations=sum(s['objective_evaluations'] for s in stages),stop_reason=LABELS.get(reason,reason),stage_termination=stages)
