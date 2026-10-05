"""Real FCC Al: existing ANNEAL (three seeds) -> existing AMOEBA."""
import argparse
import csv
import json
from dataclasses import replace,asdict
from pathlib import Path
import numpy as np
from qcbed import compile_config_from_sources
from qcbed.depth_refinement import enable_requested_depth,evaluate_depth_point
from qcbed.anneal_local_refinement import run_optimizer_stage,case_settings,SEEDS,select_global_best,boundary_hit
from qcbed.crystallography_refinement import apply_crystallographic_parameters
from qcbed.config_serialization import save_yaml_config
from qcbed.objective_diagnostics import evaluate_preflight_objectives
from qcbed.data_preparation import prepared_dataset_manifest
if __package__:
    from .b7a_refinement_layout import real_config
    from .a9_real_variance_closure import prepare_reference_variance
else:
    from b7a_refinement_layout import real_config
    from a9_real_variance_closure import prepare_reference_variance


def write_json(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path)
    p.add_argument('--resume',action='store_true',help='Replay compatible stage journals in an existing output directory')
    p.add_argument('--depth-window',nargs=2,type=float,default=(1000.,1100.))
    args=p.parse_args(argv);out=args.output.resolve()
    window=tuple(args.depth_window);initial_depth=sum(window)/2
    if args.resume:
        if not out.is_dir():raise ValueError('resume requires an existing output directory')
        if not (out/'checkpoints'/f'anneal_{SEEDS[0]}'/'checkpoint_latest.json').is_file():
            raise ValueError('resume requires a first-stage checkpoint to verify configuration/data identity')
    else:out.mkdir(parents=True,exist_ok=False)
    c=compile_config_from_sources(enable_requested_depth(real_config()));data,_=prepare_reference_variance()
    settings=case_settings(c.source_configuration)
    if args.resume:
        previous=json.loads((out/'pre_profile.json').read_text())
        if previous['settings']!=asdict(settings) or tuple(previous.get('window',(1000.,1100.)))!=window:raise ValueError('profile settings mismatch')
        baseline,sparse=previous['baseline'],previous['profile']
    else:
        baseline=evaluate_depth_point(c,data,840.,.86)
        sparse=[evaluate_depth_point(c,data,d,.86) for d in np.linspace(*window,5)]
    assert (baseline['ynorm'],baseline['objective'])==(14870390.0181119,.939675432747781)
    write_json(out/'pre_profile.json',dict(baseline=baseline,profile=sparse,settings=asdict(settings),window=window))
    print('Baseline/profile PASS',flush=True)
    runs=[]
    for seed in SEEDS:
        cp=run_optimizer_stage(c,data,out/'checkpoints'/f'anneal_{seed}',optimizer='ANNEAL',
            initial=[initial_depth,.86],seed=seed,settings=settings,window=window,
            resume=args.resume and (out/'checkpoints'/f'anneal_{seed}').exists())
        runs.append(cp);write_json(out/'anneal_summary.json',runs)
        print('ANNEAL',seed,'best',cp['best']['requested_depth'],cp['best']['B'],cp['best']['objective'],flush=True)
    best=select_global_best(runs);boundary=boundary_hit(best,window)
    def snapshot(name,row):
        changed=apply_crystallographic_parameters(c,np.array([row['requested_depth'],row['B']]))
        save_yaml_config(replace(changed.source_configuration,title=f'B7d {name} diagnostic'),
            out/f'{name}_snapshot.yml',base_directory=out,overwrite=args.resume)
        write_json(out/f'{name}_state.json',dict(state=row,settings=asdict(settings),case_window=window))
    snapshot('anneal_best',best)
    local=None
    if not boundary:
        local=run_optimizer_stage(c,data,out/'checkpoints'/'local',optimizer='AMOEBA',
            initial=[best['requested_depth'],best['B']],seed=2026,settings=settings,window=window,
            resume=args.resume and (out/'checkpoints'/'local').exists())
        snapshot('local_best',local['best'])
        print('AMOEBA best',local['best']['requested_depth'],local['best']['B'],local['best']['objective'],flush=True)
    else:print('Boundary hit: local stage deliberately deferred; revisit search window',flush=True)
    replay={}
    for name,opt,initial,reference in [('anneal','ANNEAL',[initial_depth,.86],runs[0]),
        *([] if local is None else [('local','AMOEBA',[best['requested_depth'],best['B']],local)])]:
        path=out/'resume_validation'/name
        if not path.exists():
            partial=run_optimizer_stage(c,data,path,optimizer=opt,initial=initial,seed=2026,settings=settings,window=window,stop_after=12)
            assert partial['status']=='interrupted'
        resumed=run_optimizer_stage(c,data,path,optimizer=opt,initial=initial,seed=2026,settings=settings,window=window,resume=True)
        for key in ('history','best','last_evaluated','retained_final','optimizer_result'):
            if resumed[key]!=reference[key]:raise ValueError(f'{name} replay differs: {key}')
        replay[name]='exact';print(name,'real replay PASS',flush=True)
    beam_reports={}
    for label,row in [('baseline',baseline),('anneal',best),*([] if local is None else [('local',local['best'])])]:
        current=apply_crystallographic_parameters(c,np.array([row['requested_depth'],row['B']]))
        r,artifacts=evaluate_preflight_objectives(current,{'conventional':data['conventional']},
            baseline_manifests={'conventional':prepared_dataset_manifest(data['conventional'])},
            input_kind='confirmed_real',physical_input_confirmation='Confirmed FCC Al B7d',verify_repeat=False)
        mode=r['modes']['conventional']
        assert mode['objective']==row['objective'] and mode['ynorm']==row['ynorm']
        beam_reports[label]=mode['decomposition']['beams']
        from qcbed.masked_pattern_diagnostics import save_masked_pattern_match
        stem={'baseline':'baseline','anneal':'anneal_best','local':'local_best'}[label]
        save_masked_pattern_match(artifacts['conventional'],c.state.scan.scan_rows,c.state.scan.scan_columns,
            [b['stable_id'] for b in beam_reports[label]],out/'figures'/f'{stem}_masked_pattern_match.png',
            title=stem,ynorm=mode['ynorm'])
    write_json(out/'per_beam_residuals.json',beam_reports)
    residual_rows=[]
    for i,b in enumerate(beam_reports['baseline']):
        row=dict(beam_id=b['stable_id'],h=b['original_h'],k=b['original_k'],l=b['original_l'])
        for stage,beams in beam_reports.items():
            row[stage]=beams[i]['squared_objective_contribution']
            row[stage+'_delta']=row[stage]-b['squared_objective_contribution']
        residual_rows.append(row)
    with (out/'per_beam_residuals.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(residual_rows[0]));w.writeheader();w.writerows(residual_rows)
    summary=dict(status='PASS',window=window,initial_depth=initial_depth,baseline=baseline,settings=asdict(settings),seeds=SEEDS,sparse_profile=sparse,
        anneal_runs=[dict(seed=s,best=r['best'],last_evaluated=r['last_evaluated'],retained_final=r['retained_final'],
            result=r['optimizer_result'],boundary_hit=boundary_hit(r['best'],window)) for s,r in zip(SEEDS,runs)],
        global_best=best,selection='lowest finite valid evaluated ANNEAL objective across seeds',boundary_hit=boundary,
        local=local,resume=replay,per_beam=beam_reports)
    write_json(out/'refinement_summary.json',summary)
    for name,records in [('anneal_history',[dict(row,seed=seed) for seed,r in zip(SEEDS,runs) for row in r['history']]),
                         ('local_history',[] if local is None else local['history'])]:
        if records:
            keys=list(dict.fromkeys(k for row in records for k in row))
            with (out/f'{name}.csv').open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(records)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    (out/'figures').mkdir(exist_ok=True)
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for seed,r in zip(SEEDS,runs):
        valid=[row for row in r['history'] if row['objective'] is not None]
        axes[0].scatter([x['requested_depth'] for x in valid],[x['objective'] for x in valid],s=10,label=str(seed))
        axes[1].scatter([x['B'] for x in valid],[x['objective'] for x in valid],s=10)
        accepted=[x for x in valid if x['best_so_far']]
        axes[0].scatter([x['requested_depth'] for x in accepted],[x['objective'] for x in accepted],marker='x',s=30)
    axes[0].set_xlabel('Requested depth (Å)');axes[1].set_xlabel('B (Å²)')
    for ax in axes:ax.set_ylabel('Recorded objective')
    axes[0].legend(title='Seed; crosses = best-so-far')
    fig.tight_layout();fig.savefig(out/'figures'/'landscape.png',dpi=150);plt.close(fig)
    lines=['# B7d ANNEAL → AMOEBA','', 'PASS for workflow integration; no physical convergence claim.',
        f'Canonical baseline unchanged. Case depth window {window} Å; legal depth bounds [600,1200] Å.',
        'ANNEAL uses unchanged source relative-normal/cooling rules, not Metropolis acceptance.',
        f'Initial [{initial_depth},0.86], seeds2026/2027/2028; explicit settings:', '```json',json.dumps(asdict(settings),indent=2),'```',
        'Local AMOEBA chosen for plateaus; both parameters remain free. Budget MXF40, existing D1/D2 and PLAMDA.',
        'Bounds reject, never clip. Bound-accepted is not an optimizer acceptance claim; engine does not expose each decision.',
        '', '| State | Requested Å | Effective Å | B Å² | YNORM | Objective |','|---|---:|---:|---:|---:|---:|']
    for label,row in [('Baseline',baseline),*[(f'ANNEAL seed{s}',r['best']) for s,r in zip(SEEDS,runs)],
                      *([] if local is None else [('Local best',local['best']),('Local retained final',local['retained_final'])])]:
        lines.append(f"| {label} | {row['requested_depth']:.17g} | {row['effective_depth']:.17g} | {row['B']:.17g} | {row['ynorm']:.17g} | {row['objective']:.17g} |")
    lines+=['',f'Global boundary hit: {boundary}. Exact real replay: {replay}.',
        'Per-beam squared-objective contributions are in per_beam_residuals.json; no per-beam fitted weighting.',
        'Binary support7974, weight1, frozen experimental variance/signal/order/corners verified at each physical evaluation.',
        'Checkpoint replay reconstructs RNG and optimizer state from initial state plus completed answers; does not restart from best.',
        'No lattice/occupancy/absorption/GUI/packaging work.']
    (out/'refinement_summary.md').write_text('\n'.join(lines)+'\n')
    print('B7d PASS',flush=True)


if __name__=='__main__':main()
