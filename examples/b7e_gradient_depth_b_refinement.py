"""Controlled gradient-only depth/B ANNEAL -> AMOEBA; existing engines only."""
import argparse
import csv
import json
from dataclasses import asdict, replace
from pathlib import Path
import numpy as np
from qcbed import compile_config_from_sources
from qcbed.anneal_local_refinement import (run_optimizer_stage, case_settings, SEEDS,
    select_global_best, boundary_hit)
from qcbed.depth_refinement import enable_requested_depth, evaluate_depth_point
from qcbed.crystallography_refinement import apply_crystallographic_parameters, physical_target_layout
from qcbed.amoeba import build_legacy_simplex
from qcbed.refinement import _physical_vector
from qcbed.config_serialization import save_yaml_config
from qcbed.data_preparation import prepared_dataset_manifest
from qcbed.objective_diagnostics import evaluate_preflight_objectives
from qcbed.masked_pattern_diagnostics import save_masked_pattern_match
from examples.b7a_refinement_layout import real_config
from examples.a9_real_variance_closure import prepare_reference_variance

WINDOW = (900., 1100.)
LOCAL_BUDGET = 100


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def write_csv(path, rows):
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader(); writer.writerows(rows)


def gradient_config(body_path):
    c = enable_requested_depth(real_config())
    return replace(c, experiment=replace(c.experiment, path=str(body_path), gradient_option=1))


def local_config(config, best):
    # SCAL remains anchored at canonical [840,.86], D1/D2=(1.05,.95).
    # PLAMDA magnitudes .05/.25 give 4.2 Å / .0215 Å² physical steps.
    # Choose direction toward bound midpoint; no clipping or input adjustment.
    steps = {'depth': -.05 if best['requested_depth'] < 1000. else .05,
             'debye_waller': -.25 if best['B'] < 1. else .25}
    return replace(config, refinement=replace(config.refinement, parameters=tuple(
        replace(p, plamda=steps[p.target.category.value]) for p in config.refinement.parameters)))


def physical_simplex(compiled, best):
    initial = np.array([best['requested_depth'], best['B']])
    return np.array([_physical_vector(compiled.layout, row, 'legacy_simplex')
        for row in build_legacy_simplex(compiled.state, initial, np.asarray(compiled.plamda))])


def boundary_status(row):
    depth = 'interior'
    if boundary_hit(row, WINDOW):
        depth = 'near 900 Å' if row['requested_depth'] < 1000. else 'near 1100 Å'
    # Explicit diagnostic tolerance: one nominal local B simplex step, not a bound change.
    b = 'near 0.5 Å²' if row['B'] <= .5215 else 'near 1.5 Å²' if row['B'] >= 1.4785 else 'interior'
    return dict(depth=depth, B=b, B_nearness_tolerance=.0215)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--resume', action='store_true')
    args = p.parse_args(argv); out = args.output.resolve()
    if args.resume:
        if not (out/'checkpoints'/f'anneal_{SEEDS[0]}'/'checkpoint_latest.json').is_file():
            raise ValueError('resume requires an existing first-stage checkpoint')
    else:
        out.mkdir(parents=True, exist_ok=False)
    all_data, _ = prepare_reference_variance(); data = {'gradient': all_data['gradient']}
    d = data['gradient']; manifest = prepared_dataset_manifest(d)
    body = out/'gradient_body.dat'
    if not args.resume:
        with body.open('x') as f:
            for a in (d.intensities, d.variances, d.weights):
                np.savetxt(f, a.T.reshape(-1, d.scan_geometry.scan_columns), fmt='%.17g')
    config = gradient_config(body); c = compile_config_from_sources(config)
    settings = case_settings(config)
    evaluate = lambda compiled, depth, b: evaluate_depth_point(compiled, data, depth, b, signal_mode='gradient')
    baseline = evaluate(c, 840., .86)
    np.testing.assert_allclose([baseline['ynorm'], baseline['objective']],
        [13447206.106343817, 1.1896911365294378], rtol=1e-14, atol=0)
    print('Gradient baseline PASS', baseline, flush=True)
    sparse = [evaluate(c, depth, .86) for depth in (900., 950., 1000., 1050., 1100.)]
    write_csv(out/'sparse_depth_profile.csv', sparse)
    print('Sparse gradient profile', [(r['requested_depth'],r['objective']) for r in sparse], flush=True)
    runs = []
    for seed in SEEDS:
        path = out/'checkpoints'/f'anneal_{seed}'
        cp = run_optimizer_stage(c, data, path, optimizer='ANNEAL', initial=[1000.,.86],
            seed=seed, settings=settings, window=WINDOW, signal_mode='gradient',
            resume=args.resume and path.exists())
        runs.append(cp)
        write_json(out/'anneal_summary.json', runs)
        print('ANNEAL',seed,'best',cp['best'],flush=True)
    best = select_global_best(runs)
    lc = compile_config_from_sources(local_config(config, best))
    simplex = physical_simplex(lc, best)
    simplex_rows = [evaluate(lc, *v) for v in simplex]
    if simplex_rows[0]['occurrence_count'] == simplex_rows[1]['occurrence_count']:
        raise ValueError('Local depth simplex step does not change physical geometry')
    if simplex_rows[0]['objective'] == simplex_rows[2]['objective']:
        raise ValueError('Local B simplex step does not change gradient objective')
    if any(not (900<=v[0]<=1100 and .5<=v[1]<=1.5) for v in simplex):
        raise ValueError('Local simplex outside case bounds')
    write_json(out/'simplex_audit.json',dict(plamda=lc.plamda,range_factors=lc.range_factors,
        physical_vertices=simplex.tolist(),evaluations=simplex_rows))
    print('Simplex adequacy PASS',simplex.tolist(),flush=True)
    path = out/'checkpoints'/'local'
    local = run_optimizer_stage(lc,data,path,optimizer='AMOEBA',initial=[best['requested_depth'],best['B']],
        seed=2026,settings=settings,window=WINDOW,signal_mode='gradient',local_budget=LOCAL_BUDGET,
        resume=args.resume and path.exists())
    print('AMOEBA',local['best'],local['optimizer_result'],flush=True)
    replay = {}
    for name,compiled,opt,initial,reference in (
        ('anneal',c,'ANNEAL',[1000.,.86],runs[0]),
        ('local',lc,'AMOEBA',[best['requested_depth'],best['B']],local)):
        path = out/'resume_validation'/name
        common = dict(optimizer=opt, initial=initial,seed=2026,settings=settings,window=WINDOW,
                      signal_mode='gradient',local_budget=LOCAL_BUDGET if name=='local' else 40)
        if not path.exists():
            partial = run_optimizer_stage(compiled,data,path,**common,stop_after=12)
            if partial['status']!='interrupted':raise ValueError('expected interrupted validation run')
        resumed = run_optimizer_stage(compiled,data,path,**common,resume=True)
        for key in ('history','best','last_evaluated','retained_final','optimizer_result'):
            if resumed[key]!=reference[key]:raise ValueError(f'{name} replay differs: {key}')
        replay[name]='exact completed-prefix/history/best/last/final/optimizer replay'
        print(name,'real replay PASS',flush=True)
    states = [('baseline',c,baseline),('anneal_best',c,best),('local_best',lc,local['best'])]
    beam_reports = {}; masks = []
    for label,compiled,row in states:
        changed = apply_crystallographic_parameters(compiled,np.array([row['requested_depth'],row['B']]))
        if label!='baseline':
            save_yaml_config(replace(changed.source_configuration,title=f'B7e gradient {label} diagnostic'),
                out/f'{label}_snapshot.yml',base_directory=out,overwrite=args.resume)
            write_json(out/f'{label}_state.json',dict(signal_mode='gradient',state=row,
                parameter_layout=physical_target_layout(changed),case_window=WINDOW))
        r,arrays = evaluate_preflight_objectives(changed,data,baseline_manifests={'gradient':manifest},
            input_kind='confirmed_real',physical_input_confirmation='Confirmed FCC Al B7e gradient depth/B',verify_repeat=False)
        mode = r['modes']['gradient']
        if (mode['ynorm'],mode['objective'])!=(row['ynorm'],row['objective']):
            raise ValueError('Presentation reevaluation differs')
        beam_reports[label] = mode['decomposition']['beams']
        labels = [b['stable_id'] for b in beam_reports[label]]
        dest = out/'figures'/f'gradient_{label}_masked_pattern_match.png'
        displayed = save_masked_pattern_match(arrays['gradient'],c.state.scan.scan_rows,c.state.scan.scan_columns,
            labels,dest,title=f'Gradient {label}',ynorm=mode['ynorm'])
        # Extend the diagnostic NPZ with ordered physical beam provenance (no rescaling).
        np.savez_compressed(dest.with_suffix('.npz'),**displayed,ynorm=mode['ynorm'],
            beam_ids=np.asarray(labels),lattice_coordinates=np.asarray(d.lattice_coordinates),
            conventional_hkl=np.asarray([[b['original_h'],b['original_k'],b['original_l']] for b in beam_reports[label]]))
        masks.append(displayed)
    for m in masks[1:]:
        np.testing.assert_array_equal(m['mask'],masks[0]['mask'])
        np.testing.assert_array_equal(m['experimental'],masks[0]['experimental'])
    rows = []
    for stage,beams in beam_reports.items():
        for i,b in enumerate(beams):
            delta = b['squared_objective_contribution']-beam_reports['baseline'][i]['squared_objective_contribution']
            rows.append(dict(stage=stage,stable_id=b['stable_id'],n1=d.lattice_coordinates[i][0],n2=d.lattice_coordinates[i][1],
                h=b['original_h'],k=b['original_k'],l=b['original_l'],active_count=b['active_count'],
                squared_objective_contribution=b['squared_objective_contribution'],fraction_of_total_squared_objective=b['fraction_of_total_error'],
                delta_vs_gradient_baseline=delta,change='improved' if delta<0 else 'worsened' if delta>0 else 'unchanged'))
    write_csv(out/'per_beam_residuals.csv',rows)
    write_json(out/'per_beam_residuals.json',beam_reports)
    write_csv(out/'anneal_history.csv',[dict(row,seed=seed) for seed,run in zip(SEEDS,runs) for row in run['history']])
    write_csv(out/'local_history.csv',local['history'])
    if prepared_dataset_manifest(d)!=manifest:raise ValueError('experimental dataset changed')
    if c.source_configuration!=config:raise ValueError('canonical source mutated')
    summary = dict(status='PASS',signal_mode='gradient',parameter_layout=physical_target_layout(c),
        baseline=baseline,sparse_profile=sparse,settings=asdict(settings),seeds=SEEDS,window=WINDOW,
        legal_depth_bounds=[600.,1200.],B_bounds=[.5,1.5],initial=[1000.,.86],
        anneal_runs=[dict(seed=s,best=r['best'],last_evaluated=r['last_evaluated'],retained_final=r['retained_final'],
            optimizer_result=r['optimizer_result'],boundary=boundary_status(r['best'])) for s,r in zip(SEEDS,runs)],
        global_best=best,local=local,local_budget=LOCAL_BUDGET,local_plamda=lc.plamda,
        simplex_audit=simplex_rows,boundary=boundary_status(local['best']),resume=replay,
        dataset_manifest=manifest,experimental_unchanged=True,canonical_unchanged=True,
        per_beam_changes={s:{k:sum(r['stage']==s and r['change']==k for r in rows) for k in ('improved','worsened','unchanged')}
            for s in ('anneal_best','local_best')},
        convergence='Workflow PASS is not physical convergence; inspect boundary and optimizer stop independently.',
        historical_comparison='Conventional best depth 941.963 Å, B .74035 Å², conventional objective .803944: not directly comparable to gradient objective.',
        dilation='DEFERRED legacy vacancy coupling; no independent lattice-scale target')
    write_json(out/'refinement_summary.json',summary)
    lines=['# B7e gradient depth/B ANNEAL → AMOEBA','',summary['convergence'],'',
        'Exactly depth then Al_fcc B; a=4.049 Å, occupancy1, absorption and geometry unchanged.',
        'Gradient mask7281, active weights[1.0]; support>0.875 is binary inclusion; variance is denominator.',
        f'ANNEAL settings: {asdict(settings)}; seeds {SEEDS}; initial [1000,.86].',
        f'Local PLAMDA {lc.plamda}, D1/D2 {lc.range_factors}, MXF {LOCAL_BUDGET} (source soft limit).',
        '', '| State | Requested Å | Effective Å | Occurrences | B Å² | YNORM | Gradient objective |',
        '|---|---:|---:|---:|---:|---:|---:|']
    for label,row in [('Baseline',baseline),*[(f'Seed {s}',r['best']) for s,r in zip(SEEDS,runs)],
                      ('Local best',local['best']),('Local retained final',local['retained_final'])]:
        lines.append(f"| {label} | {row['requested_depth']:.17g} | {row['effective_depth']:.17g} | {row['occurrence_count']} | {row['B']:.17g} | {row['ynorm']:.17g} | {row['objective']:.17g} |")
    lines += ['',f"Local stop: {local['optimizer_result']['stop_reason']}; boundaries: {summary['boundary']}.",
        f"Per-beam changes: {summary['per_beam_changes']}.",f'Replay: {replay}.',
        'Histories distinguish bound acceptance from unexposed optimizer acceptance; no acceptance decisions invented.',
        'Snapshots are diagnostic state/configuration with the exported gradient body; no GUI/production-run enablement is implied.',
        summary['historical_comparison'],summary['dilation']]
    (out/'refinement_summary.md').write_text('\n'.join(lines)+'\n')
    print('B7e gradient depth/B PASS',flush=True)


if __name__=='__main__':main()
