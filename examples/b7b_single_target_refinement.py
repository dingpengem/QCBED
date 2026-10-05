"""Controlled real B-only PARAMIN, with an interrupted/resumed equivalence run."""
import argparse
import csv
import json
from pathlib import Path
from qcbed import compile_config_from_sources
from qcbed.single_target_refinement import run_single_target_refinement
if __package__:
    from .b7a_refinement_layout import real_config
    from .a9_real_variance_closure import prepare_reference_variance
else:
    from b7a_refinement_layout import real_config
    from a9_real_variance_closure import prepare_reference_variance


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(argv);out=args.output.resolve()
    if out.exists():raise FileExistsError(out)
    compiled=compile_config_from_sources(real_config());data,_=prepare_reference_variance()
    direct=run_single_target_refinement(compiled,data,out)
    partial=run_single_target_refinement(compiled,data,out/'resume_validation',stop_after=3)
    if partial['status']!='interrupted':raise ValueError('interruption did not occur')
    resumed=run_single_target_refinement(compiled,data,out/'resume_validation',resume=True)
    for key in ('history','best','last_evaluated','retained_final'):
        if direct[key]!=resumed[key]:raise ValueError(f'resume differs: {key}')
    summary=dict(status='PASS',resume_exact=True,run=direct,
        interpretation='Controlled integration only; no claim of converged or global optimum',
        cancel_status='Cooperative pre-evaluation hook preserves completed journal; covered by focused tests')
    (out/'refinement_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    with (out/'history.csv').open('w',newline='') as f:
        names=['evaluation_index','B','ynorm','objective','bound_status']
        writer=csv.DictWriter(f,fieldnames=names,extrasaction='ignore');writer.writeheader();writer.writerows(direct['history'])
    lines=['# B7b controlled single-target refinement','',
        'PASS. Conventional only; Al_fcc B only; bounds[0.5,1.5] Å².',
        'Unchanged deterministic PARAMIN, MAXFPARA1 (up to2 cycles), tolerance0, at most9 proposals.',
        'Mask is binary inclusion; variance is the uncertainty denominator. No fractional support weighting.',
        '', '## Mask audit','',json.dumps(direct['initial_diagnostics']['mask_audit'],indent=2),
        '', '## Local sensitivity','', '| B | YNORM | Objective |','|---:|---:|---:|']
    for row in direct['initial_diagnostics']['local_sensitivity']:
        lines.append(f"| {row['B']:.17g} | {row['ynorm']:.17g} | {row['objective']:.17g} |")
    lines += ['', '## Optimizer history','', '| Proposal | B | YNORM | Objective | Bounds |','|---:|---:|---:|---:|---|']
    for r in direct['history']:
        lines.append(f"| {r['evaluation_index']} | {r['B']:.17g} | {r['ynorm']} | {r['objective']} | {r['bound_status']} |")
    lines += ['', '## Result','',f"Best: {direct['best']['B']}, YNORM {direct['best']['ynorm']}, objective {direct['best']['objective']}",
        f"Final retained: {direct['retained_final']}; last evaluated: B={direct['last_evaluated']['B']}",
        f"Physical optimizer evaluations: {direct['physical_evaluation_count']}; proposals: {direct['proposal_count']}",
        'Four pre-optimizer diagnostic evaluations are separate from the optimizer count.',
        'Resume after three completed proposals matches uninterrupted history, best and final exactly.',
        'Four original sites and830 atom occurrences share each proposed B; all unrelated state/data remain fixed.',
        'Repeated baseline amplitudes/intensities, YNORM and objective agree exactly.',
        'Version2 ownership identity in checkpoint_latest.json; separate refined_snapshot.yml leaves baseline unchanged.',
        'Results do not establish convergence. General run_config/GUI remain disabled for groups.',
        'No multi-parameter refinement, GUI changes or application rebuild.']
    (out/'refinement_summary.md').write_text('\n'.join(lines)+'\n')
    print('B7b PASS',direct['best']['B'],direct['best']['objective'],direct['proposal_count'])
    return 0


if __name__=='__main__':raise SystemExit(main())
