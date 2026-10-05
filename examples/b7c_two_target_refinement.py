"""Depth audit, sensitivity grid, and bounded depth-candidate/B refinement."""
import argparse
import csv
import json
from pathlib import Path
from qcbed import compile_config_from_sources
from qcbed.depth_refinement import enable_requested_depth,run_depth_profile
if __package__:
    from .b7a_refinement_layout import real_config
    from .a9_real_variance_closure import prepare_reference_variance
    from .b7c_depth_audit import depth_audit
else:
    from b7a_refinement_layout import real_config
    from a9_real_variance_closure import prepare_reference_variance
    from b7c_depth_audit import depth_audit


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(argv);out=args.output.resolve()
    if out.exists():raise FileExistsError(out)
    audit=depth_audit() # Gate is evaluated before activating the depth binding.
    compiled=compile_config_from_sources(enable_requested_depth(real_config()))
    data,_=prepare_reference_variance()
    direct=run_depth_profile(compiled,data,audit,out)
    partial=run_depth_profile(compiled,data,audit,out/'resume_validation',stop_after=7)
    if partial['status']!='interrupted':raise ValueError('interruption gate failed')
    resumed=run_depth_profile(compiled,data,audit,out/'resume_validation',resume=True)
    for key in ('history','best','last_evaluated','retained_final'):
        if direct[key]!=resumed[key]:raise ValueError(f'B7c resume differs: {key}')
    summary=dict(status='PASS',depth='READY',run=direct,resume_exact=True,
        scope='three explicit depth geometries; one PARAMIN cycle for B at each; not global convergence')
    (out/'refinement_summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    for name,rows in [('depth_audit',audit['rows']),('sensitivity_grid',direct['sensitivity_grid']),('history',direct['history'])]:
        with (out/f'{name}.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    lines=['# B7c controlled depth/B refinement','',
        'PASS; requested depth READY for explicit candidate search only. Canonical baseline unchanged.',
        'Vector order: [sample.requested_depth, site_group.Al_fcc.debye_waller_b]. Bounds[600,1200] Å and[0.5,1.5] Å².',
        'Strict-< whole-slice sequence; no fractional slices or depth interpolation.',
        '', '## Depth audit','', '| Requested Å | Effective Å | Occurrences | YNORM | Objective |','|---:|---:|---:|---:|---:|']
    for r in audit['rows']:lines.append(f"| {r['requested_depth']:.17g} | {r['effective_depth']:.17g} | {r['occurrence_count']} | {r['ynorm']:.17g} | {r['objective']:.17g} |")
    lines += ['', '## Sensitivity grid','', '| B | Requested Å | Effective Å | Objective |','|---:|---:|---:|---:|']
    for r in direct['sensitivity_grid']:lines.append(f"| {r['B']} | {r['requested_depth']:.17g} | {r['effective_depth']:.17g} | {r['objective']:.17g} |")
    for label in ('baseline','best','last_evaluated','retained_final'):
        lines+=['',f'## {label}','', '```json',json.dumps(direct[label],indent=2),'```']
    lines+=['','Search budget: three depth candidates × five PARAMIN proposals =15; one source cycle each.',
        'The nine-point sensitivity grid and baseline precede search. Depth audit is separate. No evaluation cache was added.',
        'Mask remains7974 active points, weight1; variance and all unrelated input state remain unchanged.',
        'Interrupted after7 proposals, resumed history/best/last/final match exactly.',
        'Outer checkpoint version1 (B7c-specific format) includes canonical ownership identity version2 and depth provenance.',
        'Keep child depth journals with the outer checkpoint. Existing B7b checkpoint format/identity unchanged.',
        'Cooperative cancellation preserves completed child journal state; no force-kill.',
        'Refined YAML and derived geometry sidecar are separate from baseline input.',
        'No global optimum/convergence claim; no lattice or additional targets, GUI or app rebuild.']
    (out/'refinement_summary.md').write_text('\n'.join(lines)+'\n')
    print('B7c PASS',direct['best'])
    return 0


if __name__=='__main__':raise SystemExit(main())
