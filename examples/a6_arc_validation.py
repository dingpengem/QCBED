"""Real-data arc diagnostics; manual 38 px prior, geometric A5, unchanged A6 fit/A7.

Run from the repository root: python examples/a6_arc_validation.py --output tmp/a6_arcs
No synthetic spacing or legacy coordinate table enters this calculation.
"""
import argparse
import json
from pathlib import Path
from qcbed.data_preparation import DiskGeometry, LatticeDetectionConfig
from qcbed.data_preparation.diagnostics import validate_pattern


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pattern',type=Path,default=Path(__file__).resolve().parents[1]/'reference_data/pattern.dat')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    rows=[]
    for tolerance in (0.,1.,2.):
        result=validate_pattern(args.pattern,args.output/f'tolerance_{tolerance:g}',geometry=DiskGeometry(38,38),
            geometry_provenance='Manual visual radius estimate 38 px; matched A5 disk-profile unavailable',
            detection_config=LatticeDetectionConfig(radius_tolerance_px=tolerance))
        rows.append(dict(tolerance_px=tolerance,status=result['status'],a6=result['a6'],a7=result['a7'],
                         diagnostics=result.get('candidate_detection'),blockers=result['blockers']))
    (args.output/'radius_study.json').write_text(json.dumps(rows,indent=2)+'\n')
    lines=['# A6 arc radius-tolerance validation','',
           'Physical semi-axes: manually estimated 38 px; center: canonical geometric localization.',
           'Gradient sigma 1 px, edge percentile 85; equal votes from both normal signs.',
           'No fit is forced to the independent visual ~85 px spacing.','',
           '| Radius tolerance | Candidate count | A6/A7 outcome |','|---:|---:|---|']
    for row in rows:
        lines.append(f"| {row['tolerance_px']:g} | {row['diagnostics']['candidate_count']} | {row['status']}: {'; '.join(row['blockers'])} |")
    lines+=['','Each tolerance directory contains candidate CSV (full coordinates, votes, edge support,',
            'angular coverage, radius residual and central distances), overview, and validation JSON',
            'with every pairwise distance. arc_validation_overview.png distinguishes initial gradient',
            'pixels from validated support; arc_validation.csv records geometric metrics, and',
            'arc_support_pixels.csv preserves exact fragment membership. No inlier/outlier',
            'classification exists when A6 fails.',
            'Minimal seeds produce bootstrap predictions and crops even without global refinement.',
            'Only sufficient independent confirmations trigger the existing A7 global fit.']
    (args.output/'radius_study.md').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines))


if __name__=='__main__':
    main()
