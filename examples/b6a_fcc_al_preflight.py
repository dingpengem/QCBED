"""Confirmed real FCC Al preflight; reuses B6a, never propagates a wave."""
import argparse
import csv
from pathlib import Path
import numpy as np
from qcbed import load_yaml_config, compile_config_from_sources
from qcbed.input_preflight import preflight_inputs
from qcbed.crystallography import (
    legacy_for007_001_basis_to_conventional_hkl,
    fcc_conventional_reflection_diagnostic,
)
if __package__:
    from .a9_real_variance_closure import prepare_reference_variance
    from .b6a_real_input_preflight import write_report
else:
    from a9_real_variance_closure import prepare_reference_variance
    from b6a_real_input_preflight import write_report

ROOT = Path(__file__).resolve().parents[1]
CONFIRMATION = 'User-confirmed FCC Al conventional cell, [001], ordered (200)/(020).'


def real_inputs():
    config = load_yaml_config(ROOT/'examples/fcc_al_real_preflight.yml')
    compiled = compile_config_from_sources(config)
    datasets, closure = prepare_reference_variance()
    report = preflight_inputs(compiled, datasets,
        baseline_manifests={m: r['dataset'] for m, r in closure['modes'].items()},
        input_kind='confirmed_real', physical_input_confirmation=CONFIRMATION)
    for row in report['beam_table']:
        hkl = tuple(row[f'original_{v}'] for v in 'hkl')
        if hkl != (2*row['n1'], 2*row['n2'], 0):
            raise ValueError('confirmed conventional (200)/(020) assignment changed')
        row['fcc_diagnostic'] = fcc_conventional_reflection_diagnostic(hkl)
    report['physical_provenance'] = {
        'structure': 'examples/fcc_al_conventional.extxyz; user-confirmed four-site FCC motif, full occupancy',
        'a0_angstrom': 4.049,
        'lattice_source': 'tests/fixtures/legacy_reference/FOR007.DAT record 4: 2 * layer thickness 2.0245; confirmed FCC [001] half-period layers',
        'thermal_source': 'tests/fixtures/legacy_reference/FOR007.DAT B=0.86 angstrom squared',
        'occupancy': 'Four fully occupied conventional sites. Legacy projected half weights are not partial chemical occupancy.',
        'other_reference_settings': 'FOR007: 200 kV, mesh16x16, depth840 A, DSAP2.5 A^-1, absorption (0,1,1)',
        'experimental_source': 'pattern.dat and sigma.dat via unchanged A9 variance closure; archived A7/A8/A10',
    }
    report['legacy_mapping'] = [dict(p=p, q=q, conventional_hkl=legacy_for007_001_basis_to_conventional_hkl(p,q))
                                for p,q in ((1,1),(1,-1),(1,0),(0,1))]
    return report, datasets


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv)
    report, datasets = real_inputs()
    write_report(report, args.output)
    with (args.output/'legacy_for007_mapping.csv').open('w', newline='') as f:
        writer = csv.writer(f); writer.writerow(('legacy_p','legacy_q','conventional_h','conventional_k','conventional_l'))
        for r in report['legacy_mapping']:
            writer.writerow((r['p'],r['q'],*r['conventional_hkl']))
    # Output-only serialization: beam-major blocks, public arrays (tilt,beam).
    # Useful future input, but no production runner is invoked here.
    data = datasets['conventional']
    with (args.output/'conventional_body.dat').open('w') as f:
        for values in (data.intensities, data.variances, data.weights):
            np.savetxt(f, values.T.reshape(-1, data.scan_geometry.scan_columns), fmt='%.17g')
    with (args.output/'preflight.md').open('a') as f:
        f.write('\n## Confirmed physical provenance\n\n')
        for key,value in report['physical_provenance'].items():
            f.write(f'- {key}: {value}\n')
        f.write('\n## FCC diagnostic (no beam filtering)\n\n')
        for r in report['beam_table']:
            f.write(f"- {r['stable_id']}: {r['fcc_diagnostic']}\n")
        f.write('\n## Legacy internal compatibility representation → canonical conventional hkl\n\n')
        for r in report['legacy_mapping']:
            f.write(f"- ({r['p']},{r['q']}) → {r['conventional_hkl']}\n")
    print(f"B6a {report['status']}: {args.output/'preflight.md'}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
