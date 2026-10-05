"""B7e audit-only gate and fixed gradient baseline; never starts optimization."""
import argparse
import json
from pathlib import Path
import numpy as np
from qcbed import compile_config_from_sources
from qcbed.depth_refinement import enable_requested_depth
from qcbed.objective_diagnostics import evaluate_preflight_objectives
from examples.b7a_refinement_layout import real_config
from examples.a9_real_variance_closure import prepare_reference_variance
from examples.b6a_fcc_al_preflight import CONFIRMATION


def run_audit():
    config = enable_requested_depth(real_config())
    compiled = compile_config_from_sources(config)
    datasets, closure = prepare_reference_variance()
    report, arrays = evaluate_preflight_objectives(compiled,
        {'gradient': datasets['gradient']},
        baseline_manifests={'gradient': closure['modes']['gradient']['dataset']},
        input_kind='confirmed_real', physical_input_confirmation=CONFIRMATION,
        verify_repeat=False)
    baseline = report['modes']['gradient']
    weights = arrays['gradient']['weights']
    active_weights = np.unique(weights[weights != 0]).tolist()
    if baseline['active_count'] != 7281 or active_weights != [1.0]:
        raise ValueError('Gradient binary-mask baseline changed')
    expected = np.array([13447206.106343817, 1.1896911365294378])
    actual = np.array([baseline['ynorm'], baseline['objective']])
    np.testing.assert_allclose(actual, expected, rtol=1e-14, atol=0)
    if compiled.source_configuration != config:
        raise ValueError('Canonical configuration changed')
    return dict(status='PASS (audit only)', dilation_status='DEFERRED',
        signal_mode='gradient', optimizer_executed=False,
        legacy_semantics=dict(source='legacy_fortran/oneit.f:893-916',
            units='dimensionless vacancy-coupling coefficient',
            coefficient_domain=[0., 1.], baseline_coefficient=compiled.state.dilation,
            applied_scale='(1 - VAC * DILATION) ** (1/3) when 0 < VAC*DILATION <= 0.5 and RTYPE=S; otherwise 1',
            current_zero_vacancy_scale=1.,
            mapping='abs; reciprocal if greater than 1 (not clipping or hard rejection)'),
        global_lattice_scale_bounds=None,
        gate_reason='Legacy vacancy coefficient is not an independent original-cell lattice scale. No validated global lattice-scale bounds were found. At zero vacancy the legacy coefficient has no physical sensitivity.',
        rebuild_status='Not implemented or claimed: existing physical application supports group B and requested depth only. No third target enabled.',
        baseline=baseline, active_weights=active_weights,
        baseline_absolute_errors=np.abs(actual-expected).tolist(),
        canonical_configuration_unchanged=True,
        next_step='Define a separate original-cell lattice-scale/strain target with explicit units, baseline and scientifically justified bounds before any three-target optimization.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    report = run_audit()
    (args.output/'audit_summary.json').write_text(json.dumps(report, indent=2)+'\n')
    b = report['baseline']
    (args.output/'audit_summary.md').write_text(
        '# B7e gradient audit: PASS (audit only)\n\n'
        'Dilation **DEFERRED**. '+report['gate_reason']+'\n\n'
        'Source: `legacy_fortran/oneit.f:893–916`; dimensionless coefficient domain 0–1; '
        'baseline coefficient 0, applied scale 1. This domain is not a global lattice-scale bound.\n\n'
        f"Gradient baseline: {b['active_count']} active points, weights [1.0], "
        f"YNORM {b['ynorm']!r}, objective {b['objective']!r}. "
        f"Absolute regression errors: {report['baseline_absolute_errors']}.\n\n"
        'No ANNEAL/AMOEBA, perturbation sweep, refined snapshots, checkpoints or best-state figures were generated. '
        'Existing checkpoints and canonical inputs are unchanged. Conventional scores are not gradient references.\n\n'
        +report['next_step']+'\n')
    print(json.dumps({k: report[k] for k in ('status','dilation_status','baseline_absolute_errors')}))


if __name__ == '__main__':
    main()
