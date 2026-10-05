"""B7a physical group/bounds bring-up: fixed baseline and one B perturbation, no fit."""
from dataclasses import replace, asdict
from pathlib import Path
import argparse
import json
import numpy as np
from qcbed import load_yaml_config, compile_config_from_sources
from qcbed.configuration import (AtomicSiteGroup, AbsorptionConfig, ParameterBinding,
    BoundsConfig, TargetBound)
from qcbed.crystallography.site_groups import group_target
from qcbed.crystallography_refinement import (apply_crystallographic_parameters,
    physical_target_layout, site_occurrence_provenance)
from qcbed.objective_diagnostics import evaluate_preflight_objectives
from qcbed.canonical_run import canonical_checkpoint_identity
if __package__:
    from .b6a_fcc_al_preflight import ROOT, CONFIRMATION
    from .a9_real_variance_closure import prepare_reference_variance
else:
    from b6a_fcc_al_preflight import ROOT, CONFIRMATION
    from a9_real_variance_closure import prepare_reference_variance

BOUNDS_SOURCE = 'examples/reference_configuration.py: Debye-Waller category absolute bounds 0.5–1.5 Å²; B starting0.86, PLAMDA0.0002, D1/D2(1.05,0.95). Project bring-up limits, not universal Al limits.'


TARGET_AUDIT = (
    ('Debye–Waller B', 'READY', 'site-group', 'Existing species B scattering; physical group alias and rebuild validated. Project bounds0.5–1.5 Å².'),
    ('Global lattice/dilation', 'DEFERRED', 'global original lattice', 'Legacy folding/derived dilation is not an original-cell lattice rebuild; no new physical lattice bounds justified.'),
    ('Requested sample depth', 'DEFERRED', 'sample/region', 'Existing strict-< cyclic whole-slice rule, not continuous effective depth. Reference600–1200 Å bounds exist; intentionally excluded from first layout.'),
    ('Occupancy/vacancy', 'DEFERRED', 'site-group', 'Full chemical occupancy metadata supported; legacy vacancy mapping is distinct and nonzero vacancy trials unsupported.'),
    ('Absorption a/b/c', 'DEFERRED', 'site-group to species', 'Fixed independent (0,1,1) is explicit. Shared species coefficients work; B7a group parameter application enables B only. Linked C remains unsupported in B4b.'),
    ('Accelerating voltage', 'DEFERRED', 'global microscope', 'Existing wavelength binding; voltage changes with absorption rejected by legacy trial semantics. No reviewed voltage bounds.'),
    ('Scan A/B/C', 'DEFERRED', 'global scan', 'Legacy targets work; preserve validated A10 geometry in first bring-up.'),
    ('Independent atom x/y', 'NOT APPLICABLE', 'future symmetry modes', 'Independent simulation-site edits violate FCC group symmetry; future symmetry-aware structural modes required.'),
    ('Reflection multipliers', 'NOT APPLICABLE', 'original reflection records', 'Current real case has no explicit records, unity multipliers; not used to absorb physical model mismatch.'),
)


def real_config():
    baseline=load_yaml_config(ROOT/'examples/fcc_al_real_preflight.yml')
    target=group_target('Al_fcc')
    return replace(baseline,crystallography=replace(baseline.crystallography,
        site_groups=(AtomicSiteGroup('Al_fcc','Al',(0,1,2,3),.86,1.,(0.,1.,1.)),),
        debye_waller=(),absorption=AbsorptionConfig()),
        refinement=replace(baseline.refinement,
            parameters=(ParameterBinding(target,(1.05,.95),.0002),),active_categories=None,
            bounds=BoundsConfig('reject',(TargetBound(target,.5,1.5),))))


def run_diagnostics():
    config=real_config();compiled=compile_config_from_sources(config)
    datasets,closure=prepare_reference_variance()
    data={'conventional':datasets['conventional']}
    manifests={'conventional':closure['modes']['conventional']['dataset']}
    def evaluate(c):
        r,a=evaluate_preflight_objectives(c,data,baseline_manifests=manifests,
            input_kind='confirmed_real',physical_input_confirmation=CONFIRMATION,verify_repeat=False)
        return r['modes']['conventional'],a['conventional']
    initial,initial_arrays=evaluate(compiled)
    vector=compiled.packed_starting_vector
    reapplied=apply_crystallographic_parameters(compiled,vector)
    repeated,repeated_arrays=evaluate(reapplied)
    for name in ('complex_amplitudes','calculated_intensities','objective_intensities'):
        if not np.array_equal(initial_arrays[name],repeated_arrays[name]):
            raise ValueError('Identical parameter-vector reapplication drifted')
    if (initial['ynorm'],initial['objective'])!=(repeated['ynorm'],repeated['objective']):
        raise ValueError('Round-trip objective/YNORM drift')
    if not np.allclose([initial['ynorm'],initial['objective']],
        [14870390.0181119,.939675432747781],rtol=1e-14,atol=0):
        raise ValueError('B6b baseline regression')
    trial=vector.copy();trial[0]+=.0002
    perturbed=apply_crystallographic_parameters(compiled,trial)
    after,after_arrays=evaluate(perturbed)
    before_sample=compiled.crystallography.simulation_sample
    after_sample=perturbed.crystallography.simulation_sample
    grating_error=max(float(np.max(np.abs(a.phase_grating-b.phase_grating)))
        for a,b in zip(before_sample.prepared.slice_types,after_sample.prepared.slice_types))
    if grating_error==0:raise ValueError('B perturbation did not change scattering')
    for site in after_sample.depth_sample.original_structure.sites+after_sample.depth_sample.recast_structure.sites:
        if site.physical_group_id!='Al_fcc' or site.debye_waller_b!=trial[0]:
            raise ValueError('Shared B ownership lost')
    occurrences=sum(len(after_sample.depth_sample.slices[i].sites)
                    for i in after_sample.parameters.slice_sequence_zero_based)
    report=dict(status='PASS',signal_mode='conventional',source_configuration=asdict(config),
        original_structure_snapshot=dict(lattice=compiled.crystallography.original_structure.direct_lattice.vectors.tolist(),
            sites=[asdict(s) for s in compiled.crystallography.original_structure.sites]),
        bounds_rationale=BOUNDS_SOURCE,target_catalogue_audit=TARGET_AUDIT,layout=physical_target_layout(compiled),
        parameter_count=compiled.parameter_count,parameter_vector=vector.tolist(),
        d1_d2=compiled.range_factors,plamda=compiled.plamda,
        provenance=site_occurrence_provenance(compiled),
        baseline=initial,identical_vector_reapplication=repeated,
        round_trip_exact=True,b6b_baseline_match=True,
        checkpoint_identity=canonical_checkpoint_identity(compiled),
        checkpoint_status='Version2 ownership-aware identity; legacy version1 unchanged. Group optimizer/resume explicitly blocked pending B7b integration.',
        perturbation=dict(target_id='site_group.Al_fcc.debye_waller_b',scope='site-group',
            initial_value=float(vector[0]),perturbed_value=float(trial[0]),
            original_sites_affected=4,stored_sites_affected=4,simulation_atom_occurrences_affected=occurrences,
            max_phase_grating_change=grating_error,
            max_amplitude_change=float(np.max(np.abs(initial_arrays['complex_amplitudes']-after_arrays['complex_amplitudes']))),
            result=after),optimizer_executed=False)
    return report,compiled,reapplied,perturbed


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args(argv)
    args.output=args.output.resolve()
    if args.output.exists():raise FileExistsError(args.output)
    report,compiled,_,_=run_diagnostics();args.output.mkdir(parents=True)
    from qcbed.config_serialization import save_yaml_config
    save_yaml_config(compiled.source_configuration,args.output/"configuration.yml",base_directory=ROOT)
    (args.output/'layout_summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    lines=['# B7a physical site-group bring-up','',
        'PASS; conventional mode; no optimizer. One baseline, one identical-vector reapplication, one transient perturbation.',
        '', '## Layout','',
        '| Target | Scope | Initial | Lower | Upper | Units |',
        '|---|---|---:|---:|---:|---|']
    for row in report['layout']:
        lines.append(f"| {row['target_id']} | {row['scope']} | {row['initial_value']} | {row['lower']} | {row['upper']} | {row['units']} |")
    lines += ['',BOUNDS_SOURCE,'',
        'Al_fcc owns original sites0,1,2,3; full occupancy1; shared B0.86 Å²; independent absorption(0,1,1).',
        f"Parameter count: {report['parameter_count']}; vector: {report['parameter_vector']}",
        f"Baseline YNORM: {report['baseline']['ynorm']:.17g}; objective: {report['baseline']['objective']:.17g}",
        'B6b baseline and vector round trip: exact; no physical baseline changes.',
        '', '## Perturbation','', '```json',json.dumps(report['perturbation'],indent=2),'```',
        '', '## Target-family audit','', '| Family | Status | Ownership | Binding / bounds rationale |','|---|---|---|---|']
    lines += ['| '+' | '.join(row)+' |' for row in TARGET_AUDIT]
    lines += ['', '## Snapshot/checkpoint boundary','',report['checkpoint_status'],
        'Exact source configuration and compact occurrence provenance are in layout_summary.json.',
        'YAML site_groups loading/saving is supported; independent same-element groups, linked-C and GUI editing are deferred.',
        'Recommended B7b candidate: shared Al_fcc B only, after controlled-optimizer integration is separately authorized.']
    (args.output/'layout_summary.md').write_text('\n'.join(lines)+'\n')
    print('B7a PASS',report['baseline']['objective'],report['perturbation']['result']['objective'])
    return 0


if __name__=='__main__':raise SystemExit(main())
