"""Assemble integer-center Patch A products; independently smoke-test synthetic data.

No optimization, noise inference, physical indexing inference or geometry fitting.
Real objective readiness remains blocked until variance and physical mapping exist.
"""
import argparse
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import platform
import sys
from importlib.metadata import version

import numpy as np
from qcbed.data_preparation import (
    ArcBeamCatalogue, BootstrapSite, DiskGeometry, ExtractionConfig, SignalMode,
    WeightMode, prepare_experimental_qcbed, scan_from_catalogue,
    assemble_prepared_dataset, adapt_prepared_dataset, prepared_dataset_manifest)
from qcbed.data_preparation.extraction_diagnostics import (
    load_a7_extraction_catalogue, write_extraction_diagnostics)


LIMITATIONS = [
    'A7 global fit unavailable: only three independent arc-confirmed observations.',
    'A4 legacy noise-mask/gain provenance unresolved.',
    'sigma.dat is a 256x256 uncertainty-related map; sigma versus variance unresolved.',
    'Noise-characterization backend implemented; real legacy validation pending provenance.',
    'Historical DM corners differ by up to 0.012662801232897314; not matched-input equivalence.',
    'Subpixel refinement is experimental/opt-in; canonical closure uses integer centers.',
    'Fractional extraction origins unsupported pending verified DM interpolation/rounding.',
    'Real objective blocked: validated variance and physical HK/scan-basis mapping unavailable.',
]


def synthetic_input(mode):
    """Declared synthetic geometry/variance/HK; never a substitution for real input."""
    from examples.reference_configuration import reference_configuration
    from qcbed.config_compiler import compile_config
    config = reference_configuration()
    indices = ((0,0),(-1,0),(1,0),(0,-1),(0,1),(-1,-1),(1,1),(-1,1),(1,-1))
    sites = tuple(BootstrapSite(n, (48.+24*n[0],48.+24*n[1]), None, None,
        'reference' if n==(0,0) else 'defining', 'declared synthetic geometry',
        beam_id=f'synthetic_{i}') for i,n in enumerate(indices))
    catalogue = ArcBeamCatalogue(sites, 'known synthetic geometry', False)
    geometry = DiskGeometry(8,8)
    cfg = ExtractionConfig.from_disk_geometry(geometry, final_binning=2)
    y,x = np.indices((96,96), dtype=float)
    pattern = 10+x*x+2*y*y+.001*x*y*y
    result = prepare_experimental_qcbed(pattern,catalogue,geometry,
        selected_beam_ids=[s.beam_id for s in sites], config=cfg, signal_mode=mode,
        weight_mode=WeightMode.BEAM_MASK, variance_map=np.full(pattern.shape,4.),
        uncertainty_provenance='Declared synthetic independent pixel variance 4 count^2')
    conversion = scan_from_catalogue(catalogue,central_beam_id='synthetic_0',
        canonical_basis_beam_ids=('synthetic_2','synthetic_4'),zone_axis_origin_px=(48.,48.),
        extraction_config=cfg,basis_transform=((1,0),(0,1)))
    dataset = assemble_prepared_dataset(result,conversion,intensity_units='synthetic counts',
        physical_geometry=geometry,provenance=('Synthetic integration test; not the real pattern.',))
    # Explicit simulation channels for this artificial test, no inference from n1/n2.
    assignments = tuple(zip((b.id for b in dataset.beams),config.experiment.beam_indices_hk))
    adapted = adapt_prepared_dataset(dataset,beam_hk_assignments=assignments,
        basis_1=config.scan.basis_1,basis_2=config.scan.basis_2,objective_mode='masked')
    state = replace(compile_config(config).state,scan=adapted.scan)
    return dataset, adapted, state


def synthetic_smoke(mode):
    from qcbed.refinement import pack_refinement_parameters
    from qcbed.refinement_trial import evaluate_refinement_trial
    from qcbed.execution import ForwardExecutionConfig
    dataset, adapted, state = synthetic_input(mode)
    trial = evaluate_refinement_trial(state,pack_refinement_parameters(state),
        adapted.experimental,objective_mode=adapted.objective_mode,
        execution=ForwardExecutionConfig('batched',64))
    if not np.isfinite(trial.objective) or not np.isfinite(trial.scaling):
        raise ValueError('Synthetic smoke test produced nonfinite objective/scaling')
    return dict(status='passed',objective=float(trial.objective),ynorm=float(trial.scaling),
        evaluations=1,optimization=False,source='declared synthetic input; existing Al reference configuration',
        dataset=prepared_dataset_manifest(dataset),physical_beam_mapping=adapted.beam_mapping)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    root=Path(__file__).resolve().parents[1]
    p.add_argument('--pattern',type=Path,default=root/'reference_data/pattern.dat')
    p.add_argument('--catalogue',type=Path,default=root/'tmp/a7_boundary_catalogue_r38_final/beam_catalogue.csv')
    p.add_argument('--a7-validation',type=Path,default=root/'tmp/a7_boundary_catalogue_r38_final/validation_summary.json')
    p.add_argument('--a8-measurement',type=Path,default=root/'tmp/a8_manual_selection/measurement.json')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    sources={name:dict(path=str(getattr(a,name).resolve()),
        sha256=hashlib.sha256(getattr(a,name).read_bytes()).hexdigest())
        for name in ('pattern','catalogue','a7_validation','a8_measurement')}
    validation=json.loads(a.a7_validation.read_text())
    if validation['sha256']!=sources['pattern']['sha256']:
        raise ValueError('A7 validation belongs to a different pattern')
    catalogue=load_a7_extraction_catalogue(a.catalogue)
    # Verify diagnostic CSV remains the archived A7 scientific catalogue.
    archived=validation['catalogue']['sites']
    actual=[(s.beam_id,s.integer_indices,s.predicted_center_px) for s in catalogue.sites]
    expected=[(s['beam_id'],tuple(s['integer_indices']),tuple(s['predicted_center_px'])) for s in archived]
    if actual!=expected: raise ValueError('Catalogue differs from validated A7 archive')
    measurement=json.loads(a.a8_measurement.read_text())
    image=np.loadtxt(a.pattern)
    geometry=DiskGeometry(38,38)
    cfg=ExtractionConfig.from_disk_geometry(geometry,final_binning=2)
    conversion=scan_from_catalogue(catalogue,central_beam_id='lattice_0_0',
        canonical_basis_beam_ids=('lattice_1_0','lattice_0_1'),
        zone_axis_origin_px=measurement['zone_axis_origin_px'],extraction_config=cfg,
        basis_transform=((1,0),(0,-1)))
    if tuple(measurement['central_disk_center_px'])!=conversion.reference_point_px:
        raise ValueError('A8 and A7 central geometry differ')
    if measurement['geometry']!=asdict(geometry): raise ValueError('A8 physical geometry differs')
    a.output.mkdir(parents=True,exist_ok=False)
    modes={}
    for mode in SignalMode:
        result=prepare_experimental_qcbed(image,catalogue,geometry,
            selected_beam_ids=[s.beam_id for s in catalogue.sites if s.status!='rejected'],
            config=cfg,signal_mode=mode,weight_mode=WeightMode.BEAM_MASK)
        dataset=assemble_prepared_dataset(result,conversion,intensity_units='input pattern units',
            physical_geometry=geometry,provenance=(
                'Manual physical radius 38 px; matched A5 disk-profile image unavailable.',
                'Integer-center canonical geometry; subpixel disabled.',
                'Pattern used as supplied; no A3/A4 provenance inferred.',
                f'Pattern SHA256={sources["pattern"]["sha256"]}'))
        expected_indices=((0,0),(-1,0),(1,0),(0,-1),(0,1),(-1,-1),(1,1),(-1,1),(1,-1))
        if dataset.lattice_coordinates!=expected_indices or dataset.intensities.shape!=(38*38,9):
            raise ValueError('Reference closure geometry/order differs from validated nine-beam selection')
        folder=a.output/mode.value
        write_extraction_diagnostics(folder,image,result,geometry,provenance=sources)
        for name in ('intensities','variances','weights'):
            array=getattr(dataset,name)
            if array is not None: np.save(folder/f'prepared_{name}.npy',array)
        modes[mode.value]=dict(dataset=prepared_dataset_manifest(dataset),
            real_objective=dict(status='blocked',objective=None,ynorm=None,
                reason='Validated scientific variance and physical HK/scan basis unavailable; no placeholders used.'),
            synthetic_smoke=synthetic_smoke(mode))
    report=dict(inputs=sources,software=dict(python=sys.version,platform=platform.platform(),
        packages={n:version(n) for n in ('numpy','scipy','qcbed')}),
        physical_geometry=asdict(geometry),canonical_center_policy='integer/grid; subpixel disabled',
        a5_a6_a7=validation,a8=measurement,a10=asdict(conversion),modes=modes,
        adapter_policy='Explicit ID/HK assignments and physical scan basis; validated variance required; no placeholders',
        limitations=LIMITATIONS)
    (a.output/'prepared_dataset_provenance.json').write_text(json.dumps(report,indent=2)+'\n')
    lines=['# Patch A prepared dataset closure','',
        'Core assembly and synthetic engine integration complete. Real objective is blocked, not evaluated.',
        '', 'Integer centers; radius 38/38 px; physical 76x76; binning 2; final 38x38; nine beams.',
        f'Zone origin: {conversion.zone_axis_origin_px}',f'Scan: {asdict(conversion.scan_geometry)}','']
    for mode,info in modes.items():
        ds=info['dataset'];smoke=info['synthetic_smoke']
        lines += [f'## {mode} + beam_mask','',
            f'Arrays: {ds["arrays"]}; raw pixels/beam={ds["raw_pixels_per_beam"]}; final pixels/beam={ds["final_pixels_per_beam"]}.',
            f'Nonzero weights: {[b["nonzero_weight_pixels"] for b in ds["beam_mapping"]]}',
            f'Real objective: BLOCKED. Synthetic one-evaluation smoke: objective={smoke["objective"]}, YNORM={smoke["ynorm"]}.','']
    lines += ['## Ordered beam mapping','', '| Column | Stable ID | (n1,n2) | Original extraction block |',
              '|---|---|---|---|']
    for entry in modes['conventional']['dataset']['beam_mapping']:
        lines.append(f"| {entry['refinement_column']} | {entry['stable_beam_id']} | {entry['lattice_coordinates']} | {entry['original_extraction_block']} |")
    lines += ['', '## Known limitations','']+['- '+v for v in LIMITATIONS]
    (a.output/'prepared_dataset_summary.md').write_text('\n'.join(lines)+'\n')
    print(a.output/'prepared_dataset_summary.md')


if __name__=='__main__': main()
