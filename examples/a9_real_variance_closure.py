"""Reference sigma replay through unchanged A9/A10; no propagation or objective.

Future production obtains sigma from A4. The archived geometry is independent of
sigma, contains no physical hkl assignment, and is never refitted here.
"""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import numpy as np
from qcbed.data_preparation import (DiskGeometry, ExtractionConfig, SignalMode,
    WeightMode, load_legacy_sigma_map, variance_map_from_sigma_map,
    prepare_experimental_qcbed, scan_from_catalogue, assemble_prepared_dataset,
    prepared_dataset_manifest)
from qcbed.data_preparation.extraction_diagnostics import load_a7_extraction_catalogue

ROOT=Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def variance_statistics(dataset):
    values=dataset.variances
    active=dataset.weights>0
    v=values[active]
    if not v.size or not np.isfinite(v).all() or np.any(v<=0):
        raise ValueError('Active variance must be positive and finite; investigate sigma/support, no floor applied')
    with np.errstate(divide='raise',invalid='raise',over='raise'):
        np.divide(dataset.weights,values,out=np.zeros_like(values),where=active)
    return dict(shape=list(values.shape),finite_count=int(np.isfinite(values).sum()),
        nan_count=int(np.isnan(values).sum()),inf_count=int(np.isinf(values).sum()),
        active_count=int(active.sum()),active_min=float(v.min()),active_max=float(v.max()),
        active_median=float(np.median(v)),active_nonpositive_count=int((v<=0).sum()),
        active_nan_count=int(np.isnan(v).sum()),active_inf_count=int(np.isinf(v).sum()),
        support_per_beam=active.sum(axis=0).tolist())


def prepare_reference_variance(*, pattern_path=ROOT/'reference_data/pattern.dat',
        sigma_path=ROOT/'reference_data/sigma.dat',
        geometry_path=ROOT/'tests/data/real_variance/geometry.json',
        catalogue_path=ROOT/'tests/data/real_variance/beam_catalogue.csv'):
    context=json.loads(Path(geometry_path).read_text())
    if digest(pattern_path)!=context['pattern_sha256'] or digest(sigma_path)!=context['sigma_sha256']:
        raise ValueError('Reference pattern/sigma hash differs from archived validation inputs')
    if digest(catalogue_path)!=context['provenance']['catalogue']['sha256']:
        raise ValueError('Catalogue differs from archived A7 geometry')
    image=np.loadtxt(pattern_path)
    sigma=load_legacy_sigma_map(sigma_path,expected_shape=image.shape)
    variance=variance_map_from_sigma_map(sigma)
    geometry=DiskGeometry(**context['geometry'])
    cfg=ExtractionConfig.from_disk_geometry(geometry,final_binning=2)
    catalogue=load_a7_extraction_catalogue(catalogue_path)
    conversion=scan_from_catalogue(catalogue,central_beam_id='lattice_0_0',
        canonical_basis_beam_ids=('lattice_1_0','lattice_0_1'),
        zone_axis_origin_px=context['zone_axis_origin_px'],extraction_config=cfg,
        basis_transform=((1,0),(0,-1)))
    if conversion.reference_point_px!=tuple(context['central_disk_center_px']):
        raise ValueError('A7/A8 central geometry mismatch')
    if json.loads(json.dumps(asdict(conversion.scan_geometry)))!=context['a10']['scan_geometry']:
        raise ValueError('A10 corners/dimensions differ from archived closure')
    provenance=('legacy_saved_sigma_map; variance=sigma**2; reference compatibility only; '
        f'source={Path(sigma_path).resolve()}; sha256={digest(sigma_path)}; '
        'exact historical save statement not recovered; source semantics and numerical validation support interpretation')
    datasets={};reports={}
    for mode in SignalMode:
        options=dict(selected_beam_ids=[s.beam_id for s in catalogue.sites if s.status!='rejected'],
            config=cfg,signal_mode=mode,weight_mode=WeightMode.BEAM_MASK)
        baseline=prepare_experimental_qcbed(image,catalogue,geometry,**options)
        extracted=prepare_experimental_qcbed(image,catalogue,geometry,
            variance_map=variance,uncertainty_provenance=provenance,**options)
        for old,new in zip(baseline.beams,extracted.beams,strict=True):
            for name in ('signal','support_fraction','support_weight','final_weight'):
                if not np.array_equal(getattr(old,name),getattr(new,name)):
                    raise ValueError(f'{mode.value}: adding variance changed {name}')
        dataset=assemble_prepared_dataset(extracted,conversion,intensity_units='input pattern units',
            physical_geometry=geometry,provenance=(provenance,'beam_mask support; no fabricated A4 noise-valid mask'))
        manifest=prepared_dataset_manifest(dataset)
        previous=context['baseline_modes'][mode.value]['dataset']
        for name in ('intensities','weights'):
            if manifest['arrays'][name]!=previous['arrays'][name]:
                raise ValueError(f'{mode.value}: {name} differs from archived variance=None closure')
        datasets[mode.value]=dataset
        reports[mode.value]=dict(statistics=variance_statistics(dataset),dataset=manifest,
            signal_support_regression='exact equality against variance=None and archived hashes',
            a10_regression='exact equality with archived corners/dimensions')
    summary=dict(source='legacy_saved_sigma_map',sigma_shape=list(sigma.shape),
        sigma_values=int(sigma.size),sigma_format='whitespace-delimited text, rows/columns',
        sigma_min=float(sigma.min()),sigma_max=float(sigma.max()),
        variance_formula='sigma**2',provenance=provenance,
        geometry=asdict(geometry),extraction_config=asdict(cfg),modes=reports,
        refinement_adapter='not run for real data: physical hkl/basis still unconfirmed; synthetic adapter tested separately',
        objective='not evaluated',a4_refit=False,
        canonical_production='corrected CBED -> A4 -> sigma map -> variance map -> A9',
        geometry_source=dict(path=str(Path(geometry_path).resolve()),sha256=digest(geometry_path)))
    return datasets,summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    datasets,summary=prepare_reference_variance()
    args.output.mkdir(parents=True,exist_ok=False)
    for mode,data in datasets.items():
        folder=args.output/mode;folder.mkdir()
        for name in ('intensities','variances','weights'):
            np.save(folder/f'{name}.npy',getattr(data,name))
    (args.output/'variance_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    lines=['# A9 real variance closure','',summary['provenance'],'',
        'Saved sigma replay only. No A4 refit, propagation, objective or optimization.','']
    for mode,report in summary['modes'].items():
        lines += [f'## {mode}', '',json.dumps(report['statistics'],indent=2),'',
            'Beam order: '+repr(datasets[mode].lattice_coordinates),'',
            'Scan corners: '+repr(asdict(datasets[mode].scan_geometry)),'']
    (args.output/'variance_summary.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({k:v['statistics'] for k,v in summary['modes'].items()},indent=2))


if __name__=='__main__': main()
