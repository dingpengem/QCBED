"""B6a real-input gate; writes PASS/FAIL diagnostics without evaluating physics.

A confirmation string records the user's source for the real physical assignment;
it is not an indexing search. Development/demo structures are never selected by
default. Missing real inputs produce a nonzero exit and a readable FAIL report.
"""
import argparse
import json
from pathlib import Path
import sys

from qcbed import load_yaml_config,compile_config_from_sources
from qcbed.input_preflight import PreflightError,require_physical_config,preflight_inputs


def write_report(report,output):
    output=Path(output)
    output.mkdir(parents=True,exist_ok=False)
    (output/'preflight.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    lines=['# B6a input preflight', '', f"Status: **{report['status']}**",'',
           'Propagation, YNORM, objective and optimization: not executed.','']
    if report['status']!='PASS':
        lines += [report['code'], '',report['message'],'']
    else:
        lines += ['Input kind: '+report['input_kind'],'',
                  '## Sample','', '```json',json.dumps(report['sample'],indent=2),'```','',
                  '## Maximum errors','', '```json',json.dumps(report['max_errors'],indent=2),'```','',
                  '## Beam indexing','',
                  '| ID | (n1,n2) | Original hkl | Internal HK | Zone law | g error | Beam-parallel g |',
                  '|---|---|---|---|---:|---:|---:|']
        for b in report['beam_table']:
            lines.append(f"| {b['stable_id']} | ({b['n1']},{b['n2']}) | ({b['original_h']},{b['original_k']},{b['original_l']}) | ({b['internal_H']},{b['internal_K']}) | {b['zone_law_value']} | {b['reciprocal_vector_error']:.17g} | {b['beam_parallel_component']:.17g} |")
        for mode,result in report['modes'].items():
            lines += ['',f'## {mode}','',f"Shape: {result['signal_shape']}; active samples: {result['active_count']}",'',
                '| Beam | Active samples | Variance min | Variance max | Variance median |',
                '|---|---:|---:|---:|---:|']
            for b in result['beams']:
                v=b['variance']
                lines.append(f"| {b['stable_id']} | {b['active_samples']} | {v['min']:.17g} | {v['max']:.17g} | {v['median']:.17g} |")
            lines += ['', 'Corners/dimensions: '+repr(result['scan_geometry']),'']
        lines += ['',report['interpretation'],'']
    (output/'preflight.md').write_text('\n'.join(lines)+'\n')


def synthetic_inputs():
    """Explicit demonstration only; no reference pattern or real indexing used."""
    from dataclasses import replace
    import numpy as np
    from qcbed.configuration import ExperimentalBeamConfig,CrystalScanConfig
    from qcbed.data_preparation import (ArcBeamCatalogue,BootstrapSite,DiskGeometry,
        ExtractionConfig,SignalMode,prepare_experimental_qcbed,scan_from_catalogue,
        assemble_prepared_dataset,prepared_dataset_manifest)
    config=load_yaml_config(Path(__file__).with_name('crystallography_config.yml'))
    indices=((0,0),(-1,0),(1,0),(0,-1),(0,1),(-1,-1),(1,1),(-1,1),(1,-1))
    center=np.array([127.,126.]); g1=np.array([81.,9.]);g2=np.array([-10.,81.])
    sites=tuple(BootstrapSite(n,tuple(center+n[0]*g1+n[1]*g2),None,None,
        'reference' if n==(0,0) else 'missing','synthetic',beam_id=f'synthetic-{i}')
        for i,n in enumerate(indices))
    catalogue=ArcBeamCatalogue(sites,'synthetic',False)
    geometry=DiskGeometry(38,38)
    extraction=ExtractionConfig.from_disk_geometry(geometry,final_binning=2)
    conversion=scan_from_catalogue(catalogue,central_beam_id='synthetic-0',
        canonical_basis_beam_ids=('synthetic-2','synthetic-4'),zone_axis_origin_px=tuple(center),
        extraction_config=extraction,basis_transform=((1,0),(0,-1)))
    y,x=np.indices((256,256),dtype=float);image=10+x*x+2*y*y
    datasets={}
    for mode in SignalMode:
        result=prepare_experimental_qcbed(image,catalogue,geometry,
            selected_beam_ids=[s.beam_id for s in sites],config=extraction,signal_mode=mode,
            variance_map=np.full(image.shape,4.),uncertainty_provenance='declared synthetic independent pixel variance')
        datasets[mode.value]=assemble_prepared_dataset(result,conversion,intensity_units='synthetic',physical_geometry=geometry)
    sg=conversion.scan_geometry
    config=replace(config,crystallography=replace(config.crystallography,
        beams=tuple(ExperimentalBeamConfig(s.beam_id,s.integer_indices) for s in sites)),
        scan=CrystalScanConfig(sg.scan_columns,sg.scan_rows,sg.corner_a,sg.corner_b,sg.corner_c))
    return compile_config_from_sources(config),datasets,{m:prepared_dataset_manifest(d) for m,d in datasets.items()}


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--config',type=Path)
    mode.add_argument('--synthetic',action='store_true',help='Explicit synthetic demonstration; never a real PASS')
    parser.add_argument('--physical-input-confirmation',help='Explicit provenance confirming this structure/zone/hkl assignment for the real pattern')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(argv)
    try:
        if args.synthetic:
            compiled,datasets,baselines=synthetic_inputs()
            report=preflight_inputs(compiled,datasets,baseline_manifests=baselines,input_kind='synthetic')
        else:
            config=load_yaml_config(args.config)
            require_physical_config(config)
            if not args.physical_input_confirmation or not args.physical_input_confirmation.strip():
                raise PreflightError('unconfirmed_physical_inputs','Provide confirmed real original structure, zone and g1_hkl/g2_hkl provenance; synthetic examples are not confirmation.')
            compiled=compile_config_from_sources(config)
            # Import only after the physical-input gate; no synthetic replacement.
            if __package__:
                from .a9_real_variance_closure import prepare_reference_variance
            else:
                from a9_real_variance_closure import prepare_reference_variance
            datasets,closure=prepare_reference_variance()
            report=preflight_inputs(compiled,datasets,
                baseline_manifests={m:r['dataset'] for m,r in closure['modes'].items()},
                input_kind='confirmed_real',physical_input_confirmation=args.physical_input_confirmation)
            report['experimental_reference']=dict(pattern='reference_data/pattern.dat',
                sigma='reference_data/sigma.dat',variance_source='legacy_saved_sigma_map; sigma**2',
                geometry_source=closure['geometry_source'])
    except (ValueError,TypeError,OSError) as exc:
        report=dict(status='FAIL',code=getattr(exc,'code','configuration_or_input_error'),message=str(exc),
            propagation_executed=False,objective_executed=False)
    write_report(report,args.output)
    print(f"B6a {report['status']}: {args.output/'preflight.md'}")
    return 0 if report['status']=='PASS' else 1


if __name__=='__main__':sys.exit(main())
