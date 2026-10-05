"""One fixed real evaluation per signal mode, plus deterministic repeats. No fit."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from qcbed import load_yaml_config, compile_config_from_sources
from qcbed.objective_diagnostics import evaluate_preflight_objectives
from qcbed.crystallography import fcc_conventional_reflection_diagnostic
if __package__:
    from .a9_real_variance_closure import prepare_reference_variance
    from .b6a_fcc_al_preflight import ROOT, CONFIRMATION
else:
    from a9_real_variance_closure import prepare_reference_variance
    from b6a_fcc_al_preflight import ROOT, CONFIRMATION


def evaluate_real():
    path = ROOT/'examples/fcc_al_real_preflight.yml'
    compiled = compile_config_from_sources(load_yaml_config(path))
    datasets, closure = prepare_reference_variance()
    report, arrays = evaluate_preflight_objectives(compiled,datasets,
        baseline_manifests={m:r['dataset'] for m,r in closure['modes'].items()},
        input_kind='confirmed_real',physical_input_confirmation=CONFIRMATION)
    report['input_config_path'] = str(path)
    report['experimental_provenance'] = closure
    for r in report['preflight']['beam_table']:
        r['fcc_diagnostic'] = fcc_conventional_reflection_diagnostic(tuple(r[f'original_{c}'] for c in 'hkl'))
    return report,arrays


def write_diagnostics(report,artifacts,output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    output = Path(output); output.mkdir(parents=True,exist_ok=False)
    (output/'objective_summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    lines=['# B6b fixed-input objective evaluation','',f"Status: {report['status']}",'',
        'One primary evaluation per mode plus one exact-repeat evaluation each. No optimization.',
        'Exact attached B5 sample; direct gratings, unshifted reciprocal ordering.',
        'Global YNORM only. Objective = sqrt(weighted squared error / weighted experimental norm).',
        'Per-beam squared contributions add to objective squared; allocated root contributions add to objective.',
        '', '## Physical/sample provenance','',
        'Source microscope and crystallographic inputs (including thermal/absorption settings):',
        '```json',json.dumps({k:report['source_configuration'][k] for k in ('microscope','crystallography')},indent=2),'```',
        'Resolved physical sample (including occupancies, positions and multiplier records):',
        '```json',json.dumps(report['resolved_configuration']['sample'],indent=2),'```',
        f"Slice occurrences: {len(report['preflight']['sample']['sequence'])}",
        '```json',json.dumps(report['preflight']['sample'],indent=2),'```',
        '',f"Active parameters: {report['active_parameter_count']}; all starting values fixed.",
        'Complete source/resolved configurations and applied sample parameters are in JSON.',
        'No matched legacy firstout.dat comparison: conventional cell/motif, channels and experimental inputs differ.','']
    rows=report['source_configuration']['scan']['scan_rows']; cols=report['source_configuration']['scan']['scan_columns']
    csv_rows=[]
    for mode,arrays in artifacts.items():
        directory=output/mode; directory.mkdir()
        for name,a in arrays.items():np.save(directory/f'{name}.npy',a)
        stats=report['modes'][mode]
        lines += [f'## {mode}','',f"YNORM: {stats['ynorm']:.17g}; objective: {stats['objective']:.17g}",
            f"Active points: {stats['active_count']}; exact repeat errors: {stats['repeat_max_absolute_errors']}",
            f"Raw intensity statistics: {stats['calculated_intensities']}",f"Amplitude statistics: {stats['amplitudes']}",
            f"Experimental variance/support: {report['preflight']['modes'][mode]}",'',
            '| Beam | Original hkl | HK | Active | Weighted squared error | Fraction | Allocated root error |',
            '|---|---|---|---:|---:|---:|---:|']
        fig,axes=plt.subplots(len(stats['decomposition']['beams']),3,figsize=(10,25),squeeze=False)
        for b,r in enumerate(stats['decomposition']['beams']):
            hkl=tuple(r[f'original_{v}'] for v in 'hkl'); hk=(r['internal_H'],r['internal_K'])
            lines.append(f"| {r['stable_id']} | {hkl} | {hk} | {r['active_count']} | {r['weighted_squared_error']:.17g} | {r['fraction_of_total_error']:.8g} | {r['allocated_root_objective']:.8g} |")
            flat={k:v for k,v in r.items() if not isinstance(v,dict)}; flat['mode']=mode; csv_rows.append(flat)
            beamdir=directory/r['stable_id'];beamdir.mkdir()
            for name in ('experimental','calculated_intensities','objective_intensities','scaled_calculated','residual','weighted_squared_residual','weights','variances'):
                np.save(beamdir/f'{name}.npy',arrays[name][:,b].reshape(rows,cols))
            images=[arrays[n][:,b].reshape(rows,cols) for n in ('experimental','scaled_calculated','residual')]
            low=min(images[0].min(),images[1].min()); high=max(images[0].max(),images[1].max())
            extent=float(np.abs(images[2]).max())
            for j,(im,title) in enumerate(zip(images,('Observed','Calculated × global YNORM','Calculated − observed'))):
                artist=axes[b,j].imshow(im,origin='upper',vmin=low if j<2 else -extent,
                    vmax=high if j<2 else extent,cmap='viridis' if j<2 else 'coolwarm')
                axes[b,j].set_title(f'{hkl}: {title}',fontsize=8);fig.colorbar(artist,ax=axes[b,j],shrink=.7)
        fig.tight_layout();fig.savefig(directory/'comparison.png',dpi=130);plt.close(fig)
    with (output/'beam_residuals.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(csv_rows[0]));writer.writeheader();writer.writerows(csv_rows)
    lines += ['', '## Interpretation and limits','',
        'PASS means finite, consistent integration, not a fitted model or calibrated experimental indexing.',
        'Figures use a shared linear observed/calculated range per beam; residuals use symmetric linear ranges.',
        'No display normalization changes scientific arrays. Support is saved unchanged; unweighted residuals may be nonzero outside support.',
        'Weighted squared residuals are exactly zero there. Signed gradient signals are retained.',
        'B7 needs a separately reviewed refinement selection/bounds and authorization; current layout is empty.']
    (output/'objective_summary.md').write_text('\n'.join(lines)+'\n')


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(argv)
    if args.output.exists():raise FileExistsError(args.output)
    report,arrays=evaluate_real();write_diagnostics(report,arrays,args.output)
    for mode,r in report['modes'].items():print(mode,'YNORM',r['ynorm'],'objective',r['objective'])
    return 0


if __name__=='__main__':raise SystemExit(main())
