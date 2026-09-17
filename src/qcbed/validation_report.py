"""Self-contained HTML validation reports with embedded Matplotlib PNGs."""

import base64
from html import escape
from io import BytesIO
from pathlib import Path
import numpy as np
from .display_mask import experimental_mask_statistics

from .visualization import (plot_cbed_comparison, plot_refinement_convergence,
                            plot_parameter_trajectories, plot_weight_diagnostics, plot_local_progress)


def _png(fig, alt):
    stream=BytesIO()
    fig.savefig(stream,format="png",dpi=100)
    fig.clear()
    encoded=base64.b64encode(stream.getvalue()).decode("ascii")
    return f'<img alt="{escape(alt,quote=True)}" src="data:image/png;base64,{encoded}">'


def generate_refinement_validation_report(path, validation, *, layout, config,
                                          reference_input_name, limitations=(), local_run=None):
    """Write portable HTML; values/maps are never rescaled for presentation."""
    result=validation.refinement
    c=validation.comparison
    initial,final=c.initial_simulated,c.final_simulated
    mode=final.objective_mode
    masks=experimental_mask_statistics(c.experimental)
    if len(layout.parameters)!=len(result.initial_vector):
        raise ValueError("layout must match validation parameter count")
    def table(headers,rows):
        return '<table><tr>'+''.join(f'<th>{escape(str(h))}</th>' for h in headers)+'</tr>'+''.join(
            '<tr>'+''.join(f'<td>{escape(str(v))}</td>' for v in row)+'</tr>' for row in rows)+'</table>'
    parts=['<!doctype html><html lang="en"><meta charset="utf-8"><title>QCBED validation</title>',
           '<style>body{font-family:sans-serif;margin:2em;max-width:1600px}img{max-width:100%}td,th{padding:.4em;border:1px solid #ccc}table{border-collapse:collapse}section{margin-top:2em}</style>',
           '<h1>QCBED source-derived validation</h1>',
           f'<h2>Objective mode: {escape(mode)}</h2>',
           '<p><strong>Non-legacy masked objective: zero-weight points are excluded before variance division.</strong> '
           'No historical Fortran executable equivalence is claimed.</p>' if mode=="masked" else
           '<p>Legacy objective: source-faithful arithmetic, including zero-weight divisions.</p>',
           '<p>Evaluation-only diagnostics; no optimizer stages were run. No refinement or improvement is claimed. '
           'Undefined objective values are retained as NaN, never repaired by masking or altering variance.</p>' if result.loop_count==0 else
           '<p>Python refinement validation; this is not a historical production reproduction.</p>',
           table(("Run field","Value"),[("Reference input",reference_input_name),("Mode",config.method),
               ("Execution backend",config.execution.execution_backend),
               ("Batch size",config.execution.batch_size),("Process workers","not used"),
               ("LOOPS",config.loops),("Seed",config.seed),("Stage sequence",config.stages),
               ("Included points",masks['n_total']-masks['n_zero_weight']),
               ("Masked points",masks['n_zero_weight']),
               ("Zero-variance masked points",masks['n_zero_variance_and_zero_weight']),
               ("Included zero-variance points",masks['n_zero_variance_and_nonzero_weight']),
               ("MAXFSIMP / MAXFPARA",f"{config.maxfsimp} / {config.maxfpara}"),
               ("Active parameters",len(layout.parameters)),("Scan positions",len(initial.tilts)),
               ("Beams",len(final.beam_indices_hk)),("Driver evaluations",result.evaluations),
               ("Initial objective",initial.objective),("Final objective",final.objective),
               ("Initial YNORM",initial.ynorm),("Final YNORM",final.ynorm),("GRADIENT option",final.gradient_option)]),
           '<p>Residual = final YNORM-scaled, optionally GRADIENT-processed simulation − experimental file values. '
           'Raw |amplitude|² is available separately in the comparison data. No display normalization is applied. '
           'Each beam shares one scale across its three data panels; residuals use a symmetric zero-centered scale. '
           'Axes are scan indices with fractional reciprocal H/K step vectors, not assumed Cartesian H/K axes.</p>',
           '<h2>Parameters: physical starting and retained values</h2>',
           table(("Parameter","Initial","Final"),[(f"{i+1}: {p.category} {p.target}",result.initial_vector[i],result.final_vector[i])
                 for i,p in enumerate(layout.parameters)]),
           '<h2>Stage history</h2>',table(("Loop","Optimizer","Input objective","Retained objective","Optimizer evaluations","Seed"),
               [(h.loop_index,h.optimizer,h.input_objective,h.retained_objective,h.optimizer_result.evaluations,h.seed)
                for h in result.stage_history]),
           '<p>History shows stage endpoints only. Validation reevaluations are not optimizer evaluations.</p>',
           _png(plot_refinement_convergence(result,objective_mode=mode),"Stage-level convergence")]
    if local_run is not None and local_run.progress:
        parts.append(_png(plot_local_progress(local_run),"Observed evaluations and stage boundaries"))
    for start in range(0,len(layout.parameters),6):
        parts.append(_png(plot_parameter_trajectories(result,layout,parameter_indices=range(start,min(start+6,len(layout.parameters)))),
                          f"Parameter trajectories {start+1} onward"))
    parts+=['<h2>Per-beam diagnostics</h2>',
            '<p>Numerator Σw·residual²/variance and denominator Σw·experimental²/variance contribute to the existing '
            f'{mode} objective. '+('All points are retained. ' if mode=="legacy" else 'Only nonzero-weight points enter these contributions. ')+
            'Unweighted RMS is a convenience diagnostic only: '
            'the default excludes weight-zero display pixels; all-point and masked-only RMS are also listed. '
            'NaN means undefined (including an empty selected subset). Summing beam contributions may '
            'differ from the original sequential total by floating-point summation roundoff.</p>',
            table(("Beam (H,K)","Numerator","Denominator","RMS: valid display","RMS: all points","RMS: masked only"),
                [(beam,c.numerator_by_beam[i],c.denominator_by_beam[i],c.rms_by_beam[i],
                  c.rms_all_by_beam[i],c.rms_masked_by_beam[i]) for i,beam in enumerate(final.beam_indices_hk)]),
            '<h2>Experimental weights and display masks</h2>',
            '<p>Display masking is separate from legacy objective arithmetic. All four comparison panels use '
            'weight != 0; grey pixels are masked, not zero intensity. Their values do not set color limits. '
            + ('No source-derived objective contribution is skipped. Zero variance causes undefined divisions '
               'even when weight is zero; masking the display does not repair a NaN objective.' if mode=="legacy" else
               'Separately, the explicitly selected masked numerical policy excludes these same points before '
               'variance division in both YNORM and objective; included variances must be positive.')+'</p>',
            table(("Statistic","Value"),[(key,value) for key,value in masks.items()
                  if key not in ("per_beam","visible_zero_variance")]),
            table(("Beam (H,K)","All points","Valid display","Masked","Zero variance","Overlap","Visible zero variance"),
                  [(m['beam_hk'],m['n_total'],m['n_total']-m['n_zero_weight'],m['n_zero_weight'],
                    m['n_zero_variance'],m['n_zero_variance_and_zero_weight'],m['n_zero_variance_and_nonzero_weight'])
                   for m in masks['per_beam']])]
    if masks['n_zero_variance'] and not masks['n_zero_variance_and_nonzero_weight']:
        parts.append('<p>In this dataset, zero variance coincides with zero display weight/mask.</p>')
    if not np.all(np.isfinite(c.numerator_by_beam)) or not np.all(np.isfinite(c.denominator_by_beam)):
        parts.append('<p>Legacy objective contributions are nonfinite; inspect zero variances below. '
                     'No variance repair or zero-weight skipping was applied.</p>')
    if masks['visible_zero_variance']:
        parts.append('<h3>Zero variance at visible points (zero-based indices)</h3>'+table(
            tuple(masks['visible_zero_variance'][0]),[tuple(row.values()) for row in masks['visible_zero_variance']]))
    for i,beam in enumerate(final.beam_indices_hk):
        parts.append(f'<section><h2>Beam {escape(str(beam))}</h2>'+_png(plot_cbed_comparison(c,beam_index=i),f"Beam {beam} comparison")+'</section>')
        parts.append(_png(plot_weight_diagnostics(c.experimental,initial.scan,beam_index=i),
                          f"Beam {beam} weight, zero variance and overlap maps"))
    if local_run is not None:
        from .validation import CbedComparison
        parts.append('<h2>Snapshot gallery</h2><p>Initial and final are mandatory. Intermediate snapshots '
                     'are evaluated trials, not necessarily accepted states. Evaluation cadence includes '
                     'pre-stage and initialization calls. Final diagnostic reevaluation is outside optimizer history.</p>')
        for snap in local_run.snapshots:
            if snap.label=="final":
                continue  # final comparisons are already shown above
            parts.append(f'<h3>{escape(snap.label)}: evaluation {snap.progress.evaluation_index}, '
                         f'{escape(snap.progress.stage_name)}, objective {snap.progress.objective}, '
                         f'YNORM {snap.progress.ynorm}</h3>')
            view=CbedComparison(c.experimental,initial,snap.simulated)
            for i,beam in enumerate(final.beam_indices_hk):
                fig=plot_cbed_comparison(view,beam_index=i)
                fig.axes[2].set_title("Initial simulation" if snap.label=="initial" else "Checkpoint trial")
                fig.axes[3].set_title("Checkpoint − experimental")
                parts.append(_png(fig,f"{snap.label} beam {beam}"))
    notes=("Source-derived format and Python validation do not establish historical executable equivalence.",
           "Historical compiler/RNG/ABI and optimizer limitations remain applicable. No per-evaluation trajectory is invented.",
           "Experimental coordinates use the initial supplied scan; fitted scan coordinates may differ, while acquisition indices remain aligned.")+tuple(limitations)
    parts+=['<h2>Limitations</h2><ul>']+[f'<li>{escape(str(note))}</li>' for note in notes]+['</ul></html>']
    target=Path(path)
    target.write_text('\n'.join(parts),encoding="utf-8")
    return target
