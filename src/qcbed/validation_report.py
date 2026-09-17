"""Self-contained HTML validation reports with embedded Matplotlib PNGs."""

import base64
from html import escape
from io import BytesIO
from pathlib import Path

from .visualization import (plot_cbed_comparison, plot_refinement_convergence,
                            plot_parameter_trajectories)


def _png(fig, alt):
    stream=BytesIO()
    fig.savefig(stream,format="png",dpi=100)
    fig.clear()
    encoded=base64.b64encode(stream.getvalue()).decode("ascii")
    return f'<img alt="{escape(alt,quote=True)}" src="data:image/png;base64,{encoded}">'


def generate_refinement_validation_report(path, validation, *, layout, config,
                                          reference_input_name, limitations=()):
    """Write portable HTML; values/maps are never rescaled for presentation."""
    result=validation.refinement
    c=validation.comparison
    initial,final=c.initial_simulated,c.final_simulated
    if len(layout.parameters)!=len(result.initial_vector):
        raise ValueError("layout must match validation parameter count")
    def table(headers,rows):
        return '<table><tr>'+''.join(f'<th>{escape(str(h))}</th>' for h in headers)+'</tr>'+''.join(
            '<tr>'+''.join(f'<td>{escape(str(v))}</td>' for v in row)+'</tr>' for row in rows)+'</table>'
    parts=['<!doctype html><html lang="en"><meta charset="utf-8"><title>QCBED validation</title>',
           '<style>body{font-family:sans-serif;margin:2em;max-width:1600px}img{max-width:100%}td,th{padding:.4em;border:1px solid #ccc}table{border-collapse:collapse}section{margin-top:2em}</style>',
           '<h1>QCBED source-derived validation</h1>',
           '<p>Evaluation-only diagnostics; optimizer driver was not run. No refinement or improvement is claimed. '
           'Undefined objective values are retained as NaN, never repaired by masking or altering variance.</p>' if result.loop_count==0 else
           '<p>Python refinement validation; this is not a historical production reproduction.</p>',
           table(("Run field","Value"),[("Reference input",reference_input_name),("Mode",config.method),
               ("LOOPS",config.loops),("Seed",config.seed),("MAXFSIMP / MAXFPARA",f"{config.maxfsimp} / {config.maxfpara}"),
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
           _png(plot_refinement_convergence(result),"Stage-level convergence")]
    for start in range(0,len(layout.parameters),6):
        parts.append(_png(plot_parameter_trajectories(result,layout,parameter_indices=range(start,min(start+6,len(layout.parameters)))),
                          f"Parameter trajectories {start+1} onward"))
    parts+=['<h2>Per-beam diagnostics</h2>',
            '<p>Numerator Σw·residual²/variance and denominator Σw·experimental²/variance contribute to the existing '
            'legacy objective. Unweighted RMS is a convenience diagnostic only. Summing beam contributions may '
            'differ from the original sequential total by floating-point summation roundoff.</p>',
            table(("Beam (H,K)","Numerator","Denominator","Unweighted RMS"),
                [(beam,c.numerator_by_beam[i],c.denominator_by_beam[i],c.rms_by_beam[i]) for i,beam in enumerate(final.beam_indices_hk)])]
    for i,beam in enumerate(final.beam_indices_hk):
        parts.append(f'<section><h2>Beam {escape(str(beam))}</h2>'+_png(plot_cbed_comparison(c,beam_index=i),f"Beam {beam} comparison")+'</section>')
    notes=("Source-derived format and Python validation do not establish historical executable equivalence.",
           "Historical compiler/RNG/ABI and optimizer limitations remain applicable. No per-evaluation trajectory is invented.",
           "Experimental coordinates use the initial supplied scan; fitted scan coordinates may differ, while acquisition indices remain aligned.")+tuple(limitations)
    parts+=['<h2>Limitations</h2><ul>']+[f'<li>{escape(str(note))}</li>' for note in notes]+['</ul></html>']
    target=Path(path)
    target.write_text('\n'.join(parts),encoding="utf-8")
    return target
