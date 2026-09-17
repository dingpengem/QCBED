"""Matplotlib diagnostics only; never normalize or mutate scientific data."""

import numpy as np

from .output_io import reshape_scan_map


def plot_cbed_comparison(comparison, *, beam_index=0):
    """One readable four-panel figure per beam, origin lower, common data scale.

    Pixel axes are scan column/row, not incorrectly assumed Cartesian H/K.
    Each panel labels its own increasing-column/row fractional H,K vectors.
    """
    from matplotlib.figure import Figure
    c=comparison
    if isinstance(beam_index,bool) or not isinstance(beam_index,(int,np.integer)) or not 0<=beam_index<len(c.experimental.beam_indices_hk):
        raise ValueError("beam_index must identify a measured beam")
    b=beam_index
    initial,final=c.initial_simulated,c.final_simulated
    values=[c.experimental.intensities[:,b],initial.scaled_intensity[:,b],
            final.scaled_intensity[:,b],c.residual[:,b]]
    lower=min(float(np.min(x)) for x in values[:3])
    upper=max(float(np.max(x)) for x in values[:3])
    if lower==upper:
        margin=max(abs(lower)*.01,1.)
        lower,upper=lower-margin,upper+margin
    extent=float(np.max(np.abs(values[3]))) or 1.
    fig=Figure(figsize=(16,4.8),layout="constrained")
    axes=fig.subplots(1,4)
    for i,(ax,title,v) in enumerate(zip(axes,("Experimental","Initial simulation","Final simulation","Final − experimental"),values)):
        sim=initial if i<2 else final
        image=ax.imshow(reshape_scan_map(v,sim.scan),origin="lower",interpolation="nearest",
            aspect="auto",cmap="RdBu_r" if i==3 else "viridis",
            vmin=-extent if i==3 else lower,vmax=extent if i==3 else upper)
        tilts=sim.tilts
        dc=tilts[1]-tilts[0]
        dr=tilts[sim.scan.scan_columns]-tilts[0]
        ax.set_xlabel(f"Scan column (fast)\nΔ(H,K)=({dc[0]:.3g},{dc[1]:.3g})")
        ax.set_ylabel(f"Scan row (slow)\nΔ(H,K)=({dr[0]:.3g},{dr[1]:.3g})")
        ax.set_title(title)
        if i==2:
            fig.colorbar(image,ax=list(axes[:3]),label="Objective-space value (no display normalization)",shrink=.8)
        if i==3:
            fig.colorbar(image,ax=ax,label="Signed residual",shrink=.8)
    fig.suptitle(f"Beam {c.experimental.beam_indices_hk[b]} | objective {initial.objective:.6g} → {final.objective:.6g}\n"
                 f"YNORM {initial.ynorm:.6g} → {final.ynorm:.6g}; GRADIENT={final.gradient_option}")
    return fig


def plot_refinement_convergence(result):
    """Stage-level retained scores; no fabricated per-evaluation history."""
    from matplotlib.figure import Figure
    fig=Figure(figsize=(10,4),layout="constrained")
    ax=fig.subplots()
    values=[result.initial_objective]+[h.retained_objective for h in result.stage_history]
    ax.plot(range(len(values)),values,"o-",label="Retained objective")
    if not np.all(np.isfinite(values)):
        ax.text(.5,.5,"Undefined objective (NaN/Inf); no finite convergence value",
                transform=ax.transAxes,ha="center",wrap=True)
    for i,h in enumerate(result.stage_history,1):
        ax.scatter(i-.15,h.input_objective,marker="x",color="black",label="Stage input" if i==1 else None)
        ax.axvline(i-.5,color="grey",linestyle=":",linewidth=.6)
    labels=["Initial"]+[f"{h.optimizer}\nloop {h.loop_index}" for h in result.stage_history]
    ax.set_xticks(range(len(values)),labels,rotation=30,ha="right")
    ax.set_ylabel("Legacy normalized weighted root error")
    ax.set_xlabel("Completed optimizer stage (not individual evaluations)")
    ax.legend()
    return fig


def plot_parameter_trajectories(result, layout, *, parameter_indices=None):
    """One labeled axis per selected parameter; physical stage endpoints only."""
    from matplotlib.figure import Figure
    indices=tuple(range(len(layout.parameters))) if parameter_indices is None else tuple(parameter_indices)
    if not indices or any(isinstance(i,bool) or not isinstance(i,(int,np.integer)) or not 0<=i<len(layout.parameters) for i in indices):
        raise ValueError("parameter_indices must select at least one active parameter")
    fig=Figure(figsize=(10,2.3*len(indices)),layout="constrained")
    axes=fig.subplots(len(indices),1,squeeze=False)[:,0]
    points=[result.initial_vector]+[h.retained_vector for h in result.stage_history]
    for ax,i in zip(axes,indices):
        spec=layout.parameters[i]
        ax.plot(range(len(points)),[p[i] for p in points],"o-")
        ax.set_ylabel(f"{i+1}: {spec.category}\n{spec.target}")
        ax.set_xlabel("Completed stage; physical parameter units")
        ax.set_xticks(range(len(points)))
    return fig
