"""Matplotlib diagnostics only; never normalize or mutate scientific data."""

import numpy as np

from .output_io import reshape_scan_map


def plot_cbed_comparison(comparison, *, beam_index=0, mask_zero_weight=True):
    """Four panels, lower origin, shared visible-value scale and grey mask.

    Pixel axes are scan column/row, not incorrectly assumed Cartesian H/K.
    Each panel labels its own increasing-column/row fractional H,K vectors.
    mask_zero_weight=False exposes unmasked maps for debugging. Stored arrays
    and source objective arithmetic are never altered by this display choice.
    """
    from matplotlib.figure import Figure
    c=comparison
    if isinstance(beam_index,bool) or not isinstance(beam_index,(int,np.integer)) or not 0<=beam_index<len(c.experimental.beam_indices_hk):
        raise ValueError("beam_index must identify a measured beam")
    b=beam_index
    initial,final=c.initial_simulated,c.final_simulated
    values=[c.experimental.intensities[:,b],initial.scaled_intensity[:,b],
            final.scaled_intensity[:,b],c.residual[:,b]]
    valid=c.valid_display[:,b] if mask_zero_weight else np.ones(len(values[0]),dtype=bool)
    lower=min(float(np.min(x[valid])) for x in values[:3]) if np.any(valid) else 0.
    upper=max(float(np.max(x[valid])) for x in values[:3]) if np.any(valid) else 1.
    if lower==upper:
        margin=max(abs(lower)*.01,1.)
        lower,upper=lower-margin,upper+margin
    extent=(float(np.max(np.abs(values[3][valid]))) or 1.) if np.any(valid) else 1.
    fig=Figure(figsize=(16,4.8),layout="constrained")
    axes=fig.subplots(1,4)
    for i,(ax,title,v) in enumerate(zip(axes,("Experimental","Initial simulation","Final simulation","Final − experimental"),values)):
        sim=initial if i<2 else final
        from matplotlib import colormaps
        cmap=colormaps["RdBu_r" if i==3 else "viridis"].with_extremes(bad="#bdbdbd")
        mapped=np.ma.array(reshape_scan_map(v,sim.scan),
                           mask=(~valid).reshape(sim.scan.scan_rows,sim.scan.scan_columns))
        image=ax.imshow(mapped,origin="lower",interpolation="nearest",
            aspect="auto",cmap=cmap,
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
                 f"YNORM {initial.ynorm:.6g} → {final.ynorm:.6g}; objective mode={final.objective_mode}; GRADIENT={final.gradient_option}; "
                 + ("grey = weight-zero display mask" if mask_zero_weight else "unmasked debugging view"))
    return fig


def plot_local_progress(run):
    """Observed trial objectives plus stage-retained states, with boundaries."""
    from matplotlib.figure import Figure
    fig=Figure(figsize=(11,4),layout="constrained")
    ax=fig.subplots()
    trials=[r for r in run.progress if r.kind=="evaluation"]
    retained=[r for r in run.progress if r.kind=="retained"]
    ax.plot([r.evaluation_index for r in trials],[r.objective for r in trials],".-",label="Evaluated trials (may be rejected)")
    ax.scatter([r.evaluation_index for r in retained],[r.objective for r in retained],
               marker="*",s=80,label="Stage-retained state")
    seen=set()
    for r in trials:
        if r.stage_index not in seen:
            seen.add(r.stage_index)
            ax.axvline(r.evaluation_index-.5,color="grey",linestyle=":")
            ax.text(r.evaluation_index,1,f"{r.stage_name}, loop {r.loop_index}",
                    transform=ax.get_xaxis_transform(),va="top",fontsize=8)
    ax.set_xlabel("Driver objective evaluation (not optimizer iteration)")
    ax.set_ylabel(f"{run.objective_mode} objective")
    ax.legend()
    return fig


def plot_weight_diagnostics(experimental, scan, *, beam_index=0):
    """Unmasked weight/zero-variance/overlap maps in column-fast acquisition order."""
    from matplotlib.figure import Figure
    b=beam_index
    if isinstance(b,bool) or not isinstance(b,(int,np.integer)) or not 0<=b<len(experimental.beam_indices_hk):
        raise ValueError("beam_index must identify a measured beam")
    if (scan.scan_rows,scan.scan_columns)!=(experimental.scan_rows,experimental.scan_columns):
        raise ValueError("scan dimensions must match experimental data")
    w=experimental.weights[:,b]
    z=experimental.variances[:,b]==0
    fig=Figure(figsize=(12,4),layout="constrained")
    axes=fig.subplots(1,3)
    for i,(ax,v,label) in enumerate(zip(axes,(w,z,(w==0)&z),
                                        ("Weight (unmodified)","Variance == 0","Weight == 0 and variance == 0"))):
        image=ax.imshow(reshape_scan_map(v.astype(np.float64),scan),origin="lower",
                        interpolation="nearest",aspect="auto",cmap="viridis",
                        vmin=min(0.,float(w.min())) if i==0 else 0,
                        vmax=max(1.,float(w.max())) if i==0 else 1)
        ax.set_title(label)
        ax.set_xlabel("Scan column (fast)")
        ax.set_ylabel("Scan row (slow)")
        fig.colorbar(image,ax=ax,shrink=.8)
    fig.suptitle(f"Beam {experimental.beam_indices_hk[b]}: display-mask diagnostics")
    return fig


def plot_refinement_convergence(result, *, objective_mode="legacy"):
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
    ax.set_ylabel(f"{objective_mode} normalized weighted root error")
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
