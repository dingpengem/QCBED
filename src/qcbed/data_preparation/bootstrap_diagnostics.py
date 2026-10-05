"""Read-only bootstrap predictions/evidence, including scientifically incomplete runs."""
import csv
from pathlib import Path
import numpy as np


def write_bootstrap_diagnostics(output,image,result,geometry,*,catalogue=False):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse
    from .diagnostics import crop_bounds, beam_filename
    output=Path(output);folder=output/('beams' if catalogue else 'bootstrap_beams');folder.mkdir()
    colors=dict(reference='white',defining='cyan',confirmed='lime',weak_evidence='orange',missing='red',edge_clipped='gray',observed='lime',weak='orange',rejected='red')
    def draw(ax,site):
        color=colors[site.status]
        ax.scatter(*site.predicted_center_px,c=color,marker='x',label='Predicted')
        ax.add_patch(Ellipse(site.predicted_center_px,width=2*geometry.radius_x_px,height=2*geometry.radius_y_px,
                             fill=False,edgecolor=color,alpha=.6))
        if site.observed_center_px is not None:
            ax.scatter(*site.observed_center_px,edgecolors=color,facecolors='none',marker='o',label='Observed')
        if site.candidate is not None and site.candidate.supporting_pixels:
            points=np.array(site.candidate.supporting_pixels)
            ax.scatter(points[:,0],points[:,1],c=color,s=5,label='Arc evidence')
        ax.annotate(f'{site.integer_indices} {site.status}',site.predicted_center_px,color=color,fontsize=7)
    low,high=np.percentile(image,[1,99.5])
    fig,ax=plt.subplots(figsize=(9,9));ax.imshow(image,cmap='gray',vmin=low,vmax=high,origin='upper')
    for site in result.sites: draw(ax,site)
    ax.set_xlim(-.5,image.shape[1]-.5);ax.set_ylim(image.shape[0]-.5,-.5)
    ax.set_title(('Arc beam catalogue\n' if catalogue else 'Bootstrap predictions: defining ≠ independently confirmed\n')+result.status)
    overview=output/('beam_catalogue_overview.png' if catalogue else 'bootstrap_overview.png');fig.savefig(overview,dpi=150);plt.close(fig)
    records=[]
    for i,site in enumerate(result.sites):
        x0,x1,y0,y1=crop_bounds(image.shape,site.predicted_center_px,geometry)[0]
        fig,ax=plt.subplots(figsize=(6,6))
        ax.imshow(image[y0:y1,x0:x1],cmap='gray',vmin=low,vmax=high,origin='upper',extent=(x0-.5,x1-.5,y1-.5,y0-.5))
        draw(ax,site);ax.set_xlim(x0-.5,x1-.5);ax.set_ylim(y1-.5,y0-.5)
        ax.set_title(f'{site.integer_indices}: {site.status}\n{site.source}\n{site.reason}');ax.legend(fontsize=7)
        path=folder/beam_filename(i,site.integer_indices);fig.savefig(path,dpi=130);plt.close(fig)
        c=site.candidate
        records.append(dict(beam_id=site.beam_id or f'lattice_{site.integer_indices[0]}_{site.integer_indices[1]}',
            initial_indices=site.initial_integer_indices,lattice_residual_px=site.lattice_residual_px,
            boundary_fraction=site.boundary_fraction,usable_boundary_fraction=site.usable_boundary_fraction,
            boundary_limited=site.boundary_limited,
            physical_radius_x=geometry.radius_x_px,physical_radius_y=geometry.radius_y_px,
            footprint_x_min=site.predicted_center_px[0]-geometry.radius_x_px,
            footprint_x_max=site.predicted_center_px[0]+geometry.radius_x_px,
            footprint_y_min=site.predicted_center_px[1]-geometry.radius_y_px,
            footprint_y_max=site.predicted_center_px[1]+geometry.radius_y_px,
            local_hypothesis_shift_px=None if c is None else float(np.linalg.norm(np.array(c.center_px)-site.predicted_center_px)),
            candidate_id=None if c is None else c.id,
            n1=site.integer_indices[0],n2=site.integer_indices[1],
            predicted_x=site.predicted_center_px[0],predicted_y=site.predicted_center_px[1],
            observed_x=None if site.observed_center_px is None else site.observed_center_px[0],
            observed_y=None if site.observed_center_px is None else site.observed_center_px[1],
            shift_px=site.shift_px,status=site.status,source=site.source,reason=site.reason,
            local_hypothesis_x=None if c is None else c.center_px[0],local_hypothesis_y=None if c is None else c.center_px[1],
            support_pixels=None if c is None else c.supporting_edge_pixels,
            angular_coverage_deg=None if c is None else 360*c.angular_coverage,
            radius_rms_px=None if c is None else c.radius_residual_rms_px,
            normal_residual_deg=None if c is None else c.mean_normal_residual_deg,
            confidence=None if c is None else c.confidence,crop=str(path)))
    table=output/('beam_catalogue.csv' if catalogue else 'bootstrap_sites.csv')
    with table.open('x',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
    if catalogue: return dict(beam_overview=str(overview),beam_table=str(table),beam_crops=str(folder))
    return dict(bootstrap_overview=str(overview),bootstrap_table=str(table),bootstrap_crops=str(folder))
