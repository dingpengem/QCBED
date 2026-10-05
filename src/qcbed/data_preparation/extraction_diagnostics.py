"""Optional A9 artifact I/O; never used by scientific extraction."""
import csv
import json
from dataclasses import asdict
from pathlib import Path

import numpy as np

from .arc_catalogue import ArcBeamCatalogue
from .lattice_bootstrap import BootstrapSite


def load_a7_extraction_catalogue(path):
    """Read ordered A7 diagnostic records without relabeling or re-detection."""
    with Path(path).open() as stream:
        rows = list(csv.DictReader(stream))
    sites = []
    for row in rows:
        observed = None if not row['observed_x'] else (float(row['observed_x']),float(row['observed_y']))
        sites.append(BootstrapSite((int(row['n1']),int(row['n2'])),
            (float(row['predicted_x']),float(row['predicted_y'])),observed,None,
            row['status'],row['source'],beam_id=row['beam_id']))
    return ArcBeamCatalogue(tuple(sites),'imported A7 diagnostic catalogue',False)


def write_extraction_diagnostics(output, pattern, result, physical_geometry, *, provenance):
    """Write numerical NPY arrays, review figures and summaries into a NEW folder."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse, Rectangle

    output = Path(output)
    output.mkdir(parents=True,exist_ok=False)
    def save(name, a):
        np.save(output / f'{name}.npy', a)
        plt.imsave(output / f'{name}.png', a, cmap='gray',
                   **({'vmin': 0, 'vmax': 1} if a.dtype == np.bool_ else {}))
    save('beam_geometry_mask',result.fields.beam_geometry_mask)
    save('signal_full',result.fields.signal)
    save('support_fraction_full',result.fields.support_fraction)
    if result.fields.variance is not None:
        save('variance_full',result.fields.variance)
    fig,ax=plt.subplots()
    ax.imshow(pattern,cmap='gray')
    for b in result.beams:
        ax.add_patch(Ellipse(b.beam.center_px,2*b.roi.geometry.radius_x_px,
            2*b.roi.geometry.radius_y_px,fill=False,color='cyan' if b.status=='extracted' else 'red'))
        ax.annotate(b.beam.id,b.beam.center_px,color='yellow',fontsize=6)
    ax.set_xlim(min(-.5,*(b.beam.center_px[0]-b.roi.geometry.radius_x_px-3 for b in result.beams)),
                max(pattern.shape[1]-.5,*(b.beam.center_px[0]+b.roi.geometry.radius_x_px+3 for b in result.beams)))
    ax.set_ylim(max(pattern.shape[0]-.5,*(b.beam.center_px[1]+b.roi.geometry.radius_y_px+3 for b in result.beams)),
                min(-.5,*(b.beam.center_px[1]-b.roi.geometry.radius_y_px-3 for b in result.beams)))
    ax.add_patch(Rectangle((-.5,-.5),pattern.shape[1],pattern.shape[0],fill=False,color='white',linestyle='--'))
    ax.set_title(f'A9 {result.config.beam_radius_factor:g}R: cyan retained / red rejected')
    fig.savefig(output/'beam_roi_overview.png',dpi=180); plt.close(fig)
    rows=[]
    for index,b in enumerate(result.beams):
        prefix=f'beam_{index:03d}'  # Never interpret IDs as filesystem paths.
        for field in ('signal','roi_geometry_mask','support_fraction','support_weight','valid_mask','detector_in_frame_mask','final_weight','variance'):
            a=getattr(b,field)
            if a is not None: save(f'{prefix}_{field}',a)
        if b.window_frame_coverage < 1 and b.status == 'extracted':
            fig, axes = plt.subplots(1, 3, figsize=(10,3))
            measured = b.roi_geometry_mask & b.detector_in_frame_mask
            missing = b.roi_geometry_mask & ~b.detector_in_frame_mask
            axes[0].imshow(measured.astype(int)+2*missing, vmin=0, vmax=2, cmap='viridis')
            axes[0].set_title('0 outside / 1 measured / 2 missing')
            axes[1].imshow(b.signal, cmap='gray'); axes[1].set_title('Signal')
            axes[2].imshow(b.final_weight, cmap='gray', vmin=0, vmax=1); axes[2].set_title('Final weight')
            fig.suptitle(b.beam.id); fig.tight_layout()
            fig.savefig(output/f'{prefix}_edge_support.png',dpi=160); plt.close(fig)
        rows.append(dict(beam_id=b.beam.id,n1=b.integer_indices[0],n2=b.integer_indices[1],
            center_x=b.beam.center_px[0],center_y=b.beam.center_px[1],center_source=b.center_source,
            observed_center=b.observed_center_px,catalogue_status=b.catalogue_status,catalogue_source=b.catalogue_source,
            a7_status=b.catalogue_status,full_beam_roi_pixels=b.expected_pixels,
            in_frame_beam_roi_pixels=b.in_frame_pixels,beam_roi_frame_coverage=b.coverage_fraction,
            window_frame_coverage=b.window_frame_coverage,retained_by_a9=b.status=='extracted',
            rejection_reason='' if b.status=='extracted' else 'ellipse frame coverage below configured threshold',
            nonzero_final_weight_pixels=int(np.count_nonzero(b.final_weight)),
            status=b.status,physical_radius_x=physical_geometry.radius_x_px,physical_radius_y=physical_geometry.radius_y_px,
            extraction_radius_x=b.roi.geometry.radius_x_px,extraction_radius_y=b.roi.geometry.radius_y_px,
            roi_bounds=b.bounds_xyxy,roi_expected_pixels=b.expected_pixels,roi_in_frame_pixels=b.in_frame_pixels,
            coverage_fraction=b.coverage_fraction,noise_valid_fraction=b.noise_valid_fraction,
            support_valid_pixels=int(b.valid_mask.sum()),signal_mode=result.fields.signal_mode.value,
            weight_mode=result.fields.weight_mode.value,signal_shape=b.signal.shape,
            signal_sum=float(b.signal.sum()),signal_mean=float(b.signal.mean()),
            variance_mean=None if b.variance is None else float(b.variance.mean()),
            final_nonzero_weight_count=int(np.count_nonzero(b.final_weight)),provenance=b.provenance,
            artifact_prefix=prefix))
    with (output/'beam_summary.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    summary=dict(config=asdict(result.config),provenance=provenance,
        dimensional_flow=dict(physical_selection_shape_px=result.config.physical_selection_shape_px,
            local_extraction_shape_px=(result.config.selection_height,result.config.selection_width),
            diagnostic_window_shape_px=result.config.diagnostic_window_shape_px,
            gradient_step_px=result.config.gradient_step_px,
            final_binning=result.config.final_binning,
            final_pattern_matching_shape_px=result.config.final_pattern_matching_shape_px),
        selected_beam_count=len(rows),extracted_beam_count=sum(b.status=='extracted' for b in result.beams),
        rejected_for_coverage=sum(b.status=='rejected_coverage' for b in result.beams),beams=rows,
        uncertainty=result.fields.uncertainty_provenance or 'Unavailable; beam-mask support only',
        order='canonical catalogue subset; beam then row then column; abstract n1,n2, no HKL',
        limitations='No A10 conversion. No scan corners. No inferred variance. DM fractional Warp unverified.')
    (output/'validation_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary
