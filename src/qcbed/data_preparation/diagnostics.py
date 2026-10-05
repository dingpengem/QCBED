"""Read-only A5–A7 validation outputs; plotting is imported only when requested."""
import argparse
import csv
from dataclasses import asdict
import hashlib
import json
from pathlib import Path

import numpy as np

from .models import DiskGeometry
from .disk_detection import estimate_disk_geometry
from .central_arc import locate_central_disk, write_central_arc_localization, CentralDiskResult, CentralArcConfig
from .lattice_detection import (LatticeDetectionConfig, detect_lattice_candidates,
                                infer_lattice_seed, arc_center_votes)
from .lattice_refinement import LatticeRefinementConfig, refine_lattice_and_generate_catalogue, generate_beam_catalogue
from .lattice_bootstrap import bootstrap_lattice
from .bootstrap_diagnostics import write_bootstrap_diagnostics


class _BootstrapIncomplete(Exception):
    """Predictions were saved, but no refined lattice may be claimed."""



def load_pattern(path):
    """Headerless rectangular whitespace text: physical rows -> y, columns -> x."""
    image = np.loadtxt(path, dtype=np.float64, ndmin=2)
    if image.ndim != 2 or not image.size or not np.isfinite(image).all():
        raise ValueError('Expected finite nonempty rectangular numeric text')
    return image


def beam_filename(index, indices):
    def signed(n):
        return f'm{abs(n)}' if n < 0 else str(n)
    return f'beam_{index:03d}_n{signed(indices[0])}_{signed(indices[1])}.png'


def crop_bounds(shape, center, geometry, scale=1.25):
    """Inclusive pixel-center bounds rounded outward; half-open clipped slices."""
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError('crop scale must be positive and finite')
    x,y=center
    bounds=(int(np.floor(x-scale*geometry.radius_x_px)),
            int(np.ceil(x+scale*geometry.radius_x_px))+1,
            int(np.floor(y-scale*geometry.radius_y_px)),
            int(np.ceil(y+scale*geometry.radius_y_px))+1)
    x0,x1,y0,y1=bounds
    clipped=(max(0,x0),min(shape[1],x1),max(0,y0),min(shape[0],y1))
    return clipped, clipped != bounds


def ellipse_overlay(center, geometry):
    """Matplotlib ellipse center/width/height in full-frame x/y pixel coordinates."""
    return dict(xy=tuple(center), width=2*geometry.radius_x_px, height=2*geometry.radius_y_px)


def beam_rows(catalogue):
    assignments={a.integer_indices:a.residual_px for a in catalogue.lattice.assignments if a.inlier}
    rows=[]
    for b in catalogue.beams:
        observed=b.observed_center_px
        shift=None if observed is None else float(np.linalg.norm(np.array(observed)-b.predicted_center_px))
        rows.append(dict(beam_id=b.observation.id,n1=b.lattice_indices[0],n2=b.lattice_indices[1],
            predicted_x_px=b.predicted_center_px[0],predicted_y_px=b.predicted_center_px[1],
            observed_x_px=None if observed is None else observed[0],
            observed_y_px=None if observed is None else observed[1],
            predicted_to_observed_shift_px=shift,status=b.status.value,response_score=b.response_score,
            lattice_residual_px=assignments.get(b.lattice_indices)))
    return rows


def write_beam_table(path,catalogue):
    rows=beam_rows(catalogue)
    if not rows:
        raise ValueError('Cannot serialize an empty catalogue')
    with Path(path).open('x',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    return rows


def write_figures(output,image,a5,seed,catalogue,*,crop_scale=1.25):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse
    output=Path(output); crops=output/'beams'; crops.mkdir()
    fig,ax=plt.subplots(figsize=(10,9))
    ax.imshow(image,origin='upper',cmap='gray')
    ax.scatter(*a5.central_disk_center_px,marker='*',s=130,c='red',label='A5 central disk')
    points=np.array([a.candidate.center_px for a in seed.assignments])
    if len(points):
        ax.scatter(points[:,0],points[:,1],marker='+',c='yellow',label='A6 candidates')
    predicted=np.array([b.predicted_center_px for b in catalogue.beams])
    ax.scatter(predicted[:,0],predicted[:,1],marker='x',c='cyan',label='A7 predictions')
    observed=np.array([b.observed_center_px for b in catalogue.beams if b.observed_center_px is not None])
    if len(observed):
        ax.scatter(observed[:,0],observed[:,1],marker='o',facecolors='none',edgecolors='lime',label='Reliable local observations')
    if len(catalogue.beams)<=150:
        for b in catalogue.beams:
            ax.annotate(str(b.lattice_indices),b.predicted_center_px,fontsize=6,color='cyan')
    ax.legend(); ax.set_title('Geometry diagnostics; no detected-edge claim')
    fig.savefig(output/'beam_catalogue_overview.png',dpi=150); plt.close(fig)
    records=[]
    for i,b in enumerate(catalogue.beams):
        bounds,clipped=crop_bounds(image.shape,b.predicted_center_px,catalogue.geometry,crop_scale)
        x0,x1,y0,y1=bounds
        filename=beam_filename(i,b.lattice_indices)
        record=dict(beam_id=b.observation.id,filename=f'beams/{filename}',bounds=bounds,clipped=clipped)
        if x1<=x0 or y1<=y0:
            record.update(filename=None,reason='No image intersection'); records.append(record); continue
        fig,ax=plt.subplots(figsize=(5,5))
        ax.imshow(image[y0:y1,x0:x1],origin='upper',cmap='gray',extent=(x0-.5,x1-.5,y1-.5,y0-.5))
        ax.add_patch(Ellipse(**ellipse_overlay(b.predicted_center_px,catalogue.geometry),
                             fill=False,edgecolor='cyan',label='Expected footprint (not detected edge)'))
        ax.scatter(*b.predicted_center_px,marker='x',c='cyan',label='Predicted')
        if b.observed_center_px is not None:
            ax.scatter(*b.observed_center_px,marker='+',c='lime',label='Observed')
            ax.plot([b.predicted_center_px[0],b.observed_center_px[0]],
                    [b.predicted_center_px[1],b.observed_center_px[1]],color='lime')
        ax.set_title(f'{b.lattice_indices} | {b.status.value} | response={b.response_score}\nclipped={clipped}')
        ax.legend(fontsize=6); fig.savefig(crops/filename,dpi=120); plt.close(fig)
        records.append(record)
    return records


def write_arc_validation(output,image,geometry,candidates,arrays):
    """Review evidence, not a ground-truth internal/boundary classification."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse
    output=Path(output)
    retained=np.zeros(image.shape,dtype=bool)
    for c in candidates:
        if c.supporting_pixels:
            points=np.array(c.supporting_pixels)
            retained[points[:,1],points[:,0]]=True
    initial=arrays['strong_gradient_pixels']
    rejected=initial & ~retained
    fig,axes=plt.subplots(1,3,figsize=(16,5))
    axes[0].imshow(initial,cmap='gray',origin='upper')
    axes[0].set_title('All initial strong-gradient pixels')
    rgb=np.zeros((*image.shape,3))
    rgb[rejected]=(0.7,0.15,0.15);rgb[retained]=(0.1,1,0.2)
    axes[1].imshow(rgb,origin='upper')
    axes[1].set_title('Green: validated support; red: not retained\nRed does not prove an internal fringe')
    low,high=np.percentile(image,[1,99.5])
    axes[2].imshow(image,cmap='gray',origin='upper',vmin=low,vmax=high)
    for i,c in enumerate(candidates):
        color=plt.get_cmap('tab20')(i%20)
        points=np.array(c.supporting_pixels)
        if len(points): axes[2].scatter(points[:,0],points[:,1],s=4,color=color)
        axes[2].scatter(*c.center_px,color=color,marker='+')
        axes[2].add_patch(Ellipse(xy=c.center_px,width=2*(geometry.radius_x_px+c.radius_offset_px),
            height=2*(geometry.radius_y_px+c.radius_offset_px),fill=False,color=color,alpha=.5))
        axes[2].annotate(f'{c.id}: radius RMS {c.radius_residual_rms_px:.2f}px',c.center_px,fontsize=6,color=color)
    axes[2].set_xlim(-.5,image.shape[1]-.5);axes[2].set_ylim(image.shape[0]-.5,-.5)
    axes[2].set_title(f'Validated fragments by center: {len(candidates)}\nEllipses are radius-prior hypotheses')
    overview=output/'arc_validation_overview.png'
    fig.tight_layout();fig.savefig(overview,dpi=150);plt.close(fig)
    table=output/'arc_validation.csv'
    with table.open('x',newline='') as stream:
        writer=csv.writer(stream)
        writer.writerow(['candidate_id','center_x','center_y','support_pixel_count','arc_length_px',
            'angular_coverage_deg','fragment_count','center_vote_rms_px','mean_radius_px','radius_std_px',
            'mean_normal_residual_deg','radius_residual_rms_px','radius_offset_px','confidence'])
        for c in candidates:
            writer.writerow([c.id,*c.center_px,c.supporting_edge_pixels,c.arc_length_px,360*c.angular_coverage,c.fragment_count,
                c.center_vote_rms_px,c.mean_radius_px,c.radius_std_px,c.mean_normal_residual_deg,
                c.radius_residual_rms_px,c.radius_offset_px,c.confidence])
    # Preserve exact membership for audit without reconstructing scientific classifications from PNGs.
    membership=output/'arc_support_pixels.csv'
    with membership.open('x',newline='') as stream:
        writer=csv.writer(stream);writer.writerow(['candidate_id','x_px','y_px'])
        for c in candidates:
            writer.writerows((c.id,x,y) for x,y in c.supporting_pixels)
    return dict(arc_overview=str(overview),arc_table=str(table),arc_support_pixels=str(membership),
                retained_support_pixels=int(retained.sum()),not_retained_gradient_pixels=int(rejected.sum()))


def write_candidate_diagnostics(output,image,geometry,center,candidates,config):
    """Visualize A6's geometric search support, including a failed seed attempt.

    The mask is diagnostic only; selection remains owned by A6. No fitted beam
    catalogue is inferred from this display or from an empty candidate list.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse
    output=Path(output)
    votes,arrays,support = arc_center_votes(image,geometry,center,config=config)
    support['candidates_after_suppression']=len(candidates)
    support.update(write_arc_validation(output,image,geometry,candidates,arrays))
    fig,axes=plt.subplots(2,2,figsize=(12,10))
    vmin,vmax=np.percentile(image,[1,99.5])
    axes[0,0].imshow(image,origin='upper',cmap='gray',vmin=vmin,vmax=vmax)
    axes[0,0].scatter(*center,c='red',marker='*',label='A5 center')
    for c in candidates:
        axes[0,0].scatter(*c.center_px,c='lime',marker='+')
        axes[0,0].add_patch(Ellipse(**ellipse_overlay(c.center_px,geometry),fill=False,edgecolor='cyan',alpha=.5))
    axes[0,0].set_title('Preliminary arc centers; physical ellipse prior')
    axes[0,1].imshow(arrays['gradient_magnitude'],origin='upper',cmap='magma')
    axes[0,1].set_title('Gradient magnitude (edge-only smoothing)')
    axes[1,0].imshow(arrays['selected_edges'],origin='upper',cmap='gray')
    axes[1,0].set_title('Selected edge pixels; central footprint excluded')
    axes[1,1].imshow(votes,origin='upper',cmap='viridis')
    axes[1,1].set_title('Bilinear normal-direction center votes')
    path=output/'candidate_detection_overview.png'
    fig.savefig(path,dpi=150); plt.close(fig)
    table=output/'a6_candidates.csv'
    with table.open('x',newline='') as stream:
        writer=csv.writer(stream); writer.writerow(['candidate_id','x_px','y_px','vote_strength',
            'supporting_edge_pixels','angular_coverage','radius_residual_rms_px','distance_from_central_px'])
        for c in candidates:
            writer.writerow([c.id,*c.center_px,c.vote_strength,c.supporting_edge_pixels,c.angular_coverage,
                c.radius_residual_rms_px,float(np.linalg.norm(np.array(c.center_px)-center))])
    support['pairwise_distances_px']=[float(np.linalg.norm(np.array(a.center_px)-b.center_px))
        for i,a in enumerate(candidates) for b in candidates[i+1:]]
    return dict(**support,candidate_count=len(candidates),overview=str(path),candidate_table=str(table),
                display='1–99.5 percentile contrast only; numerical intensity unchanged')


def write_seed_diagnostics(output, image, seed):
    """Retain A6 classification even if subsequent A7 cannot complete."""
    import matplotlib.pyplot as plt
    output=Path(output)
    fig,ax=plt.subplots()
    ax.imshow(image,origin='upper',cmap='gray')
    for flag,color,label in ((True,'lime','A6 inlier'),(False,'red','A6 outlier')):
        points=np.array([a.candidate.center_px for a in seed.assignments if a.inlier==flag])
        if len(points): ax.scatter(points[:,0],points[:,1],c=color,marker='+',label=label)
    ax.legend(); ax.set_title('A6 lattice consistency; preliminary observations')
    path=output/'a6_lattice_assignments.png'
    fig.savefig(path,dpi=150); plt.close(fig)
    return str(path)


def _json_default(value):
    if isinstance(value,np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def validate_pattern(path,output,*,expected_shape=None,profile=None,radii=None,geometry=None,
                     geometry_provenance=None,crop_scale=1.25,detection_config=None):
    """First/default geometry run; no A2–A4, tuning or fallback. New output directory only."""
    detection_config=detection_config if detection_config is not None else LatticeDetectionConfig()
    if not isinstance(detection_config,LatticeDetectionConfig):
        raise TypeError('detection_config must be LatticeDetectionConfig')
    output=Path(output).resolve(); output.mkdir(parents=True,exist_ok=False)
    summary=dict(status='blocked',input_file=str(Path(path).resolve()),
        preprocessing='Already legacy PSF-corrected/noise-prepared, per user provenance; A2–A4 not run',
        expected_shape=None if expected_shape is None else list(expected_shape),configuration=dict(a5_search=asdict(CentralArcConfig()),
            a6=asdict(detection_config),a7=asdict(LatticeRefinementConfig()),crop_scale=crop_scale),
        blockers=[],a5=None,a6=None,a7=None,outputs={},
        stability_assessment='Not established: no completed lattice fit')
    stage='input_validation'
    try:
        image=load_pattern(path)
        summary.update(input_shape=list(image.shape),value_count=image.size,dtype=str(image.dtype),
                       sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),
                       layout='Whitespace text, physical rows=y, columns=x; no reshape/resampling')
        if expected_shape is not None and image.shape != tuple(expected_shape):
            summary['blockers'].append(f'Expected {tuple(expected_shape)}, found {image.shape}; no upsampling/reshape performed')
        if sum(v is not None for v in (profile,radii,geometry)) > 1:
            raise ValueError('Supply only one of profile, radii or DiskGeometry')
        if geometry is not None and not isinstance(geometry,DiskGeometry):
            raise TypeError('geometry must be DiskGeometry')
        if profile is None and radii is None and geometry is None:
            summary['blockers'].append('Missing matched A5 disk-profile input or validated radii in pattern pixel units')
        elif not geometry_provenance:
            summary['blockers'].append('Explicit geometry provenance is required')
        summary['geometry_input']=dict(profile=None if profile is None else str(Path(profile).resolve()),
                                      radii=(geometry.radius_x_px,geometry.radius_y_px) if geometry is not None else radii,
                                      source='profile_estimation' if profile is not None else 'manually_supplied',
                                      provenance=geometry_provenance)
        if not summary['blockers']:
            geometry=geometry if geometry is not None else (estimate_disk_geometry(load_pattern(profile)).geometry if profile else DiskGeometry(*radii))
            summary['a5']=dict(geometry=asdict(geometry),geometry_estimated_by_a5=profile is not None,
                               central_center_method='known-radius partial disk-boundary geometry',
                               central_disk_center_px=None)
            stage='a5_central_localization'
            a5=locate_central_disk(image,geometry)
            summary['a5']['central_disk_center_px']=a5.central_disk_center_px
            if isinstance(a5,CentralDiskResult):
                summary['a5']['localization']=asdict(a5)
                summary['outputs']['central_localization']=write_central_arc_localization(output,image,a5)
            stage='a6_candidate_detection'
            candidates=detect_lattice_candidates(image,geometry,a5.central_disk_center_px,config=detection_config)
            summary['candidate_detection']=write_candidate_diagnostics(output,image,geometry,
                a5.central_disk_center_px,candidates,detection_config)
            summary['outputs'].update(candidate_overview=summary['candidate_detection']['overview'],
                                      candidate_table=summary['candidate_detection']['candidate_table'],
                                      arc_overview=summary['candidate_detection']['arc_overview'],
                                      arc_table=summary['candidate_detection']['arc_table'],
                                      arc_support_pixels=summary['candidate_detection']['arc_support_pixels'])
            stage='a6_lattice_seed'
            seed=infer_lattice_seed(candidates,a5.central_disk_center_px,geometry,config=detection_config)
            summary['a6']=asdict(seed)
            summary['a6'].update(candidate_count=len(seed.assignments),
                inlier_count=sum(a.inlier for a in seed.assignments),
                outlier_count=sum(not a.inlier for a in seed.assignments))
            summary['outputs']['a6_assignments']=write_seed_diagnostics(output,image,seed)
            stage='a7_lattice_and_catalogue'
            if seed.seed_method=='two_validated_arc_centers':
                stage='a7_bootstrap'
                bootstrap=bootstrap_lattice(image,seed,geometry)
                summary['bootstrap']=asdict(bootstrap)
                summary['outputs'].update(write_bootstrap_diagnostics(output,image,bootstrap,geometry))
                from .arc_catalogue import generate_arc_catalogue, refinement_comparison
                arc_catalogue=generate_arc_catalogue(bootstrap,image.shape,geometry)
                summary['arc_catalogue']=asdict(arc_catalogue)
                summary['refinement_comparison']=refinement_comparison(bootstrap)
                summary['counts']=arc_catalogue.status_counts
                summary['outputs'].update(write_bootstrap_diagnostics(output,image,arc_catalogue,geometry,catalogue=True))
                if bootstrap.refined_lattice is None:
                    raise _BootstrapIncomplete(bootstrap.message)
                stage='a7_lattice_and_catalogue'
                summary['a7']=asdict(bootstrap.refined_lattice)
                summary.update(status='completed',beam_count=len(arc_catalogue.sites),
                    stability_assessment='Global fit of independently arc-supported observations; review residuals')
            else:
                catalogue=refine_lattice_and_generate_catalogue(image,seed,geometry)
                summary['a7']=asdict(catalogue.lattice)
                rows=write_beam_table(output/'beam_catalogue.csv',catalogue)
                crops=write_figures(output,image,a5,seed,catalogue,crop_scale=crop_scale)
                shifts=[r['predicted_to_observed_shift_px'] for r in rows if r['predicted_to_observed_shift_px'] is not None]
                summary.update(status='completed',stability_assessment='Single completed run; review residuals and support, no radius-sensitivity claim',beam_count=len(rows),counts=catalogue.status_counts,
                    max_shift_px=max(shifts) if shifts else None,median_shift_px=float(np.median(shifts)) if shifts else None,
                    largest_shifts=sorted([r for r in rows if r['predicted_to_observed_shift_px'] is not None],
                                         key=lambda r:-r['predicted_to_observed_shift_px'])[:10],
                    lowest_response_observed=sorted([r for r in rows if r['status']=='observed'],
                                                    key=lambda r:r['response_score'])[:10],
                    missing_sites=[(r['n1'],r['n2']) for r in rows if r['status']=='missing'],crops=crops,
                    outputs=dict(**summary['outputs'],beam_table=str(output/'beam_catalogue.csv'),overview=str(output/'beam_catalogue_overview.png'),
                                 crops=str(output/'beams')))
    except _BootstrapIncomplete as exc:
        summary['status']='bootstrap_only';summary['blockers'].append(str(exc))
        summary['stability_assessment']='Predictions and local evidence only; insufficient/rejected global refinement'
    except (ValueError,TypeError,OSError,RuntimeError,np.linalg.LinAlgError) as exc:
        summary['status']='failed'; summary['failed_stage']=stage; summary['blockers'].append(f'{type(exc).__name__}: {exc}')
    text=json.dumps(summary,indent=2,default=_json_default,allow_nan=False)
    (output/'validation_summary.json').write_text(text+'\n')
    lines=['# A7 reference validation', '', 'Status: '+summary['status'], '',
           f"Input: {summary['input_file']}", f"Shape: {summary.get('input_shape', 'unavailable')}", '',
           *['- '+b for b in summary['blockers']], '',
           'Geometry input: '+str(summary.get('geometry_input')), '',
           'A5: '+str(summary['a5']), '',
           'Candidate detection: '+str(summary.get('candidate_detection')), '',
           'A6: '+str(None if summary['a6'] is None else {k:summary['a6'][k] for k in
               ('seed_method','primitive_basis_certain','initial_basis_vectors_px','basis_angle_deg',
                'basis_condition_number','basis_vector_1_px','basis_vector_2_px','candidate_count','inlier_count','outlier_count',
                'support_fraction','fit_rms_residual_px')}), '',
           'A7: '+str(None if summary['a7'] is None else {k:summary['a7'][k] for k in
               ('origin_px','basis_vector_1_px','basis_vector_2_px','initial_rms_residual_px','fit_rms_residual_px')}), '',
           'Stability: '+summary['stability_assessment'], '',
           'Beam counts: '+str(summary.get('counts','unavailable')), '',
           'Shift median/max: '+str((summary.get('median_shift_px'),summary.get('max_shift_px'))), '',
           'Outputs: '+str(summary['outputs']), '',
           'Full configuration, provenance, assignments and review diagnostics: validation_summary.json',
           'No automatic scientific conclusion. No A2–A4 or A8 operations performed.']
    (output/'validation_summary.md').write_text('\n'.join(lines)+'\n')
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pattern',type=Path,default=Path('reference_data/pattern.dat'))
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--expected-shape',type=int,nargs=2,default=None,metavar=('NY','NX'))
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--profile',type=Path)
    group.add_argument('--radii',type=float,nargs=2,metavar=('RX','RY'))
    parser.add_argument('--geometry-provenance')
    parser.add_argument('--edge-smoothing-sigma-px',type=float,default=1.)
    parser.add_argument('--edge-percentile',type=float,default=85.)
    parser.add_argument('--radius-tolerance-px',type=float,default=2.)
    parser.add_argument('--crop-scale',type=float,default=1.25)
    args=parser.parse_args()
    result=validate_pattern(args.pattern,args.output,expected_shape=args.expected_shape,profile=args.profile,
        radii=args.radii,geometry_provenance=args.geometry_provenance,crop_scale=args.crop_scale,
        detection_config=LatticeDetectionConfig(edge_smoothing_sigma_px=args.edge_smoothing_sigma_px,
            edge_percentile=args.edge_percentile,radius_tolerance_px=args.radius_tolerance_px))
    print(result['status']); print('\n'.join(result['blockers']))
    print(args.output/'validation_summary.json')
    return 0 if result['status']=='completed' else 1


if __name__=='__main__':
    raise SystemExit(main())
