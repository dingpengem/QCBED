"""Fixed-radius subpixel experiment; A10 comparison is downstream only."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import numpy as np
from qcbed.data_preparation import (DiskGeometry,SubpixelArcConfig,locate_central_disk,
    detect_lattice_candidates,infer_lattice_seed,bootstrap_lattice,generate_arc_catalogue,
    ExtractionConfig,SignalMode,prepare_experimental_qcbed,scan_from_catalogue)
from qcbed.data_preparation.arc_geometry import arc_edge_support
from qcbed.data_preparation.subpixel_arc import refine_validated_candidate
from qcbed.data_preparation.bootstrap_diagnostics import write_bootstrap_diagnostics


def main():
    p=argparse.ArgumentParser(description=__doc__)
    root=Path(__file__).resolve().parents[1]
    p.add_argument('--pattern',type=Path,default=root/'reference_data/pattern.dat')
    p.add_argument('--a8-measurement',type=Path,default=root/'tmp/a8_manual_selection/measurement.json')
    p.add_argument('--accepted-validation',type=Path,default=root/'tmp/a7_boundary_catalogue_r38_final/validation_summary.json')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--bound',type=float,default=1.)
    a=p.parse_args(); a.output.mkdir(parents=True,exist_ok=False)
    image=np.loadtxt(a.pattern); geometry=DiskGeometry(38,38); cfg=SubpixelArcConfig(bound_px=a.bound)
    original=locate_central_disk(image,geometry)
    central=locate_central_disk(image,geometry,subpixel_config=cfg)
    # Keep the previously reviewed/revalidated candidate selection, not a fresh
    # automatic proposal list (the original real seed was explicitly selected).
    from qcbed.data_preparation import LatticeCandidate,LatticeDetectionConfig,BootstrapConfig
    accepted=json.loads(a.accepted_validation.read_text())
    import hashlib
    if hashlib.sha256(a.pattern.read_bytes()).hexdigest()!=accepted['sha256']:
        raise ValueError('Accepted arc evidence belongs to a different input pattern')
    detected=tuple(LatticeCandidate(**v['candidate']) for v in accepted['seed']['assignments'])
    arc=LatticeDetectionConfig(**accepted['seed']['config'])
    arrays,_=arc_edge_support(image,geometry,original.central_disk_center_px,config=arc)
    refined=tuple(refine_validated_candidate(c,arrays,geometry,arc,config=cfg) for c in detected)
    oldseed=infer_lattice_seed(detected,original.central_disk_center_px,geometry,config=arc)
    seed=infer_lattice_seed(refined,central.central_disk_center_px,geometry,config=arc)
    bootstrap=bootstrap_lattice(image,seed,geometry,config=BootstrapConfig(**accepted['bootstrap']['config']))
    catalogue=generate_arc_catalogue(bootstrap,image.shape,geometry)
    write_bootstrap_diagnostics(a.output,image,catalogue,geometry,catalogue=True)
    # Source unchanged A9; exact fractional-origin failure must not be bypassed.
    extraction=ExtractionConfig.from_disk_geometry(geometry,final_binning=2,diagnostic_window_shape_px=(64,64))
    ids=[s.beam_id for s in catalogue.sites if s.status!='rejected']
    from qcbed.data_preparation.extraction import _roi_counts
    from qcbed.data_preparation import ReflectionROI
    coverage=[]
    for s in catalogue.sites:
        n,mask=_roi_counts(image.shape,ReflectionROI(s.beam_id,s.beam_id,s.predicted_center_px,DiskGeometry(34.2,34.2)))
        coverage.append(dict(beam_id=s.beam_id,indices=s.integer_indices,coverage=float(mask.sum()/n),
            retained=bool(mask.sum()/n>=extraction.minimum_beam_frame_coverage),status=s.status,
            predicted_center=s.predicted_center_px,observed_center=s.observed_center_px))
    a9={}
    for mode in SignalMode:
        try:
            result=prepare_experimental_qcbed(image,catalogue,geometry,selected_beam_ids=ids,config=extraction,signal_mode=mode)
            a9[mode.value]=dict(status='extracted',retained=sum(b.status=='extracted' for b in result.beams))
        except ValueError as exc:
            if 'fractional selection origin' not in str(exc): raise
            a9[mode.value]=dict(status='blocked',reason=str(exc))
    origin=json.loads(a.a8_measurement.read_text())['zone_axis_origin_px']
    scan=scan_from_catalogue(catalogue,central_beam_id='lattice_0_0',canonical_basis_beam_ids=('lattice_1_0','lattice_0_1'),
        zone_axis_origin_px=origin,extraction_config=extraction,basis_transform=((1,0),(0,-1)))
    old=np.array([oldseed.basis_vector_1_px,oldseed.basis_vector_2_px]);new=np.array([seed.basis_vector_1_px,seed.basis_vector_2_px])
    def angle(matrix): return float(np.rad2deg(np.arccos(np.clip(matrix[0]@matrix[1]/np.prod(np.linalg.norm(matrix,axis=1)),-1,1))))
    # Historical values enter ONLY after optimization, indexing, A7 and A10.
    historical=np.array([[-.508156,.407641],[.407641,.508156],[-.407641,-.508156]])
    corners=np.array([scan.scan_geometry.corner_a,scan.scan_geometry.corner_b,scan.scan_geometry.corner_c])
    stages={}
    for label,path in [('A_diagnostic64',root/'tmp/a10_current_binning2/scan_geometry.json'),
                       ('B_integer76',root/'tmp/a10_physical76_binning2/scan_geometry.json')]:
        previous=json.loads(path.read_text()); stages[label]=dict(scan_geometry=previous['scan_geometry'],
            maximum_absolute_difference=previous['maximum_absolute_difference'])
    stages['C_subpixel76']=dict(scan_geometry=asdict(scan.scan_geometry),differences=(corners-historical).tolist(),
        absolute_differences=abs(corners-historical).tolist(),maximum_absolute_difference=float(abs(corners-historical).max()))
    report=dict(centers=[asdict(central.arc_support.subpixel_refinement)]+[asdict(c.subpixel_refinement) for c in refined],
        physical_radii=asdict(geometry),old_basis=old.tolist(),new_basis=new.tolist(),component_changes=(new-old).tolist(),
        length_changes=(np.linalg.norm(new,axis=1)-np.linalg.norm(old,axis=1)).tolist(),angle_change_deg=angle(new)-angle(old),
        a7_status=bootstrap.status,a7_message=bootstrap.message,catalogue=coverage,a9=a9,a10=asdict(scan),comparison=stages,
        a8='Retained independently validated direct-beam symmetry origin; no origin rerun or geometric-center substitution',
        protocol='One fixed objective/configuration run; no historical targets enter geometric fit')
    (a.output/'subpixel_comparison.json').write_text(json.dumps(report,indent=2)+'\n')
    (a.output/'subpixel_comparison.md').write_text('# Subpixel geometric validation\n\n```json\n'+json.dumps(report,indent=2)+'\n```\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse
    for index,c in enumerate((central.arc_support,)+refined):
        fit=c.subpixel_refinement
        fig,ax=plt.subplots(figsize=(7,7));ax.imshow(image,cmap='gray')
        ax.scatter(*fit.initial_center_px,marker='+',s=120,c='yellow',label='Grid center')
        ax.scatter(*fit.refined_center_px,marker='x',s=90,c='red',label='Subpixel center')
        ax.add_patch(Ellipse(fit.refined_center_px,76,76,fill=False,color='cyan'))
        support=np.array(fit.support_points_px);ax.scatter(support[:,0],support[:,1],s=7,c='lime',label='Arc support')
        x,y=fit.initial_center_px;ax.set_xlim(x-42,x+42);ax.set_ylim(y+42,y-42)
        ax.set_title(f'{c.id}: dx={fit.delta_px[0]:+.5f}, dy={fit.delta_px[1]:+.5f} px')
        ax.legend(loc='upper left', fontsize=8)
        inset=ax.inset_axes([.65,.04,.31,.31])
        inset.imshow(image,cmap='gray',interpolation='nearest')
        inset.plot([x,fit.refined_center_px[0]],[y,fit.refined_center_px[1]],color='cyan')
        inset.scatter(x,y,marker='+',s=120,c='yellow')
        inset.scatter(*fit.refined_center_px,marker='x',s=90,c='red')
        inset.set_xlim(x-1.4,x+1.4);inset.set_ylim(y+1.4,y-1.4)
        inset.set_title('Center detail (px)',fontsize=8);inset.tick_params(labelsize=6)
        fig.savefig(a.output/f'center_{index}.png',dpi=160);plt.close(fig)
    print(json.dumps({k:report[k] for k in ('old_basis','new_basis','length_changes','angle_change_deg','a7_status','a9','comparison')},indent=2))


if __name__=='__main__': main()
