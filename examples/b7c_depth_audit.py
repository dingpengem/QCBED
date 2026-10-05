"""Audit requested depth before enabling a depth refinement binding."""
from dataclasses import replace, asdict
import numpy as np
from qcbed.depth_geometry import requested_depth_to_sample_geometry
from qcbed import compile_config_from_sources
from qcbed.objective_diagnostics import evaluate_preflight_objectives
from qcbed.single_target_refinement import audit_binary_masks
if __package__:
    from .b7a_refinement_layout import real_config
    from .a9_real_variance_closure import prepare_reference_variance
else:
    from b7a_refinement_layout import real_config
    from a9_real_variance_closure import prepare_reference_variance


def depth_audit():
    baseline=real_config(); compiled=compile_config_from_sources(baseline)
    data,closure=prepare_reference_variance()
    thicknesses=tuple(s.thickness_angstrom for s in compiled.crystallography.simulation_sample.depth_sample.slices)
    extent=requested_depth_to_sample_geometry(840.,thicknesses)
    lo,hi=extent.plateau_lower_exclusive,extent.effective_depth
    requests=sorted(set((838.,839.,840.,841.,842.,lo,np.nextafter(lo,np.inf),hi,np.nextafter(hi,np.inf))))
    rows=[]
    for requested in requests:
        cfg=replace(baseline,crystallography=replace(baseline.crystallography,depth_angstrom=float(requested)))
        current=compile_config_from_sources(cfg)
        r,_=evaluate_preflight_objectives(current,{'conventional':data['conventional']},
            baseline_manifests={'conventional':closure['modes']['conventional']['dataset']},
            input_kind='confirmed_real',physical_input_confirmation='Confirmed FCC Al depth audit',verify_repeat=False)
        geometry=requested_depth_to_sample_geometry(requested,thicknesses)
        if tuple(r['preflight']['sample']['sequence'])!=geometry.sequence:raise ValueError('depth mapping differs')
        result=r['modes']['conventional']
        rows.append(dict(requested_depth=float(requested),effective_depth=geometry.effective_depth,
            occurrence_count=geometry.occurrence_count,stored_slice_count=geometry.stored_slice_count,
            plateau_lower_exclusive=geometry.plateau_lower_exclusive,
            ynorm=result['ynorm'],objective=result['objective']))
    for a in rows:
        for b in rows:
            if a['occurrence_count']==b['occurrence_count']:
                if (a['effective_depth'],a['ynorm'],a['objective'])!=(b['effective_depth'],b['ynorm'],b['objective']):
                    raise ValueError('same sequence did not reproduce exact plateau')
    initial=next(r for r in rows if r['requested_depth']==840.)
    if (initial['ynorm'],initial['objective'])!=(14870390.0181119,.939675432747781):
        raise ValueError('baseline drift')
    return dict(status='PASS',depth_gate='READY',mask_audit=audit_binary_masks(data),rows=rows,
        candidate_depths=[lo,840.,float(np.nextafter(hi,np.inf))],bounds=[600.,1200.],
        bounds_source='examples/reference_configuration.py requested region-depth bounds; same strict-< whole-slice expansion',
        strategy='Explicit three-depth candidate search; unchanged 1D PARAMIN for B within each candidate')
