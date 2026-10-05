"""A7 catalogue of recorded arc observations; no intensity-based confirmations."""
from dataclasses import dataclass, replace
import numpy as np
from .boundary_sites import enumerate_boundary_sites, boundary_fraction
from .lattice_bootstrap import BootstrapSite


@dataclass(frozen=True)
class ArcBeamCatalogue:
    sites: tuple[BootstrapSite,...]
    status: str
    refined: bool

    @property
    def status_counts(self):
        return {name:sum(s.status==name for s in self.sites)
                for name in ('reference','observed','weak','missing','edge_clipped','rejected')}


def generate_arc_catalogue(bootstrap, image_shape, geometry):
    """Reproject recorded observations using existing fit labels; never invent data.

    Without a fit this is explicitly a partial seed catalogue. With a fit, all
    footprint-intersecting predictions are regenerated. Sites newly entering the
    frame remain unobserved. Candidate IDs and old labels retain provenance.
    """
    fit=bootstrap.refined_lattice;seed=bootstrap.initial_seed
    origin=seed.reference_center_px if fit is None else fit.origin_px
    matrix=np.column_stack((seed.basis_vector_1_px,seed.basis_vector_2_px) if fit is None
                           else (fit.basis_vector_1_px,fit.basis_vector_2_px))
    original_matrix=np.column_stack((seed.basis_vector_1_px,seed.basis_vector_2_px))
    assignments={} if fit is None else {a.candidate.id:a for a in fit.assignments}
    # Recover the integer basis transform from accepted fit label correspondences,
    # not from rounded physical distances (which can change near a tie).
    transform=np.eye(2,dtype=int)
    if fit is not None:
        a=np.array([r.initial_integer_indices for r in fit.assignments if not r.reassigned])
        b=np.array([r.integer_indices for r in fit.assignments if not r.reassigned])
        if len(a)>=2 and np.linalg.matrix_rank(a)==2:
            transform=np.rint(np.linalg.lstsq(a,b,rcond=None)[0]).astype(int)
            if not np.array_equal(a@transform,b): raise ValueError('Inconsistent canonical label mapping')
        else:
            # No reliable mapping for unobserved old sites: accepted candidates
            # still use the fitter's exact labels, predictions remain unobserved.
            transform=None
    mapped={}
    for site in bootstrap.sites:
        assignment=None if site.candidate is None else assignments.get(site.candidate.id)
        index=assignment.integer_indices if assignment is not None else (
            tuple(int(v) for v in np.array(site.integer_indices)@transform) if transform is not None else None)
        if index is None: continue
        if index in mapped: raise ValueError('Ambiguous catalogue label mapping')
        mapped[index]=(site,assignment)
    limit=fit.config.max_enumeration_sites if fit is not None else 100000
    sites=[]
    for index,point in enumerate_boundary_sites(origin,matrix,image_shape,geometry,max_sites=limit):
        previous,assignment=mapped.get(index,(None,None))
        usable=boundary_fraction(point,image_shape,geometry,margin=int(np.ceil(4*seed.config.edge_smoothing_sigma_px))+1)
        fraction=boundary_fraction(point,image_shape,geometry)
        observed=None;candidate=None;residual=None;source='lattice_prediction';reason='No recorded arc confirmation'
        status='edge_clipped' if usable<bootstrap.config.minimum_boundary_fraction else 'missing'
        boundary=False
        if previous is not None:
            candidate=previous.candidate;source=previous.source;reason=previous.reason;boundary=previous.boundary_limited
            status={'defining':'observed','confirmed':'observed','weak_evidence':'weak'}.get(previous.status,previous.status)
            observed=previous.observed_center_px
        if assignment is not None:
            residual=assignment.residual_px
            if not assignment.inlier: status='rejected';observed=None;reason='Rejected by global lattice fit'
        if index==(0,0): status='reference';observed=None
        sites.append(BootstrapSite(index,point,observed,
            None if observed is None else float(np.linalg.norm(np.array(observed)-point)),status,source,candidate,reason,
            fraction,usable,boundary,f'lattice_{index[0]}_{index[1]}',
            None if previous is None else previous.integer_indices,residual))
    return ArcBeamCatalogue(tuple(sites),'refined arc catalogue' if fit is not None else 'partial seed catalogue — global fit unavailable',fit is not None)


def refinement_comparison(bootstrap):
    """Measured seed/fit differences, or explicit absent quantities without a fit."""
    fit=bootstrap.refined_lattice;seed=bootstrap.initial_seed
    if fit is None:
        return dict(origin_refined=None,g1_refined=None,g2_refined=None,rms_before=None,rms_after=None,
                    origin_change=None,g1_change=None,g2_change=None,length_changes=None,angle_change_deg=None)
    old=np.array([seed.basis_vector_1_px,seed.basis_vector_2_px])
    new=np.array([fit.basis_vector_1_px,fit.basis_vector_2_px])
    def angle(v):
        return float(np.degrees(np.arccos(np.clip(np.dot(*v)/np.prod(np.linalg.norm(v,axis=1)),-1,1))))
    return dict(origin_refined=fit.origin_px,g1_refined=fit.basis_vector_1_px,g2_refined=fit.basis_vector_2_px,
        rms_before=fit.initial_rms_residual_px,rms_after=fit.fit_rms_residual_px,
        origin_change=tuple(np.array(fit.origin_px)-seed.reference_center_px),
        g1_change=tuple(new[0]-old[0]),g2_change=tuple(new[1]-old[1]),
        length_changes=tuple(np.linalg.norm(new,axis=1)-np.linalg.norm(old,axis=1)),
        angle_change_deg=angle(new)-angle(old))
