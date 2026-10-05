"""Prediction is not observation: locally validate new sites before the A7 fit."""
from dataclasses import dataclass, replace
import numpy as np
from ..detector.models import _real
from .noise import _pattern
from .boundary_sites import boundary_fraction, enumerate_boundary_sites
from .lattice_detection import (LatticeCandidate, LatticeAssignment, ReciprocalLatticeSeed)
from .arc_geometry import arc_edge_support, _validated_center
from .lattice_refinement import (LatticeRefinementConfig, RefinedReciprocalLattice,
                                 enumerate_seed_sites, refine_reciprocal_lattice)


@dataclass(frozen=True)
class BootstrapConfig:
    search_radius_px: float = 4.0
    search_step_px: float = 1.0
    minimum_boundary_fraction: float = .25

    def __post_init__(self):
        for name in ('search_radius_px','search_step_px'):
            value=_real(getattr(self,name),name)
            if value<=0: raise ValueError(f'{name} must be positive')
            object.__setattr__(self,name,value)
        fraction=_real(self.minimum_boundary_fraction,'minimum_boundary_fraction')
        if not 0<fraction<=1: raise ValueError('minimum_boundary_fraction must be in (0,1]')
        object.__setattr__(self,'minimum_boundary_fraction',fraction)
        if self.search_step_px>self.search_radius_px:
            raise ValueError('search_step_px must not exceed search_radius_px')
        if self.search_radius_px/self.search_step_px>100:
            raise ValueError('Local search exceeds 100 steps per radius')


@dataclass(frozen=True)
class BootstrapSite:
    integer_indices: tuple[int,int]
    predicted_center_px: tuple[float,float]
    observed_center_px: tuple[float,float] | None
    shift_px: float | None
    status: str
    source: str
    candidate: LatticeCandidate | None = None
    reason: str = ''
    boundary_fraction: float = 1.0
    usable_boundary_fraction: float = 1.0
    boundary_limited: bool = False
    beam_id: str = ""
    initial_integer_indices: tuple[int,int] | None = None
    lattice_residual_px: float | None = None


@dataclass(frozen=True)
class BootstrapResult:
    initial_seed: ReciprocalLatticeSeed
    sites: tuple[BootstrapSite,...]
    expanded_seed: ReciprocalLatticeSeed
    newly_confirmed_count: int
    required_observations: int
    config: BootstrapConfig
    refined_lattice: RefinedReciprocalLattice | None
    status: str
    message: str


def bootstrap_lattice(corrected_pattern, seed, geometry, *, config=None, refinement_config=None):
    """Search seed-predicted sites, then call the unchanged fit only with redundancy.

    No whole-image candidate detection, intensity maxima, or invented observations.
    Original A6 assignments are retained once. The known central reference is
    displayed but is not counted as a new independent diffracted observation.
    """
    config=BootstrapConfig() if config is None else config
    refinement_config=LatticeRefinementConfig() if refinement_config is None else refinement_config
    if not isinstance(config,BootstrapConfig) or not isinstance(refinement_config,LatticeRefinementConfig):
        raise TypeError('Expected bootstrap and refinement configurations')
    if not isinstance(seed,ReciprocalLatticeSeed): raise TypeError('Expected ReciprocalLatticeSeed')
    image=_pattern(corrected_pattern)
    matrix=np.column_stack((seed.basis_vector_1_px,seed.basis_vector_2_px))
    sites=enumerate_boundary_sites(seed.reference_center_px,matrix,image.shape,geometry,
        max_sites=refinement_config.max_enumeration_sites)
    # Keep local windows disjoint for the canonical reduced basis. No silent clipping.
    if config.search_radius_px>=.25*min(np.linalg.norm(matrix,axis=0)):
        raise ValueError('Bootstrap search radius must be < one quarter of the seed basis lengths')
    arrays,_=arc_edge_support(image,geometry,seed.reference_center_px,config=seed.config)
    original={a.integer_indices:a for a in seed.assignments if a.inlier}
    if len(original)<2:
        raise ValueError('Bootstrap requires at least two defining observations')
    if (len(original)!=sum(a.inlier for a in seed.assignments) or
        len({a.candidate.id for a in seed.assignments})!=len(seed.assignments) or
        len({a.candidate.center_px for a in original.values()})!=len(original)):
        raise ValueError('Duplicate defining observations or lattice sites')
    assignments=list(seed.assignments)
    used_ids={a.candidate.id for a in assignments}
    used_pixels=set(p for a in assignments if a.inlier for p in a.candidate.supporting_pixels)
    records=[];new_count=0
    n=int(np.floor(config.search_radius_px/config.search_step_px))
    offsets=[(i*config.search_step_px,j*config.search_step_px) for j in range(-n,n+1) for i in range(-n,n+1)
             if np.hypot(i*config.search_step_px,j*config.search_step_px)<=config.search_radius_px+1e-12]
    margin=int(np.ceil(4*seed.config.edge_smoothing_sigma_px))+1
    for indices,predicted in sites:
        fraction=boundary_fraction(predicted,image.shape,geometry)
        usable=boundary_fraction(predicted,image.shape,geometry,margin=margin)
        visibility=dict(boundary_fraction=fraction,usable_boundary_fraction=usable)
        if indices==(0,0):
            records.append(BootstrapSite(indices,predicted,None,None,'reference','central_reference',reason='Known A5 center; not independent validation',**visibility))
            continue
        if indices in original:
            c=original[indices].candidate
            records.append(BootstrapSite(indices,predicted,c.center_px,float(np.linalg.norm(np.array(c.center_px)-predicted)),
                                         'defining','a6_defining',c,'Defines seed; not new validation',**visibility))
            continue
        best=None;clipped=False;searched=False
        for offset in offsets:
            center=np.array(predicted)+offset
            if boundary_fraction(center,image.shape,geometry,margin=margin)<config.minimum_boundary_fraction:
                clipped=True;continue
            searched=True
            metrics=_validated_center(center,arrays,geometry,seed.config)
            if metrics is None: continue
            key=(metrics['confidence'],-np.linalg.norm(offset),-center[0],-center[1])
            if best is None or key>best[0]: best=(key,tuple(center),metrics)
        if best is None:
            records.append(BootstrapSite(indices,predicted,None,None,'edge_clipped' if not searched else 'missing',
                                         'a7_bootstrap_prediction',reason='Insufficient usable circumference' if not searched else 'No validated local arc support',**visibility))
            continue
        _,center,metrics=best
        shift=float(np.linalg.norm(np.array(center)-predicted))
        boundary=shift>=config.search_radius_px-config.search_step_px*.5
        overlap=bool(used_pixels.intersection(metrics['supporting_pixels']))
        identifier=f'bootstrap_{indices[0]}_{indices[1]}'
        if identifier in used_ids: raise ValueError('Bootstrap candidate ID collision')
        # Candidate is local geometric evidence; only confirmed records have an observed center.
        candidate=LatticeCandidate(identifier,center,metrics['confidence'],**metrics)
        if boundary or overlap:
            records.append(BootstrapSite(indices,predicted,None,None,'weak_evidence','a7_bootstrap_prediction',candidate,
                'Best support at search boundary' if boundary else 'Support overlaps an existing observation',
                boundary_limited=boundary,**visibility))
            continue
        assignments.append(LatticeAssignment(candidate,indices,shift,True))
        used_ids.add(identifier);used_pixels.update(candidate.supporting_pixels);new_count+=1
        records.append(BootstrapSite(indices,predicted,center,shift,'confirmed','a7_bootstrap_confirmed',candidate,**visibility))
    inlier_residuals=[a.residual_px for a in assignments if a.inlier]
    expanded=replace(seed,assignments=tuple(assignments),
        support_fraction=float(np.mean([a.inlier for a in assignments])),
        fit_rms_residual_px=float(np.sqrt(np.mean(np.square(inlier_residuals)))) if inlier_residuals else 0.)
    count=sum(a.inlier for a in assignments)
    if count<refinement_config.min_inliers:
        return BootstrapResult(seed,tuple(records),expanded,new_count,refinement_config.min_inliers,config,None,
            'insufficient_observations','Bootstrap lattice generated, but insufficient independent observations for global refinement')
    try:
        refined=refine_reciprocal_lattice(expanded,config=refinement_config)
    except ValueError as exc:
        return BootstrapResult(seed,tuple(records),expanded,new_count,refinement_config.min_inliers,config,None,
                               'global_fit_rejected',str(exc))
    return BootstrapResult(seed,tuple(records),expanded,new_count,refinement_config.min_inliers,config,refined,
                           'refined','Existing global fit completed after bootstrap confirmation')
