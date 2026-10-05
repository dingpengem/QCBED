"""Bounded center-only refinement of frozen, already validated arc evidence.

No intensities, beam indexing, scan corners or historical targets enter the fit.
Ellipse radius errors and normal signs follow the existing arc validator.
"""
from dataclasses import dataclass, replace
import numpy as np
from scipy.optimize import least_squares
from .models import DiskGeometry, _pair
from ..detector.models import _real


@dataclass(frozen=True)
class SubpixelArcConfig:
    bound_px: float = 1.
    normal_scale_deg: float = 12.
    vote_scale_px: float = 2.5

    def __post_init__(self):
        for name in ('bound_px','normal_scale_deg','vote_scale_px'):
            value=_real(getattr(self,name),name)
            if value<=0: raise ValueError(f'{name} must be positive')
            object.__setattr__(self,name,value)
        if self.normal_scale_deg>=90: raise ValueError('normal_scale_deg must be <90')


@dataclass(frozen=True)
class SubpixelArcResult:
    initial_center_px: tuple[float,float]
    refined_center_px: tuple[float,float]
    delta_px: tuple[float,float]
    delta_magnitude: float
    objective_before: float
    objective_after: float
    accepted: bool
    reason: str
    config: SubpixelArcConfig
    support_count: int
    support_points_px: tuple[tuple[float,float], ...] = ()
    bound_active: bool = False
    method: str = 'deterministic bounded least_squares; fixed radius and support membership'


def refine_arc_center(initial_center_px, points_px, normals, geometry, *, config=None):
    """Float64 center-only least squares, frozen support and fixed normal polarity.

    Per-pixel residuals: signed radial-ray distance in pixels; normal cross product
    scaled by sin(normal_scale_deg); two center-vote components / vote_scale_px.
    The sum of squared residuals divided by support count is reported.
    """
    cfg=SubpixelArcConfig() if config is None else config
    if not isinstance(cfg,SubpixelArcConfig) or not isinstance(geometry,DiskGeometry):
        raise TypeError('Expected SubpixelArcConfig and DiskGeometry')
    initial=np.array(_pair(initial_center_px,'initial_center_px'))
    points=np.array(points_px,dtype=float,copy=True); normal=np.array(normals,dtype=float,copy=True)
    if points.ndim!=2 or points.shape[1]!=2 or points.shape!=normal.shape or len(points)<3:
        raise ValueError('Need matching Nx2 support points/normals with N>=3')
    if not np.isfinite(points).all() or not np.isfinite(normal).all(): raise ValueError('Support must be finite')
    norms=np.linalg.norm(normal,axis=1)
    if np.any(norms==0): raise ValueError('Normals must be nonzero')
    normal/=norms[:,None]
    axes=np.array([geometry.radius_x_px,geometry.radius_y_px])
    polarity=np.where(np.sum(normal*(points-initial),axis=1)>=0,1.,-1.)
    implied=points-polarity[:,None]*normal*axes**2/np.linalg.norm(normal*axes,axis=1)[:,None]
    def residual(center):
        d=points-center
        distance=np.linalg.norm(d,axis=1); radial=np.linalg.norm(d/axes,axis=1)
        radius=distance-distance/np.maximum(radial,np.finfo(float).tiny)
        expected=d/axes**2
        expected/=np.maximum(np.linalg.norm(expected,axis=1)[:,None],np.finfo(float).tiny)
        cross=normal[:,0]*expected[:,1]-normal[:,1]*expected[:,0]
        return np.concatenate((radius,cross/np.sin(np.deg2rad(cfg.normal_scale_deg)),
                               ((implied-center)/cfg.vote_scale_px).ravel()))
    before=float(np.sum(residual(initial)**2)/len(points))
    def result(center,after,accepted,reason,bound=False):
        delta=center-initial
        return SubpixelArcResult(tuple(initial),tuple(center),tuple(delta),float(np.linalg.norm(delta)),
            before,after,accepted,reason,cfg,len(points),tuple(tuple(float(v) for v in p) for p in points),bound)
    try:
        fit=least_squares(residual,initial,bounds=(initial-cfg.bound_px,initial+cfg.bound_px),
            method='trf',ftol=1e-12,xtol=1e-12,gtol=1e-12,max_nfev=200)
        after=float(np.sum(residual(fit.x)**2)/len(points))
        condition=float(np.linalg.cond(fit.jac))
        if (not fit.success or not np.isfinite(after) or not np.isfinite(condition) or condition>1e8
                or np.any(abs(fit.x-initial)>cfg.bound_px) or after>before):
            return result(initial,before,False,'Unstable, out-of-bound or non-improving fit; initial center retained')
        return result(fit.x,after,True,'Geometric objective non-increasing',bool(np.any(fit.active_mask)))
    except (ValueError,FloatingPointError,np.linalg.LinAlgError):
        return result(initial,before,False,'Numerical failure; initial center retained')


def refine_validated_candidate(candidate, arrays, geometry, arc_config, *, config=None):
    """Freeze candidate support; reapply unchanged acceptance criteria once after fit.

    No new edge pixels enter. If the final existing validator fails on frozen
    evidence, keep the original candidate and report the rejected proposal.
    """
    from .arc_geometry import _validated_center
    indices={(int(x),int(y)):i for i,(x,y) in enumerate(zip(arrays['edge_x'],arrays['edge_y']))}
    if not candidate.supporting_pixels: raise ValueError('Candidate has no validated support')
    try: ids=np.array([indices[p] for p in candidate.supporting_pixels])
    except KeyError as exc: raise ValueError('Validated support absent from supplied edge evidence') from exc
    points=np.column_stack((arrays['edge_x'][ids],arrays['edge_y'][ids]))
    normals=np.column_stack((arrays['normal_x'][ids],arrays['normal_y'][ids]))
    result=refine_arc_center(candidate.center_px,points,normals,geometry,config=config)
    frozen={k:arrays[k][ids] for k in ('edge_x','edge_y','normal_x','normal_y')}
    metrics=_validated_center(np.array(result.refined_center_px),frozen,geometry,arc_config) if result.accepted else None
    if result.accepted and metrics is None:
        result=replace(result,refined_center_px=result.initial_center_px,delta_px=(0.,0.),delta_magnitude=0.,
            objective_after=result.objective_before,accepted=False,reason='Existing arc validation failed; original center retained')
    if not result.accepted: return replace(candidate,subpixel_refinement=result)
    return replace(candidate,center_px=result.refined_center_px,subpixel_refinement=result,**metrics)
