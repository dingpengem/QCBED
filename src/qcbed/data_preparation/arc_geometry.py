"""Shared known-radius geometry for central localization, A6 and A7 confirmation."""
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from .subpixel_arc import SubpixelArcResult

from dataclasses import dataclass
from itertools import combinations

import numpy as np
from scipy.ndimage import maximum_filter, gaussian_filter, sobel, map_coordinates
from scipy.signal import fftconvolve
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from .models import DiskGeometry, _pair
from .noise import _pattern
from .disk_detection import _center
from ..detector.models import _integer, _real


@dataclass(frozen=True)
class LatticeDetectionConfig:
    edge_smoothing_sigma_px: float = 1.0
    edge_percentile: float = 85.0
    radius_tolerance_px: float = 2.0
    radius_samples: int = 5
    vote_smoothing_sigma_px: float = 1.0
    minimum_vote_strength: float = 0.2
    relative_vote_threshold: float = .3
    candidate_merge_radius_px: float = 6.0
    minimum_support_pixels: int = 12
    minimum_arc_length_px: float = 8.0
    center_vote_rms_scale_px: float = 2.5
    maximum_normal_residual_deg: float = 12.0
    maximum_neighbor_normal_change_deg: float = 25.0
    minimum_angular_coverage: float = .15
    central_exclusion_scale: float = 1.2
    minimum_lattice_vector_scale: float = 1.5
    max_candidates: int = 64
    displacement_tolerance_px: float = 1.5
    lattice_residual_tolerance_px: float = 2.
    min_family_votes: int = 3
    max_families: int = 32
    min_inliers: int = 6
    min_inlier_fraction: float = .6
    min_basis_sine: float = .2

    def __post_init__(self):
        value = _real(self.edge_percentile, 'edge_percentile')
        if not 0 < value < 100:
            raise ValueError('edge_percentile must be between 0 and 100')
        object.__setattr__(self, 'edge_percentile', value)
        tolerance = _real(self.radius_tolerance_px, 'radius_tolerance_px')
        if tolerance < 0:
            raise ValueError('radius_tolerance_px must be nonnegative')
        object.__setattr__(self, 'radius_tolerance_px', tolerance)
        for name in ('radius_samples', 'minimum_support_pixels', 'max_candidates', 'min_family_votes', 'max_families', 'min_inliers'):
            value = _integer(getattr(self, name), name)
            if value < 1:
                raise ValueError(f'{name} must be positive')
            object.__setattr__(self, name, value)
        if self.min_inliers < 3 or self.max_candidates < self.min_inliers or self.max_families < 2:
            raise ValueError('Need >=3 inliers, enough candidates, and >=2 vector families')
        if self.radius_samples < 3 or self.radius_samples % 2 != 1:
            raise ValueError('radius_samples must be odd and >=3')
        for name in ('edge_smoothing_sigma_px', 'vote_smoothing_sigma_px',
                     'relative_vote_threshold', 'minimum_vote_strength', 'candidate_merge_radius_px',
                     'minimum_arc_length_px', 'center_vote_rms_scale_px', 'maximum_normal_residual_deg',
                     'maximum_neighbor_normal_change_deg', 'minimum_angular_coverage', 'central_exclusion_scale', 'minimum_lattice_vector_scale',
                     'displacement_tolerance_px', 'lattice_residual_tolerance_px',
                     'min_inlier_fraction', 'min_basis_sine'):
            value = _real(getattr(self, name), name)
            if value <= 0:
                raise ValueError(f'{name} must be positive')
            if name in ('relative_vote_threshold', 'minimum_angular_coverage', 'min_inlier_fraction', 'min_basis_sine') and value > 1:
                raise ValueError(f'{name} must be <=1')
            if name in ('maximum_normal_residual_deg', 'maximum_neighbor_normal_change_deg') and value >= 90:
                raise ValueError(f'{name} must be <90 degrees')
            object.__setattr__(self, name, value)


@dataclass(frozen=True)
class LatticeCandidate:
    id: str
    center_px: tuple[float, float]
    response_score: float  # compatibility name: now center-vote strength, not intensity
    supporting_edge_pixels: int = 0
    angular_coverage: float = 0.0  # union of supported phase intervals / (2*pi)
    radius_residual_rms_px: float = 0.0

    fragment_count: int = 0
    arc_length_px: float = 0.0
    center_vote_rms_px: float = 0.0
    mean_radius_px: float = 0.0  # circularized radius for ellipse: sqrt(a*b)*elliptic radius
    radius_std_px: float = 0.0
    mean_normal_residual_deg: float = 0.0
    confidence: float = 0.0  # geometric ranking score, not a probability
    supporting_pixels: tuple[tuple[int, int], ...] = ()  # (x,y), immutable review evidence
    radius_offset_px: float = 0.0
    subpixel_refinement: "SubpixelArcResult | None" = None

    @property
    def vote_strength(self):
        return self.response_score

    def __post_init__(self):
        if not isinstance(self.id, str) or not self.id:
            raise ValueError('candidate id must be nonempty text')
        object.__setattr__(self, 'center_px', _pair(self.center_px, 'center_px'))
        object.__setattr__(self, 'response_score', _real(self.response_score, 'response_score'))
        fragments = _integer(self.fragment_count, 'fragment_count')
        if fragments < 0:
            raise ValueError('fragment_count must be nonnegative')
        object.__setattr__(self, 'fragment_count', fragments)
        count = _integer(self.supporting_edge_pixels, 'supporting_edge_pixels')
        coverage = _real(self.angular_coverage, 'angular_coverage')
        residual = _real(self.radius_residual_rms_px, 'radius_residual_rms_px')
        if count < 0 or not 0 <= coverage <= 1 or residual < 0:
            raise ValueError('Invalid arc support diagnostics')
        object.__setattr__(self, 'supporting_edge_pixels', count)
        object.__setattr__(self, 'angular_coverage', coverage)
        object.__setattr__(self, 'radius_residual_rms_px', residual)
        for name in ('arc_length_px', 'center_vote_rms_px', 'mean_radius_px', 'radius_std_px',
                     'mean_normal_residual_deg', 'confidence'):
            value = _real(getattr(self, name), name)
            if value < 0:
                raise ValueError(f'{name} must be nonnegative')
            object.__setattr__(self, name, value)
        object.__setattr__(self, 'radius_offset_px', _real(self.radius_offset_px, 'radius_offset_px'))
        pixels = tuple(tuple(_integer(v, 'supporting pixel') for v in pair) for pair in self.supporting_pixels)
        if any(len(pair) != 2 for pair in pixels):
            raise ValueError('supporting_pixels must be (x,y) pairs')
        object.__setattr__(self, 'supporting_pixels', pixels)


def _config(config):
    if config is None:
        return LatticeDetectionConfig()
    if not isinstance(config, LatticeDetectionConfig):
        raise TypeError('config must be LatticeDetectionConfig')
    return config


def arc_edge_support(corrected_pattern, geometry: DiskGeometry, central_disk_center_px, *, config=None,
                     exclude_central=True, edge_mask=None):
    """Shared edge preparation only: never runs candidate detection or center voting."""
    config = _config(config)
    image = _pattern(corrected_pattern)
    if not isinstance(geometry, DiskGeometry):
        raise TypeError('geometry must be DiskGeometry')
    cx, cy = _center(central_disk_center_px, image.shape)
    rx, ry = geometry.radius_x_px, geometry.radius_y_px
    if config.radius_tolerance_px >= min(rx, ry):
        raise ValueError('radius tolerance must be smaller than both physical semi-axes')
    smooth = gaussian_filter(image, config.edge_smoothing_sigma_px, mode='nearest')
    gx = sobel(smooth, axis=1, mode='nearest')/8
    gy = sobel(smooth, axis=0, mode='nearest')/8
    magnitude = np.hypot(gx, gy)
    yy, xx = np.indices(image.shape)
    margin = int(np.ceil(4*config.edge_smoothing_sigma_px))+1
    valid = (xx>=margin)&(xx<image.shape[1]-margin)&(yy>=margin)&(yy<image.shape[0]-margin)
    if not isinstance(exclude_central, bool):
        raise TypeError('exclude_central must be Boolean')
    # Defaults preserve the A6 diffracted-disk path exactly.
    if exclude_central:
        valid &= ((xx-cx)/rx)**2+((yy-cy)/ry)**2 > config.central_exclusion_scale**2
    if edge_mask is not None:
        if not isinstance(edge_mask,np.ndarray) or edge_mask.dtype!=np.bool_ or edge_mask.shape!=image.shape:
            raise ValueError('edge_mask must be Boolean and match the image shape')
        valid &= edge_mask
    values = magnitude[valid]
    threshold = float(np.percentile(values, config.edge_percentile)) if values.size else 0.
    strong = valid & (magnitude > threshold + 100*np.finfo(float).eps*np.max(abs(image)))
    # Keep derivative ridges, not weak smooth tails of an unrelated small feature.
    ux=np.divide(gx,magnitude,out=np.zeros_like(gx),where=magnitude>0)
    uy=np.divide(gy,magnitude,out=np.zeros_like(gy),where=magnitude>0)
    forward=map_coordinates(magnitude,[yy+uy,xx+ux],order=1,mode='nearest')
    backward=map_coordinates(magnitude,[yy-uy,xx-ux],order=1,mode='nearest')
    selected=strong & (magnitude>=forward) & (magnitude>=backward)
    y, x = np.nonzero(selected)
    nx, ny = gx[selected]/magnitude[selected], gy[selected]/magnitude[selected]
    eligible = ((xx-cx)/rx)**2+((yy-cy)/ry)**2 > config.central_exclusion_scale**2
    return dict(gradient_magnitude=magnitude,strong_gradient_pixels=strong,selected_edges=selected,
        eligible=eligible,edge_x=x,edge_y=y,normal_x=nx,normal_y=ny), dict(
        edge_threshold=threshold,strong_gradient_pixels=int(strong.sum()),selected_edge_pixels=len(x),
        eligible_center_pixels=int(eligible.sum()),physical_disk_radii=(rx,ry),
        central_exclusion_radii=(rx*config.central_exclusion_scale,ry*config.central_exclusion_scale))


def arc_center_votes(corrected_pattern, geometry, central_disk_center_px, *, config=None):
    """Known-radius normal-direction votes with both signs and bilinear accumulation."""
    config=_config(config)
    arrays,diagnostics=arc_edge_support(corrected_pattern,geometry,central_disk_center_px,config=config)
    image=np.asarray(corrected_pattern)
    rx,ry=geometry.radius_x_px,geometry.radius_y_px
    x,y,nx,ny=(arrays[k] for k in ('edge_x','edge_y','normal_x','normal_y'))
    votes = np.zeros(image.shape)
    offsets = np.linspace(-config.radius_tolerance_px, config.radius_tolerance_px, config.radius_samples)
    for offset in offsets:
        a, b = rx+offset, ry+offset
        denominator = np.hypot(a*nx, b*ny)
        dx, dy = a*a*nx/denominator, b*b*ny/denominator
        for sign in (-1, 1):
            vx, vy = x+sign*dx, y+sign*dy
            inside = (vx>=0)&(vx<=image.shape[1]-1)&(vy>=0)&(vy<=image.shape[0]-1)
            vx, vy = vx[inside], vy[inside]
            ix, iy = np.floor(vx).astype(int), np.floor(vy).astype(int)
            fx, fy = vx-ix, vy-iy
            for ox, oy, weight in ((0,0,(1-fx)*(1-fy)),(1,0,fx*(1-fy)),(0,1,(1-fx)*fy),(1,1,fx*fy)):
                ok = (ix+ox<image.shape[1])&(iy+oy<image.shape[0])
                np.add.at(votes, (iy[ok]+oy,ix[ok]+ox),weight[ok]/config.radius_samples)
    votes = gaussian_filter(votes, config.vote_smoothing_sigma_px, mode='constant')
    eligible=arrays['eligible']
    vote_threshold = max(config.minimum_vote_strength, config.relative_vote_threshold*float(votes[eligible].max(initial=0)))
    proposal_peaks=eligible & (votes>=vote_threshold) & (votes==maximum_filter(votes,size=3,mode='constant'))
    diagnostics.update(proposal_count=int(proposal_peaks.sum()),vote_threshold=vote_threshold,radius_offsets_px=offsets.tolist())
    arrays['proposal_peaks']=proposal_peaks
    return votes, arrays, diagnostics


def _coherent_fragments(points, normals, config):
    """8-neighbor fragments, cutting abrupt sign-invariant normal changes."""
    if not len(points):
        return ()
    pairs = cKDTree(points).query_pairs(np.sqrt(2)+1e-8, output_type='ndarray')
    if len(pairs):
        smooth = abs(np.sum(normals[pairs[:,0]]*normals[pairs[:,1]],axis=1)) >= np.cos(
            np.deg2rad(config.maximum_neighbor_normal_change_deg))
        pairs = pairs[smooth]
    graph = coo_matrix((np.ones(len(pairs)),(pairs[:,0],pairs[:,1])),shape=(len(points),len(points)))
    count, labels = connected_components(graph,directed=False)
    return tuple(np.flatnonzero(labels==i) for i in range(count))


def _angular_support(displacement, axes):
    """Union of projected pixel footprints; missing angular gaps contribute nothing."""
    phase=np.mod(np.arctan2(displacement[:,1]/axes[1],displacement[:,0]/axes[0]),2*np.pi)
    speed=np.hypot(axes[0]*np.sin(phase),axes[1]*np.cos(phase))
    # Projection of a square pixel onto the local unit tangent, in phase radians.
    half_width=.5*(abs(axes[0]*np.sin(phase))+abs(axes[1]*np.cos(phase)))/speed**2
    intervals=[]
    for angle,width in zip(phase,half_width):
        lo,hi=angle-width,angle+width
        for shift in (-2*np.pi,0,2*np.pi):
            left,right=max(0.,lo+shift),min(2*np.pi,hi+shift)
            if left<right: intervals.append((left,right))
    merged=[]
    for lo,hi in sorted(intervals):
        if merged and lo<=merged[-1][1]:
            merged[-1]=(merged[-1][0],max(hi,merged[-1][1]))
        else: merged.append((lo,hi))
    span=sum(hi-lo for lo,hi in merged)
    length=0.
    for lo,hi in merged:
        t=np.linspace(lo,hi,16)
        velocity=np.hypot(axes[0]*np.sin(t),axes[1]*np.cos(t))
        length+=float(np.sum((velocity[:-1]+velocity[1:])*.5*np.diff(t)))
    return span/(2*np.pi),length


def _validated_center(center, arrays, geometry, config):
    """Gather all nominal-radius-compatible ridges around one center hypothesis.

    Connectivity describes independent fragments; it is never an acceptance gate.
    The radial-distance error is exactly abs(distance-R) for circles. For an
    ellipse it is distance to the nominal ellipse along the same radial ray.
    No fitted radius offset or extra support tolerance broadens the physical prior.
    """
    points=np.column_stack((arrays['edge_x'],arrays['edge_y']))
    normals=np.column_stack((arrays['normal_x'],arrays['normal_y']))
    displacement=points-center
    axes=np.array([geometry.radius_x_px,geometry.radius_y_px])
    radial=np.linalg.norm(displacement/axes,axis=1)
    distance=np.linalg.norm(displacement,axis=1)
    boundary_distance=np.divide(distance,radial,out=np.full_like(distance,np.inf),where=radial>0)
    radius_error=abs(distance-boundary_distance)
    expected=displacement/axes**2
    expected/=np.maximum(np.linalg.norm(expected,axis=1)[:,None],np.finfo(float).tiny)
    alignment=np.clip(abs(np.sum(normals*expected,axis=1)),0,1)
    normal_angles=np.rad2deg(np.arccos(alignment))
    eligible=(radius_error<=config.radius_tolerance_px)&(normal_angles<=config.maximum_normal_residual_deg)
    ids=np.flatnonzero(eligible)
    if len(ids)<config.minimum_support_pixels:
        return None
    coverage,length=_angular_support(displacement[ids],axes)
    if coverage<config.minimum_angular_coverage or length<config.minimum_arc_length_px:
        return None
    fragments=len(_coherent_fragments(points[ids],normals[ids],config))
    shift=normals[ids]*axes**2/np.linalg.norm(normals[ids]*axes,axis=1)[:,None]
    sign=np.where(np.sum(normals[ids]*displacement[ids],axis=1)>=0,1.,-1.)
    implied=points[ids]-sign[:,None]*shift
    # Relative to this hypothesis, not an independently fitted fragment center.
    rms=float(np.sqrt(np.mean(np.sum((implied-center)**2,axis=1))))
    radius_rms=float(np.sqrt(np.mean(radius_error[ids]**2)))
    angle=float(normal_angles[ids].mean())
    confidence=float(len(ids)*coverage*(1+.1*np.log1p(fragments-1))/
        (1+rms/config.center_vote_rms_scale_px+radius_rms+angle/config.maximum_normal_residual_deg))
    radii=np.sqrt(np.prod(axes))*radial[ids]
    return dict(supporting_edge_pixels=len(ids),angular_coverage=float(coverage),fragment_count=fragments,
        radius_residual_rms_px=radius_rms,arc_length_px=length,center_vote_rms_px=rms,
        mean_radius_px=float(radii.mean()),radius_std_px=float(radii.std()),
        mean_normal_residual_deg=angle,confidence=confidence,
        supporting_pixels=tuple(tuple(int(v) for v in point) for point in points[ids]),radius_offset_px=0.)


