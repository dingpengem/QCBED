"""A8 explicit beam-local centrosymmetry measurement; no catalogue-order selection.

DM N:594–623: y-major unit-step upper-exclusive search, half-open crops,
first minimum of sum(abs(Cut-Cut[::-1,::-1])), and the legacy origin relation.
Fractional DM subarea rounding is unverified: reject nonintegral crop bounds.
"""
from dataclasses import dataclass
import numpy as np
from .models import DiskGeometry, _pair
from .noise import _pattern
from ..detector.models import _real, _integer, _text


@dataclass(frozen=True)
class SymmetryReference:
    beam_id: str
    beam_center_px: tuple[float,float]
    feature_offset_px: tuple[float,float]  # mandatory, relative to the beam
    lattice_indices: tuple[int,int] | None = None

    def __post_init__(self):
        _text(self.beam_id,'beam_id')
        for name in ('beam_center_px','feature_offset_px'):
            object.__setattr__(self,name,_pair(getattr(self,name),name))
        if self.lattice_indices is not None:
            values=tuple(_integer(v,'lattice index') for v in self.lattice_indices)
            if len(values)!=2: raise ValueError('Expected two lattice indices, not HKL')
            object.__setattr__(self,'lattice_indices',values)


@dataclass(frozen=True)
class SymmetrySearchConfig:
    search_half_size_px: tuple[float,float] | None = None
    crop_half_size_px: tuple[float,float] | None = None

    def __post_init__(self):
        for name in ('search_half_size_px','crop_half_size_px'):
            value=getattr(self,name)
            if value is not None:
                value=_pair(value,name)
                if min(value)<=0: raise ValueError(f'{name} must be positive')
                object.__setattr__(self,name,value)


@dataclass(frozen=True)
class SymmetryMeasurement:
    reference: SymmetryReference
    geometry: DiskGeometry
    central_disk_center_px: tuple[float,float]
    feature_seed_px: tuple[float,float]
    symmetry_center_px: tuple[float,float]  # DM trial coordinate, not crop midpoint
    crop_pixel_midpoint_px: tuple[float,float]
    zone_axis_origin_px: tuple[float,float]
    search_half_size_px: tuple[float,float]
    crop_half_size_px: tuple[float,float]
    crop_bounds_xyxy: tuple[int,int,int,int]
    mismatch: float
    trial_count: int
    tied_minimum_count: int
    boundary_limited: bool
    crop_has_contrast: bool
    trials: tuple[tuple[float,float,float], ...]
    reference_kind: str = "diffracted"


@dataclass(frozen=True)
class OriginConsistency:
    measurements: tuple[SymmetryMeasurement,...]
    tolerance_px: float
    maximum_pairwise_distance_px: float
    consistent: bool
    combined_origin_px: tuple[float,float] | None


def rotate180(crop):
    """Exact row/column reversal for any even/odd rectangular crop."""
    return _pattern(crop)[::-1,::-1].copy()


def legacy_zone_axis_origin(central,reference,symmetry):
    central=np.array(_pair(central,'central'))
    reference=np.array(_pair(reference,'reference'))
    symmetry=np.array(_pair(symmetry,'symmetry'))
    return tuple(float(v) for v in symmetry+.5*(central-reference))


def select_symmetry_reference(catalogue, *, feature_offset_px, beam_id=None, lattice_indices=None):
    """Resolve an explicitly selected direct or observed diffracted A7 beam by ID/indices.

    The geometric center is the catalogue lattice prediction, NOT the local
    observed arc offset. Both selectors, if supplied, must identify the same beam.
    """
    if beam_id is None and lattice_indices is None:
        raise ValueError('The symmetry reference beam and local feature must be specified')
    if lattice_indices is not None:
        lattice_indices=tuple(_integer(v,'lattice index') for v in lattice_indices)
        if len(lattice_indices)!=2: raise ValueError('Expected two lattice indices')
    records=[b for b in catalogue.sites if (beam_id is None or b.beam_id==beam_id)
             and (lattice_indices is None or b.integer_indices==lattice_indices)]
    if len(records)!=1: raise ValueError('Symmetry reference selection is missing, ambiguous, or conflicting')
    beam=records[0]
    if beam.status!='observed' and not (beam.status=='reference' and beam.integer_indices==(0,0)):
        raise ValueError('Symmetry reference must be a reliable direct or observed diffracted beam')
    return SymmetryReference(beam.beam_id,beam.predicted_center_px,feature_offset_px,beam.integer_indices)


def measure_symmetry_origin(pattern, geometry, central_disk_center_px, reference=None, *, config=None):
    if reference is None: raise ValueError('The symmetry reference beam and local feature must be specified')
    if not isinstance(reference,SymmetryReference): raise TypeError('Expected SymmetryReference')
    if not isinstance(geometry,DiskGeometry): raise TypeError('Expected DiskGeometry')
    config=SymmetrySearchConfig() if config is None else config
    if not isinstance(config,SymmetrySearchConfig): raise TypeError('Expected SymmetrySearchConfig')
    image=_pattern(pattern);central=_pair(central_disk_center_px,'central_disk_center_px')
    seed=np.array(reference.beam_center_px)+reference.feature_offset_px
    axes=np.array([geometry.radius_x_px,geometry.radius_y_px])
    if np.sum((np.array(reference.feature_offset_px)/axes)**2)>1:
        raise ValueError('Feature seed must lie inside the selected physical disk')
    search=config.search_half_size_px or tuple(axes/4)  # width/8
    half=config.crop_half_size_px or tuple(axes/2)  # width/4
    nx=int(np.ceil(2*search[0]));ny=int(np.ceil(2*search[1]))
    if nx*ny>100000: raise ValueError('Symmetry search exceeds 100000 trials')
    trials=[];best=None;ties=0
    for y in np.arange(seed[1]-search[1],seed[1]+search[1],1.):
        for x in np.arange(seed[0]-search[0],seed[0]+search[0],1.):
            bounds=(x-half[0],y-half[1],x+half[0],y+half[1])
            if any(not float(v).is_integer() for v in bounds):
                raise ValueError('Fractional DM crop bounds are unverified; supply explicit search/crop sizes yielding integer bounds')
            x0,y0,x1,y1=map(int,bounds)
            if x0<0 or y0<0 or x1>image.shape[1] or y1>image.shape[0] or x1<=x0 or y1<=y0:
                raise ValueError('Entire symmetry search/crop footprint must fit the frame; no clipping or padding')
            crop=image[y0:y1,x0:x1]
            with np.errstate(over='raise',invalid='raise'):
                mismatch=float(np.sum(np.abs(crop-crop[::-1,::-1]),dtype=np.float64))
            trials.append((float(x),float(y),mismatch))
            if best is None or mismatch<best[0]: best=(mismatch,float(x),float(y),(x0,y0,x1,y1));ties=1
            elif mismatch==best[0]: ties+=1
    mismatch,x,y,bounds=best;x0,y0,x1,y1=bounds
    boundary=(x==seed[0]-search[0] or y==seed[1]-search[1] or
              x==trials[nx-1][0] or y==trials[-1][1])
    crop=image[y0:y1,x0:x1]
    return SymmetryMeasurement(reference,geometry,central,tuple(seed),(x,y),
        ((x0+x1-1)/2,(y0+y1-1)/2),legacy_zone_axis_origin(central,reference.beam_center_px,(x,y)),
        tuple(search),tuple(half),bounds,mismatch,len(trials),ties,bool(boundary),bool(np.ptp(crop)>0),tuple(trials),
        reference_kind="direct" if reference.lattice_indices==(0,0) or reference.beam_center_px==central else "diffracted")


def combine_origin_measurements(measurements, *, tolerance_px, combine=False):
    """Explicit Euclidean pairwise tolerance; equal-weight mean only on request."""
    measurements=tuple(measurements)
    if not measurements: raise ValueError('At least one symmetry reference measurement is required')
    if not all(isinstance(m,SymmetryMeasurement) for m in measurements): raise TypeError('Expected symmetry measurements')
    tolerance=_real(tolerance_px,'tolerance_px')
    if tolerance<0: raise ValueError('tolerance_px must be nonnegative')
    if type(combine) is not bool: raise TypeError('combine must be bool')
    positions=np.array([m.zone_axis_origin_px for m in measurements])
    distance=float(np.max(np.linalg.norm(positions[:,None]-positions[None,:],axis=2)))
    consistent=distance<=tolerance
    origin=tuple(float(v) for v in positions.mean(axis=0)) if consistent and combine else None
    return OriginConsistency(measurements,tolerance,distance,consistent,origin)


def write_symmetry_diagnostics(output, pattern, measurement):
    """Fresh directory; unchanged intensities displayed with explicit x/y extents."""
    import json
    from dataclasses import asdict
    from pathlib import Path
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse,Rectangle
    image=_pattern(pattern);m=measurement;output=Path(output);output.mkdir(parents=True,exist_ok=False)
    x0,y0,x1,y1=m.crop_bounds_xyxy;crop=image[y0:y1,x0:x1]
    rotated=rotate180(crop);difference=np.abs(crop-rotated)
    fig,axes=plt.subplots(2,3,figsize=(15,9))
    bx,by=m.reference.beam_center_px;rx,ry=m.geometry.radius_x_px,m.geometry.radius_y_px
    for ax in axes[0,:2]:
        ax.imshow(image,origin='upper',cmap='gray')
        ax.add_patch(Ellipse((bx,by),2*rx,2*ry,fill=False,color='cyan'))
        ax.scatter(bx,by,c='cyan',marker='+',label='Beam geometry')
        ax.scatter(*m.feature_seed_px,c='orange',marker='x',label='Feature seed')
        ax.scatter(*m.symmetry_center_px,c='lime',marker='+',label='DM trial symmetry center')
        ax.scatter(*m.zone_axis_origin_px,c='magenta',marker='*',label='Derived zone origin')
        ax.add_patch(Rectangle((x0-.5,y0-.5),x1-x0,y1-y0,fill=False,color='lime'))
        sx,sy=m.search_half_size_px
        ax.add_patch(Rectangle((m.feature_seed_px[0]-sx,m.feature_seed_px[1]-sy),2*sx,2*sy,fill=False,color='orange'))
    axes[0,0].legend(fontsize=7);axes[0,0].set_title('Detector / derived origin')
    axes[0,1].set(xlim=(max(-.5,bx-rx-2),min(image.shape[1]-.5,bx+rx+2)),
                  ylim=(min(image.shape[0]-.5,by+ry+2),max(-.5,by-ry-2)),title='Selected disk and local feature')
    for ax,data,title in zip(axes[1],(crop,rotated,difference),('Optimal crop','Exact 180° rotation','Absolute difference')):
        ax.imshow(data,origin='upper',cmap='gray',extent=(x0-.5,x1-.5,y1-.5,y0-.5));ax.set_title(title)
    axes[0,2].axis('off');axes[0,2].text(0,1,f'{m.reference.beam_id} ({m.reference_kind})\nFeature offset: {m.reference.feature_offset_px}\nSymmetry: {m.symmetry_center_px}\nCrop midpoint: {m.crop_pixel_midpoint_px}\nOrigin: {m.zone_axis_origin_px}\nMismatch: {m.mismatch}\nTies: {m.tied_minimum_count}\nBoundary: {m.boundary_limited}\nContrast: {m.crop_has_contrast}',va='top')
    fig.tight_layout();fig.savefig(output/'symmetry_diagnostic.png',dpi=150);plt.close(fig)
    (output/'measurement.json').write_text(json.dumps(asdict(m),indent=2,allow_nan=False)+'\n')
    return output/'symmetry_diagnostic.png'
