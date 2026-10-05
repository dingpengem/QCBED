"""Canonical A5 localization from physical disk-boundary geometry."""
from dataclasses import dataclass
import numpy as np
from ..detector.models import _real
from .models import DiskGeometry
from .noise import _pattern
from .disk_detection import _center
from .arc_geometry import LatticeDetectionConfig, LatticeCandidate, arc_edge_support, _validated_center


@dataclass(frozen=True)
class CentralArcConfig:
    search_radius_px: float = 8.0
    search_step_px: float = .5

    def __post_init__(self):
        for name in ('search_radius_px','search_step_px'):
            value=_real(getattr(self,name),name)
            if value<=0: raise ValueError(f'{name} must be positive')
            object.__setattr__(self,name,value)
        if self.search_step_px>self.search_radius_px or self.search_radius_px/self.search_step_px>100:
            raise ValueError('Require step <= radius and at most 100 steps per radius')


@dataclass(frozen=True)
class CentralDiskCenterComparison:
    intensity_center_px: tuple[float,float]
    arc_center_px: tuple[float,float] | None
    difference_vector_px: tuple[float,float] | None  # arc minus intensity, x/y pixels
    difference_magnitude_px: float | None
    arc_support: LatticeCandidate | None
    geometry: DiskGeometry
    search_config: CentralArcConfig
    arc_config: LatticeDetectionConfig
    status: str
    intensity_seed_source: str
    geometry_provenance: str
    canonical_center_method: str = 'known_radius_partial_arcs'


def _central_edges(image,geometry,seed,search,arc):
    y,x=np.indices(image.shape)
    # Local threshold distribution; derivatives remain identical to A6's shared helper.
    reach=np.array([geometry.radius_x_px,geometry.radius_y_px])+search.search_radius_px+arc.radius_tolerance_px+1
    mask=(abs(x-seed[0])<=reach[0])&(abs(y-seed[1])<=reach[1])
    return arc_edge_support(image,geometry,seed,config=arc,exclude_central=False,edge_mask=mask)


@dataclass(frozen=True)
class _ArcSearchResult:
    search_prior_px: tuple[float, float]
    center_px: tuple[float, float] | None
    support: LatticeCandidate | None
    geometry: DiskGeometry
    search_config: CentralArcConfig
    arc_config: LatticeDetectionConfig
    status: str


@dataclass(frozen=True)
class CentralDiskResult:
    geometry: DiskGeometry
    central_disk_center_px: tuple[float, float]
    search_prior_px: tuple[float, float]
    arc_support: LatticeCandidate
    search_config: CentralArcConfig
    arc_config: LatticeDetectionConfig
    status: str
    method: str = 'known_radius_partial_arcs'


def locate_central_disk(corrected_pattern, geometry, *, initial_center_px=None,
                        search_config=None, arc_config=None, subpixel_config=None) -> CentralDiskResult:
    """Locate the transmitted disk from combined known-radius boundary evidence.

    Default prior is the pixel-frame midpoint ((nx-1)/2, (ny-1)/2). An explicit
    initial_center_px is only a search prior, never an accepted preset center.
    The physical semi-axes are fixed. No intensity localization or fallback runs.
    Failure raises ValueError; a boundary optimum carries search_boundary status.
    """
    image = _pattern(corrected_pattern)
    seed = ((image.shape[1]-1)/2, (image.shape[0]-1)/2) if initial_center_px is None else initial_center_px
    result = _estimate(image, geometry, seed, search_config, arc_config)
    if result.center_px is None:
        raise ValueError('No validated central disk boundary in the configured geometric search region')
    support=result.support
    if subpixel_config is not None:
        from .subpixel_arc import refine_validated_candidate
        arrays,_=_central_edges(image,geometry,result.search_prior_px,result.search_config,result.arc_config)
        support=refine_validated_candidate(support,arrays,geometry,result.arc_config,config=subpixel_config)
    return CentralDiskResult(geometry, support.center_px, result.search_prior_px,
                            support, result.search_config, result.arc_config, result.status)


def compare_central_disk_centers(corrected_pattern, geometry, *, intensity_center_px=None,
        search_config=None, arc_config=None, geometry_provenance='Supplied physical DiskGeometry'):
    """Compatibility-only comparison seeded by the historical intensity algorithm.

    Never called by canonical localization or the preparation workflow.
    """
    from .legacy_disk_localization import locate_central_disk_intensity
    seed = (locate_central_disk_intensity(corrected_pattern, geometry).central_disk_center_px
            if intensity_center_px is None else intensity_center_px)
    if not isinstance(geometry_provenance,str) or not geometry_provenance.strip():
        raise ValueError('geometry_provenance must be nonempty text')
    result = _estimate(corrected_pattern, geometry, seed, search_config, arc_config)
    shift = None if result.center_px is None else tuple(float(v) for v in np.array(result.center_px)-result.search_prior_px)
    distance = None if shift is None else float(np.linalg.norm(shift))
    source = 'legacy_integrated_intensity' if intensity_center_px is None else 'supplied_intensity_seed'
    return CentralDiskCenterComparison(result.search_prior_px, result.center_px, shift, distance,
        result.support, geometry, result.search_config, result.arc_config, result.status, source, geometry_provenance)


def _estimate(corrected_pattern, geometry, seed, search_config=None, arc_config=None):
    image=_pattern(corrected_pattern)
    if not isinstance(geometry,DiskGeometry): raise TypeError('geometry must be DiskGeometry')
    search=CentralArcConfig() if search_config is None else search_config
    arc=LatticeDetectionConfig() if arc_config is None else arc_config
    if not isinstance(search,CentralArcConfig) or not isinstance(arc,LatticeDetectionConfig):
        raise TypeError('Invalid central arc configuration')
    if search.search_radius_px>min(geometry.radius_x_px,geometry.radius_y_px)/2:
        raise ValueError('Central search radius must not exceed half the smaller physical radius')
    seed=tuple(float(v) for v in _center(seed,image.shape))
    arrays,_=_central_edges(image,geometry,seed,search,arc)
    lower=np.array([geometry.radius_x_px,geometry.radius_y_px]);upper=np.array(image.shape[::-1])-1-lower
    n=int(np.floor(search.search_radius_px/search.search_step_px));best=None
    for j in range(-n,n+1):
        for i in range(-n,n+1):
            shift=np.array([i,j])*search.search_step_px
            distance=float(np.linalg.norm(shift))
            if distance>search.search_radius_px+1e-12: continue
            center=np.array(seed)+shift
            if np.any(center<lower) or np.any(center>upper): continue
            metrics=_validated_center(center,arrays,geometry,arc)
            if metrics is None: continue
            key=(metrics['confidence'],-distance,-center[0],-center[1])
            if best is None or key>best[0]: best=(key,tuple(float(v) for v in center),tuple(float(v) for v in shift),distance,metrics)
    if best is None:
        return _ArcSearchResult(seed,None,None,geometry,search,arc,'no_validated_arc')
    _,center,shift,distance,metrics=best
    candidate=LatticeCandidate('central_boundary',center,metrics['confidence'],**metrics)
    status='search_boundary' if distance>=search.search_radius_px-search.search_step_px*.5 else 'estimated'
    return _ArcSearchResult(seed,center,candidate,geometry,search,arc,status)


def write_central_center_comparison(output, corrected_pattern, result):
    """Fresh diagnostic directory with original intensities, separate centers and JSON."""
    import json
    from pathlib import Path
    from dataclasses import asdict
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse
    image=_pattern(corrected_pattern);output=Path(output);output.mkdir(parents=True,exist_ok=False)
    arrays,_=_central_edges(image,result.geometry,result.intensity_center_px,result.search_config,result.arc_config)
    fig,axes=plt.subplots(1,2,figsize=(13,6));low,high=np.percentile(image,[1,99.5])
    for ax in axes:
        ax.imshow(image,cmap='gray',origin='upper',vmin=low,vmax=high)
        ax.scatter(*result.intensity_center_px,c='cyan',marker='x',label='Legacy intensity (comparison only)')
        centers=[(result.intensity_center_px,'cyan')]
        retained=set() if result.arc_support is None else set(result.arc_support.supporting_pixels)
        ridges=list(zip(arrays['edge_x'],arrays['edge_y']))
        rejected=np.array([p for p in ridges if p not in retained])
        if len(rejected): ax.scatter(rejected[:,0],rejected[:,1],c='orange',s=2,alpha=.4,label='Unretained ridges')
        if result.arc_center_px is not None:
            ax.scatter(*result.arc_center_px,c='magenta',marker='+',label='Geometric estimate')
            centers.append((result.arc_center_px,'magenta'))
            points=np.array(result.arc_support.supporting_pixels)
            ax.scatter(points[:,0],points[:,1],c='lime',s=5,label='Boundary support')
        for center,color in centers:
            ax.add_patch(Ellipse(center,width=2*result.geometry.radius_x_px,height=2*result.geometry.radius_y_px,
                                 fill=False,edgecolor=color))
    axes[0].set_xlim(-.5,image.shape[1]-.5);axes[0].set_ylim(image.shape[0]-.5,-.5)
    cx,cy=result.intensity_center_px;rx,ry=result.geometry.radius_x_px,result.geometry.radius_y_px
    pad=result.search_config.search_radius_px+5
    axes[1].set_xlim(max(-.5,cx-rx-pad),min(image.shape[1]-.5,cx+rx+pad))
    axes[1].set_ylim(min(image.shape[0]-.5,cy+ry+pad),max(-.5,cy-ry-pad))
    axes[1].legend(fontsize=7);fig.suptitle(f'Central center comparison: {result.status}; delta={result.difference_vector_px} px')
    path=output/'central_disk_center_comparison.png';fig.savefig(path,dpi=150);plt.close(fig)
    (output/'central_disk_center_comparison.json').write_text(json.dumps(asdict(result),indent=2,allow_nan=False)+'\n')
    return path


def write_central_arc_localization(output, corrected_pattern, result):
    """Write canonical boundary evidence to an existing diagnostic directory."""
    import json
    from pathlib import Path
    from dataclasses import asdict
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse
    output=Path(output)
    path=output/'central_disk_arc_localization.png'
    metadata=output/'central_disk_arc_localization.json'
    if path.exists() or metadata.exists():
        raise FileExistsError('Central localization diagnostics already exist')
    image=_pattern(corrected_pattern)
    arrays,_=_central_edges(image,result.geometry,result.search_prior_px,result.search_config,result.arc_config)
    retained=set(result.arc_support.supporting_pixels)
    rejected=np.array([p for p in zip(arrays['edge_x'],arrays['edge_y']) if p not in retained])
    points=np.array(result.arc_support.supporting_pixels)
    fig,ax=plt.subplots(figsize=(7,7))
    low,high=np.percentile(image,[1,99.5])
    ax.imshow(image,cmap='gray',origin='upper',vmin=low,vmax=high)
    if len(rejected): ax.scatter(rejected[:,0],rejected[:,1],c='orange',s=3,alpha=.5,label='Non-supporting ridges')
    ax.scatter(points[:,0],points[:,1],c='lime',s=6,label='Combined boundary support')
    ax.scatter(*result.central_disk_center_px,c='magenta',marker='+',label='Canonical geometric center')
    ax.add_patch(Ellipse(result.central_disk_center_px,2*result.geometry.radius_x_px,
                         2*result.geometry.radius_y_px,fill=False,edgecolor='magenta'))
    ax.set(xlabel='x (pixel)',ylabel='y (pixel)',title=f'Known-radius central localization: {result.central_disk_center_px}')
    ax.legend(fontsize=8)
    fig.savefig(path,dpi=150);plt.close(fig)
    metadata.write_text(json.dumps(asdict(result),indent=2,allow_nan=False)+'\n')
    return str(path)
