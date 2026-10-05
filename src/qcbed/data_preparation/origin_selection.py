"""Frontend-independent A8 coarse selection/provenance; numerical science in origin.py."""
from dataclasses import dataclass,asdict
from pathlib import Path
import hashlib
import json
import csv
import numpy as np
from .models import DiskGeometry,_pair
from .origin import SymmetryReference,SymmetrySearchConfig,measure_symmetry_origin,write_symmetry_diagnostics
from .diagnostics import load_pattern
from ..detector.models import _real


@dataclass(frozen=True)
class SelectableBeam:
    beam_id: str
    lattice_indices: tuple[int,int]
    center_px: tuple[float,float]
    status: str


def load_selection_catalogue(path):
    with Path(path).open() as stream: rows=list(csv.DictReader(stream))
    beams=tuple(SelectableBeam(r['beam_id'],(int(r['n1']),int(r['n2'])),
        _pair((float(r['predicted_x']),float(r['predicted_y'])),'beam center'),r['status']) for r in rows)
    if len({b.beam_id for b in beams})!=len(beams): raise ValueError('Duplicate catalogue beam IDs')
    return beams


def nearest_beam(beams,click_px,*,selection_radius_px,allow_unreliable=False):
    click=np.array(_pair(click_px,'click'))
    radius=_real(selection_radius_px,'selection_radius_px')
    if radius<=0: raise ValueError('Selection radius must be positive')
    if type(allow_unreliable) is not bool: raise TypeError('allow_unreliable must be bool')
    allowed={'observed'} | ({'weak','missing','edge_clipped'} if allow_unreliable else set())
    candidates=[(float(np.linalg.norm(np.array(b.center_px)-click)),b.beam_id,b) for b in beams
                if b.status in allowed or (b.status=='reference' and b.lattice_indices==(0,0))]
    if not candidates: raise ValueError('No selectable symmetry-reference beams')
    distance,_,beam=min(candidates,key=lambda item:(item[0],item[1]))
    if distance>radius: raise ValueError(f'No selectable beam within {radius:g} px; nearest is {distance:g} px away')
    return beam,distance


def roi_bounds(bounds):
    values=tuple(_real(v,'ROI bound') for v in bounds)
    if len(values)!=4: raise ValueError('ROI must be (x0,y0,x1,y1)')
    x0,y0,x1,y1=values
    if x0>=x1 or y0>=y1: raise ValueError('ROI must have positive width and height')
    return values


def detector_to_local_roi(bounds,beam_center):
    return tuple(np.array(roi_bounds(bounds))-np.tile(_pair(beam_center,'beam center'),2))


def local_to_detector_roi(bounds,beam_center):
    return tuple(np.array(roi_bounds(bounds))+np.tile(_pair(beam_center,'beam center'),2))


def snap_roi_to_grid(bounds):
    """Explicit user action only: nearest integer coordinate, ties to even."""
    return roi_bounds(tuple(float(v) for v in np.rint(roi_bounds(bounds))))


def prepare_roi(bounds,beam,geometry,image_shape,*,search_half_size_px=None,crop_half_size_px=None):
    """Validate the entire coarse ROI and all sampled crops before numerical work."""
    x0,y0,x1,y1=roi_bounds(bounds)
    center=np.array([(x0+x1)/2,(y0+y1)/2]);axes=np.array([geometry.radius_x_px,geometry.radius_y_px])
    if x0<0 or y0<0 or x1>image_shape[1]-1 or y1>image_shape[0]-1:
        raise ValueError('ROI must lie inside the detector pixel-center frame')
    corners=np.array([(x0,y0),(x0,y1),(x1,y0),(x1,y1)])
    if np.any(np.sum(((corners-beam.center_px)/axes)**2,axis=1)>1):
        raise ValueError('ROI corners must lie inside the selected physical beam ellipse')
    config=SymmetrySearchConfig(search_half_size_px or ((x1-x0)/2,(y1-y0)/2),crop_half_size_px)
    search=np.array(config.search_half_size_px);crop=np.array(config.crop_half_size_px or tuple(axes/2))
    lower=center-search;upper=center+search
    if np.any(lower<[x0,y0]) or np.any(upper>[x1,y1]):
        raise ValueError('Advanced search must stay within the selected coarse ROI')
    counts=np.ceil(2*search).astype(int)
    if np.prod(counts.astype(object))>100000: raise ValueError('Symmetry search exceeds 100000 trials')
    last=lower+counts-1
    if np.any((lower-crop)!=np.rint(lower-crop)) or np.any((lower+crop)!=np.rint(lower+crop)):
        raise ValueError('Fractional crop bounds: explicitly snap ROI to grid or adjust advanced crop/search sizes')
    if np.any(lower-crop<0) or np.any(last+crop>np.array(image_shape[::-1])):
        raise ValueError('Insufficient detector room for all trial crops; reselect a smaller search ROI')
    reference=SymmetryReference(beam.beam_id,beam.center_px,tuple(center-beam.center_px),beam.lattice_indices)
    return reference,config


def _hash(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def make_selection(pattern_path,catalogue_path,beam,click_px,roi_detector,geometry,central_center,*,
                   selection_radius_px,allow_unreliable=False,raw_roi_detector=None,
                   search_half_size_px=None,crop_half_size_px=None):
    """Create an explicit JSON-safe record from an accepted preview, without running A8."""
    pattern_path=Path(pattern_path).resolve();catalogue_path=Path(catalogue_path).resolve()
    image=load_pattern(pattern_path)
    selected,distance=nearest_beam(load_selection_catalogue(catalogue_path),click_px,
        selection_radius_px=selection_radius_px,allow_unreliable=allow_unreliable)
    if selected!=beam: raise ValueError('Selected beam no longer matches the catalogue/click')
    reference,config=prepare_roi(roi_detector,beam,geometry,image.shape,
        search_half_size_px=search_half_size_px,crop_half_size_px=crop_half_size_px)
    return dict(reference_kind="direct" if beam.lattice_indices==(0,0) else "diffracted",schema_version=1,pattern_path=str(pattern_path),pattern_sha256=_hash(pattern_path),
        catalogue_path=str(catalogue_path),catalogue_sha256=_hash(catalogue_path),
        beam=asdict(beam),mouse_click_px=list(_pair(click_px,'click')),click_distance_px=distance,
        selection_radius_px=selection_radius_px,allow_unreliable=allow_unreliable,
        raw_roi_detector=None if raw_roi_detector is None else list(roi_bounds(raw_roi_detector)),
        roi_detector=list(roi_bounds(roi_detector)),roi_local=list(detector_to_local_roi(roi_detector,beam.center_px)),
        roi_grid_policy='Explicit nearest-integer conversion only when requested; ties to even',
        geometry=asdict(geometry),central_disk_center_px=list(_pair(central_center,'central center')),
        search_config=asdict(config),resolved_crop_half_size_px=list(config.crop_half_size_px or (geometry.radius_x_px/2,geometry.radius_y_px/2)))


def replay_selection(selection_path,output):
    """Rerun the stored selection without a UI, verifying input identities and geometry."""
    def invalid(value): raise ValueError('Nonfinite JSON value: '+value)
    record=json.loads(Path(selection_path).read_text(),parse_constant=invalid)
    return execute_selection(record,output)


def execute_selection(record,output):
    if record.get('schema_version')!=1: raise ValueError('Unsupported centrosymmetry selection schema')
    for kind in ('pattern','catalogue'):
        path=Path(record[kind+'_path'])
        if not path.is_absolute(): raise ValueError('Selection input paths must be absolute')
        if _hash(path)!=record[kind+'_sha256']: raise ValueError(f'{kind} content changed since selection')
    beam,distance=nearest_beam(load_selection_catalogue(record['catalogue_path']),record['mouse_click_px'],
        selection_radius_px=record['selection_radius_px'],allow_unreliable=record['allow_unreliable'])
    stored=record['beam']
    if (beam.beam_id!=stored['beam_id'] or beam.center_px!=tuple(stored['center_px']) or
        beam.lattice_indices!=tuple(stored['lattice_indices']) or beam.status!=stored['status']):
        raise ValueError('Stored beam identity/geometry disagrees with catalogue')
    local=detector_to_local_roi(record['roi_detector'],beam.center_px)
    if tuple(record['roi_local'])!=local: raise ValueError('Detector/local ROI provenance disagrees')
    geometry=DiskGeometry(**record['geometry']);image=load_pattern(record['pattern_path'])
    config=SymmetrySearchConfig(**record['search_config'])
    reference,config=prepare_roi(record['roi_detector'],beam,geometry,image.shape,
        search_half_size_px=config.search_half_size_px,crop_half_size_px=config.crop_half_size_px)
    result=measure_symmetry_origin(image,geometry,record['central_disk_center_px'],reference,config=config)
    write_symmetry_diagnostics(output,image,result)
    _write_selection_context(output,image,record,result)
    saved=dict(record,reference_kind=result.reference_kind,result=asdict(result))
    Path(output,'a8_selection.json').write_text(json.dumps(saved,indent=2,allow_nan=False)+'\n')
    return result


def _write_selection_context(output,image,record,result):
    """Display-only user ROI provenance alongside final A8 geometry, also on replay."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Ellipse,Rectangle
    fig,axes=plt.subplots(1,2,figsize=(12,6))
    bx,by=record['beam']['center_px'];rx,ry=result.geometry.radius_x_px,result.geometry.radius_y_px
    x0,y0,x1,y1=record['roi_detector'];low,high=np.percentile(image,[1,99.5])
    for ax in axes:
        ax.imshow(image,origin='upper',cmap='gray',vmin=low,vmax=high)
        ax.add_patch(Ellipse((bx,by),2*rx,2*ry,fill=False,color='cyan'))
        ax.add_patch(Rectangle((x0,y0),x1-x0,y1-y0,fill=False,color='orange',label='Accepted coarse ROI'))
        for point,color,marker,label in ((result.reference.beam_center_px,'cyan','+','Beam geometry'),
            (result.feature_seed_px,'orange','x','Coarse feature seed'),(result.symmetry_center_px,'lime','+','Refined Symmetry Center'),
            (result.zone_axis_origin_px,'magenta','*','Derived origin'),(result.central_disk_center_px,'white','+','Central disk')):
            ax.scatter(*point,c=color,marker=marker,label=label)
    axes[0].legend(fontsize=7)
    axes[1].set(xlim=(max(-.5,bx-rx-2),min(image.shape[1]-.5,bx+rx+2)),
        ylim=(min(image.shape[0]-.5,by+ry+2),max(-.5,by-ry-2)))
    fig.suptitle(f"{record['beam']['beam_id']} {record['beam']['lattice_indices']} | ROI local={record['roi_local']}\n"
                 f"symmetry={result.symmetry_center_px}, mismatch={result.mismatch:g}, origin={result.zone_axis_origin_px}")
    fig.tight_layout();fig.savefig(Path(output)/'selection_context.png',dpi=150);plt.close(fig)
