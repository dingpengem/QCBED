"""Preparation document I/O, independent of Qt and refinement execution."""
from dataclasses import asdict, replace
from pathlib import Path
import hashlib
import json
import numpy as np
from qcbed.data_preparation.live import prepare_from_raw_cbed
from qcbed.data_preparation.dataset_diagnostics import prepared_dataset_manifest
from qcbed.config_serialization import _plain


def load_image(path):
    path=Path(path)
    if path.suffix.lower()=='.npy':image=np.load(path,allow_pickle=False)
    elif path.suffix.lower() in ('.tif','.tiff','.png'):
        from PIL import Image
        with Image.open(path) as source:image=np.asarray(source)
    else:image=np.loadtxt(path,ndmin=2)
    if image.ndim!=2 or image.dtype.kind not in 'iuf' or not image.size or not np.isfinite(image).all():
        raise ValueError('Expected a finite, single-channel two-dimensional numeric image; no color conversion or resampling is performed.')
    return np.array(image,dtype=float,copy=True)


def checksum(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare_files(pattern, selection, output, *, calibration_aperture=None, disk_profile=None, progress=lambda s:None):
    inputs={k:dict(path=str(Path(v).resolve()),sha256=checksum(v)) for k,v in
            (('pattern',pattern),('calibration_aperture',calibration_aperture),('disk_profile',disk_profile)) if v}
    result=prepare_from_raw_cbed(load_image(pattern),selection,
        calibration_aperture=load_image(calibration_aperture) if calibration_aperture else None,
        disk_profile=load_image(disk_profile) if disk_profile else None,progress=progress)
    # Check again before publishing: no provenance for inputs modified in flight.
    if any(checksum(v['path'])!=v['sha256'] for v in inputs.values()):raise ValueError('Preparation source changed during processing; rerun with stable inputs.')
    return publish_preparation(result,output,inputs)


def publish_preparation(result,output,inputs):
    """Exclusive cache output; all numerical arrays retain their values/order."""
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    data=result.dataset
    body=out/'experimental.dat'
    with body.open('x') as stream:
        for array in (data.intensities,data.variances,data.weights):np.savetxt(stream,array.T.ravel(),fmt='%.17g')
    np.save(out/'corrected_pattern.npy',result.corrected_pattern)
    np.savez(out/'uncertainty.npz',sigma=result.sigma,variance=result.variance,validity=result.validity)
    # Existing overview reader consumes numeric text without altering image values.
    np.savetxt(out/'overview.dat',result.corrected_pattern,fmt='%.17g')
    manifest=prepared_dataset_manifest(data)
    manifest.update(pattern_sha256=checksum(out/'overview.dat'),inputs=inputs,live_preparation=_plain(result.provenance))
    manifest_path=out/'patch_a.json';manifest_path.write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
    overview=out/'overview.dat'
    selection=dict(method='live symmetry-center search; awaiting visual confirmation',confirmed=False,
        zone_axis_origin_px=list(result.symmetry.zone_axis_origin_px),measurement=_plain(asdict(result.symmetry)))
    metadata=dict(format='qcbed-manual-mask-v1',source=dict(path=str(body),header='absent',sha256=checksum(body)),
        data_sha256=checksum(body),excluded_indices=[],zone_axis_selection=selection,
        overview=dict(pattern=str(overview),manifest=str(manifest_path),pattern_sha256=checksum(overview),manifest_sha256=checksum(manifest_path)))
    side=body.with_suffix('.dat.prep.json');side.write_text(json.dumps(metadata,indent=2,allow_nan=False)+'\n')
    descriptor=dict(format='qcbed-prepared-input-v1',experiment=dict(path=str(body),gradient_option=int(data.mode.value=='gradient'),header='absent',format='legacy_intensity_blocks'),
        scan=asdict(data.scan_geometry),beams=[dict(id=b.id,lattice_coordinates=list(index)) for b,index in zip(data.beams,data.lattice_coordinates)],
        scan_basis_transform=[[1,0],[0,1]],preparation_metadata=dict(path=str(side),sha256=checksum(side)))
    path=out/'prepared.json';path.write_text(json.dumps(descriptor,indent=2,allow_nan=False)+'\n')
    return str(path)


def apply_preparation(config,path):
    """Replace only A inputs. B lattice meaning must be revalidated explicitly."""
    from qcbed.configuration import CrystalScanConfig, CrystalExperimentalDataConfig, ExperimentalBeamConfig
    descriptor=json.loads(Path(path).read_text())
    if descriptor.get('format')!='qcbed-prepared-input-v1':raise ValueError('Unsupported saved preparation format')
    if config.crystallography is None:raise ValueError('Live preparation requires a modern Atomic Model document; legacy replay remains unchanged.')
    base=Path(path).resolve().parent
    exp=dict(descriptor['experiment']);exp['path']=str((base/exp['path']).resolve())
    metadata=descriptor.get('preparation_metadata')
    if metadata:
        side=(base/metadata['path']).resolve()
        if side!=Path(exp['path']+'.prep.json') or checksum(side)!=metadata['sha256']:
            raise ValueError('Prepared mask/origin metadata reference or checksum mismatch')
    c=replace(config,scan=CrystalScanConfig(**{k:tuple(v) if isinstance(v,list) else v for k,v in descriptor['scan'].items()}),
        experiment=CrystalExperimentalDataConfig(**exp),crystallography=replace(config.crystallography,
            beams=tuple(ExperimentalBeamConfig(b['id'],tuple(b['lattice_coordinates'])) for b in descriptor['beams']),
            scan_basis_transform=tuple(map(tuple,descriptor['scan_basis_transform']))))
    return c


def new_preparation_document(path):
    """An explicitly incomplete B/C draft, never a borrowed reference problem."""
    from qcbed.configuration import (QcbedConfig,MicroscopeConfig,CrystallographyConfig,StructureFileConfig,
        CbedReflectionBasisConfig,CrystalScanConfig,CrystalExperimentalDataConfig,RefinementConfig,OptimizerConfig,RunConfig,BoundsConfig)
    from qcbed.modern_yaml import OPTIMIZER_DEFAULTS
    from qcbed.models import MeshConfig
    from qcbed.execution import ForwardExecutionConfig
    base=Path(path).resolve().parent
    # These visible starting controls are not a validated Atomic Model; its file
    # is intentionally absent. Validate B must succeed before execution is possible.
    crystal=CrystallographyConfig(StructureFileConfig(str(base/'SELECT_STRUCTURE_FILE')),(0,0,1),
        CbedReflectionBasisConfig((1,0,0),(0,1,0)),MeshConfig(16,16),1000.,(0.,1.),(),2.5,(),((1,0),(0,1)))
    descriptor=json.loads(Path(path).read_text());scan=descriptor['scan'];exp=descriptor['experiment']
    config=QcbedConfig('New preparation — configure Atomic Model and QCBED controls',MicroscopeConfig(200.,(0,0,1)),None,
        CrystalScanConfig(**{k:tuple(v) if isinstance(v,list) else v for k,v in scan.items()}),CrystalExperimentalDataConfig(**exp),
        RefinementConfig(OptimizerConfig(method='SIMP',**OPTIMIZER_DEFAULTS),(),-1.,objective_mode='masked',invalid_trial_policy='reject',bounds=BoundsConfig(policy='reject')),
        ForwardExecutionConfig('batched'),RunConfig(str(base/'run'),str(base/'run/final_out.dat')),crystallography=crystal)
    return apply_preparation(config,path)


def save_raw_project(document,path,*,base_directory,overwrite=False):
    """Atomic frontend draft. It has no invented prepared arrays or scan corners."""
    import os,tempfile,yaml,copy
    base=Path(base_directory)
    if not base.is_absolute():raise ValueError('An explicit absolute project base is required')
    destination=(base/path).resolve()
    prep=document['data_prep']
    from qcbed.modern_yaml import _mapping
    _mapping(prep,('pattern','aperture','disk_profile'),('pattern',),path='data_prep')
    copied=copy.deepcopy(document)
    copied['data_prep']={k:str((base/v).resolve()) for k,v in prep.items()}
    atomic=copied.get('atomic_model',{})
    if 'structure' in atomic:
        source=atomic['structure']
        if isinstance(source,str):atomic['structure']=str((base/source).resolve())
        else:source['path']=str((base/source['path']).resolve())
    output=copied.get('qcbed',{}).get('output',{})
    for key in ('directory','filename'):
        if key in output:output[key]=str((base/output[key]).resolve())
    if destination in [Path(v) for v in copied['data_prep'].values()]:raise ValueError('Cannot overwrite a raw scientific input')
    if destination.exists() and not overwrite:raise FileExistsError(destination)
    fd,temp=tempfile.mkstemp(prefix='.'+destination.name,dir=destination.parent)
    try:
        with os.fdopen(fd,'w') as stream:
            yaml.safe_dump(copied,stream,sort_keys=False);stream.flush();os.fsync(stream.fileno())
        if overwrite:os.replace(temp,destination)
        else:os.link(temp,destination)
    finally:Path(temp).unlink(missing_ok=True)


def derive_preparation_dimensions(geometry, roi, binning):
    """Expand array extents to bin multiples; never alter measured geometry."""
    from qcbed.data_preparation.extraction import ExtractionConfig
    from qcbed.data_preparation.origin_selection import roi_bounds
    from math import ceil
    from qcbed.data_preparation.models import DiskGeometry
    from qcbed.detector.models import _integer
    if not isinstance(geometry,DiskGeometry):raise TypeError('Expected DiskGeometry')
    binning=_integer(binning,'binning')
    if binning<=0:raise ValueError('binning must be positive')
    width=binning*ceil(2*geometry.radius_x_px/binning)
    height=binning*ceil(2*geometry.radius_y_px/binning)
    cfg=ExtractionConfig(width,height,final_binning=binning,
        physical_selection_shape_px=(height,width))
    x0,y0,x1,y1=finalize_symmetry_roi(roi)
    half=((x1-x0)/2,(y1-y0)/2)
    return cfg.selection_width,cfg.selection_height,tuple(int(v) for v in half)


def finalize_symmetry_roi(bounds, image_shape=None):
    """Outward even-grid preview, shifted only to keep all backend trial crops in frame."""
    from math import floor,ceil
    from qcbed.data_preparation.origin_selection import roi_bounds
    x0,y0,x1,y1=roi_bounds(bounds)
    result=[]
    for axis,(lower,upper) in enumerate(((x0,x1),(y0,y1))):
        size=2*ceil((ceil(upper)-floor(lower))/2)
        start=min(floor(lower),max(ceil(upper)-size,floor((lower+upper-size)/2+.5)))
        if image_shape is not None:
            length=image_shape[1-axis]
            minimum=size//2
            maximum=min(length-1-size,length+1-size-size//2)
            if minimum>maximum:
                raise ValueError('The selected region is too close to the image boundary. Please select a slightly smaller region around the centrosymmetric feature.')
            start=min(maximum,max(minimum,start))
        result.append((start,start+size))
    return (result[0][0],result[1][0],result[0][1],result[1][1])
