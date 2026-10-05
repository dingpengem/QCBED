"""Frontend-independent experimental inspection and non-destructive mask snapshots.

Numerical inputs remain ordinary canonical experimental bodies. Companion JSON is
editor provenance only; runners consume the body through their existing reader.
"""
from dataclasses import replace
from pathlib import Path
import hashlib
import io
import json
import numpy as np
from qcbed import compile_config,resolve_config_paths
from qcbed.experimental_config import load_experimental_data


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_data(config,base):
    config=resolve_config_paths(config,base_directory=base)
    if config.crystallography is not None:
        # A inspection needs only ordered geometric labels, not a valid B model.
        # These labels never enter a runner; production still compiles physical HK.
        from qcbed.experimental import _parse_experimental_body
        from qcbed.legacy_io import _Records
        spec=config.experiment;scan=config.scan;records=_Records(Path(spec.path))
        if spec.header=='validate':
            record=records.read('experimental header')
            actual=tuple(record.integer(n) for n in ('NUMP','NUML','GRADIENTOPT'))
            if actual!=(scan.scan_columns,scan.scan_rows,spec.gradient_option):raise ValueError('Experimental header differs from canonical scan metadata')
        return _parse_experimental_body(records,scan.scan_columns,scan.scan_rows,spec.gradient_option,
            tuple(b.lattice_coordinates for b in config.crystallography.beams))
    return load_experimental_data(compile_config(config),base_directory=base)


def provenance(config):
    path=Path(config.experiment.path)
    if not path.is_absolute():raise ValueError('Resolve experimental path against the explicit document base first')
    side=path.with_suffix(path.suffix+'.prep.json')
    if not side.exists():return {}
    meta=json.loads(side.read_text())
    if not isinstance(meta,dict) or meta.get('format')!='qcbed-manual-mask-v1' or meta.get('data_sha256')!=digest(path):
        raise ValueError('Manual-mask provenance does not match experimental file')
    return meta


def edit_mask(config,base,*,beam=None,box=None,restore=False,reset=False,overview=None,zone_axis_selection=None):
    """Half-open integer (x0,y0,x1,y1), one ordered beam. Never enable source-invalid pixels."""
    config=resolve_config_paths(config,base_directory=base)
    data=read_data(config,base);meta=provenance(config)
    if not np.all((data.weights==0)|(data.weights==1)):
        raise ValueError('Manual masking requires binary input weights; legacy nonbinary weights are unchanged')
    source=meta.get('source',dict(path=str(Path(config.experiment.path).resolve()),
        header=config.experiment.header,sha256=digest(config.experiment.path)))
    if digest(source['path'])!=source['sha256']:raise ValueError('Original mask source changed; restore is unsafe')
    original=read_data(replace(config,experiment=replace(config.experiment,path=source['path'],header=source['header'])),base)
    weights=data.weights.copy()
    if reset:weights[:]=original.weights
    if box is not None:
        x0,y0,x1,y1=box
        if any(type(v) is not int for v in (*box,beam)) or not (0<=beam<weights.shape[1] and 0<=x0<x1<=data.scan_columns and 0<=y0<y1<=data.scan_rows):
            raise ValueError('Select an integer pixel or half-open region inside this beam')
        indices=(np.arange(y0,y1)[:,None]*data.scan_columns+np.arange(x0,x1)).ravel()
        weights[indices,beam]=original.weights[indices,beam] if restore else 0.
    metadata=dict(format='qcbed-manual-mask-v1',source=source,
        excluded_indices=np.argwhere((original.weights!=0)&(weights==0)).tolist(),
        overview=overview if overview is not None else meta.get('overview'),
        zone_axis_selection=zone_axis_selection if zone_axis_selection is not None else meta.get('zone_axis_selection'))
    stream=io.StringIO()
    for a in (data.intensities,data.variances,weights):np.savetxt(stream,a.T.ravel(),fmt='%.17g')
    content=stream.getvalue().encode();metadata['data_sha256']=hashlib.sha256(content).hexdigest()
    identity=hashlib.sha256(content+json.dumps(metadata,sort_keys=True).encode()).hexdigest()
    folder=Path(base)/'.qcbed-data';folder.mkdir(exist_ok=True)
    path=folder/(identity+'.dat');side=path.with_suffix('.dat.prep.json')
    # Exclusive publication: existing snapshots are immutable, including during active runs.
    for target,payload in ((path,content),(side,(json.dumps(metadata,indent=2)+'\n').encode())):
        try:
            with target.open('xb') as f:f.write(payload)
        except FileExistsError:
            if target.read_bytes()!=payload:raise ValueError('Snapshot collision; refusing overwrite')
    return replace(config,experiment=replace(config.experiment,path=str(path),header='absent'),
        legacy_export=replace(config.legacy_export,experimental_filename_raw=None))


def load_overview(config,base,pattern_path,manifest_path):
    """Use an existing Patch A dataset manifest, verified against the selected signal.

    Accept a dataset manifest or closure/variance report containing modes.*.dataset.
    No catalogue/origin is inferred from file order or from physical HK.
    """
    from qcbed.data_preparation.diagnostics import load_pattern
    data=read_data(config,base);raw=json.loads(Path(manifest_path).read_text())
    mode='gradient' if data.gradient_option else 'conventional'
    d=raw['modes'][mode]['dataset'] if 'modes' in raw else raw
    if d.get('signal_mode',mode)!=mode:raise ValueError('Overview signal mode differs from configured input')
    scan=d.get('scan_geometry')
    if scan and (scan['scan_rows'],scan['scan_columns'])!=(data.scan_rows,data.scan_columns):raise ValueError('Overview scan dimensions differ from configured input')
    expected=hashlib.sha256(data.intensities.tobytes(order='C')).hexdigest()
    if d['arrays']['intensities']['sha256']!=expected:raise ValueError('Saved preparation overview does not match the experimental data/order')
    mapping=d['beam_mapping']
    if [b['refinement_column'] for b in mapping]!=list(range(data.intensities.shape[1])):raise ValueError('Invalid overview beam ordering')
    image=load_pattern(pattern_path)
    origin=np.asarray(d['reciprocal_frame']['zone_axis_origin_px'],float)
    if origin.shape!=(2,) or not np.isfinite(origin).all():raise ValueError('Invalid recorded zone-axis origin')
    rois=d['rois'];ids={b['stable_beam_id'] for b in mapping}
    if len(ids)!=len(mapping):raise ValueError('Duplicate overview beam stable IDs')
    from qcbed.data_preparation.models import DiskGeometry
    for roi in rois:
        DiskGeometry(**roi['geometry'])
        center=np.asarray(roi['center_px'],float)
        if center.shape!=(2,) or not np.isfinite(center).all():raise ValueError('Invalid ROI center')
    if any(r['beam_id'] not in ids for r in rois):raise ValueError('Unknown beam in saved preparation regions')
    pattern_record=raw.get('inputs',{}).get('pattern',{})
    expected_pattern=raw.get('pattern_sha256',pattern_record.get('sha256'))
    if expected_pattern and digest(pattern_path)!=expected_pattern:raise ValueError('Full pattern differs from the saved preparation source')
    return image,d


def attach_overview(config,base,pattern,manifest):
    pattern=Path(pattern).resolve();manifest=Path(manifest).resolve()
    load_overview(config,base,pattern,manifest)
    return edit_mask(config,base,zone_axis_selection={},overview=dict(pattern=str(pattern),manifest=str(manifest),
        pattern_sha256=digest(pattern),manifest_sha256=digest(manifest)))


def restored_overview(config,base):
    info=provenance(config).get('overview')
    if not info:return None
    if any(digest(info[k])!=info[k+'_sha256'] for k in ('pattern','manifest')):raise ValueError('Attached overview source changed')
    image,manifest=load_overview(config,base,info['pattern'],info['manifest'])
    selection=provenance(config).get('zone_axis_selection')
    if selection:
        manifest=dict(manifest,reciprocal_frame=dict(manifest['reciprocal_frame'],zone_axis_origin_px=selection['zone_axis_origin_px']))
    return image,manifest

def select_zone_axis(config,base,*,beam,feature=None):
    """Manual symmetry-center pick, or confirmation of the recorded symmetry-center origin.

    Selection is never automatic detection. New positions use unchanged symmetry-center/scan-coordinate
    equations and explicit existing extraction geometry, not inferred dimensions.
    """
    from qcbed.data_preparation.origin import legacy_zone_axis_origin
    from qcbed.data_preparation.scan_conversion import convert_scan_geometry
    image,m=restored_overview(config,base)
    frame=m['reciprocal_frame']
    central=frame.get('central_disk_center_px')
    if feature is None:
        origin=frame['zone_axis_origin_px']
        record=dict(provenance(config).get('zone_axis_selection') or {},method='explicit confirmation of recorded symmetry-center origin',
                    zone_axis_origin_px=list(origin),confirmed=True)
    else:
        feature=np.asarray(feature,float)
        if feature.shape!=(2,) or not np.isfinite(feature).all() or not (0<=feature[0]<image.shape[1] and 0<=feature[1]<image.shape[0]):
            raise ValueError('Select a feature inside the full pattern')
        reference=m['beam_mapping'][beam]
        if central is None:raise ValueError('Central disk center metadata required for symmetry-center relation')
        origin=legacy_zone_axis_origin(central,reference['center_px'],feature.tolist())
        extraction=m.get('extraction_config')
        if extraction is None:raise ValueError('Explicit A9 selection/binning metadata required to update A10 scan geometry')
        result=convert_scan_geometry(reference_point_px=central,zone_axis_origin_px=origin,
            reference_point_1_px=tuple(np.asarray(origin)+frame['basis_vector_1_px']),
            reference_point_2_px=tuple(np.asarray(origin)+frame['basis_vector_2_px']),
            selection_width=extraction['selection_width'],selection_height=extraction['selection_height'],
            final_binning=extraction['final_binning'],expansion_factor=extraction.get('expansion_factor'))
        scan=result.scan_geometry
        config=replace(config,scan=replace(config.scan,corner_a=scan.corner_a,corner_b=scan.corner_b,corner_c=scan.corner_c))
        record=dict(method='manual symmetry-center selection; existing symmetry-center/scan-coordinate equations',
            beam_id=reference['stable_beam_id'],lattice_coordinates=reference.get('lattice_coordinates'),
            reference_beam_center_px=reference['center_px'],symmetry_center_px=feature.tolist(),
            zone_axis_origin_px=list(origin))
    return edit_mask(config,base,zone_axis_selection=record)
