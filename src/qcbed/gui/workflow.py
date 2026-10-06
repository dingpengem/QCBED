"""Editor-only ordered validation receipts; never an alternative numerical model."""
from dataclasses import dataclass, field, replace
from pathlib import Path
from .data_model import digest, provenance, read_data, restored_overview
import numpy as np


def keys(config):
    def identity(path):
        p=Path(path)
        return (str(p),digest(p) if p.is_file() else None)
    exp=identity(config.experiment.path)
    side=identity(str(config.experiment.path)+'.prep.json')
    try:
        meta=provenance(config);overview=meta.get('overview') or {}
        sources=tuple(identity(overview[k]) for k in ('pattern','manifest') if k in overview)
    except (OSError,ValueError):sources=('invalid provenance',)
    a=(config.experiment,config.scan,exp,side,sources,
       (config.crystallography.beams,config.crystallography.scan_basis_transform) if config.crystallography else None)
    crystal=config.crystallography
    b=(config.sample,crystal,config.microscope.zone_axis if hasattr(config.microscope,'zone_axis') else config.microscope,
       identity(crystal.structure.path) if crystal else None)
    if crystal is not None:
        b=(crystal.structure,identity(crystal.structure.path),crystal.zone_axis,
           crystal.cbed_basis,crystal.site_groups,crystal.debye_waller,crystal.absorption,
           config.microscope.zone_axis)
    c=(config.microscope,config.sample,config.scan,config.experiment,
       config.crystallography,config.refinement) if crystal is not None else config
    return a,b,c


def refinement_key(config):
    r=config.refinement
    return (r.parameters,r.bounds,r.active_categories,
            config.crystallography.projected_elastic if config.crystallography else ())


def simulation_key(config):
    """Physical C inputs only; execution/binding edits do not change this receipt."""
    crystal=replace(config.crystallography,projected_elastic=()) if config.crystallography else None
    return (keys(config)[:2], crystal, config.microscope, config.scan,
            config.experiment, config.sample)


@dataclass
class Workflow:
    receipts: list = field(default_factory=lambda:[None,None,None])
    simulation_receipt: object = None
    refinement_receipt: object = None
    snapshot: object = None
    preparation_descriptor: object = None

    def invalidate(self,index):
        if index < 2:
            self.simulation_receipt=None
            self.refinement_receipt=None
        self.receipts[index:]=[None]*(3-index)
        if index==0:
            self.snapshot=None
            self.preparation_descriptor=None

    def sync(self,config):
        if self.simulation_receipt != simulation_key(config):
            self.simulation_receipt=None
            self.refinement_receipt=None
        if self.refinement_receipt != refinement_key(config):
            self.refinement_receipt=None
        current=keys(config)
        for i,receipt in enumerate(self.receipts):
            if receipt is not None and receipt!=current[i]:self.invalidate(i);break

    def accept(self,index,config,*,base_directory=None):
        self.sync(config)
        if any(v is None for v in self.receipts[:index]):raise ValueError('Validate previous stages first')
        if index==0 and config.crystallography and base_directory is not None:
            from qcbed.modern_yaml import publish_prepared_descriptor
            self.preparation_descriptor=publish_prepared_descriptor(config,base_directory)
        self.receipts[index]=keys(config)[index]
        if index==0:self.snapshot=config

    @property
    def ready(self):return all(v is not None for v in self.receipts)


def validate_data(config,base):
    data=read_data(config,base);overview=restored_overview(config,base)
    if overview is None:raise ValueError('Use New Data Preparation or Load Previous Preparation before Validate A')
    meta=provenance(config)
    if not meta.get('zone_axis_selection'):raise ValueError('Select or explicitly confirm the zone-axis position before Validate A')
    if meta['zone_axis_selection'].get('confirmed') is False:raise ValueError('Inspect and explicitly confirm the refined zone-axis position before Validate A')
    mapping=overview[1]['beam_mapping'];rois=overview[1]['rois']
    if any(sum(r['beam_id']==b['stable_beam_id'] for r in rois)!=1 for b in mapping):raise ValueError('Each beam needs exactly one ROI')
    if not all(np.isfinite(a).all() for a in (data.intensities,data.variances,data.weights)):raise ValueError('Nonfinite experimental arrays')
    if not np.all((data.weights==0)|(data.weights==1)):raise ValueError('Prepared inclusion weights must be binary')
    if not np.any(data.weights) or np.any(data.variances[data.weights!=0]<=0):raise ValueError('Active pixels require positive variance')
    manifest=overview[1]
    live=manifest.get('live_preparation',{})
    if live:
        # The source manifest owns physical labels; B supplies their HKL meaning.
        expected=tuple((b['stable_beam_id'],tuple(b['lattice_coordinates'])) for b in mapping)
        actual=tuple((b.id,b.lattice_coordinates) for b in config.crystallography.beams)
        if actual!=expected:raise ValueError('Live preparation beam identity/order differs from A snapshot')
        from qcbed.data_preparation.scan_conversion import convert_scan_geometry
        frame=manifest['reciprocal_frame'];origin=np.asarray(frame['zone_axis_origin_px'])
        cfg=manifest['extraction_config']
        height,width=cfg['physical_selection_shape_px']
        scan=convert_scan_geometry(reference_point_px=frame['central_disk_center_px'],zone_axis_origin_px=tuple(origin),
            reference_point_1_px=tuple(origin+frame['basis_vector_1_px']),reference_point_2_px=tuple(origin+frame['basis_vector_2_px']),
            selection_width=width,selection_height=height,
            final_binning=cfg['final_binning'],expansion_factor=cfg.get('expansion_factor')).scan_geometry
        if any(getattr(config.scan,k)!=getattr(scan,k) for k in ('scan_rows','scan_columns','corner_a','corner_b','corner_c')):
            raise ValueError('Scan geometry differs from live A10 preparation')
        if config.crystallography.scan_basis_transform!=((1,0),(0,1)):
            raise ValueError('Live preparation requires its recorded identity scan basis transform')
    detail=f"; sampling: {live.get('sampling_policy','recorded integer extraction')}; disk geometry: {manifest.get('physical_disk_geometry')}; noise: {live.get('noise_source','recorded prepared input')}"
    return detail+f"{len(mapping)} beams; {np.count_nonzero(data.weights)} active pixels; {len(meta.get('excluded_indices',()))} manually excluded; zone-axis position: {meta['zone_axis_selection']['zone_axis_origin_px']}"


def validate_atomic_model(config):
    """B-only source/indexing/metadata checks; no simulation or optimizer compile."""
    from qcbed.crystallography_config import load_config_structure
    from qcbed.crystallography import ZoneAxis, CrystalOrientation
    from qcbed.crystallography.reflection_basis import OriginalReflectionBasis
    from qcbed.crystallography.site_groups import group_compiler_inputs
    c=config.crystallography
    structure=load_config_structure(config)
    zone=ZoneAxis(c.zone_axis)
    if ZoneAxis(config.microscope.zone_axis)!=zone:
        raise ValueError('Atomic Model zone axis must match the configured orientation.')
    CrystalOrientation(structure.direct_lattice,zone)
    OriginalReflectionBasis(c.cbed_basis.g1_hkl,c.cbed_basis.g2_hkl,zone)
    if c.site_groups:
        for g in c.site_groups:
            if g.debye_waller_b is None:
                raise ValueError(f'Debye–Waller B is required for site group {g.id}.')
        group_compiler_inputs(config,structure)
    else:
        symbols={s.symbol for s in structure.sites}
        thermal={t.species:t.value_angstrom_squared for t in c.debye_waller}
        if len(thermal)!=len(c.debye_waller) or set(thermal)!=symbols:
            raise ValueError('Define physical site groups and Debye–Waller B (Å²).')
        if any(s.debye_waller_b is not None and s.debye_waller_b!=thermal[s.symbol] for s in structure.sites):
            raise ValueError('Source and configured Debye–Waller factors conflict.')
        coefficients={v.species for v in c.absorption.coefficients}
        if c.absorption.enabled and (coefficients!=symbols or len(coefficients)!=len(c.absorption.coefficients)):
            raise ValueError('Absorption coefficients must cover each element exactly once.')
        if not c.absorption.enabled and c.absorption.coefficients:
            raise ValueError('Disabled absorption must have no coefficients.')
    return structure
