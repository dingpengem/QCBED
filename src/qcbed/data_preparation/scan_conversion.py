"""A10: DM point1/2/3 reciprocal conversion, with no image processing or I/O.

Detector coordinates are (x,y). Basis vectors are MATRIX COLUMNS. The DM output
order is point1 -> corner_a, point2 -> corner_b, point3 -> corner_c; the engine
advances columns along b-a and rows along c-a. No physical HKL is inferred.
"""
from dataclasses import dataclass
import numpy as np
from .models import ScanGeometry, ReciprocalFrame, PreparedExperimentalDataset, _pair
from ..detector.models import _integer, _real


@dataclass(frozen=True)
class ScanConversionResult:
    scan_geometry: ScanGeometry
    reference_point_px: tuple[float,float]
    zone_axis_origin_px: tuple[float,float]
    reference_point_1_px: tuple[float,float]
    reference_point_2_px: tuple[float,float]
    vec1: tuple[float,float]
    vec2: tuple[float,float]
    point4: tuple[float,float]
    points: tuple[tuple[float,float], ...]
    determinant: float
    condition_number: float
    selection_width: int
    selection_height: int
    final_binning: int
    expansion_factor: int | None
    canonical_to_scan_basis: tuple[tuple[int,int],tuple[int,int]] | None
    provenance: tuple[str,...]

    @property
    def handedness(self):
        return 'positive determinant' if self.determinant > 0 else 'negative determinant'


def convert_scan_geometry(*, reference_point_px, zone_axis_origin_px,
                          reference_point_1_px, reference_point_2_px,
                          selection_width, selection_height, final_binning,
                          expansion_factor=None, maximum_condition=1e10):
    """Literal DM point construction, followed by an equivalent 2x2 solve.

    expansion_factor is recorded only: it cancels out of the source corner rule.
    None means unavailable historical metadata, not an assumed magnification.
    """
    reference, origin, r1, r2 = [np.array(_pair(v,n)) for v,n in (
        (reference_point_px,'reference_point_px'),(zone_axis_origin_px,'zone_axis_origin_px'),
        (reference_point_1_px,'reference_point_1_px'),(reference_point_2_px,'reference_point_2_px'))]
    w,h,q = [_integer(v,n) for v,n in ((selection_width,'selection_width'),
        (selection_height,'selection_height'),(final_binning,'final_binning'))]
    if min(w,h,q)<=0 or w%q or h%q:
        raise ValueError('Positive selection dimensions must be divisible by final_binning')
    if expansion_factor is not None:
        expansion_factor=_integer(expansion_factor,'expansion_factor')
        if expansion_factor<=0: raise ValueError('expansion_factor must be positive')
    maximum_condition=_real(maximum_condition,'maximum_condition')
    if maximum_condition<1: raise ValueError('maximum_condition must be >=1')
    matrix=np.column_stack((r1-origin,r2-origin))
    if not np.isfinite(matrix).all(): raise ValueError('Non-finite reciprocal basis')
    determinant=float(np.linalg.det(matrix)); condition=float(np.linalg.cond(matrix))
    if not np.isfinite(determinant) or determinant==0 or not np.isfinite(condition) or condition>maximum_condition:
        raise ValueError(f'Singular/pathological reciprocal basis: determinant={determinant}, condition={condition}')
    point4=reference-origin
    half_x,half_y=(w-q)/2,(h-q)/2
    points=point4+np.array([[-half_x,-half_y],[half_x,-half_y],[-half_x,half_y]])
    corners=np.linalg.solve(matrix,points.T).T
    if not np.isfinite(corners).all(): raise ValueError('Non-finite scan corners')
    pair=lambda a: tuple(float(v) for v in a)
    return ScanConversionResult(ScanGeometry(*(pair(c) for c in corners),scan_rows=h//q,scan_columns=w//q),
        pair(reference),pair(origin),pair(r1),pair(r2),pair(matrix[:,0]),pair(matrix[:,1]),pair(point4),
        tuple(pair(p) for p in points),determinant,condition,w,h,q,expansion_factor,None,
        ('DM N:810-841 / E:198-231: explicit reference points minus A8 origin',
         'point1=a; point2=b; point3=c; x=column,y=row',
         'expansion does not enter source corner equation',
         f'maximum accepted basis condition number={maximum_condition:g}'))


def scan_from_catalogue(catalogue, *, central_beam_id, canonical_basis_beam_ids,
                        zone_axis_origin_px, extraction_config, basis_transform):
    """Resolve IDs, not row positions. Explicit integer columns map canonical -> scan.

    Reference endpoint translation follows N:625-630: preserve displacements
    while reanchoring endpoints at the measured zone origin. This is not selecting
    observed centers or subtracting a displaced origin from unshifted beam centers.
    """
    from dataclasses import replace
    def lookup(id):
        matches=[s for s in catalogue.sites if s.beam_id==id]
        if len(matches)!=1: raise ValueError(f'Ambiguous or missing beam ID: {id}')
        return matches[0]
    central=lookup(central_beam_id)
    if central.integer_indices!=(0,0): raise ValueError('Central beam must be explicitly identified as (0,0)')
    if len(canonical_basis_beam_ids)!=2: raise ValueError('Two ordered basis beam IDs required')
    basis_sites=[lookup(id) for id in canonical_basis_beam_ids]
    if [s.integer_indices for s in basis_sites]!=[(1,0),(0,1)]:
        raise ValueError('Canonical basis IDs must explicitly identify (1,0) then (0,1)')
    transform=np.asarray(basis_transform)
    if transform.shape!=(2,2) or transform.dtype.kind not in 'iu':
        raise ValueError('basis_transform must be a 2x2 integer matrix')
    if abs(int(transform[0,0])*int(transform[1,1])-int(transform[0,1])*int(transform[1,0]))!=1:
        raise ValueError('basis_transform must be unimodular (determinant +/-1)')
    origin=np.array(_pair(zone_axis_origin_px,'zone_axis_origin_px'))
    canonical=np.column_stack([np.array(s.predicted_center_px)-central.predicted_center_px for s in basis_sites])
    scan=canonical@transform
    cfg=extraction_config
    if cfg.physical_selection_shape_px is None:
        raise ValueError('A10 requires an explicit physical selection; diagnostic windows are not scan geometry')
    physical_height,physical_width=cfg.physical_selection_shape_px
    result=convert_scan_geometry(reference_point_px=central.predicted_center_px,
        zone_axis_origin_px=tuple(origin),reference_point_1_px=tuple(origin+scan[:,0]),reference_point_2_px=tuple(origin+scan[:,1]),
        selection_width=physical_width,selection_height=physical_height,
        final_binning=cfg.final_binning,expansion_factor=cfg.expansion_factor)
    return replace(result,canonical_to_scan_basis=tuple(tuple(int(v) for v in row) for row in transform),
        provenance=result.provenance+(f'point4 beam ID={central_beam_id}; canonical basis IDs={tuple(canonical_basis_beam_ids)}',
        'reference endpoints=A8 origin + canonical displacement matrix @ explicit integer transform'))


def assemble_prepared_dataset(extraction, conversion, *, intensity_units, physical_geometry=None, provenance=()):
    """Attach A10 geometry to retained A9 products, without recomputing any values.

    Missing variance stays absent; engine readiness is a separate adapter concern.
    Rejected A9 diagnostic slots are not active beams.
    """
    cfg=extraction.config
    if (cfg.selection_width,cfg.selection_height,cfg.final_binning,cfg.expansion_factor)!=(
        conversion.selection_width,conversion.selection_height,conversion.final_binning,conversion.expansion_factor):
        raise ValueError('A9 extraction geometry does not match A10 conversion')
    beams=tuple(b for b in extraction.beams if b.status=='extracted')
    if not beams: raise ValueError('No retained A9 beams')
    frame=ReciprocalFrame(conversion.reference_point_px,conversion.zone_axis_origin_px,conversion.vec1,conversion.vec2)
    arrays=[None if name == 'variance' and extraction.fields.variance is None else
            np.column_stack([getattr(b,name).ravel(order='C') for b in beams])
            for name in ('signal','variance','final_weight')]
    return PreparedExperimentalDataset(*arrays,tuple(b.beam for b in beams),tuple(b.roi for b in beams),
        frame,conversion.scan_geometry,extraction.fields.signal_mode,intensity_units,
        provenance=conversion.provenance+tuple(provenance)+(
            extraction.fields.uncertainty_provenance or 'Scientific variance unavailable; no uncertainty inferred',
            'A9 values copied unchanged; no physical HKL assignment',
            *(f'{b.beam.id}: {item}' for b in beams for item in b.provenance)),
        physical_disk_geometry=physical_geometry, extraction_config=cfg,
        weight_mode=extraction.fields.weight_mode,
        lattice_coordinates=tuple(b.integer_indices for b in beams),
        center_sources=tuple(b.center_source for b in beams),
        extraction_block_indices=tuple(i for i,b in enumerate(extraction.beams) if b.status=='extracted'))
