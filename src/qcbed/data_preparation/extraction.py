"""A9 detector-grid extraction, independent of scan conversion and refinement.

Coordinates are (x,y), arrays (row,column). Fractional Warp interpolation is not
established: this first source-supported path requires integer selection origins.
Out-of-frame samples are zero extended and marked unsupported, never wrapped.
"""
from dataclasses import dataclass, fields
from enum import Enum

import numpy as np

from .models import BeamObservation, DiskGeometry, ReflectionROI, SignalMode
from .noise import _pattern
from .variance import variance_map_from_sigma_map
from ..detector.models import _integer, _real


class WeightMode(str, Enum):
    BEAM_MASK = 'beam_mask'
    NOISE_CHARACTERIZATION = 'noise_characterization'


@dataclass(frozen=True)
class ExtractionConfig:
    selection_width: int
    selection_height: int
    final_binning: int = 1
    expansion_factor: int = 1
    beam_radius_factor: float = .9
    minimum_weight_coverage: float = .875
    minimum_beam_frame_coverage: float = .90
    conventional_minimum_subtraction: bool = True
    physical_selection_shape_px: tuple[int,int] | None = None  # (rows,columns)
    diagnostic_window_shape_px: tuple[int,int] | None = None  # presentation only
    sampling_policy: str = 'strict_integer'

    @property
    def gradient_step_px(self):
        """Legacy binding: original-grid gradient step equals final binning."""
        return self.final_binning

    @property
    def final_pattern_matching_shape_px(self):
        return (self.selection_height//self.final_binning, self.selection_width//self.final_binning)

    @classmethod
    def from_disk_geometry(cls, geometry, **kwargs):
        """Explicit physical-diameter selection policy; never silently round pixels."""
        if not isinstance(geometry, DiskGeometry): raise TypeError('geometry must be DiskGeometry')
        width,height=2*geometry.radius_x_px,2*geometry.radius_y_px
        if not width.is_integer() or not height.is_integer():
            raise ValueError('Physical diameter must be integral; provide an explicit selection instead of rounding')
        return cls(int(width),int(height),physical_selection_shape_px=(int(height),int(width)),**kwargs)

    def __post_init__(self):
        if self.sampling_policy not in ('strict_integer','nearest_detector_grid'):
            raise ValueError('Unknown extraction sampling policy')
        for name in ('selection_width', 'selection_height', 'final_binning', 'expansion_factor'):
            value = _integer(getattr(self, name), name)
            if value <= 0:
                raise ValueError(f'{name} must be positive')
        if self.selection_width % self.final_binning or self.selection_height % self.final_binning:
            raise ValueError('Selection dimensions must be divisible by final_binning')
        for name in ('beam_radius_factor', 'minimum_weight_coverage', 'minimum_beam_frame_coverage'):
            value = _real(getattr(self, name), name)
            if not 0 < value <= 1:
                raise ValueError(f'{name} must be in (0,1]')
        if type(self.conventional_minimum_subtraction) is not bool:
            raise TypeError('conventional_minimum_subtraction must be bool')
        for name in ('physical_selection_shape_px','diagnostic_window_shape_px'):
            shape=getattr(self,name)
            if shape is not None:
                if not isinstance(shape,(tuple,list)) or len(shape)!=2:
                    raise ValueError(f'{name} must contain (rows,columns)')
                shape=tuple(_integer(v,name) for v in shape)
                if min(shape)<=0: raise ValueError(f'{name} must be positive')
                object.__setattr__(self,name,shape)
        if self.physical_selection_shape_px is not None and self.physical_selection_shape_px!=(self.selection_height,self.selection_width):
            raise ValueError('Local numerical extraction must match the declared physical selection')
        # The source Round((e-1)/2) tie rule is unverified. Avoid guessing.
        if self.expansion_factor % 2 == 0:
            raise ValueError('Use odd expansion_factor: DM Round half-tie semantics are unverified')


class _FrozenArrays:
    def __post_init__(self):
        for f in fields(self):
            a = getattr(self, f.name)
            if isinstance(a, np.ndarray):
                object.__setattr__(self, f.name,
                    np.frombuffer(a.tobytes(), dtype=a.dtype).reshape(a.shape))


@dataclass(frozen=True, eq=False)
class SignalFields(_FrozenArrays):
    signal: np.ndarray
    variance: np.ndarray | None
    support_fraction: np.ndarray
    valid: np.ndarray
    beam_geometry_mask: np.ndarray
    signal_mode: SignalMode
    weight_mode: WeightMode
    uncertainty_provenance: str


@dataclass(frozen=True, eq=False)
class BeamExtraction(_FrozenArrays):
    beam: BeamObservation
    roi: ReflectionROI
    integer_indices: tuple[int, int]  # Abstract A7 labels, NOT physical HKL.
    predicted_center_px: tuple[float, float]
    observed_center_px: tuple[float, float] | None
    center_source: str
    catalogue_status: str
    catalogue_source: str
    bounds_xyxy: tuple[int, int, int, int]  # Half-open detector coordinates.
    roi_geometry_mask: np.ndarray
    signal: np.ndarray
    variance: np.ndarray | None
    support_fraction: np.ndarray
    support_weight: np.ndarray
    noise_weight: np.ndarray
    valid_mask: np.ndarray
    detector_in_frame_mask: np.ndarray
    final_weight: np.ndarray
    expected_pixels: int
    in_frame_pixels: int
    noise_valid_pixels: int | None
    status: str
    provenance: tuple[str, ...]

    @property
    def coverage_fraction(self):
        return self.in_frame_pixels / self.expected_pixels if self.expected_pixels else 0.

    @property
    def window_frame_coverage(self):
        return float(self.detector_in_frame_mask.mean())

    @property
    def beam_roi_frame_coverage(self):
        return self.coverage_fraction

    @property
    def noise_valid_fraction(self):
        if self.noise_valid_pixels is None:
            return None
        return self.noise_valid_pixels / self.in_frame_pixels if self.in_frame_pixels else 0.

    @property
    def sigma(self):
        return None if self.variance is None else np.sqrt(self.variance)


@dataclass(frozen=True)
class ExtractionResult:
    """Intermediate A1 ROI/beam products, not a replacement scan dataset.

    No scan corners or physical HKL are inferred. Missing variance remains None;
    such products cannot yet be passed to the refinement experimental loader.
    """
    beams: tuple[BeamExtraction, ...]
    fields: SignalFields
    config: ExtractionConfig

    def blocks(self):
        """Three beam-major flattened blocks; includes rejected zero-weight beams."""
        signal = np.concatenate([b.signal.ravel(order='C') for b in self.beams])
        variance = None if self.fields.variance is None else np.concatenate(
            [b.variance.ravel(order='C') for b in self.beams])
        weight = np.concatenate([b.final_weight.ravel(order='C') for b in self.beams])
        return signal, variance, weight

    def tilt_beam_arrays(self):
        """Same values in A1/engine (local row-major pixel, beam) orientation."""
        return tuple(None if a is None else a.reshape(len(self.beams), -1).T.copy()
                     for a in self.blocks())


def _mask(value, shape, name):
    if not isinstance(value, np.ndarray) or value.dtype != np.bool_ or value.shape != shape:
        raise ValueError(f'{name} must be a Boolean array matching pattern shape')
    return value.copy()


def ellipse_mask(shape, center, geometry):
    y, x = np.indices(shape)
    return ((x-center[0])/geometry.radius_x_px)**2 + ((y-center[1])/geometry.radius_y_px)**2 < 1


def build_beam_geometry_mask(shape, rois):
    """Binary UNION: inside an enabled ellipse=1; outside all ellipses=0."""
    mask = np.zeros(shape, dtype=bool)
    for roi in rois:
        if roi.enabled:
            mask |= ellipse_mask(shape, roi.center_px, roi.geometry)
    return mask


def _shift(a, dx, dy):
    """Return a[y+dy,x+dx] with zero extension (no periodic boundary)."""
    ny, nx = a.shape
    out = np.zeros_like(a)
    x0, x1 = max(0, -dx), min(nx, nx-dx)
    y0, y1 = max(0, -dy), min(ny, ny-dy)
    if x1 > x0 and y1 > y0:
        out[y0:y1, x0:x1] = a[y0+dy:y1+dy, x0+dx:x1+dx]
    return out


def construct_signal_fields(pattern, beam_geometry_mask, *, signal_mode=SignalMode.CONVENTIONAL,
                            weight_mode=WeightMode.BEAM_MASK, final_binning=1,
                            variance_map=None, sigma_map=None, noise_valid_mask=None,
                            detector_valid_mask=None, uncertainty_provenance=''):
    """Legacy eight-neighbor mean minus center; variance is propagated separately.

    Noise mode requires explicit variance/sigma, A4 final mask and provenance.
    Its weight is support, NOT inverse variance: the objective already divides by V.
    """
    image = _pattern(pattern)
    if not isinstance(signal_mode, SignalMode) or not isinstance(weight_mode, WeightMode):
        raise TypeError('signal_mode and weight_mode must be typed enums')
    q = _integer(final_binning, 'final_binning')
    if q <= 0:
        raise ValueError('final_binning must be positive')
    geometry = _mask(beam_geometry_mask, image.shape, 'beam_geometry_mask')
    valid = np.ones(image.shape, dtype=bool)
    if detector_valid_mask is not None:
        valid &= _mask(detector_valid_mask, image.shape, 'detector_valid_mask')
    if variance_map is not None and sigma_map is not None:
        raise ValueError('Supply either sigma_map or variance_map, not both')
    variance = None
    source = sigma_map if sigma_map is not None else variance_map
    if source is not None:
        variance = _pattern(source)
        if variance.shape != image.shape or np.any(variance < 0):
            raise ValueError('Uncertainty must match pattern and be nonnegative')
        if not isinstance(uncertainty_provenance, str) or not uncertainty_provenance.strip():
            raise ValueError('Explicit uncertainty provenance is required')
        if sigma_map is not None:
            variance = variance_map_from_sigma_map(variance)
    if weight_mode == WeightMode.NOISE_CHARACTERIZATION:
        if variance is None or noise_valid_mask is None:
            raise ValueError('noise_characterization requires explicit A4 uncertainty and final noise-valid mask')
        valid &= _mask(noise_valid_mask, image.shape, 'noise_valid_mask')
    elif noise_valid_mask is not None:
        raise ValueError('noise_valid_mask requires noise_characterization; modes are not mixed')
    base = (geometry & valid).astype(float)
    signal = image.copy()
    coverage = base.copy()
    if signal_mode == SignalMode.GRADIENT:
        total = np.zeros_like(image)
        vtotal = np.zeros_like(image) if variance is not None else None
        coverage = np.zeros_like(image)
        full_valid = valid.copy()
        offsets = ((-q,0),(q,0),(0,-q),(0,q),(-q,-q),(q,q),(-q,q),(q,-q))
        for dx, dy in offsets:
            total += _shift(image, dx, dy)
            coverage += _shift(base, dx, dy)
            full_valid &= _shift(valid, dx, dy)
            if variance is not None:
                vtotal += _shift(variance, dx, dy)
        signal = total/8 - image
        coverage = (coverage + base)/9
        variance = None if variance is None else variance + vtotal/64
        valid = full_valid
    if not np.isfinite(signal).all() or (variance is not None and not np.isfinite(variance).all()):
        raise ValueError('Non-finite signal/variance produced by A9 arithmetic')
    semantic = 'sigma squared' if sigma_map is not None else 'variance supplied directly'
    provenance = f'{semantic}: {uncertainty_provenance}' if source is not None else ''
    return SignalFields(signal, variance, coverage, valid, geometry, signal_mode,
                        weight_mode, provenance)


def threshold_support(coverage, minimum=.875):
    return (np.asarray(coverage) > minimum).astype(float)


def _crop(a, bounds):
    x0, y0, x1, y1 = bounds
    out = np.zeros((y1-y0, x1-x0), dtype=a.dtype)
    left, top, right, bottom = max(0,x0), max(0,y0), min(a.shape[1],x1), min(a.shape[0],y1)
    if right > left and bottom > top:
        out[top-y0:bottom-y0, left-x0:right-x0] = a[top:bottom, left:right]
    return out


def _integrate(a, config, *, variance=False, threshold=False):
    """Literal expansion, q-by-q accumulation, then e-spaced sampling.

    Variance sum scaling e²/q⁴ is the source semi-integration convention, not an
    assertion that repeated/interpolated detector samples are independent.
    """
    e, q = config.expansion_factor, config.final_binning
    huge = np.repeat(np.repeat(a, e, axis=0), e, axis=1)
    summed = np.zeros((huge.shape[0]//q, huge.shape[1]//q))
    for y in range(q):
        for x in range(q):
            summed += huge[y::q, x::q]
    result = summed * (e*e / q**4) if variance else summed / q**2
    if threshold:
        result = threshold_support(result, config.minimum_weight_coverage)
    offset = (e-1)//2
    return result[offset::e, offset::e]


def _roi_counts(shape, roi):
    cx, cy = roi.center_px
    g = roi.geometry
    gx0, gy0 = int(np.floor(cx-g.radius_x_px)), int(np.floor(cy-g.radius_y_px))
    gx1, gy1 = int(np.ceil(cx+g.radius_x_px))+1, int(np.ceil(cy+g.radius_y_px))+1
    expected = int(ellipse_mask((gy1-gy0,gx1-gx0),(cx-gx0,cy-gy0),g).sum())
    geom = ellipse_mask(shape, (cx,cy), g)
    return expected, geom


def extract_beam_data(site, roi, fields, config, *, noise_valid_mask=None):
    """Extract at predicted center; integer origins only, no hidden snapping.

    Local radial mask follows N exactly: radius width/2, including its boundary,
    centered at (width/2,height/2). Rectangles retain this explicit x-radius rule.
    """
    w, h = config.selection_width, config.selection_height
    cx, cy = roi.center_px
    x0, y0 = cx-w/2, cy-h/2
    ideal_origin=(x0,y0)
    if config.sampling_policy=='nearest_detector_grid':
        x0,y0=map(float,np.rint((x0,y0)))
    if not x0.is_integer() or not y0.is_integer():
        raise ValueError(f'{site.beam_id}: fractional selection origin; DM interpolation is unverified; no snapping applied')
    bounds = (int(x0), int(y0), int(x0)+w, int(y0)+h)
    expected, geom = _roi_counts(fields.signal.shape, roi)
    inside = int(geom.sum())
    noise_count = None if noise_valid_mask is None else int((geom & noise_valid_mask).sum())
    status = 'extracted' if expected and inside/expected >= config.minimum_beam_frame_coverage else 'rejected_coverage'
    local = _crop(fields.signal, bounds)
    in_frame = _crop(np.ones(fields.signal.shape, dtype=bool), bounds)
    y,x = np.indices((h,w))
    radial = (x-w/2)**2 + (y-h/2)**2 <= (w/2)**2
    local_center=(cx-x0,cy-y0)
    local_geom = ellipse_mask((h,w), local_center, roi.geometry)
    local_valid = _crop(fields.valid, bounds)
    participating = local_geom & radial & local_valid
    minimum = float(local[participating].min()) if (participating.any() and
        fields.signal_mode == SignalMode.CONVENTIONAL and config.conventional_minimum_subtraction) else 0.
    signal = _integrate(np.where(in_frame, local-minimum, 0.), config)
    variance = None if fields.variance is None else _integrate(_crop(fields.variance,bounds),config,variance=True)
    support = _crop(fields.support_fraction,bounds)*radial*local_geom
    fraction = _integrate(support,config)
    support_weight = _integrate(support,config,threshold=True)
    # Hard validity is distinct from fractional geometric coverage: rejected
    # detector/noise pixels or out-of-frame stencil terms cannot acquire weight.
    usable = _integrate(_crop(fields.valid,bounds).astype(float),config) == 1
    support_weight *= usable
    # No incomplete bin/stencil is scientifically evaluated. Finite positive V
    # prevents 0/0 in the unchanged legacy objective even when its weight is zero.
    detector_valid = np.ones(fields.signal.shape, dtype=bool)
    if fields.signal_mode == SignalMode.GRADIENT:
        q = config.final_binning
        frame = detector_valid.copy()
        for dx,dy in ((-q,0),(q,0),(0,-q),(0,q),(-q,-q),(q,q),(-q,q),(q,-q)):
            detector_valid &= _shift(frame,dx,dy)
    measured = _integrate(_crop(detector_valid,bounds).astype(float),config) == 1
    signal = np.where(measured, signal, 0.)
    if variance is not None:
        variance = np.where(measured, variance, 1.)
    noise_weight = np.ones_like(support_weight)
    final = support_weight * noise_weight
    if status != 'extracted':
        final[:] = 0
    beam = BeamObservation(site.beam_id, roi.center_px)
    return BeamExtraction(beam,roi,site.integer_indices,site.predicted_center_px,
        site.observed_center_px,'predicted',site.status,site.source,bounds,
        local_geom,signal,variance,fraction,support_weight,noise_weight,usable,in_frame,final,
        expected,inside,noise_count,status,
        ('center=canonical predicted; no center search', 'local_mask=legacy x-radius circle',
         f'conventional_minimum_subtracted={minimum!r}; domain=valid local beam support',
         'unmeasured/incomplete output: signal=0, supplied variance placeholder=1, weight=0',
         f'sampling_policy={config.sampling_policy}; integer sampling; zero extension; invalid stencil excluded',
         f'ideal_origin_px={ideal_origin}; sampled_origin_px={(x0,y0)}; sampling_offset_px={(x0-ideal_origin[0],y0-ideal_origin[1])}; physical center unchanged',
         f'expansion={config.expansion_factor}; binning={config.final_binning}',
         fields.uncertainty_provenance or 'uncertainty unavailable; not fabricated'))


def prepare_experimental_qcbed(pattern, catalogue, geometry, *, selected_beam_ids,
                               config, signal_mode=SignalMode.CONVENTIONAL,
                               weight_mode=WeightMode.BEAM_MASK, **noise_inputs):
    """Preserve canonical catalogue order, irrespective of selection-list order.

    All statuses require explicit selection. Unconfirmed (missing), weak and
    edge-clipped lattice predictions are eligible; rejected/unknown sites are refused.
    A7 observation status is preserved, not promoted by geometric eligibility.
    """
    if not isinstance(config, ExtractionConfig) or not isinstance(geometry, DiskGeometry):
        raise TypeError('config/geometry must be ExtractionConfig/DiskGeometry')
    image = _pattern(pattern)
    ids = tuple(selected_beam_ids)
    sites = tuple(catalogue.sites)
    catalogue_ids = [s.beam_id for s in sites]
    if len(set(catalogue_ids)) != len(sites) or any(not i for i in catalogue_ids):
        raise ValueError('Catalogue must have unique nonempty stable beam IDs')
    if not ids or len(set(ids)) != len(ids) or not set(ids) <= set(catalogue_ids):
        raise ValueError('Select unique existing beam IDs explicitly')
    chosen = tuple(s for s in sites if s.beam_id in ids)
    allowed = {'reference', 'observed', 'defining', 'confirmed', 'weak', 'weak_evidence', 'edge_clipped', 'missing'}
    if any(s.status not in allowed for s in chosen):
        raise ValueError('Rejected or unknown catalogue sites are not extractable')
    scaled = DiskGeometry(geometry.radius_x_px*config.beam_radius_factor,
                          geometry.radius_y_px*config.beam_radius_factor)
    rois = tuple(ReflectionROI(f'roi:{s.beam_id}', s.beam_id, s.predicted_center_px, scaled) for s in chosen)
    retained = []
    for roi in rois:
        expected, geom = _roi_counts(image.shape, roi)
        if expected and geom.sum()/expected >= config.minimum_beam_frame_coverage:
            retained.append(roi)
    mask = build_beam_geometry_mask(image.shape,retained)
    result_fields = construct_signal_fields(image,mask,signal_mode=signal_mode,
        weight_mode=weight_mode,final_binning=config.final_binning,**noise_inputs)
    beams = tuple(extract_beam_data(s,r,result_fields,config,
                  noise_valid_mask=noise_inputs.get('noise_valid_mask')) for s,r in zip(chosen,rois))
    return ExtractionResult(beams,result_fields,config)
