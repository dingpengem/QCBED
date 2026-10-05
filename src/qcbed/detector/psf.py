"""AMMPSF aperture calibration, independent of image loaders and CBED correction.

Source: Nakashima AMMPSF.s (2012), documented in data_preparation_algorithm.md.
Spatial PSFs are centered at (ny//2,nx//2); transfer arrays use unshifted NumPy
frequency order. Bilinear resampling uses zero extension. These explicit Python
conventions still require comparison against exported DigitalMicrograph fixtures.
"""
from dataclasses import asdict, dataclass
import copy
import json

import numpy as np

from .models import DetectorResponseModel, PSFCalibrationMetrics, _integer, _real


SMALL_CORE_RADIUS_PX = 4.99
LARGE_CORE_RADIUS_PX = 24.99
LEGACY_MAX_DIM = 13
LEGACY_ROTATION_STEPS = 16
LEGACY_ITERATION_MISMATCH_LIMIT = 50.
LEGACY_CANDIDATE_MISMATCH_LIMIT = 25.


class PSFCalibrationError(ValueError):
    """Degenerate data, unusable numerical state, or exhausted candidate budget."""


@dataclass(frozen=True)
class PSFCalibrationConfig:
    """Legacy radii/selection defaults plus explicit finite Python safety budgets."""
    small_core_radius_px: float = SMALL_CORE_RADIUS_PX
    large_core_radius_px: float = LARGE_CORE_RADIUS_PX
    accepted_measurements: int = 10
    max_attempts: int = 100
    max_iterations: int = 100
    fourier_relative_floor: float = 1e-12

    def __post_init__(self):
        for name in ('small_core_radius_px', 'large_core_radius_px', 'fourier_relative_floor'):
            value = _real(getattr(self, name), name)
            if value <= 0:
                raise ValueError(f'{name} must be positive')
            object.__setattr__(self, name, value)
        if self.small_core_radius_px > self.large_core_radius_px:
            raise ValueError('small core radius must not exceed large core radius')
        if self.fourier_relative_floor >= 1:
            raise ValueError('fourier_relative_floor must be less than one')
        for name in ('accepted_measurements', 'max_attempts', 'max_iterations'):
            value = _integer(getattr(self, name), name)
            if value <= 0:
                raise ValueError(f'{name} must be positive')
            object.__setattr__(self, name, value)
        if self.max_attempts < self.accepted_measurements:
            raise ValueError('max_attempts must allow the requested accepted measurements')


def _finite(values, name):
    if not np.all(np.isfinite(values)):
        raise PSFCalibrationError(f'{name}: non-finite numerical state')
    return values


def _sum_normalize(values):
    total = float(np.sum(_finite(values, 'PSF')))
    if not np.isfinite(total) or total <= 0:
        raise PSFCalibrationError('PSF/aperture must have a positive finite sum')
    return _finite(values / total, 'sum normalization')


def _peak_normalize(values):
    peak = float(np.max(_finite(values, 'peak input')))
    if peak <= 0:
        raise PSFCalibrationError('No positive peak for normalization')
    return _finite(values / peak, 'peak normalization')


def _sample(image, x, y):
    """Bilinear samples at (x=column,y=row); zero beyond the pixel grid."""
    x0, y0 = np.floor(x).astype(int), np.floor(y).astype(int)
    fx, fy = x-x0, y-y0
    result = np.zeros_like(x, dtype=np.float64)
    ny, nx = image.shape
    for dx, dy, weight in ((0, 0, (1-fx)*(1-fy)), (1, 0, fx*(1-fy)),
                           (0, 1, (1-fx)*fy), (1, 1, fx*fy)):
        xx, yy = x0+dx, y0+dy
        inside = (xx >= 0) & (xx < nx) & (yy >= 0) & (yy < ny)
        result[inside] += image[yy[inside], xx[inside]]*weight[inside]
    return result


def _rotational_sum(image, start=1, stop=LEGACY_ROTATION_STEPS):
    """Successive legacy angle additions, then x/y reflections of updated sums.

    Angles are 360/2**k degrees. This is not a one-pass radial bin average.
    Reflection is about the centered pixel, not about the array's half-pixel edge.
    """
    ny, nx = image.shape
    y, x = np.indices(image.shape, dtype=float)
    cy, cx = ny//2, nx//2
    result = image.copy()
    for k in range(start, stop+1):
        angle = 2*np.pi/(2**k)
        c, s = np.cos(angle), np.sin(angle)
        rotated = _sample(result, (x-cx)*c+(y-cy)*s+cx, -(x-cx)*s+(y-cy)*c+cy)
        result = result + rotated
        result = result + _sample(result, 2*cx-x, y)
        result = result + _sample(result, x, 2*cy-y)
    return _finite(result, 'rotational averaging')


def _sharp_aperture(current, threshold, dim, measured_sum):
    """Sum shifted hard threshold masks, then match measured aperture sum."""
    y, x = np.indices(current.shape, dtype=float)
    shifts = np.array([0.]) if dim == 1 else np.linspace(.5, -.5, dim)
    sharp = np.zeros_like(current)
    for dy in shifts:
        for dx in shifts:
            sharp += _sample(current, x+dx, y+dy) > threshold
    return _sum_normalize(sharp)*measured_sum


def _fourier_divide(numerator, denominator, relative_floor):
    """Phase-preserving magnitude floor, not a Wiener/noise-spectrum filter.

    Exactly zero bins have no invertible transfer and their quotient is zero
    (null-space omission). No arbitrary phase is assigned; counts are recorded.
    """
    _finite(numerator, 'Fourier numerator')
    _finite(denominator, 'Fourier denominator')
    magnitude = np.abs(denominator)
    scale = float(np.max(magnitude))
    if scale <= 0:
        raise PSFCalibrationError('Unusable Fourier denominator: all components are zero')
    floor = max(scale*relative_floor, np.finfo(float).tiny)
    zeros = magnitude == 0
    safe = denominator.copy()
    small = (magnitude < floor) & ~zeros
    safe[small] *= floor/magnitude[small]
    safe[zeros] = 1
    quotient = numerator/safe
    quotient[zeros] = 0
    return _finite(quotient, 'Fourier division'), int(np.count_nonzero(small | zeros))


def _otf_from_psf(centered_psf):
    return np.fft.fft2(np.fft.ifftshift(centered_psf))


def _real_inverse(spectrum):
    # Full-complex equivalent of DM RealIFFT; real-valued image output is explicit.
    return _finite(np.fft.ifft2(spectrum).real, 'inverse FFT')


def _regularize_psf(rough, coverage, small_mask, large_mask, dim):
    less = _rotational_sum(rough, 1, 2)
    more = _rotational_sum(less, 3, LEGACY_ROTATION_STEPS)/coverage
    core = _peak_normalize(less)*small_mask + _peak_normalize(more)*(1-small_mask)
    # Frequency display centering must be explicit with unshifted NumPy FFTs.
    magnitude = np.fft.fftshift(np.abs(_otf_from_psf(rough)))
    averaged = _peak_normalize(_rotational_sum(magnitude)/coverage)
    psf = np.fft.fftshift(_real_inverse(np.fft.ifftshift(averaged)))
    psf = _sum_normalize(psf)
    outer = _peak_normalize(_rotational_sum(psf)/coverage)
    combined = _peak_normalize(core)*large_mask + outer*(1-large_mask)
    if dim == 1:
        combined[combined.shape[0]//2, combined.shape[1]//2] = 2.
    return _sum_normalize(combined)


def _positive_otf(otf, dim):
    """AMMPSF min-real decision, returning OTF, next dim and adjustment flag.

    Negative minimum adds twice its absolute value and divides by max real OTF;
    the subsequent +2 cancels the source's -2 (repeat current dim). Positive
    minimum advances dim by two. Exact zero / below -1e30 was ill-defined in
    source control flow and is explicitly rejected rather than looping forever.
    """
    _finite(otf, 'OTF')
    minimum = float(np.min(otf.real))
    if minimum == 0 or minimum < -1e30:
        raise PSFCalibrationError(f'Unusable OTF minimum real value: {minimum}')
    if minimum < 0:
        adjusted = otf + 2*abs(minimum)
        peak = float(np.max(adjusted.real))
        if not np.isfinite(peak) or peak <= 0:
            raise PSFCalibrationError('Unusable adjusted OTF normalization')
        return _finite(adjusted/peak, 'adjusted OTF'), dim, True
    return otf.copy(), dim+2, False


def _quality(measured, reconstructed):
    half = float(np.max(_finite(reconstructed, 'reconvolved aperture')))/2
    if half <= 0:
        raise PSFCalibrationError('No positive reconstructed aperture peak')
    weights = np.abs(1-np.abs(reconstructed-half)/half)
    total = float(np.sum(weights))
    if total <= 0 or not np.isfinite(total):
        raise PSFCalibrationError('Unusable aperture fit weights')
    gof = float(np.sum(weights*(reconstructed-measured)**2)/total)
    mean = float(np.sum(weights*measured)/total)
    if mean <= 0 or not np.isfinite(gof):
        raise PSFCalibrationError('Unusable weighted aperture mean/fit')
    return gof, float(np.sqrt(gof/mean))


def _candidate(measured, config, rng, coverage, small, large):
    measured_fft = np.fft.fft2(measured)
    current = measured.copy()
    sharp = measured.copy()
    dim, best, guarded, adjustments = 1, None, 0, 0
    for iteration in range(1, config.max_iterations+1):
        if dim > LEGACY_MAX_DIM:
            return best, guarded, adjustments
        threshold = float(np.max(sharp))/2
        if dim == 1:
            current = measured.copy()
            threshold = float(np.max(measured))*(1+rng.random())/3
        sharp = _sharp_aperture(current, threshold, dim, float(np.sum(measured)))
        quotient, count = _fourier_divide(measured_fft, np.fft.fft2(sharp), config.fourier_relative_floor)
        guarded += count
        rough = np.fft.fftshift(_sum_normalize(_real_inverse(quotient)))
        centered = _regularize_psf(rough, coverage, small, large, dim)
        otf, next_dim, adjusted = _positive_otf(_otf_from_psf(centered), dim)
        adjustments += int(adjusted)
        if adjusted:
            centered = np.fft.fftshift(_real_inverse(otf))
        quotient, count = _fourier_divide(measured_fft, otf, config.fourier_relative_floor)
        guarded += count
        deconvolved = _real_inverse(quotient)
        current = deconvolved  # measured remains separate and immutable to this routine.
        reconstructed = _real_inverse(np.fft.fft2(sharp)*otf)
        gof, mismatch = _quality(measured, reconstructed)
        if mismatch <= LEGACY_ITERATION_MISMATCH_LIMIT and (best is None or gof < best[1]):
            best = (centered.copy(), gof, mismatch, iteration)
        dim = next_dim
    if dim <= LEGACY_MAX_DIM:
        raise PSFCalibrationError('Candidate iteration budget exhausted before dim advanced beyond 13')
    return best, guarded, adjustments


def calibrate_psf_from_aperture(aperture: np.ndarray, *,
        config: PSFCalibrationConfig | None = None, seed: int | None = None,
        rng: np.random.Generator | None = None) -> DetectorResponseModel:
    """Calibrate decoded (ny,nx) aperture data; no file I/O or CBED correction.

    Return a centered unit-sum PSF, unshifted OTF=FFT(ifftshift(PSF)), and abs(OTF).
    As in AMMPSF, accepted candidates are peak-normalized before averaging; final
    unit-sum scaling is an explicit output convention (DM displays peak-normalized).
    Metrics describe the best accepted *candidate*, not a fictitious ensemble fit;
    per-attempt selections and parameters are retained as JSON provenance.
    A supplied RNG is consumed; seed and RNG are mutually exclusive.
    """
    if config is None:
        config = PSFCalibrationConfig()
    if not isinstance(config, PSFCalibrationConfig):
        raise TypeError('config must be PSFCalibrationConfig')
    if not isinstance(aperture, np.ndarray) or aperture.ndim != 2 or aperture.dtype.kind not in 'iuf':
        raise TypeError('aperture must be a 2D real numeric NumPy array (ny,nx)')
    measured = np.array(aperture, dtype=np.float64, copy=True)
    if not all(measured.shape):
        raise ValueError('aperture must be nonempty')
    _finite(measured, 'aperture')
    ny, nx = measured.shape
    # Largest centered circle must stay within pixel-center bounds on all sides.
    available = min(ny//2, ny-1-ny//2, nx//2, nx-1-nx//2)
    if config.large_core_radius_px > available:
        raise ValueError('image too small for configured core radii about (ny//2,nx//2)')
    if np.max(measured) == np.min(measured) or not np.isfinite(np.sum(measured)) or np.sum(measured) <= 0:
        raise PSFCalibrationError('Degenerate aperture: contrast and positive finite total intensity required')
    if seed is not None and rng is not None:
        raise ValueError('supply seed or rng, not both')
    if seed is not None:
        seed = _integer(seed, 'seed')
        if seed < 0:
            raise ValueError('seed must be nonnegative')
    if rng is None:
        rng = np.random.default_rng(seed)
    if not isinstance(rng, np.random.Generator):
        raise TypeError('rng must be numpy.random.Generator')
    initial_rng = copy.deepcopy(rng.bit_generator.state)
    y, x = np.indices(measured.shape)
    radius = np.hypot(x-nx//2, y-ny//2)
    small = radius <= config.small_core_radius_px
    large = radius <= config.large_core_radius_px
    coverage = _peak_normalize(_rotational_sum(np.ones(measured.shape)))
    accepted, records = [], []
    for attempt in range(1, config.max_attempts+1):
        try:
            with np.errstate(over='raise', divide='raise', invalid='raise'):
                result, guarded, adjustments = _candidate(measured, config, rng, coverage, small, large)
            good = result is not None and result[2] <= LEGACY_CANDIDATE_MISMATCH_LIMIT
            record = dict(attempt=attempt, accepted=good, guarded_fourier_bins=guarded,
                          otf_adjustments=adjustments)
            if result is not None:
                record.update(goodness_of_fit=result[1], poisson_mismatch=result[2], selected_iteration=result[3])
            else:
                record['reason'] = 'No improving candidate with Poisson mismatch <= 50'
            if good:
                accepted.append(result)
        except (PSFCalibrationError, FloatingPointError) as exc:
            record = dict(attempt=attempt, accepted=False, reason=str(exc))
        records.append(record)
        if len(accepted) == config.accepted_measurements:
            break
    if len(accepted) < config.accepted_measurements:
        raise PSFCalibrationError(f'Accepted {len(accepted)}/{config.accepted_measurements} PSFs after '
            f'{len(records)} attempts; last outcome: {records[-1]}')
    psf = _sum_normalize(sum((_peak_normalize(r[0]) for r in accepted), np.zeros_like(measured)))
    otf = _finite(_otf_from_psf(psf), 'ensemble OTF')
    best = min(accepted, key=lambda r: r[1])
    provenance = dict(algorithm='AMMPSF source-derived A2', parameters=asdict(config), seed=seed,
        initial_rng_state=initial_rng, shape=list(measured.shape),
        legacy_limits=dict(max_dim=LEGACY_MAX_DIM, rotation_steps=LEGACY_ROTATION_STEPS,
                           iteration_mismatch=LEGACY_ITERATION_MISMATCH_LIMIT,
                           candidate_mismatch=LEGACY_CANDIDATE_MISMATCH_LIMIT),
        interpolation='bilinear, zero extension', psf_center=[ny//2, nx//2],
        fft='numpy backward normalization, unshifted OTF; PSF fftshift-centered',
        output_normalization='peak-normalized candidate ensemble, then unit-sum PSF',
        metrics_scope='best accepted candidate, not ensemble fit', candidates=records)
    encoded = json.dumps(provenance, sort_keys=True, allow_nan=False,
                         default=lambda v: v.tolist() if isinstance(v, np.ndarray) else int(v))
    return DetectorResponseModel('AMMPSF aperture', psf=psf, otf=otf, mtf=np.abs(otf),
        metrics=PSFCalibrationMetrics(best[1], best[2], len(accepted)), provenance=(encoded,))
