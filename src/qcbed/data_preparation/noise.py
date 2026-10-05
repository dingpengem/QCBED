"""Legacy whole-frame noise estimation after A3 correction; no geometry/extraction.

Arrays use (ny, nx). See docs/data_preparation_algorithm.md for the literal
Fourier-crop interpretation and the legacy omitted remainder fit group.
"""
from dataclasses import dataclass

import numpy as np
from scipy.optimize import nnls

from .models import NoiseModel
from ..detector.models import _integer, _real, _text


def _pattern(value):
    if not isinstance(value, np.ndarray) or value.dtype.kind not in 'iuf':
        raise TypeError('corrected_pattern must be a real numerical NumPy array')
    if value.ndim != 2 or not all(value.shape):
        raise ValueError('corrected_pattern must be nonempty and two-dimensional (ny, nx)')
    with np.errstate(over='raise', invalid='raise'):
        result = np.array(value, dtype=np.float64, copy=True)
    if not np.isfinite(result).all():
        raise ValueError('corrected_pattern must contain finite float64 values')
    return result


@dataclass(frozen=True)
class NoiseFitConfig:
    """Source defaults; sampling_stride is N's final overwritten binning value."""
    gain_normalization_limit: float = 40000.
    blemish_sigma: float = 10.
    intensity_slices: int = 100
    sampling_stride: int = 16

    def __post_init__(self):
        for name in ('gain_normalization_limit', 'blemish_sigma'):
            value = _real(getattr(self, name), name)
            if value <= 0:
                raise ValueError(f'{name} must be positive')
            object.__setattr__(self, name, value)
        for name in ('intensity_slices', 'sampling_stride'):
            value = _integer(getattr(self, name), name)
            if value <= 0:
                raise ValueError(f'{name} must be positive')
            object.__setattr__(self, name, value)


@dataclass(frozen=True, eq=False)
class NoiseFitResult:
    """Read-only diagnostics. Fit weights are populations, never objective weights.

    slice_points/fit_points/remainder_point columns: intensity, sigma, population.
    high_frequency is the normalized residual; smooth_estimate is unnormalized.
    goodness_of_fit uses the source denominator sum(weights)*(1-3/n_fit_points).
    """
    model: NoiseModel
    config: NoiseFitConfig
    valid_mask: np.ndarray
    blemish_mask: np.ndarray
    noise_mask: np.ndarray
    smooth_estimate: np.ndarray
    high_frequency: np.ndarray
    slice_points: np.ndarray
    fit_points: np.ndarray
    remainder_point: np.ndarray
    reference_rms: float
    high_frequency_normalization: float
    goodness_of_fit: float

    def __post_init__(self):
        for name in ('valid_mask', 'blemish_mask', 'noise_mask', 'smooth_estimate',
                     'high_frequency', 'slice_points', 'fit_points', 'remainder_point'):
            array = getattr(self, name)
            frozen = np.frombuffer(array.tobytes(), dtype=array.dtype).reshape(array.shape)
            object.__setattr__(self, name, frozen)


def _valid_mask(pattern, supplied, limit):
    # N applies Warp(mask,x-4,y-4), then multiplies by its (+4,+4) shift.
    # Zero extension makes the four-pixel outer border unavailable.
    if supplied is None:
        base = np.ones(pattern.shape, dtype=bool)
    else:
        if not isinstance(supplied, np.ndarray) or supplied.dtype != np.bool_:
            raise TypeError('valid_mask must be a Boolean NumPy array')
        if supplied.shape != pattern.shape:
            raise ValueError('valid_mask must match corrected_pattern shape')
        base = supplied
    shifted = np.zeros_like(base)
    shifted[4:, 4:] = base[:-4, :-4]
    mask = np.zeros_like(base)
    mask[:-4, :-4] = shifted[4:, 4:]
    mask &= shifted & (pattern > 0)
    if supplied is None:
        mask &= pattern <= limit
    return mask


def _smooth(pattern):
    """Periodic normalized 3x3 convolution and 4*S-6*S²+4*S³-S⁴."""
    starting = pattern.copy()
    smooth = np.zeros_like(pattern)
    for coefficient in (4., -6., 4., -1.):
        # Separable [1,2,1]/4 kernel = center1, axial.5, diagonal.25 /4.
        starting = (np.roll(starting, 1, 0) + 2*starting + np.roll(starting, -1, 0))/4
        starting = (np.roll(starting, 1, 1) + 2*starting + np.roll(starting, -1, 1))/4
        smooth += coefficient*starting
    return smooth


def _reference_rms(pattern, mask):
    """Literal centered half-width/half-height FFT crop, inverse real part /2."""
    ny, nx = pattern.shape
    spectrum = np.fft.fftshift(np.fft.fft2(pattern))
    crop = spectrum[ny//4:3*ny//4, nx//4:3*nx//4]
    noise = np.fft.ifft2(np.fft.ifftshift(crop)).real / 2
    mini_mask = mask[::2, ::2]
    if not mini_mask.any():
        raise ValueError('No valid pixels in reduced Fourier reference mask')
    return float(np.sqrt(np.mean(noise[mini_mask]**2)))


def _blemishes(residual, mask, threshold):
    values = residual[mask]
    if not values.size:
        raise ValueError('No valid pixels for noise characterization')
    mean = np.mean(values)
    rms = float(np.sqrt(np.mean((values - mean)**2)))
    return np.abs(residual - mean) <= threshold*rms, rms


def _slice_points(smooth, high_frequency, mask, count, stride):
    intensity = (smooth*mask)[::stride, ::stride]
    noise = (high_frequency*mask)[::stride, ::stride]
    weights = mask[::stride, ::stride]
    interval = (np.max(intensity) - np.min(intensity))/count
    if interval <= 0 or not weights.any():
        raise ValueError('No intensity range or valid samples for noise slicing')
    points = np.zeros((count, 3))
    for index in range(count):
        high = np.max(intensity) - index*interval
        low = high - interval
        selected = weights & (intensity > low) & (intensity <= high)
        population = np.count_nonzero(selected)
        if population:
            points[index] = (np.mean(intensity[selected]),
                             np.sqrt(np.mean(noise[selected]**2)), population)
    return points


def _rebin(points):
    """Accumulate >=largest slice population; keep last remainder out of fit (N)."""
    maximum = np.max(points[:, 2])
    remaining = np.sum(points[:, 2])
    if maximum <= 0:
        raise ValueError('No populated noise slices')
    groups = []
    cursor = 0

    def combine(rows):
        weight = np.sum(rows[:, 2])
        if weight <= 0:
            return np.zeros(3)
        return np.array([np.sum(rows[:, 0]*rows[:, 2])/weight,
                         np.sqrt(np.sum(rows[:, 1]**2*rows[:, 2])/weight), weight])

    while remaining > maximum:
        start = cursor
        weight = 0.
        while weight < maximum:
            weight += points[cursor, 2]
            cursor += 1
        groups.append(combine(points[start:cursor]))
        remaining -= weight
    return np.asarray(groups).reshape(-1, 3), combine(points[cursor:])


def _fit_points(points, dark_noise, intensity_units):
    """Deterministic nonnegative weighted LS in physical (a,b,c) coefficients."""
    if len(points) <= 3:
        raise ValueError('Need more than three effective fit points after legacy remainder exclusion')
    intensity, sigma, weights = points.T
    design = np.column_stack((intensity, np.sqrt(intensity), np.ones(len(points))))
    target = sigma.copy()
    if dark_noise is not None:
        design = design[:, :2]
        target -= dark_noise
    if np.linalg.matrix_rank(design) < design.shape[1]:
        raise ValueError('Noise fit is underdetermined: insufficient distinct intensities')
    weighted = design*np.sqrt(weights[:, None])
    # Column scaling changes conditioning, not the weighted objective/constraints.
    scale = np.linalg.norm(weighted, axis=0)
    coefficients = nnls(weighted/scale, target*np.sqrt(weights))[0]/scale
    a, b = coefficients[:2]
    c = coefficients[2] if dark_noise is None else dark_noise
    prediction = a*intensity + b*np.sqrt(intensity) + c
    gof = np.sum(weights*(sigma-prediction)**2)/(np.sum(weights)*(1-3/len(points)))
    return NoiseModel(a, b, c, intensity_units), float(gof)


def characterize_noise(corrected_pattern, *, valid_mask=None, dark_noise=None,
                       config=None, intensity_units='input intensity units') -> NoiseFitResult:
    """Fit A3-corrected whole CBED and return immutable intermediate diagnostics.

    Dimensions must be >=12 and divisible by four so the source half-frame FFT
    crop uses integer coordinates; no implicit resampling. Supplied Boolean masks
    replace automatic gain exclusion, as in N. None fits c; supplied nonnegative
    dark_noise fixes c (zero is an explicit Python extension of N's positive case).
    """
    pattern = _pattern(corrected_pattern)
    if any(n < 12 or n % 4 for n in pattern.shape):
        raise ValueError('Noise fitting requires dimensions >=12 divisible by four for the legacy FFT crop')
    if config is None:
        config = NoiseFitConfig()
    if not isinstance(config, NoiseFitConfig):
        raise TypeError('config must be NoiseFitConfig')
    _text(intensity_units, 'intensity_units')
    if dark_noise is not None:
        dark_noise = _real(dark_noise, 'dark_noise')
        if dark_noise < 0:
            raise ValueError('dark_noise must be nonnegative; use None to fit the constant')
    try:
        with np.errstate(over='raise', invalid='raise', divide='raise'):
            valid = _valid_mask(pattern, valid_mask, config.gain_normalization_limit)
            smooth = _smooth(pattern)
            residual = pattern - smooth
            blemishes, hf_rms = _blemishes(residual, valid, config.blemish_sigma)
            if hf_rms <= np.finfo(float).eps * max(1., float(np.max(np.abs(pattern)))):
                raise ValueError('Degenerate high-frequency signal: cannot normalize noise')
            reference = _reference_rms(pattern, valid)
            normalization = reference/hf_rms
            high_frequency = residual*normalization
            final_mask = valid & blemishes & (smooth > 0)
            points = _slice_points(smooth, high_frequency, final_mask,
                                   config.intensity_slices, config.sampling_stride)
            fit_points, remainder = _rebin(points)
            model, gof = _fit_points(fit_points, dark_noise, intensity_units)
            if not np.isfinite(gof):
                raise ValueError('Non-finite noise fit')
            return NoiseFitResult(model, config, valid, blemishes, final_mask, smooth,
                                  high_frequency, points, fit_points, remainder,
                                  reference, normalization, gof)
    except FloatingPointError as exc:
        raise ValueError(f'Non-finite noise characterization: {exc}') from exc


def fit_noise_model(corrected_pattern, *, valid_mask=None, dark_noise=None,
                    config=None, intensity_units='input intensity units') -> NoiseModel:
    """Return just the model; see characterize_noise for diagnostic arrays."""
    return characterize_noise(corrected_pattern, valid_mask=valid_mask, dark_noise=dark_noise,
                              config=config, intensity_units=intensity_units).model


def sigma_map_from_noise_model(corrected_pattern, model: NoiseModel) -> np.ndarray:
    """Evaluate existing a*abs(I)+b*sqrt(abs(I))+c, with no mask/normalization.

    For positive intensities this is a*I+b*sqrt(abs(I))+c. The established
    negative-intensity extension is retained; this does not change A4 fitting.
    """
    pattern = _pattern(corrected_pattern)
    if not isinstance(model, NoiseModel):
        raise TypeError('model must be NoiseModel')
    if min(model.a, model.b, model.c) < 0:
        raise ValueError('Legacy noise evaluation requires nonnegative coefficients')
    try:
        with np.errstate(over='raise', invalid='raise'):
            intensity = np.abs(pattern)
            sigma = model.a*intensity + model.b*np.sqrt(intensity) + model.c
            from .variance import _sigma_map
            return _sigma_map(sigma)
    except FloatingPointError as exc:
        raise ValueError(f'Non-finite sigma map: {exc}') from exc


def variance_map_from_noise_model(corrected_pattern, model: NoiseModel) -> np.ndarray:
    """Canonical A4 sigma stage followed by the shared sigma² conversion."""
    from .variance import variance_map_from_sigma_map
    return variance_map_from_sigma_map(sigma_map_from_noise_model(corrected_pattern,model))
