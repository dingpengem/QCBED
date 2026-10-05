"""Whole-frame direct OTF division, separate from aperture calibration and noise."""
import numpy as np
from numpy.typing import NDArray

from .models import DetectorResponseModel, _real
from .psf import _fourier_divide, PSFCalibrationError


def correct_cbed_with_otf(
    raw_pattern: np.ndarray,
    detector_response: DetectorResponseModel,
    *,
    fourier_relative_floor: float = 1e-12,
    imaginary_relative_tolerance: float = 1e-12,
    imaginary_absolute_tolerance: float = 1e-12,
) -> NDArray[np.float64]:
    """Return IFFT(FFT(raw)/OTF).real for a complete (ny, nx) frame.

    OTF is unshifted, as in A2's fft2(ifftshift(centered_psf)); no shifts,
    padding, intensity normalization, clipping, or ROI operations are applied.
    Small nonzero OTF magnitudes are floored relative to max(abs(OTF)), retaining
    phase; exactly zero transfer bins are omitted, as in A2. An all-zero OTF
    fails. These safeguards cannot recover information in null frequencies.

    Imaginary residuals must satisfy max(abs(imag)) <= absolute_tolerance +
    relative_tolerance * max(abs(real)), in input intensity units. Otherwise
    ValueError is raised. Neither input nor detector arrays are modified.
    """
    if not isinstance(raw_pattern, np.ndarray) or raw_pattern.dtype.kind not in 'iuf':
        raise TypeError('raw_pattern must be a real numerical NumPy array')
    if raw_pattern.ndim != 2 or not all(raw_pattern.shape):
        raise ValueError('raw_pattern must be nonempty and two-dimensional (ny, nx)')
    if not isinstance(detector_response, DetectorResponseModel):
        raise TypeError('detector_response must be DetectorResponseModel')
    otf = detector_response.otf
    if otf is None:
        raise ValueError('detector_response must provide an OTF; calibration is separate')
    if raw_pattern.shape != otf.shape:
        raise ValueError('raw_pattern and OTF must have identical (ny, nx) shapes; no resampling')
    floor = _real(fourier_relative_floor, 'fourier_relative_floor')
    rtol = _real(imaginary_relative_tolerance, 'imaginary_relative_tolerance')
    atol = _real(imaginary_absolute_tolerance, 'imaginary_absolute_tolerance')
    if not 0 < floor < 1:
        raise ValueError('fourier_relative_floor must be between zero and one')
    if rtol < 0 or atol < 0:
        raise ValueError('imaginary tolerances must be nonnegative')
    try:
        with np.errstate(over='raise', invalid='raise', divide='raise'):
            raw = np.array(raw_pattern, dtype=np.float64, copy=True)
            if not np.all(np.isfinite(raw)):
                raise ValueError('raw_pattern must contain finite float64 values')
            spectrum, _ = _fourier_divide(np.fft.fft2(raw), otf, floor)
            result = np.fft.ifft2(spectrum)
            if not np.all(np.isfinite(result)):
                raise ValueError('OTF correction produced non-finite values')
            residual = float(np.max(np.abs(result.imag)))
            limit = atol + rtol * float(np.max(np.abs(result.real)))
            if residual > limit:
                raise ValueError(f'OTF correction has significant imaginary residual {residual:g} '
                                 f'(allowed {limit:g}); check OTF phase/Hermitian symmetry')
            return result.real.copy()
    except (FloatingPointError, PSFCalibrationError) as exc:
        raise ValueError(f'Unusable OTF correction: {exc}') from exc
