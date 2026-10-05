"""A5 profile geometry and direct-disk localization, never reciprocal origin."""
from dataclasses import dataclass

import numpy as np

from .models import DiskGeometry, _pair
from .noise import _pattern
from ..detector.models import _integer


@dataclass(frozen=True)
class DiskGeometryEstimate:
    geometry: DiskGeometry
    initial_threshold: float
    final_threshold: float
    selected_mean_intensity: float
    selected_pixels: int


def estimate_disk_geometry(reference_image: np.ndarray) -> DiskGeometryEstimate:
    """Threshold a separate disk-profile image, not a detector calibration PSF.

    Strictly above mean, then above half that selection's mean (N:419–448).
    Max row population measures x diameter; max column population y diameter.
    Diameters are pixel counts, divided by two to give A1 ellipse semi-axes.
    """
    image = _pattern(reference_image)
    with np.errstate(over='raise', invalid='raise'):
        threshold = float(np.mean(image))
        selected = image > threshold
        if not selected.any():
            raise ValueError('Disk reference threshold selected no pixels')
        final_threshold = float(np.mean(image[selected])/2)
        selected = image > final_threshold
        if not selected.any() or selected.all():
            raise ValueError('Disk reference threshold does not isolate a disk')
        diameter_x = int(np.max(np.sum(selected, axis=1)))
        diameter_y = int(np.max(np.sum(selected, axis=0)))
        return DiskGeometryEstimate(DiskGeometry(diameter_x/2, diameter_y/2),
                                    threshold, final_threshold,
                                    float(np.mean(image[selected])), int(selected.sum()))


def _center(center, shape):
    x, y = _pair(center, 'central_disk_center_px')
    if not (0 <= x <= shape[1]-1 and 0 <= y <= shape[0]-1):
        raise ValueError('central_disk_center_px must lie inside the pixel-center frame')
    return x, y


