"""Compatibility-only integrated-intensity localization; never used by canonical A5."""
from dataclasses import dataclass
import numpy as np
from .models import DiskGeometry
from .noise import _pattern
from .disk_detection import _center
from ..detector.models import _integer

@dataclass(frozen=True)
class DiskSearchLevel:
    bin_factor: int
    center_px: tuple[float, float]  # coordinates on this level's grid
    integrated_intensity: float
    radius_scale: float


@dataclass(frozen=True)
class CentralDiskResult:
    geometry: DiskGeometry
    central_disk_center_px: tuple[float, float]
    mode: str
    integrated_intensity: float
    levels: tuple[DiskSearchLevel, ...]


def sum_bin_image(image: np.ndarray, factor: int) -> np.ndarray:
    """Sum nonoverlapping blocks from [0,0]; reject indivisible frames, never crop."""
    image = _pattern(image)
    factor = _integer(factor, 'factor')
    if factor < 1 or any(n % factor for n in image.shape):
        raise ValueError('Positive bin factor must divide both image dimensions exactly')
    ny, nx = image.shape
    with np.errstate(over='raise', invalid='raise'):
        return image.reshape(ny//factor, factor, nx//factor, factor).sum(axis=(1, 3))


def _score(image, center, rx, ry):
    x, y = center
    # Restrict arithmetic to the ellipse's bounding box, preserving row-major sum.
    x0, x1 = max(0, int(np.floor(x-rx))), min(image.shape[1], int(np.ceil(x+rx))+1)
    y0, y1 = max(0, int(np.floor(y-ry))), min(image.shape[0], int(np.ceil(y+ry))+1)
    yy, xx = np.ogrid[y0:y1, x0:x1]
    mask = ((xx-x)/rx)**2 + ((yy-y)/ry)**2 < 1
    with np.errstate(over='raise', invalid='raise'):
        return float(np.sum(image[y0:y1, x0:x1][mask]))


def ellipse_integrated_intensity(image, geometry: DiskGeometry, center_px) -> float:
    """Sum pixels strictly inside ellipse; off-frame ellipse support is absent."""
    image = _pattern(image)
    if not isinstance(geometry, DiskGeometry):
        raise TypeError('geometry must be DiskGeometry')
    center = _center(center_px, image.shape)
    return _score(image, center, geometry.radius_x_px, geometry.radius_y_px)


def _best(image, candidates, rx, ry):
    best = None
    score = -np.inf
    for x, y in candidates:
        if not (0 <= x <= image.shape[1]-1 and 0 <= y <= image.shape[0]-1):
            continue
        value = _score(image, (x, y), rx, ry)
        if best is None or value > score:  # first candidate wins ties
            best, score = (float(x), float(y)), value
    if best is None:
        raise ValueError('No in-frame ellipse search candidates; reduce binning_steps or disk size')
    return best, score


def locate_central_disk_intensity(corrected_pattern, geometry: DiskGeometry, *,
                        binning_steps: int = 5, preset_center_px=None) -> CentralDiskResult:
    """Search A3-corrected whole CBED, not A4 residual/variance or smoothed input.

    Coarse grid starts at scaled semi-axes, as in N. Each finer level searches
    offsets -1,0,+1 around twice the previous center. Final radii shrink by 0.9.
    Coordinates are (x,y); indexing is [y,x]. No reciprocal origin is assigned.
    Preset mode validates and returns the center unchanged without any binning.
    """
    image = _pattern(corrected_pattern)
    if not isinstance(geometry, DiskGeometry):
        raise TypeError('geometry must be DiskGeometry')
    if preset_center_px is not None:
        center = _center(preset_center_px, image.shape)
        score = _score(image, center, geometry.radius_x_px, geometry.radius_y_px)
        return CentralDiskResult(geometry, center, 'preset', score, ())
    steps = _integer(binning_steps, 'binning_steps')
    if steps < 1 or steps > min(image.shape).bit_length()-1:
        raise ValueError('binning_steps must be >=1 and fit the frame')
    factor = 2**steps
    coarse = sum_bin_image(image, factor)
    rx, ry = geometry.radius_x_px/factor, geometry.radius_y_px/factor
    if rx <= coarse.shape[1]-1 and ry <= coarse.shape[0]-1:
        if _score(np.ones_like(coarse), (rx, ry), rx, ry) == 0:
            raise ValueError('Coarse ellipse contains no sampled pixels; reduce binning_steps')
    candidates = ((x, y) for y in np.arange(ry, coarse.shape[0], 1)
                  for x in np.arange(rx, coarse.shape[1], 1))
    center, score = _best(coarse, candidates, rx, ry)
    levels = [DiskSearchLevel(factor, center, score, 1.)]
    for power in range(steps-1, -1, -1):
        factor = 2**power
        target = sum_bin_image(image, factor)
        scale = .9 if factor == 1 else 1.
        candidates = ((2*center[0]+dx, 2*center[1]+dy)
                      for dx in (-1, 0, 1) for dy in (-1, 0, 1))
        center, score = _best(target, candidates,
                              scale*geometry.radius_x_px/factor,
                              scale*geometry.radius_y_px/factor)
        levels.append(DiskSearchLevel(factor, center, score, scale))
    return CentralDiskResult(geometry, center, 'automatic', score, tuple(levels))
