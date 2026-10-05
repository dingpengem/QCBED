"""A7 global pixel lattice and predicted catalogue; no physical HKL or final origin."""
from dataclasses import dataclass
from enum import Enum
from itertools import product

import numpy as np
from scipy.ndimage import map_coordinates

from .models import BeamObservation, DiskGeometry
from .noise import _pattern
from .lattice_detection import ReciprocalLatticeSeed, LatticeCandidate, _canonical, disk_response_map
from ..detector.models import _integer, _real


@dataclass(frozen=True)
class LatticeRefinementConfig:
    residual_threshold_px: float = 1.5
    reassignment_improvement_px: float = .1
    max_iterations: int = 20
    min_inliers: int = 4
    local_radius_scale: float = .5
    local_step_px: float = .25
    observed_response_threshold: float = 1.
    weak_response_threshold: float = .1
    minimum_mask_coverage: float = .95
    max_enumeration_sites: int = 100000

    def __post_init__(self):
        for name in ('max_iterations', 'min_inliers', 'max_enumeration_sites'):
            value = _integer(getattr(self, name), name)
            if value < 1:
                raise ValueError(f'{name} must be positive')
            object.__setattr__(self, name, value)
        if self.min_inliers < 3:
            raise ValueError('min_inliers must be >=3')
        for name in ('residual_threshold_px', 'reassignment_improvement_px',
                     'local_radius_scale', 'local_step_px', 'observed_response_threshold',
                     'weak_response_threshold', 'minimum_mask_coverage'):
            value = _real(getattr(self, name), name)
            if value <= 0:
                raise ValueError(f'{name} must be positive')
            object.__setattr__(self, name, value)
        if self.minimum_mask_coverage > 1 or self.local_radius_scale > 1:
            raise ValueError('coverage and local_radius_scale must be <=1')
        if self.weak_response_threshold > self.observed_response_threshold:
            raise ValueError('weak threshold must not exceed observed threshold')


@dataclass(frozen=True)
class RefinedCandidateAssignment:
    candidate: LatticeCandidate
    initial_integer_indices: tuple[int, int]
    integer_indices: tuple[int, int]
    residual_px: float
    inlier: bool
    reassigned: bool


@dataclass(frozen=True)
class RefinedReciprocalLattice:
    initial_seed: ReciprocalLatticeSeed
    origin_px: tuple[float, float]
    basis_vector_1_px: tuple[float, float]
    basis_vector_2_px: tuple[float, float]
    initial_rms_residual_px: float
    fit_rms_residual_px: float
    assignments: tuple[RefinedCandidateAssignment, ...]
    iterations: int
    config: LatticeRefinementConfig


class BeamStatus(str, Enum):
    OBSERVED = 'observed'
    WEAK = 'weak'
    MISSING = 'missing'
    REJECTED = 'rejected'


@dataclass(frozen=True)
class PredictedBeam:
    observation: BeamObservation  # center is always the lattice prediction; indices=None
    lattice_indices: tuple[int, int]
    predicted_center_px: tuple[float, float]
    observed_center_px: tuple[float, float] | None
    is_central: bool
    status: BeamStatus
    response_score: float | None
    mask_coverage: float
    local_search_radius_px: float


@dataclass(frozen=True)
class LatticeBeamCatalogue:
    lattice: RefinedReciprocalLattice
    image_shape: tuple[int, int]
    geometry: DiskGeometry
    beams: tuple[PredictedBeam, ...]

    @property
    def observations(self) -> tuple[BeamObservation, ...]:
        return tuple(b.observation for b in self.beams)

    @property
    def status_counts(self) -> dict[str, int]:
        return {status.value: sum(b.status == status for b in self.beams) for status in BeamStatus}


def _config(config):
    if config is None:
        return LatticeRefinementConfig()
    if not isinstance(config, LatticeRefinementConfig):
        raise TypeError('config must be LatticeRefinementConfig')
    return config


def _matrix(lattice):
    matrix = np.column_stack((lattice.basis_vector_1_px, lattice.basis_vector_2_px)).astype(float)
    if not np.isfinite(matrix).all() or np.linalg.cond(matrix) > 1e6:
        raise ValueError('Lattice basis is singular or ill-conditioned')
    return matrix


def _unique(mask, indices, residual, candidates):
    result = mask.copy()
    used = set()
    for i in sorted(range(len(candidates)), key=lambda i: (residual[i], -candidates[i].response_score, candidates[i].id)):
        key = tuple(indices[i])
        if result[i]:
            if key in used:
                result[i] = False
            else:
                used.add(key)
    return result


def refine_reciprocal_lattice(seed: ReciprocalLatticeSeed, *, config=None) -> RefinedReciprocalLattice:
    """Fit origin and both vectors to fixed integer labels, trim and reassign.

    Strictly improving reassignments use a pixel margin. Repeated discrete states
    or an exhausted iteration budget fail rather than return an oscillating fit.
    All candidates remain in diagnostics, including rejected/duplicate observations.
    """
    config = _config(config)
    if not isinstance(seed, ReciprocalLatticeSeed):
        raise TypeError('seed must be ReciprocalLatticeSeed')
    matrix = _matrix(seed)
    origin = np.array(seed.reference_center_px, dtype=float)
    if origin.shape != (2,) or not np.isfinite(origin).all():
        raise ValueError('Invalid seed reference center')
    records = tuple(sorted(seed.assignments, key=lambda a: a.candidate.id))
    candidates = tuple(a.candidate for a in records)
    if len({c.id for c in candidates}) != len(candidates):
        raise ValueError('Duplicate candidate IDs')
    if len(records) < config.min_inliers:
        raise ValueError('Too few candidates for global lattice fit')
    positions = np.array([c.center_px for c in candidates])
    initial_indices = np.array([a.integer_indices for a in records])
    if initial_indices.shape != (len(records), 2) or initial_indices.dtype.kind not in 'iu':
        raise ValueError('Provisional coordinates must be integer pairs')
    indices = initial_indices.copy()
    residual = np.linalg.norm(positions-origin-indices@matrix.T, axis=1)
    active = _unique(np.array([a.inlier for a in records]), indices, residual, candidates)
    if active.sum() < config.min_inliers:
        raise ValueError('Too few distinct initial inlier sites')
    initial_rms = float(np.sqrt(np.mean(residual[active]**2)))
    seen = set()
    reassigned = np.zeros(len(records), dtype=bool)
    for iteration in range(1, config.max_iterations+1):
        state = (indices.tobytes(), active.tobytes())
        if state in seen:
            raise ValueError('Global lattice assignments oscillate; inspect seed/threshold')
        seen.add(state)
        design = np.column_stack((np.ones(len(indices)), indices))
        if active.sum() < config.min_inliers or np.linalg.matrix_rank(design[active]) < 3:
            raise ValueError('Too few non-collinear inlier sites for affine lattice fit')
        coefficients = np.linalg.lstsq(design[active], positions[active], rcond=None)[0]
        origin, matrix = coefficients[0], coefficients[1:].T
        if not np.isfinite(coefficients).all() or np.linalg.cond(matrix) > 1e6:
            raise ValueError('Refined lattice is singular or ill-conditioned')
        old_residual = np.linalg.norm(positions-origin-indices@matrix.T, axis=1)
        proposed = np.rint(np.linalg.solve(matrix, (positions-origin).T).T).astype(int)
        new_residual = np.linalg.norm(positions-origin-proposed@matrix.T, axis=1)
        improve = new_residual+config.reassignment_improvement_px < old_residual
        updated = np.where(improve[:, None], proposed, indices)
        reassigned |= np.any(updated != indices, axis=1)
        residual = np.linalg.norm(positions-origin-updated@matrix.T, axis=1)
        selected = _unique(residual <= config.residual_threshold_px, updated, residual, candidates)
        if np.array_equal(updated, indices) and np.array_equal(selected, active):
            break
        indices, active = updated, selected
    else:
        raise ValueError('Global lattice refinement did not converge within max_iterations')
    # Reuse A6 convention and transform labels exactly, never silently relabel HKL.
    canonical = _canonical(matrix)
    transformed = np.rint(np.linalg.solve(canonical, (indices@matrix.T).T).T).astype(int)
    residual = np.linalg.norm(positions-origin-transformed@canonical.T, axis=1)
    assignments = tuple(RefinedCandidateAssignment(c, tuple(int(n) for n in old),
                        tuple(int(n) for n in new), float(r), bool(ok), bool(changed))
        for c, old, new, r, ok, changed in zip(candidates, initial_indices, transformed, residual, active, reassigned))
    return RefinedReciprocalLattice(seed, tuple(origin), tuple(canonical[:, 0]), tuple(canonical[:, 1]),
        initial_rms, float(np.sqrt(np.mean(residual[active]**2))), assignments, iteration, config)


def enumerate_lattice_sites(lattice: RefinedReciprocalLattice, image_shape, geometry: DiskGeometry):
    """All footprint-contained integer sites; central first then radius and n1,n2."""
    if not isinstance(lattice, RefinedReciprocalLattice) or not isinstance(geometry, DiskGeometry):
        raise TypeError('Expected RefinedReciprocalLattice and DiskGeometry')
    return _enumerate_sites(lattice.origin_px, _matrix(lattice), image_shape, geometry,
                            lattice.config.max_enumeration_sites)


def enumerate_seed_sites(seed: ReciprocalLatticeSeed, image_shape, geometry: DiskGeometry, *, config=None):
    """Predictions from an unrefined seed; same ordering/coverage as the A7 catalogue."""
    if not isinstance(seed, ReciprocalLatticeSeed) or not isinstance(geometry, DiskGeometry):
        raise TypeError('Expected ReciprocalLatticeSeed and DiskGeometry')
    return _enumerate_sites(seed.reference_center_px, _matrix(seed), image_shape, geometry,
                            _config(config).max_enumeration_sites)


def _enumerate_sites(origin, matrix, image_shape, geometry, max_enumeration_sites):
    shape = tuple(_integer(n, 'image dimension') for n in image_shape)
    if len(shape) != 2 or min(shape) < 1:
        raise ValueError('image_shape must be positive (ny,nx)')
    origin = np.array(origin, dtype=float)
    lower = np.array([geometry.radius_x_px, geometry.radius_y_px])
    upper = np.array([shape[1]-1, shape[0]-1])-lower
    edge_tolerance = 64*np.finfo(float).eps*max(shape)
    if np.any(lower > upper) or np.any(origin < lower-edge_tolerance) or np.any(origin > upper+edge_tolerance):
        raise ValueError('Central predicted disk footprint must fit the frame')
    corners = np.array(list(product((lower[0], upper[0]), (lower[1], upper[1]))))
    coordinates = np.linalg.solve(matrix, (corners-origin).T).T
    lo = np.floor(coordinates.min(axis=0)).astype(int)
    hi = np.ceil(coordinates.max(axis=0)).astype(int)
    if int(hi[0]-lo[0]+1)*int(hi[1]-lo[1]+1) > max_enumeration_sites:
        raise ValueError('Enumeration exceeds max_enumeration_sites; no catalogue was truncated')
    sites = []
    for n1 in range(lo[0], hi[0]+1):
        for n2 in range(lo[1], hi[1]+1):
            point = origin+matrix@[n1,n2]
            if np.all(point >= lower-edge_tolerance) and np.all(point <= upper+edge_tolerance):
                sites.append(((n1,n2), tuple(point)))
    sites.sort(key=lambda item: (item[0] != (0,0), float(np.sum((np.array(item[1])-origin)**2)), item[0]))
    return tuple(sites)


def generate_beam_catalogue(corrected_pattern, lattice: RefinedReciprocalLattice,
                            geometry: DiskGeometry, *, valid_mask=None) -> LatticeBeamCatalogue:
    """Keep geometric predictions separate from bounded local response evidence.

    Absolute response thresholds are in input intensity units. Weak/missing sites
    retain predictions but no reliable observed center. Inadequate mask/background
    support or a local maximum on the search boundary is marked rejected.
    """
    image = _pattern(corrected_pattern)
    sites = enumerate_lattice_sites(lattice, image.shape, geometry)
    config = lattice.config
    if valid_mask is None:
        valid_mask = np.ones(image.shape, bool)
    if not isinstance(valid_mask, np.ndarray) or valid_mask.dtype != np.bool_ or valid_mask.shape != image.shape:
        raise ValueError('valid_mask must be Boolean and match image shape')
    response = disk_response_map(image, geometry)
    matrix = _matrix(lattice)
    radius = min(config.local_radius_scale*min(geometry.radius_x_px, geometry.radius_y_px),
                 .2*min(np.linalg.norm(matrix, axis=0)))
    count = max(1, int(np.ceil(radius/config.local_step_px)))
    offsets = [(float(x),float(y)) for x in np.linspace(-radius,radius,2*count+1)
               for y in np.linspace(-radius,radius,2*count+1) if x*x+y*y <= radius*radius*(1+1e-12)]
    offsets.sort(key=lambda p: (p[0]*p[0]+p[1]*p[1], p))
    beams = []
    for indices, predicted in sites:
        px,py = predicted; rx,ry = geometry.radius_x_px,geometry.radius_y_px
        # Require full annular response/search support; no edge-padding evidence.
        marginx,marginy = int(np.ceil(2*rx+radius)),int(np.ceil(2*ry+radius))
        x0,x1 = int(np.floor(px))-marginx, int(np.ceil(px))+marginx+1
        y0,y1 = int(np.floor(py))-marginy, int(np.ceil(py))+marginy+1
        support = x0 >= 0 and y0 >= 0 and x1 <= image.shape[1] and y1 <= image.shape[0]
        coverage = float(np.mean(valid_mask[max(0,y0):min(image.shape[0],y1),max(0,x0):min(image.shape[1],x1)]))
        status, observed, score = BeamStatus.REJECTED, None, None
        if support and coverage >= config.minimum_mask_coverage:
            points = np.array(offsets)+predicted
            values = map_coordinates(response, [points[:,1],points[:,0]], order=1, mode='constant', cval=0.)
            best = int(np.argmax(values)); score = float(values[best])
            if score < config.weak_response_threshold:
                status = BeamStatus.MISSING
            elif score < config.observed_response_threshold:
                status = BeamStatus.WEAK
            elif np.linalg.norm(offsets[best]) >= radius-1e-10:
                status = BeamStatus.REJECTED
            else:
                status, observed = BeamStatus.OBSERVED, tuple(points[best])
        observation = BeamObservation(f'lattice_{indices[0]}_{indices[1]}', predicted,
                                      label=f'lattice {indices[0]},{indices[1]}')
        beams.append(PredictedBeam(observation, indices, predicted, observed, indices == (0,0),
                                   status, score, coverage, radius))
    return LatticeBeamCatalogue(lattice, tuple(image.shape), geometry, tuple(beams))


def refine_lattice_and_generate_catalogue(corrected_pattern, seed, geometry, *, config=None, valid_mask=None):
    lattice = refine_reciprocal_lattice(seed, config=config)
    return generate_beam_catalogue(corrected_pattern, lattice, geometry, valid_mask=valid_mask)
