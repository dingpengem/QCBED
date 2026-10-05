"""Input-table-free initial pixel lattice inference, not crystallographic indexing."""
from dataclasses import dataclass
from itertools import combinations

import numpy as np
from scipy.signal import fftconvolve

from .models import DiskGeometry, _pair
from .noise import _pattern
from .disk_detection import _center
from .arc_geometry import (LatticeDetectionConfig, LatticeCandidate, _config,
                           arc_edge_support, arc_center_votes, _validated_center)
from ..detector.models import _integer, _real


@dataclass(frozen=True)
class LatticeAssignment:
    candidate: LatticeCandidate
    integer_indices: tuple[int, int]  # abstract n1,n2, never HKL
    residual_px: float
    inlier: bool


@dataclass(frozen=True)
class DisplacementFamily:
    vector_px: tuple[float, float]
    votes: int


@dataclass(frozen=True)
class ReciprocalLatticeSeed:
    reference_center_px: tuple[float, float]
    basis_vector_1_px: tuple[float, float]
    basis_vector_2_px: tuple[float, float]
    assignments: tuple[LatticeAssignment, ...]
    displacement_families: tuple[DisplacementFamily, ...]
    support_fraction: float
    fit_rms_residual_px: float
    config: LatticeDetectionConfig
    seed_method: str = 'multi_candidate_displacement_voting'
    source_candidates: tuple[LatticeCandidate, ...] = ()
    initial_basis_vectors_px: tuple[tuple[float, float], ...] = ()
    basis_angle_deg: float | None = None
    basis_condition_number: float | None = None
    primitive_basis_certain: bool = False


def disk_response_map(corrected_pattern, geometry: DiskGeometry) -> np.ndarray:
    """A7 local-evidence response; not the canonical A6 candidate detector.

    Inner weight (1-r²)² for r<1; annulus sin²(pi*(r-1)) for 1<r<2.
    Each lobe sums to one: DC/linear background cancels in the full-support interior.
    Zero-padded correlation is used; candidate selection excludes truncated windows.
    """
    image = _pattern(corrected_pattern)
    if not isinstance(geometry, DiskGeometry):
        raise TypeError('geometry must be DiskGeometry')
    rx, ry = geometry.radius_x_px, geometry.radius_y_px
    if 4*rx+1 > image.shape[1] or 4*ry+1 > image.shape[0]:
        raise ValueError('Disk response footprint and background annulus must fit the frame')
    y, x = np.mgrid[-int(np.ceil(2*ry)):int(np.ceil(2*ry))+1,
                    -int(np.ceil(2*rx)):int(np.ceil(2*rx))+1]
    radius = np.sqrt((x/rx)**2+(y/ry)**2)
    inner = np.where(radius < 1, (1-radius**2)**2, 0.)
    outer = np.where((radius > 1) & (radius < 2), np.sin(np.pi*(radius-1))**2, 0.)
    if outer.sum() <= 0:
        raise ValueError('Disk footprint too small to sample background annulus')
    response = fftconvolve(image, inner/inner.sum()-outer/outer.sum(), mode='same')
    if not np.isfinite(response).all():
        raise ValueError('Non-finite disk response')
    return response


def detect_lattice_candidates(corrected_pattern, geometry: DiskGeometry,
                              central_disk_center_px, *, config=None, subpixel_config=None) -> tuple[LatticeCandidate, ...]:
    """Validate combined known-radius support before reciprocal-lattice inference."""
    config = _config(config)
    votes, arrays, diagnostics = arc_center_votes(corrected_pattern, geometry, central_disk_center_px, config=config)
    peaks = np.argwhere(arrays['proposal_peaks'])
    validated=[]
    for iy,ix in peaks:
        center=(float(ix),float(iy))
        metrics=_validated_center(np.array(center),arrays,geometry,config)
        if metrics is not None:
            validated.append(LatticeCandidate(f'arc_{ix}_{iy}',center,float(votes[iy,ix]),**metrics))
    # Confidence uses geometry alone. Intensity/vote strength cannot rescue a rejected arc.
    validated.sort(key=lambda c:(-c.confidence,c.center_px,c.id))
    selected=[]
    for candidate in validated:
        if any(np.linalg.norm(np.array(candidate.center_px)-c.center_px)<config.candidate_merge_radius_px for c in selected):
            continue
        selected.append(candidate)
        if len(selected)==config.max_candidates:
            break
    if subpixel_config is not None:
        from .subpixel_arc import refine_validated_candidate
        selected=[refine_validated_candidate(c,arrays,geometry,config,config=subpixel_config) for c in selected]
    return tuple(sorted(selected,key=lambda c:(c.center_px,c.id)))


def _positive(vector):
    return -vector if vector[0] < 0 or (vector[0] == 0 and vector[1] < 0) else vector


def _families(positions, tolerance, min_length):
    vectors = [_positive(b-a) for a, b in combinations(positions, 2) if np.linalg.norm(b-a) >= min_length]
    vectors.sort(key=lambda v: (float(np.linalg.norm(v)), *v))
    groups = []
    for vector in vectors:
        matches = [(min(np.linalg.norm(vector-center), np.linalg.norm(vector+center)), i)
                   for i, (center, count) in enumerate(groups)]
        distance, index = min(matches, default=(np.inf, -1))
        if distance <= tolerance:
            center, count = groups[index]
            aligned = vector if vector@center >= 0 else -vector
            groups[index] = ((center*count+aligned)/(count+1), count+1)
        else:
            groups.append((vector.copy(), 1))
    return groups


def _canonical(matrix):
    """Gauss-reduced short basis; sign of first leading x/y positive, det positive."""
    a, b = matrix.T.copy()
    for _ in range(64):
        if np.dot(b, b) < np.dot(a, a):
            a, b = b, a
        multiple = int(np.rint(np.dot(a, b)/np.dot(a, a)))
        if multiple == 0:
            break
        b = b-multiple*a
    a = _positive(a)
    if np.linalg.det(np.column_stack((a, b))) < 0:
        b = -b
    return np.column_stack((a, b))


def _assign(matrix, relative, tolerance):
    indices = np.rint(np.linalg.solve(matrix, relative.T).T)
    residual = np.linalg.norm(relative-indices@matrix.T, axis=1)
    inliers = (residual <= tolerance) & np.any(indices != 0, axis=1)
    # Two candidate maxima may not claim the same lattice site as independent support.
    used = set()
    for i in np.argsort(residual, kind='stable'):
        key = tuple(indices[i])
        if key in used:
            inliers[i] = False
        elif inliers[i]:
            used.add(key)
    return indices, residual, inliers


def _basis_diagnostics(matrix):
    lengths=np.linalg.norm(matrix,axis=0)
    angle=float(np.rad2deg(np.arccos(np.clip(matrix[:,0]@matrix[:,1]/np.prod(lengths),-1,1))))
    return angle,float(np.linalg.cond(matrix))


def _minimal_arc_seed(candidates, center, geometry, config):
    """Two real arc observations initialize geometry, never prove primitiveness."""
    for c in candidates:
        if (c.supporting_edge_pixels<config.minimum_support_pixels or
            len(set(c.supporting_pixels))!=c.supporting_edge_pixels or
            c.angular_coverage<config.minimum_angular_coverage or
            c.arc_length_px<config.minimum_arc_length_px or c.fragment_count<1 or
            c.confidence<=0 or c.mean_normal_residual_deg>config.maximum_normal_residual_deg or
            c.radius_residual_rms_px>config.radius_tolerance_px):
            raise ValueError('Minimal seed requires two geometrically validated arc candidates with support evidence')
        displacement=np.array(c.supporting_pixels)-c.center_px
        radial=np.linalg.norm(displacement/np.array([geometry.radius_x_px,geometry.radius_y_px]),axis=1)
        distance=np.linalg.norm(displacement,axis=1)
        expected=np.divide(distance,radial,out=np.full_like(distance,np.inf),where=radial>0)
        if np.any(abs(distance-expected)>config.radius_tolerance_px+1e-12):
            raise ValueError('Minimal seed arc support violates the physical radius prior')
    initial=np.column_stack([np.array(c.center_px)-center for c in candidates])
    lengths=np.linalg.norm(initial,axis=0)
    minimum=config.minimum_lattice_vector_scale*min(geometry.radius_x_px,geometry.radius_y_px)
    if np.any(lengths<minimum):
        raise ValueError('Minimal seed displacement is too short for the physical disk geometry')
    sine=abs(np.linalg.det(initial))/np.prod(lengths)
    if sine<config.min_basis_sine:
        raise ValueError('Minimal seed requires non-collinear displacement vectors')
    if not np.isfinite(initial).all() or np.linalg.cond(initial)>1e6:
        raise ValueError('Minimal seed basis is singular or ill-conditioned')
    matrix=_canonical(initial)
    if min(np.linalg.norm(matrix,axis=0))<minimum or np.linalg.cond(matrix)>1e6:
        raise ValueError('Minimal reduced basis is incompatible with disk geometry or conditioning')
    indices,residual,inliers=_assign(matrix,initial.T,config.lattice_residual_tolerance_px)
    if not inliers.all():
        raise ValueError('Minimal seed canonicalization did not preserve both observed sites')
    assignments=tuple(LatticeAssignment(c,tuple(int(n) for n in index),float(r),True)
                      for c,index,r in zip(candidates,indices,residual))
    angle,condition=_basis_diagnostics(matrix)
    return ReciprocalLatticeSeed(tuple(center),tuple(matrix[:,0]),tuple(matrix[:,1]),assignments,(),
        1.,float(np.sqrt(np.mean(residual**2))),config,
        seed_method='two_validated_arc_centers',source_candidates=candidates,
        initial_basis_vectors_px=tuple(tuple(v) for v in initial.T),
        basis_angle_deg=angle,basis_condition_number=condition,primitive_basis_certain=False)


def infer_lattice_seed(candidates, central_disk_center_px, geometry: DiskGeometry, *, config=None) -> ReciprocalLatticeSeed:
    """Fit observed candidates only; no legacy tables, origin refinement or HKL."""
    config = _config(config)
    center = _pair(central_disk_center_px, 'central_disk_center_px')
    if not isinstance(geometry, DiskGeometry):
        raise TypeError('geometry must be DiskGeometry')
    candidates = tuple(candidates)
    if any(not isinstance(c, LatticeCandidate) for c in candidates):
        raise TypeError('candidates must contain LatticeCandidate objects')
    candidates = tuple(sorted(candidates, key=lambda c: (c.center_px, c.id)))
    if len({c.id for c in candidates}) != len(candidates):
        raise ValueError('Duplicate candidate IDs')
    if len(candidates)==2:
        return _minimal_arc_seed(candidates,np.array(center),geometry,config)
    if len(candidates) < config.min_inliers:
        raise ValueError('Too few diffracted candidates for lattice support')
    relative = np.array([c.center_px for c in candidates])-center
    minimum = config.minimum_lattice_vector_scale*min(geometry.radius_x_px, geometry.radius_y_px)
    families = _families(np.vstack((np.zeros(2), relative)), config.displacement_tolerance_px, minimum)
    families = [(v, n) for v, n in families if n >= config.min_family_votes]
    families.sort(key=lambda item: (-item[1], np.linalg.norm(item[0]), *item[0]))
    families = families[:config.max_families]
    best = None
    for (a, va), (b, vb) in combinations(families, 2):
        matrix = np.column_stack((a, b))
        initial_hypothesis = matrix.copy()
        if abs(np.linalg.det(matrix))/(np.linalg.norm(a)*np.linalg.norm(b)) < config.min_basis_sine:
            continue
        matrix = _canonical(matrix)
        for _ in range(3):
            if not np.isfinite(matrix).all() or np.linalg.cond(matrix) > 1e6:
                break
            indices, residual, inliers = _assign(matrix, relative, config.lattice_residual_tolerance_px)
            if inliers.sum() < config.min_inliers or np.linalg.matrix_rank(indices[inliers]) < 2:
                break
            matrix = np.linalg.lstsq(indices[inliers], relative[inliers], rcond=None)[0].T
        if not np.isfinite(matrix).all() or np.linalg.cond(matrix) > 1e6:
            continue
        matrix = _canonical(matrix)
        indices, residual, inliers = _assign(matrix, relative, config.lattice_residual_tolerance_px)
        count = int(inliers.sum())
        if count < config.min_inliers or count/len(candidates) < config.min_inlier_fraction:
            continue
        lengths = np.linalg.norm(matrix, axis=0)
        if min(lengths) < minimum or abs(np.linalg.det(matrix))/np.prod(lengths) < config.min_basis_sine:
            continue
        rms = float(np.sqrt(np.mean(residual[inliers]**2)))
        key = (count, -float(sum(lengths)), va+vb, -rms, *(-matrix.ravel()))
        if best is None or key > best[0]:
            best = (key, matrix, indices, residual, inliers, rms, initial_hypothesis)
    if best is None:
        raise ValueError('No stable non-collinear lattice basis explains enough candidates')
    _, matrix, indices, residual, inliers, rms, initial_hypothesis = best
    assignments = tuple(LatticeAssignment(c, tuple(int(n) for n in index), float(r), bool(ok))
                        for c, index, r, ok in zip(candidates, indices, residual, inliers))
    return ReciprocalLatticeSeed(center, tuple(matrix[:, 0]), tuple(matrix[:, 1]), assignments,
        tuple(DisplacementFamily(tuple(_positive(v)), n) for v, n in families),
        float(np.mean(inliers)), rms, config,
        source_candidates=candidates,initial_basis_vectors_px=tuple(tuple(v) for v in initial_hypothesis.T),
        basis_angle_deg=_basis_diagnostics(matrix)[0],basis_condition_number=_basis_diagnostics(matrix)[1])


def detect_reciprocal_lattice_seed(corrected_pattern, geometry: DiskGeometry,
                                   central_disk_center_px, *, config=None, subpixel_config=None) -> ReciprocalLatticeSeed:
    """A3 whole pattern + A5 footprint/center -> initial reciprocal pixel seed."""
    config = _config(config)
    candidates = detect_lattice_candidates(corrected_pattern, geometry, central_disk_center_px, config=config, subpixel_config=subpixel_config)
    return infer_lattice_seed(candidates, central_disk_center_px, geometry, config=config)
