"""Right-handed 3-D row-vector lattices; Å and Å^-1, without 2π."""
from dataclasses import dataclass
from numbers import Integral
import numpy as np


def _readonly(values):
    a = np.asarray(values, dtype=np.float64)
    return np.frombuffer(a.tobytes(), dtype=np.float64).reshape(a.shape)


def _cell(values):
    raw = np.asarray(values)
    if raw.shape != (3, 3):
        raise ValueError('lattice vectors must have shape (3,3), with vectors as rows')
    if raw.dtype.kind not in 'iuf':
        raise TypeError('lattice vectors must be real numbers')
    a = _readonly(raw)
    if not np.isfinite(a).all():
        raise ValueError('lattice vectors must be finite')
    volume = float(np.linalg.det(a))
    if not np.isfinite(volume) or volume <= 0:
        raise ValueError('lattice must be nonsingular and right-handed with positive finite volume')
    try:
        inverse = np.linalg.inv(a)
    except np.linalg.LinAlgError as exc:
        raise ValueError('lattice must be nonsingular') from exc
    if not np.isfinite(inverse).all():
        raise ValueError('lattice inverse must be finite')
    return a


@dataclass(frozen=True, eq=False)
class DirectLattice:
    """Rows a,b,c in Å. No rotation, cell reduction or handedness repair."""
    vectors: np.ndarray

    def __post_init__(self):
        object.__setattr__(self, 'vectors', _cell(self.vectors))

    @property
    def volume(self) -> float:
        """Signed (validated positive) volume in Å^3."""
        return float(np.linalg.det(self.vectors))

    @property
    def metric_tensor(self) -> np.ndarray:
        return _readonly(self.vectors @ self.vectors.T)

    @property
    def lengths(self) -> np.ndarray:
        """Lengths (a,b,c) in Å."""
        return _readonly(np.linalg.norm(self.vectors, axis=1))

    @property
    def angles(self) -> np.ndarray:
        """(alpha=b∠c, beta=a∠c, gamma=a∠b), in degrees."""
        unit = self.vectors / self.lengths[:, None]
        cosines = [unit[i] @ unit[j] for i,j in ((1,2),(0,2),(0,1))]
        # Clip only roundoff in derived cosine, never user coordinates.
        return _readonly(np.degrees(np.arccos(np.clip(cosines, -1., 1.))))

    def reciprocal(self) -> 'ReciprocalLattice':
        return ReciprocalLattice(np.linalg.inv(self.vectors).T)


@dataclass(frozen=True, eq=False)
class ReciprocalLattice:
    """Rows a*,b*,c* in Å^-1; B @ A.T = I. No 2π."""
    vectors: np.ndarray

    def __post_init__(self):
        object.__setattr__(self, 'vectors', _cell(self.vectors))

    def g_vector(self, hkl) -> np.ndarray:
        """Cartesian h*a*+k*b*+l*c* in Å^-1, for integer physical HKL."""
        values = tuple(hkl)
        if len(values) != 3 or any(isinstance(v, (bool, np.bool_)) or
                not isinstance(v, Integral) for v in values):
            raise ValueError('hkl must contain three integers excluding bool')
        return _readonly(np.asarray(values, dtype=np.float64) @ self.vectors)

    def d_spacing(self, hkl) -> float:
        """Plane spacing in Å; the origin reflection has no finite spacing."""
        g = self.g_vector(hkl)
        norm = float(np.linalg.norm(g))
        if norm == 0:
            raise ValueError('d_spacing requires a nonzero reflection; (0,0,0) is invalid')
        return 1. / norm
