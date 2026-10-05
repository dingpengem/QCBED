"""Integer crystal axes and right-handed Cartesian beam frames; no indexing."""
from dataclasses import dataclass, field
from functools import reduce
from math import gcd
from numbers import Integral
import numpy as np
from .lattice import DirectLattice, _readonly


def _indices(values, name):
    try:
        values = tuple(values)
    except TypeError as exc:
        raise ValueError(f'{name} requires three integers') from exc
    if len(values) != 3 or any(isinstance(v, (bool, np.bool_)) or
            not isinstance(v, Integral) for v in values):
        raise ValueError(f'{name} requires three integers excluding bool')
    return tuple(int(v) for v in values)


@dataclass(frozen=True)
class ZoneAxis:
    """Unoriented integer axis: gcd reduced, first nonzero entry positive.

    Opposite incidence directions deliberately canonicalize to the same axis.
    A future signed beam-direction API must be separate.
    """
    uvw: tuple[int, int, int]

    def __post_init__(self):
        uvw = _indices(self.uvw, 'uvw')
        divisor = reduce(gcd, uvw)
        if divisor == 0:
            raise ValueError('zone axis must not be (0,0,0)')
        sign = 1 if next(v for v in uvw if v) > 0 else -1
        object.__setattr__(self, 'uvw', tuple(sign*v//divisor for v in uvw))


def zone_axis_cartesian(lattice: DirectLattice, zone_axis: ZoneAxis) -> np.ndarray:
    """Normalize u*a+v*b+w*c (direct rows in Å); return Cartesian unit vector."""
    if not isinstance(lattice, DirectLattice) or not isinstance(zone_axis, ZoneAxis):
        raise TypeError('lattice/zone_axis must be DirectLattice/ZoneAxis')
    vector = np.asarray(zone_axis.uvw, dtype=float) @ lattice.vectors
    scale = np.max(np.abs(vector))
    if not np.isfinite(vector).all() or scale == 0:
        raise ValueError('zone direction must be finite and nonzero')
    vector = vector/scale
    return _readonly(vector/np.linalg.norm(vector))


def _vectors(values):
    a = np.asarray(values)
    if a.ndim < 1 or a.shape[-1] != 3 or a.dtype.kind not in 'iuf' or not np.isfinite(a).all():
        raise ValueError('vectors must be finite real arrays of shape (...,3)')
    return np.asarray(a, dtype=np.float64)


@dataclass(frozen=True, eq=False)
class CrystalOrientation:
    """Rows R=(x_beam,y_beam,z_beam); row vectors transform as v @ R.T.

    Choose the least-parallel Cartesian axis (x,y,z priority on ties within
    32 float64 eps). x_beam=normalize(reference cross z); y_beam=z cross x.
    Tie band prevents roundoff-driven branch changes at exact symmetric axes.
    This deterministic piecewise frame is not globally continuous across selector
    boundaries, nor is it a crystallographic integer reciprocal basis.
    """
    direct_lattice: DirectLattice
    zone_axis: ZoneAxis
    beam_direction_cartesian: np.ndarray = field(init=False)
    beam_basis: np.ndarray = field(init=False)
    reference_axis_index: int = field(init=False)

    def __post_init__(self):
        z = zone_axis_cartesian(self.direct_lattice, self.zone_axis)
        parallel = np.abs(z)
        index = int(np.flatnonzero(parallel <= parallel.min()+32*np.finfo(float).eps)[0])
        reference = np.eye(3)[index]
        x = np.cross(reference, z); x /= np.linalg.norm(x)
        y = np.cross(z, x); y /= np.linalg.norm(y)
        object.__setattr__(self, 'beam_direction_cartesian', z)
        object.__setattr__(self, 'beam_basis', _readonly(np.array([x,y,z])))
        object.__setattr__(self, 'reference_axis_index', index)

    def cartesian_to_beam(self, values) -> np.ndarray:
        """Transform (...,3) Cartesian row vectors; retain their original units."""
        return _readonly(_vectors(values) @ self.beam_basis.T)

    def beam_to_cartesian(self, values) -> np.ndarray:
        return _readonly(_vectors(values) @ self.beam_basis)


def zone_law_value(hkl, zone_axis: ZoneAxis) -> int:
    """Exact integer h*u+k*v+l*w; no enumeration or tolerance-based filtering."""
    if not isinstance(zone_axis, ZoneAxis):
        raise TypeError('zone_axis must be ZoneAxis')
    return sum(h*u for h,u in zip(_indices(hkl, 'hkl'), zone_axis.uvw))
