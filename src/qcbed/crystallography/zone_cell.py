"""B3 exact unimodular direct-cell completion; not a Cartesian rotation."""
from dataclasses import dataclass, field, replace
import numpy as np
from .lattice import DirectLattice, ReciprocalLattice, _readonly
from .orientation import ZoneAxis, _vectors
from .structure import AtomicStructureModel


def _bezout(a, b):
    """Deterministic Euclidean divisions on magnitudes; return g,s,t: s*a+t*b=g."""
    r0,r1 = abs(a),abs(b)
    s0,s1,t0,t1 = 1,0,0,1
    while r1:
        q = r0//r1
        r0,r1 = r1,r0-q*r1
        s0,s1 = s1,s0-q*s1
        t0,t1 = t1,t0-q*t1
    return r0,s0*(1 if a >= 0 else -1),t0*(1 if b >= 0 else -1)


def _cross(a,b):
    return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])


def _det(rows):
    return sum(a*b for a,b in zip(rows[0],_cross(rows[1],rows[2])))


def _integer_array(rows):
    # Exact algebra precedes storage conversion. Never allow integer overflow.
    info = np.iinfo(np.int64)
    if any(v < info.min or v > info.max for row in rows for v in row):
        raise OverflowError('zone-cell transform exceeds immutable int64 storage; no truncation applied')
    a = np.array(rows,dtype=np.int64)
    return np.frombuffer(a.tobytes(),dtype=np.int64).reshape(3,3)


def _completion(zone_axis):
    u,v,w = zone_axis.uvw
    g,s,t = _bezout(u,v)
    if g == 0:  # Canonical primitive axis must be [001].
        return ((1,0,0),(0,1,0),(0,0,1))
    common,alpha,beta = _bezout(g,w)
    if common != 1:
        raise ValueError('zone axis must be primitive')
    # cross(p,q)=(s*alpha,t*alpha,beta), so det=alpha*g+beta*w=1.
    return ((t,-s,0),(beta*(u//g),beta*(v//g),-alpha),(u,v,w))


@dataclass(frozen=True, eq=False)
class ZoneAxisCell:
    """Derived immutable zone cell retaining its original lattice provenance.

    A_zone=T@A_original; det(T)=+1; T[2]=canonical [uvw].
    Constructor computes all derived fields; inconsistent injected arrays cannot
    be supplied. Original hkl remains the public physical reflection identity.
    """
    original_lattice: DirectLattice
    zone_axis: ZoneAxis
    transform: np.ndarray = field(init=False)
    inverse_transform: np.ndarray = field(init=False)
    direct_lattice: DirectLattice = field(init=False)
    reciprocal_lattice: ReciprocalLattice = field(init=False)

    def __post_init__(self):
        if not isinstance(self.original_lattice,DirectLattice) or not isinstance(self.zone_axis,ZoneAxis):
            raise TypeError('expected DirectLattice and ZoneAxis')
        rows = _completion(self.zone_axis)
        if _det(rows) != 1 or rows[2] != self.zone_axis.uvw:
            raise ValueError('invalid zone completion: require det=+1 and exact third row')
        columns = (_cross(rows[1],rows[2]),_cross(rows[2],rows[0]),_cross(rows[0],rows[1]))
        inverse = tuple(zip(*columns))  # Integer adjugate, determinant exactly 1.
        transform = _integer_array(rows)
        inverse_transform = _integer_array(inverse)
        direct = DirectLattice(transform @ self.original_lattice.vectors)
        object.__setattr__(self,'transform',transform)
        object.__setattr__(self,'inverse_transform',inverse_transform)
        object.__setattr__(self,'direct_lattice',direct)
        object.__setattr__(self,'reciprocal_lattice',direct.reciprocal())


def build_zone_axis_cell(lattice: DirectLattice, zone_axis: ZoneAxis) -> ZoneAxisCell:
    return ZoneAxisCell(lattice,zone_axis)


def original_fractional_to_zone(values, zone_cell: ZoneAxisCell) -> np.ndarray:
    """f_zone=f_original@inv(T); float64, unwrapped, shape (...,3)."""
    return _readonly(_vectors(values) @ zone_cell.inverse_transform)


def zone_fractional_to_original(values, zone_cell: ZoneAxisCell) -> np.ndarray:
    """f_original=f_zone@T; no wrapping or periodic site merging."""
    return _readonly(_vectors(values) @ zone_cell.transform)


def recast_structure_to_zone_cell(structure: AtomicStructureModel,
                                  zone_cell: ZoneAxisCell) -> AtomicStructureModel:
    """Pure basis change; preserve site order, count, metadata and Cartesian positions."""
    if not isinstance(structure,AtomicStructureModel) or not isinstance(zone_cell,ZoneAxisCell):
        raise TypeError('expected AtomicStructureModel and ZoneAxisCell')
    if not np.array_equal(structure.direct_lattice.vectors,zone_cell.original_lattice.vectors):
        raise ValueError('structure must match zone cell original lattice')
    positions = original_fractional_to_zone([s.fractional_position for s in structure.sites],zone_cell)
    return AtomicStructureModel(zone_cell.direct_lattice,tuple(
        replace(s,fractional_position=tuple(p)) for s,p in zip(structure.sites,positions,strict=True)))
