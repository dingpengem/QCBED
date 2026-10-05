"""User original-cell CBED basis primitives; no Patch A runtime association."""
from dataclasses import dataclass
from .orientation import ZoneAxis, _indices, zone_law_value
from .zone_cell import _cross
from .lattice import ReciprocalLattice, _readonly


@dataclass(frozen=True)
class OriginalReflectionBasis:
    """Ordered user-facing ORIGINAL hkl pair, independent and in the supplied zone.

    Nonprimitive bases (e.g. 002) are valid and are not silently reduced; they may
    span a sublattice. Neither basis length nor angle is orthonormalized.
    """
    hkl_1: tuple[int,int,int]
    hkl_2: tuple[int,int,int]
    zone_axis: ZoneAxis

    def __post_init__(self):
        if not isinstance(self.zone_axis,ZoneAxis):
            raise TypeError('zone_axis must be ZoneAxis')
        for name in ('hkl_1','hkl_2'):
            hkl = _indices(getattr(self,name),name)
            if hkl == (0,0,0):
                raise ValueError('basis reflections must be nonzero')
            if zone_law_value(hkl,self.zone_axis) != 0:
                raise ValueError('basis reflection is outside the selected zone')
            object.__setattr__(self,name,hkl)
        if _cross(self.hkl_1,self.hkl_2) == (0,0,0):
            raise ValueError('basis reflections must be linearly independent, not collinear')


def reflection_from_basis_coordinates(basis: OriginalReflectionBasis, n1: int, n2: int) -> tuple[int,int,int]:
    """Return original hkl=n1*hkl_1+n2*hkl_2 using exact integer algebra only."""
    if not isinstance(basis,OriginalReflectionBasis):
        raise TypeError('basis must be OriginalReflectionBasis')
    n1,n2,_ = _indices((n1,n2,0),'basis coordinates')
    return tuple(n1*h+n2*k for h,k in zip(basis.hkl_1,basis.hkl_2))


def original_reflection_basis_vectors(reciprocal_lattice: ReciprocalLattice,
                                      basis: OriginalReflectionBasis):
    """Return two immutable Cartesian reciprocal rows (Å^-1), preserving metric."""
    if not isinstance(reciprocal_lattice,ReciprocalLattice) or not isinstance(basis,OriginalReflectionBasis):
        raise TypeError('expected ReciprocalLattice and OriginalReflectionBasis')
    return _readonly([reciprocal_lattice.g_vector(basis.hkl_1),reciprocal_lattice.g_vector(basis.hkl_2)])
