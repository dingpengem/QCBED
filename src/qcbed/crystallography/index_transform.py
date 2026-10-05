"""Exact original hkl ↔ internal zone HKL conversions (row-vector conventions)."""
from dataclasses import dataclass
from .orientation import ZoneAxis, _indices, zone_law_value
from .zone_cell import ZoneAxisCell


def _multiply_transpose(indices,matrix):
    # Python ints throughout: reflection indices need not fit fixed-width arrays.
    return tuple(sum(v*int(coefficient) for v,coefficient in zip(indices,row)) for row in matrix)


def original_hkl_to_zone_HKL(hkl, zone_cell: ZoneAxisCell) -> tuple[int,int,int]:
    """HKL=hkl@T.T, derived from g=hkl@B_original=HKL@B_zone."""
    return _multiply_transpose(_indices(hkl,'original hkl'),zone_cell.transform)


def zone_HKL_to_original_hkl(HKL, zone_cell: ZoneAxisCell) -> tuple[int,int,int]:
    """hkl=HKL@inv(T).T; exact integer inverse, no floating inversion/rounding."""
    return _multiply_transpose(_indices(HKL,'zone HKL'),zone_cell.inverse_transform)


def is_zone_reflection_original(hkl, zone_axis: ZoneAxis) -> bool:
    return zone_law_value(hkl,zone_axis) == 0


def is_zolz_zone_index(HKL) -> bool:
    return _indices(HKL,'zone HKL')[2] == 0


@dataclass(frozen=True)
class ZoneReflectionHK:
    """Internal refinement shorthand ONLY: (H,K,0) in recast reciprocal rows."""
    H: int
    K: int

    def __post_init__(self):
        h,k,_ = _indices((self.H,self.K,0),'internal HK')
        object.__setattr__(self,'H',h)
        object.__setattr__(self,'K',k)

    def as_HKL(self) -> tuple[int,int,int]:
        return self.H,self.K,0


def original_hkl_to_internal_HK(hkl, zone_cell: ZoneAxisCell) -> ZoneReflectionHK:
    h,k,l = original_hkl_to_zone_HKL(hkl,zone_cell)
    if l != 0:
        raise ValueError(f'original reflection is outside the zone: recast L={l}, require L=0')
    return ZoneReflectionHK(h,k)


def internal_HK_to_original_hkl(HK: ZoneReflectionHK, zone_cell: ZoneAxisCell) -> tuple[int,int,int]:
    if not isinstance(HK,ZoneReflectionHK):
        raise TypeError('HK must be explicitly typed ZoneReflectionHK')
    return zone_HKL_to_original_hkl(HK.as_HKL(),zone_cell)
