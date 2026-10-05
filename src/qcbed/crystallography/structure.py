"""Immutable fractional atomic structure; independent of ASE and visualization."""
from dataclasses import dataclass
from numbers import Real
import numpy as np
from .lattice import DirectLattice

# IUPAC symbols, atomic numbers 1–118. No dummy atoms/isotope aliases.
_ELEMENTS = tuple(('H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn '
    'Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd '
    'Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf '
    'Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og').split())
_SYMBOLS = frozenset(_ELEMENTS)


def _real(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real) or not np.isfinite(value):
        raise ValueError(f'{name} must be a finite real number excluding bool')
    return float(value)


@dataclass(frozen=True)
class AtomicSite:
    """Unwrapped fractional position; occupancy fraction; optional isotropic B in Å².

    B is stored only: no B/U conversion or scattering behavior is implied.
    """
    symbol: str
    fractional_position: tuple[float, float, float]
    occupancy: float = 1.
    debye_waller_b: float | None = None
    physical_group_id: str | None = None
    original_site_index: int | None = None

    def __post_init__(self):
        if not isinstance(self.symbol, str) or self.symbol not in _SYMBOLS:
            raise ValueError('symbol must be a valid chemical element symbol')
        if self.physical_group_id is not None and (not isinstance(self.physical_group_id,str) or not self.physical_group_id.strip()):
            raise ValueError('physical_group_id must be a nonempty string')
        if self.original_site_index is not None:
            from numbers import Integral
            if isinstance(self.original_site_index,(bool,np.bool_)) or not isinstance(self.original_site_index,Integral) or self.original_site_index < 0:
                raise ValueError('original_site_index must be a nonnegative integer')
            object.__setattr__(self,'original_site_index',int(self.original_site_index))
        values = tuple(self.fractional_position)
        if len(values) != 3:
            raise ValueError('fractional_position must contain three coordinates')
        object.__setattr__(self, 'fractional_position', tuple(_real(v, 'fractional coordinate') for v in values))
        occupancy = _real(self.occupancy, 'occupancy')
        if not 0 < occupancy <= 1:
            raise ValueError('occupancy must satisfy 0 < occupancy <= 1')
        object.__setattr__(self, 'occupancy', occupancy)
        if self.debye_waller_b is not None:
            b = _real(self.debye_waller_b, 'debye_waller_b')
            if b < 0:
                raise ValueError('debye_waller_b must be nonnegative in Å²')
            object.__setattr__(self, 'debye_waller_b', b)


@dataclass(frozen=True, eq=False)
class AtomicStructureModel:
    """Ordered atomic sites, including duplicates; no merging or coordinate wrapping."""
    direct_lattice: DirectLattice
    sites: tuple[AtomicSite, ...]

    def __post_init__(self):
        if not isinstance(self.direct_lattice, DirectLattice):
            raise TypeError('direct_lattice must be DirectLattice')
        sites = tuple(self.sites)
        if not sites or any(not isinstance(s, AtomicSite) for s in sites):
            raise ValueError('sites must contain at least one AtomicSite')
        object.__setattr__(self, 'sites', sites)
