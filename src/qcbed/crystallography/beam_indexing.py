"""Explicit experimental integer axes -> original hkl -> internal zone HK.

No metric fitting, sign/swap search, detector calibration or reflection search.
"""
from dataclasses import dataclass
from .orientation import _indices
from .reflection_basis import OriginalReflectionBasis, reflection_from_basis_coordinates
from .zone_cell import ZoneAxisCell
from .index_transform import ZoneReflectionHK, original_hkl_to_internal_HK


def _experimental_index(values):
    try:
        values = tuple(values)
    except TypeError as exc:
        raise ValueError('experimental_index must contain two integers') from exc
    if len(values) != 2:
        raise ValueError('experimental_index must contain two integers (n1,n2)')
    return _indices((*values,0),'experimental_index')[:2]


def _stable_id(value):
    if not isinstance(value,str) or not value.strip():
        raise ValueError('stable_id must be nonempty text')
    return value


def _validate_context(basis, zone_cell):
    if not isinstance(basis,OriginalReflectionBasis) or not isinstance(zone_cell,ZoneAxisCell):
        raise TypeError('expected OriginalReflectionBasis and ZoneAxisCell')
    if basis.zone_axis != zone_cell.zone_axis:
        raise ValueError('reflection basis zone_axis must match zone_cell.zone_axis')


@dataclass(frozen=True)
class IndexedBeamReflection:
    """All identities retained; original_hkl is the user-facing reflection label."""
    stable_id: str
    experimental_index: tuple[int,int]
    original_hkl: tuple[int,int,int]
    internal_hk: ZoneReflectionHK

    def __post_init__(self):
        _stable_id(self.stable_id)
        object.__setattr__(self,'experimental_index',_experimental_index(self.experimental_index))
        object.__setattr__(self,'original_hkl',_indices(self.original_hkl,'original_hkl'))
        if not isinstance(self.internal_hk,ZoneReflectionHK):
            raise TypeError('internal_hk must be ZoneReflectionHK')


def index_experimental_beam(*, stable_id: str, experimental_index,
                            reflection_basis: OriginalReflectionBasis,
                            zone_cell: ZoneAxisCell) -> IndexedBeamReflection:
    """Authoritative g1↔hkl_1, g2↔hkl_2; no automatic sign/basis correction."""
    _validate_context(reflection_basis,zone_cell)
    indices = _experimental_index(experimental_index)
    hkl = reflection_from_basis_coordinates(reflection_basis,*indices)
    # The direct beam follows the same exact algebra; no d-spacing is requested.
    internal = original_hkl_to_internal_HK(hkl,zone_cell)
    return IndexedBeamReflection(stable_id,indices,hkl,internal)


def index_experimental_beams(beams, reflection_basis: OriginalReflectionBasis,
                             zone_cell: ZoneAxisCell) -> tuple[IndexedBeamReflection,...]:
    """Map an ordered iterable of (stable_id,(n1,n2)); preserve duplicates of indices.

    Stable IDs must be unique. Empty input produces an empty tuple. No sorting or
    filtering is performed; geometric statuses are the upstream caller's concern.
    """
    _validate_context(reflection_basis,zone_cell)
    result = tuple(index_experimental_beam(stable_id=id,experimental_index=indices,
        reflection_basis=reflection_basis,zone_cell=zone_cell) for id,indices in beams)
    if len({r.stable_id for r in result}) != len(result):
        raise ValueError('experimental beams must have unique stable IDs')
    return result
