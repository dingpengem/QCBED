"""Atomic and supplied-HKL projections in a Cartesian beam frame, not indexing."""
from dataclasses import dataclass
import numpy as np
from .lattice import ReciprocalLattice, _readonly
from .structure import AtomicStructureModel, AtomicSite, _real
from .orientation import CrystalOrientation, _indices
from .replication import replica_offsets, replicated_cartesian_positions


@dataclass(frozen=True)
class ProjectedAtomicSite:
    source_site_index: int
    cell_offset: tuple[int,int,int]
    source_site: AtomicSite
    x: float
    y: float
    depth: float

    def __post_init__(self):
        index = _indices((self.source_site_index,0,0), 'source_site_index')[0]
        if index < 0 or not isinstance(self.source_site, AtomicSite):
            raise ValueError('projected site requires nonnegative source index and AtomicSite')
        object.__setattr__(self, 'source_site_index', index)
        object.__setattr__(self, 'cell_offset', _indices(self.cell_offset, 'cell_offset'))
        for name in ('x','y','depth'):
            object.__setattr__(self, name, _real(getattr(self,name), name))

    @property
    def symbol(self):
        return self.source_site.symbol


def project_atomic_structure(structure: AtomicStructureModel, orientation: CrystalOrientation,
                             repeats=(1,1,1)) -> tuple[ProjectedAtomicSite, ...]:
    """Unmerged sites in cell-major/site-major order; x,y,depth in Å.

    Positive direct-cell repeats are not a commensurate zone-axis supercell.
    """
    if not isinstance(structure, AtomicStructureModel) or not isinstance(orientation, CrystalOrientation):
        raise TypeError('structure/orientation must be AtomicStructureModel/CrystalOrientation')
    if not np.array_equal(structure.direct_lattice.vectors, orientation.direct_lattice.vectors):
        raise ValueError('orientation must use the same direct lattice as the structure')
    offsets = replica_offsets(repeats)
    positions = orientation.cartesian_to_beam(replicated_cartesian_positions(structure, offsets))
    return tuple(ProjectedAtomicSite(i, offset, site, *positions[cell*len(structure.sites)+i])
                 for cell,offset in enumerate(offsets) for i,site in enumerate(structure.sites))


@dataclass(frozen=True, eq=False)
class ProjectedReciprocalVector:
    hkl: tuple[int,int,int]
    g_cartesian: np.ndarray
    g_beam: np.ndarray

    def __post_init__(self):
        object.__setattr__(self, 'hkl', _indices(self.hkl, 'hkl'))
        for name in ('g_cartesian','g_beam'):
            raw = np.asarray(getattr(self,name))
            if raw.shape != (3,) or raw.dtype.kind not in 'iuf' or not np.isfinite(raw).all():
                raise ValueError(f'{name} must be a finite real three-vector')
            object.__setattr__(self, name, _readonly(raw))

    @property
    def g_projected_xy(self):
        return self.g_beam[:2]

    @property
    def beam_parallel_component(self):
        return float(self.g_beam[2])


def project_hkl(reciprocal_lattice: ReciprocalLattice, orientation: CrystalOrientation,
                hkl) -> ProjectedReciprocalVector:
    """Project one explicitly supplied physical HKL; units Å^-1, no 2π.

    The supplied reciprocal lattice must be dual to the orientation's direct cell
    when interpreting the zone law; no reflections are generated or selected here.
    """
    if not isinstance(reciprocal_lattice, ReciprocalLattice) or not isinstance(orientation, CrystalOrientation):
        raise TypeError('expected ReciprocalLattice and CrystalOrientation')
    indices = _indices(hkl, 'hkl')
    g = reciprocal_lattice.g_vector(indices)
    return ProjectedReciprocalVector(indices,g,orientation.cartesian_to_beam(g))
