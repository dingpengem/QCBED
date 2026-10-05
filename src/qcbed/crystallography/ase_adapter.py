"""ASE is an input frontend only; no ASE object is retained in canonical models."""
from typing import TYPE_CHECKING
from .lattice import DirectLattice
from .structure import AtomicSite, AtomicStructureModel

if TYPE_CHECKING:
    from ase import Atoms


def atomic_structure_from_ase(atoms: 'Atoms') -> AtomicStructureModel:
    """Copy cell, symbols and unwrapped fractional positions, preserving atom order.

    A full, positive-volume cell is mandatory even for nonperiodic ASE inputs.
    PBC flags, calculators, velocities and arbitrary info/arrays are not imported.
    Occupancy defaults to 1 and B to None: no ASE metadata convention is inferred.
    Disordered CIF occupancies need explicit canonical sites, not this basic adapter.
    """
    from ase import Atoms
    if not isinstance(atoms, Atoms):
        raise TypeError('atoms must be ase.Atoms')
    lattice = DirectLattice(atoms.cell.array)
    positions = atoms.get_scaled_positions(wrap=False)
    sites = tuple(AtomicSite(symbol, tuple(position)) for symbol, position in
                  zip(atoms.get_chemical_symbols(), positions, strict=True))
    return AtomicStructureModel(lattice, sites)
