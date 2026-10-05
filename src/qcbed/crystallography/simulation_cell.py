"""Strict B5 transverse-periodicity proof and explicit physical-depth slicing.

No projected-full-period fallback, motif search or approximate commensuration.
This module does not change the legacy propagation or refinement path.
"""
from dataclasses import dataclass, field
from fractions import Fraction
import math

import numpy as np

from .lattice import DirectLattice, _readonly
from .structure import AtomicStructureModel, AtomicSite
from .zone_cell import ZoneAxisCell, _integer_array, recast_structure_to_zone_cell


@dataclass(frozen=True)
class SimulationCompatibility:
    compatible: bool
    longitudinal_offsets: tuple[Fraction, Fraction]
    reason: str


def simulation_cell_compatibility(cell: ZoneAxisCell) -> SimulationCompatibility:
    """Prove integer C shears using exact rational arithmetic on input floats.

    Conservative by design: near-integer offsets are not snapped. Other integer
    supercells or motif translations may support geometries rejected here.
    """
    original = [[Fraction.from_float(float(v)) for v in row]
                for row in cell.original_lattice.vectors]
    rows = [[sum(int(cell.transform[i, k])*original[k][j] for k in range(3))
             for j in range(3)] for i in range(3)]
    c2 = sum(v*v for v in rows[2])
    offsets = tuple(sum(a*b for a, b in zip(row, rows[2]))/c2 for row in rows[:2])
    valid = all(v.denominator == 1 for v in offsets)
    return SimulationCompatibility(valid, offsets,
        'Integer C shears establish transverse lattice translations.' if valid else
        'Unsupported depth-resolved geometry: projected transverse periods are not '
        'proven lattice translations at fixed depth (noninteger C shears). '
        'A different simulation cell may be required; full-period projection is disabled.')


@dataclass(frozen=True, eq=False)
class MultisliceSimulationCell:
    """Restricted simulation cell, distinct from the B3 indexing/recast cell.

    P=A-mC, Q=B-nC with exact integer m,n. This unimodular change preserves
    internal ZOLZ HK identities. Future larger cells require explicit HK embedding.
    """
    zone_cell: ZoneAxisCell
    transform_from_zone: np.ndarray = field(init=False)
    direct_lattice: DirectLattice = field(init=False)
    frame: np.ndarray = field(init=False)

    def __post_init__(self):
        report = simulation_cell_compatibility(self.zone_cell)
        if not report.compatible:
            raise ValueError(report.reason)
        m, n = (int(v) for v in report.longitudinal_offsets)
        transform = _integer_array(((1, 0, -m), (0, 1, -n), (0, 0, 1)))
        direct = DirectLattice(transform @ self.zone_cell.direct_lattice.vectors)
        p, _, c = direct.vectors
        x, z = p/np.linalg.norm(p), c/np.linalg.norm(c)
        object.__setattr__(self, 'transform_from_zone', transform)
        object.__setattr__(self, 'direct_lattice', direct)
        object.__setattr__(self, 'frame', _readonly([x, np.cross(z, x), z]))


@dataclass(frozen=True)
class DepthSlice:
    lower_fraction: float
    upper_fraction: float
    thickness_angstrom: float
    source_indices: tuple[int, ...]
    sites: tuple[AtomicSite, ...]


@dataclass(frozen=True)
class DepthResolvedSample:
    original_structure: AtomicStructureModel
    simulation_cell: MultisliceSimulationCell
    recast_structure: AtomicStructureModel
    slices: tuple[DepthSlice, ...]
    c_period_translations: tuple[int, ...]


def build_depth_resolved_sample(structure: AtomicStructureModel,
                                cell: ZoneAxisCell, *,
                                slice_boundaries_fractional) -> DepthResolvedSample:
    """Partition one C motif into caller-specified half-open physical depth slabs.

    Two endpoints explicitly request one whole-period slice. This never bypasses
    the transverse-periodicity proof or substitutes for unsupported geometry.
    Site x/y and metadata are retained, never wrapped/merged. Integer C shifts
    used for depth ownership are recorded by original site index. No thermal,
    absorption, mesh, beam or refinement defaults are invented here.
    """
    simulation = MultisliceSimulationCell(cell)  # Reject before constructing slices.
    boundaries = tuple(slice_boundaries_fractional)
    if (len(boundaries) < 2 or any(isinstance(v, (bool, np.bool_)) or
            not isinstance(v, (int, float, np.integer, np.floating)) or
            not np.isfinite(v) for v in boundaries) or
            boundaries[0] != 0 or boundaries[-1] != 1 or
            any(a >= b for a, b in zip(boundaries, boundaries[1:]))):
        raise ValueError('explicit increasing slice boundaries must span 0 to 1 with at least two boundary values')
    recast = recast_structure_to_zone_cell(structure, cell)
    m, n = (-int(simulation.transform_from_zone[i, 2]) for i in (0, 1))
    groups = [[] for _ in boundaries[:-1]]
    translations = []
    for index, site in enumerate(recast.sites):
        x, y, z = site.fractional_position
        depth = z + m*x + n*y
        translation = math.floor(depth)
        local = depth - translation
        translations.append(translation)
        owner = int(np.searchsorted(boundaries, local, side='right')) - 1
        groups[owner].append((index, AtomicSite(site.symbol, (x, y, local),
                                               site.occupancy, site.debye_waller_b,
                                               site.physical_group_id, site.original_site_index)))
    period = float(np.linalg.norm(simulation.direct_lattice.vectors[2]))
    slices = tuple(DepthSlice(float(lo), float(hi), (hi-lo)*period,
                    tuple(i for i, _ in group), tuple(s for _, s in group))
                   for lo, hi, group in zip(boundaries, boundaries[1:], groups))
    return DepthResolvedSample(structure, simulation, recast, slices, tuple(translations))
