"""Narrow immutable Patch A annotation boundary; science arrays are never edited.

HKs passed to the existing refinement adapter are INTERNAL zone-cell HKs. They
are meaningful only with a simulation sample expressed in that same recast basis;
B4 does not construct/recast multislice samples or authorize an objective run.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING
import numpy as np
from .beam_indexing import IndexedBeamReflection, index_experimental_beams
from .reflection_basis import OriginalReflectionBasis
from .index_transform import original_hkl_to_internal_HK
from .zone_cell import ZoneAxisCell
from .orientation import _indices
if TYPE_CHECKING:
    from ..data_preparation.models import PreparedExperimentalDataset


@dataclass(frozen=True,eq=False)
class IndexedPreparedDataset:
    """Single indexed metadata tuple alongside, not inside, the unchanged A1 model."""
    dataset: PreparedExperimentalDataset
    reflection_basis: OriginalReflectionBasis
    zone_cell: ZoneAxisCell
    indexed_beams: tuple[IndexedBeamReflection,...] = field(init=False)

    def __post_init__(self):
        from ..data_preparation.models import PreparedExperimentalDataset
        if not isinstance(self.dataset,PreparedExperimentalDataset):
            raise TypeError('dataset must be PreparedExperimentalDataset')
        if len(self.dataset.lattice_coordinates) != len(self.dataset.beams):
            raise ValueError('prepared dataset requires one explicit (n1,n2) per beam; BeamObservation.indices is not inferred')
        records = index_experimental_beams(
            ((b.id,n) for b,n in zip(self.dataset.beams,self.dataset.lattice_coordinates,strict=True)),
            self.reflection_basis,self.zone_cell)
        object.__setattr__(self,'indexed_beams',records)


def index_prepared_dataset(dataset, reflection_basis, zone_cell) -> IndexedPreparedDataset:
    return IndexedPreparedDataset(dataset,reflection_basis,zone_cell)


def index_beam_catalogue(catalogue, reflection_basis, zone_cell):
    """Annotate A7 sites in their existing order, with no status-based selection.

    Only beam_id and integer_indices are read. Observed centers/ROIs/images and
    pixel-vector magnitudes play no part in physical indexing.
    """
    return index_experimental_beams(((s.beam_id,s.integer_indices) for s in catalogue.sites),
                                   reflection_basis,zone_cell)


def internal_scan_basis(indexed: IndexedPreparedDataset, *, scan_basis_transform):
    """Derive internal scan rows from an explicitly supplied A10 basis transform.

    A10 convention: scan-basis columns = experimental-basis columns @ M.
    M is an integer unimodular 2x2 matrix, supplied explicitly, never inferred from
    pixel lengths, frame orientation, corners or a textual provenance entry.
    For the canonical A10 example M=diag(1,-1), not identity. This preserves the
    A10 coefficient meaning without changing corners or reindexing the beams.
    """
    if not isinstance(indexed,IndexedPreparedDataset):
        raise TypeError('indexed must be IndexedPreparedDataset')
    return internal_scan_basis_from_crystallography(indexed.reflection_basis, indexed.zone_cell,
                                                   scan_basis_transform=scan_basis_transform)


def internal_scan_basis_from_crystallography(basis, zone_cell, *, scan_basis_transform):
    """Same B4a scan conversion before a prepared dataset is available."""
    rows = tuple(tuple(row) for row in scan_basis_transform)
    if len(rows)!=2 or any(len(row)!=2 for row in rows):
        raise ValueError('scan_basis_transform must be an explicit 2x2 integer matrix')
    rows = tuple(_indices((*row,0),'scan_basis_transform')[:2] for row in rows)
    if abs(rows[0][0]*rows[1][1]-rows[0][1]*rows[1][0])!=1:
        raise ValueError('scan_basis_transform must be unimodular')
    hk = [original_hkl_to_internal_HK(h,zone_cell) for h in (basis.hkl_1,basis.hkl_2)]
    scan = tuple((hk[0].H*rows[0][j]+hk[1].H*rows[1][j],
                  hk[0].K*rows[0][j]+hk[1].K*rows[1][j],0) for j in range(2))
    # The existing scan object converts these to float64; refuse lost integer precision.
    if any(not np.isfinite(float(v)) or int(float(v))!=v for row in scan for v in row):
        raise ValueError('internal scan basis is not exactly representable in float64')
    return scan


def adapt_indexed_prepared_dataset(indexed: IndexedPreparedDataset, *,
                                   scan_basis_transform, objective_mode):
    """Feed derived internal HK to the unchanged experimental/refinement adapter.

    Requires existing variance readiness. User inputs are original hkl and explicit
    A10 transform only, never manually calculated internal HK. Does not evaluate
    physics or verify that an arbitrary runner sample matches the recast cell.
    """
    from ..data_preparation.refinement_adapter import adapt_prepared_dataset
    scan = internal_scan_basis(indexed,scan_basis_transform=scan_basis_transform)
    return adapt_prepared_dataset(indexed.dataset,
        beam_hk_assignments=tuple((r.stable_id,(r.internal_hk.H,r.internal_hk.K)) for r in indexed.indexed_beams),
        basis_1=scan[0],basis_2=scan[1],objective_mode=objective_mode)
