"""Explicit prepared-data boundary to the unchanged refinement engine.

Geometric lattice labels are never interpreted as physical HK. Callers supply
both ordered ID/HK assignments and the scan's physical reciprocal basis. No
variance, intensity normalization, gradient or weighting is manufactured here.
"""
from dataclasses import dataclass
import numpy as np

from .models import PreparedExperimentalDataset, SignalMode
from ..experimental import ExperimentalIntensityData
from ..scan import LegacyScanDefinition


@dataclass(frozen=True)
class PreparedRefinementInput:
    experimental: ExperimentalIntensityData
    scan: LegacyScanDefinition
    # Zero-based index is the unchanged prepared column/extracted retained block.
    beam_mapping: tuple[tuple[str, tuple[int, int] | None, tuple[int, int], int], ...]
    objective_mode: str


def adapt_prepared_dataset(dataset, *, beam_hk_assignments, basis_1, basis_2,
                           objective_mode):
    """Return native experimental data and scan without modifying scientific values.

    beam_hk_assignments is a sequence of (stable ID, physical (h,k)) pairs.
    Assignment order is irrelevant; dataset order determines numerical ordering.
    Duplicate physical HK entries are retained, duplicate IDs rejected.
    basis_1/2 must describe the physical reciprocal basis of the A10 coefficients.
    Missing scientific variance is an actionable readiness error even in beam-mask
    mode. Zero-weight/zero-variance entries are passed verbatim for masked mode;
    legacy mode requires positive variance everywhere to avoid 0/0.
    """
    if not isinstance(dataset, PreparedExperimentalDataset):
        raise TypeError('dataset must be PreparedExperimentalDataset')
    if objective_mode not in ('legacy', 'masked'):
        raise ValueError('objective_mode must be legacy or masked')
    if dataset.variances is None:
        raise ValueError('Refinement not ready: validated scientific variance is absent; support is not variance')
    if any(not b.enabled for b in dataset.beams):
        raise ValueError('Disabled beams require explicit dataset selection before refinement')
    assignments = tuple(beam_hk_assignments)
    ids = [entry[0] for entry in assignments]
    expected = [b.id for b in dataset.beams]
    if len(set(ids)) != len(ids) or set(ids) != set(expected):
        raise ValueError('Physical HK assignments must identify every beam ID exactly once')
    by_id = dict(assignments)
    ordered = tuple(by_id[id] for id in expected)
    active = dataset.weights != 0 if objective_mode == 'masked' else np.ones(dataset.weights.shape, dtype=bool)
    if np.any(dataset.variances[active] <= 0):
        raise ValueError('Refinement not ready: objective requires positive variance at all evaluated points')
    sg = dataset.scan_geometry
    scan = LegacyScanDefinition(basis_1, basis_2, sg.corner_a, sg.corner_b, sg.corner_c,
                                sg.scan_columns, sg.scan_rows)
    experimental = ExperimentalIntensityData(sg.scan_columns, sg.scan_rows,
        1 if dataset.mode is SignalMode.GRADIENT else 0, ordered,
        dataset.intensities, dataset.variances, dataset.weights)
    mapping = tuple((b.id, dataset.lattice_coordinates[i] if dataset.lattice_coordinates else None,
                     experimental.beam_indices_hk[i], i) for i, b in enumerate(dataset.beams))
    return PreparedRefinementInput(experimental, scan, mapping, objective_mode)
