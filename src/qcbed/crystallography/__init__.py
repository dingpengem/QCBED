"""B1 crystallography: row-vector Å/Å^-1 lattices, no 2π or beam indexing."""
from .lattice import DirectLattice, ReciprocalLattice
from .structure import AtomicSite, AtomicStructureModel
from .ase_adapter import atomic_structure_from_ase
from .visualization import plot_atomic_structure

__all__ = ['DirectLattice', 'ReciprocalLattice', 'AtomicSite',
           'AtomicStructureModel', 'atomic_structure_from_ase', 'plot_atomic_structure']

from .orientation import ZoneAxis, CrystalOrientation, zone_axis_cartesian, zone_law_value
from .projection import ProjectedAtomicSite, ProjectedReciprocalVector, project_atomic_structure, project_hkl
from .visualization import plot_projected_atomic_structure
__all__ += ['ZoneAxis', 'CrystalOrientation', 'zone_axis_cartesian', 'zone_law_value',
            'ProjectedAtomicSite', 'ProjectedReciprocalVector', 'project_atomic_structure',
            'project_hkl', 'plot_projected_atomic_structure']

from .zone_cell import (ZoneAxisCell, build_zone_axis_cell, original_fractional_to_zone,
                        zone_fractional_to_original, recast_structure_to_zone_cell)
from .index_transform import (original_hkl_to_zone_HKL, zone_HKL_to_original_hkl,
    is_zone_reflection_original, is_zolz_zone_index, ZoneReflectionHK,
    original_hkl_to_internal_HK, internal_HK_to_original_hkl)
from .reflection_basis import (OriginalReflectionBasis, reflection_from_basis_coordinates,
                               original_reflection_basis_vectors)
__all__ += ['ZoneAxisCell', 'build_zone_axis_cell', 'original_fractional_to_zone',
    'zone_fractional_to_original', 'recast_structure_to_zone_cell',
    'original_hkl_to_zone_HKL', 'zone_HKL_to_original_hkl', 'is_zone_reflection_original',
    'is_zolz_zone_index', 'ZoneReflectionHK', 'original_hkl_to_internal_HK',
    'internal_HK_to_original_hkl', 'OriginalReflectionBasis',
    'reflection_from_basis_coordinates', 'original_reflection_basis_vectors']

from .beam_indexing import IndexedBeamReflection, index_experimental_beam, index_experimental_beams
from .prepared_indexing import (IndexedPreparedDataset, index_prepared_dataset,
    index_beam_catalogue, internal_scan_basis, adapt_indexed_prepared_dataset)
__all__ += ['IndexedBeamReflection', 'index_experimental_beam', 'index_experimental_beams',
    'IndexedPreparedDataset', 'index_prepared_dataset', 'index_beam_catalogue',
    'internal_scan_basis', 'adapt_indexed_prepared_dataset']

from .simulation_cell import (
    SimulationCompatibility, MultisliceSimulationCell, DepthSlice,
    DepthResolvedSample, simulation_cell_compatibility, build_depth_resolved_sample,
)
__all__ += ['SimulationCompatibility', 'MultisliceSimulationCell', 'DepthSlice',
            'DepthResolvedSample', 'simulation_cell_compatibility',
            'build_depth_resolved_sample']

from .simulation_sample import ZoneAxisSimulationSample, build_zone_axis_simulation_sample
__all__ += ['ZoneAxisSimulationSample', 'build_zone_axis_simulation_sample']

from .legacy_for007_basis import (
    legacy_for007_001_basis_to_conventional_hkl,
    conventional_hkl_to_legacy_for007_001_basis,
    fcc_conventional_reflection_diagnostic,
)
__all__ += ['legacy_for007_001_basis_to_conventional_hkl',
            'conventional_hkl_to_legacy_for007_001_basis',
            'fcc_conventional_reflection_diagnostic']
