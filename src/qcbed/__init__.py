"""Foundational types and grid helpers for the QCBED rewrite."""

from .models import (
    MeshConfig,
    MicroscopeConfig,
    SimulationConfig,
    SliceGeometry,
    TiltCoordinate,
)
from .legacy_io import AtomSite, LegacyInput, LegacyInputParseError, Species, parse_for007
from .electron import electron_wavelength_angstrom
from .scattering import elastic_atomic_scattering_factor
from .structure_factors import elastic_structure_factors
from .absorption import (
    AbsorptionCoefficients,
    absorptive_atomic_scattering_factor,
    absorptive_structure_factors,
    combine_elastic_absorptive_structure_factors,
)
from .reciprocal import DirectLattice2D, ReciprocalLattice2D, ReciprocalIndexGrid, ReciprocalVectorGrid
from .phase_grating import (
    reciprocal_cutoff_mask,
    legacy_nyquist_mask,
    projected_potential_direct,
    electron_interaction_constant,
    phase_grating_from_potential,
    phase_grating_from_structure_factors,
)

__all__ = [
    "MeshConfig",
    "MicroscopeConfig",
    "SimulationConfig",
    "SliceGeometry",
    "TiltCoordinate",
    "AtomSite",
    "LegacyInput",
    "LegacyInputParseError",
    "Species",
    "parse_for007",
    "electron_wavelength_angstrom",
    "elastic_atomic_scattering_factor",
    "elastic_structure_factors",
    "AbsorptionCoefficients",
    "absorptive_atomic_scattering_factor",
    "absorptive_structure_factors",
    "combine_elastic_absorptive_structure_factors",
    "DirectLattice2D",
    "ReciprocalLattice2D",
    "ReciprocalIndexGrid",
    "ReciprocalVectorGrid",
    "reciprocal_cutoff_mask",
    "legacy_nyquist_mask",
    "projected_potential_direct",
    "electron_interaction_constant",
    "phase_grating_from_potential",
    "phase_grating_from_structure_factors",
]
