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
    projected_potential_fft,
    electron_interaction_constant,
    phase_grating_from_potential,
    phase_grating_from_structure_factors,
)
from .fft_backend import (
    centered_to_unshifted,
    unshifted_to_centered,
    reciprocal_to_real_sum,
    real_to_reciprocal_coefficients,
)
from .propagator import (
    PropagationParameters,
    propagation_phase_radians,
    fractional_slice_shift_phase_radians,
    propagation_boundary_mask,
    propagation_kernel,
)
from .wave_propagation import (
    WaveStepResult,
    apply_phase_grating,
    real_to_reciprocal_wave,
    apply_propagator,
    reciprocal_to_real_wave,
    propagate_one_step,
)
from .multislice import (
    MultisliceResult,
    legacy_slice_sequence_to_zero_based,
    unit_incident_wave,
    beam_indices_to_row_col,
    beam_indices_to_legacy_ib,
    extract_beam_amplitudes,
    beam_intensities,
    simulate_slice_sequence,
)
from .simulation import (
    SliceAppliedParameters,
    SamplePreparationParameters,
    PreparedSliceType,
    PreparedMultisliceSample,
    prepare_multislice_sample,
    build_tilt_propagators,
    simulate_one_tilt,
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
    "projected_potential_fft",
    "electron_interaction_constant",
    "phase_grating_from_potential",
    "phase_grating_from_structure_factors",
    "centered_to_unshifted",
    "unshifted_to_centered",
    "reciprocal_to_real_sum",
    "real_to_reciprocal_coefficients",
    "PropagationParameters",
    "propagation_phase_radians",
    "fractional_slice_shift_phase_radians",
    "propagation_boundary_mask",
    "propagation_kernel",
    "WaveStepResult",
    "apply_phase_grating",
    "real_to_reciprocal_wave",
    "apply_propagator",
    "reciprocal_to_real_wave",
    "propagate_one_step",
    "MultisliceResult",
    "legacy_slice_sequence_to_zero_based",
    "unit_incident_wave",
    "beam_indices_to_row_col",
    "beam_indices_to_legacy_ib",
    "extract_beam_amplitudes",
    "beam_intensities",
    "simulate_slice_sequence",
    "SliceAppliedParameters",
    "SamplePreparationParameters",
    "PreparedSliceType",
    "PreparedMultisliceSample",
    "prepare_multislice_sample",
    "build_tilt_propagators",
    "simulate_one_tilt",
]
