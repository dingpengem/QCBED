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
from .masked_objective import MaskedObjectiveResult, compute_masked_objective, compute_selected_objective
from .local_runner import (RefinementProgressRecord, RefinementSnapshot,
                           LocalOptimizationRunResult, run_local_optimization)
from .output_io import (SimulatedCbedData, LegacySimulatedOutput, reshape_scan_map,
                        write_legacy_out_dat, parse_legacy_out_dat)
from .validation import CbedComparison, RefinementValidation, validate_refinement
from .visualization import (plot_cbed_comparison, plot_refinement_convergence,
                            plot_parameter_trajectories)
from .validation_report import generate_refinement_validation_report
from .anneal import LegacyAnnealResult, run_legacy_anneal
from .amoeba import LegacyAmoebaResult, build_legacy_simplex, run_legacy_amoeba
from .paramin import LegacyParaminResult, run_legacy_paramin
from .simpara import LegacySimparaResult, run_legacy_simpara
from .refinement_driver import (LegacyRefinementConfig, LegacyRefinementResult,
                                parse_legacy_refinement, run_legacy_refinement)
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

from .tilt_series import (
    ParallelExecutionMetadata, TiltSeriesResult, simulate_tilt_series_serial, simulate_tilt_series_batched,
)
from .parallel import simulate_tilt_series_parallel
from .scan import LegacyScanDefinition, parse_for007_scan, generate_scan_tilts
from .experimental import ExperimentalIntensityData, parse_experimental_intensities
from .intensity_scaling import compute_legacy_intensity_scaling, apply_legacy_intensity_scaling
from .objective import LegacyObjectiveResult, compute_legacy_objective
from .refinement import (
    RefinementRegion, RefinementReflection, RefinementState,
    RefinementParameterSpec, RefinementParameterLayout,
    build_refinement_parameter_layout, pack_refinement_parameters, apply_refinement_parameters,
)
from .refinement_io import parse_refinement_state
from .refinement_trial import RefinementTrialResult, evaluate_refinement_trial

__all__ = [
    "MaskedObjectiveResult", "compute_masked_objective", "compute_selected_objective",
    "RefinementProgressRecord", "RefinementSnapshot", "LocalOptimizationRunResult", "run_local_optimization",
    "SimulatedCbedData", "LegacySimulatedOutput", "reshape_scan_map",
    "write_legacy_out_dat", "parse_legacy_out_dat", "CbedComparison",
    "RefinementValidation", "validate_refinement", "plot_cbed_comparison",
    "plot_refinement_convergence", "plot_parameter_trajectories",
    "generate_refinement_validation_report",
    "RefinementRegion",
    "RefinementReflection",
    "RefinementState",
    "RefinementParameterSpec",
    "RefinementParameterLayout",
    "build_refinement_parameter_layout",
    "pack_refinement_parameters",
    "apply_refinement_parameters",
    "parse_refinement_state",
    "RefinementTrialResult",
    "evaluate_refinement_trial",
    "LegacyObjectiveResult",
    "compute_legacy_objective",
    "compute_legacy_intensity_scaling",
    "apply_legacy_intensity_scaling",
    "ExperimentalIntensityData",
    "parse_experimental_intensities",
    "LegacyScanDefinition",
    "parse_for007_scan",
    "generate_scan_tilts",
    "ParallelExecutionMetadata",
    "simulate_tilt_series_parallel",
    "TiltSeriesResult",
    "simulate_tilt_series_serial",
    "simulate_tilt_series_batched",
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
    "LegacyAnnealResult",
    "run_legacy_anneal",
    "LegacyAmoebaResult",
    "build_legacy_simplex",
    "run_legacy_amoeba",
    "LegacyParaminResult",
    "run_legacy_paramin",
    "LegacySimparaResult",
    "run_legacy_simpara",
    "LegacyRefinementConfig",
    "LegacyRefinementResult",
    "parse_legacy_refinement",
    "run_legacy_refinement",
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
