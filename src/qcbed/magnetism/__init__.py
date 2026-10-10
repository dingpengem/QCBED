"""Research-stage magnetic QCBED module; NOT a validated magnetic refinement engine.

The established nonmagnetic QCBED path remains unchanged. See
``docs/magnetic_qcbed_theory.md`` and ``docs/magnetic_qcbed_development.md``.
"""
from .forward import add_magnetic_gratings, magnetic_phase_grating, project_az_midpoint
from .research import (
    ElectronBeam, PeriodicMagneticField, SensitivityResult,
    diffraction_intensity, finite_difference_sensitivity,
    gaussian_magnetization_grid, marginal_magnetic_fisher,
    nio_type_ii_sites, periodic_field_from_magnetization,
    scalar_magnetic_multislice,
)

__all__ = [
    'ElectronBeam', 'PeriodicMagneticField', 'SensitivityResult',
    'add_magnetic_gratings', 'magnetic_phase_grating', 'project_az_midpoint',
    'diffraction_intensity', 'finite_difference_sensitivity',
    'gaussian_magnetization_grid', 'marginal_magnetic_fisher',
    'nio_type_ii_sites', 'periodic_field_from_magnetization',
    'scalar_magnetic_multislice',
]
