"""Educational NiO field-free magnetic sensitivity; NOT calibrated NiO QCBED.
Run: PYTHONPATH=src python examples/niO_magnetic_sensitivity.py
"""
import numpy as np
from qcbed.magnetism import (
    ElectronBeam, diffraction_intensity, finite_difference_sensitivity,
    gaussian_magnetization_grid, nio_type_ii_sites,
    periodic_field_from_magnetization, scalar_magnetic_multislice,
)

beam = ElectronBeam(200000.)
nz, ny, nx = 17, 33, 33  # odd mesh avoids ambiguous Nyquist derivatives

def simulate(p):
    # Magnetic terms only (no NiO electrostatic potential, absorption or TDS).
    sites, moments, lengths = nio_type_ii_sites(
        moment_mu_b_per_ni=p[0], canting_angle_rad=p[1])
    m = gaussian_magnetization_grid((nz, ny, nx), lengths, sites, moments, 0.6e-10)
    dz, dy, dx = lengths[2]/nz, lengths[1]/ny, lengths[0]/nx
    az = periodic_field_from_magnetization(m, (dz, dy, dx)).a_t_m[..., 2]
    out = scalar_magnetic_multislice(
        np.ones((ny, nx), dtype=np.complex128),
        np.zeros((nz, ny, nx)), az, np.full(nz, dz), (dy, dx), beam)
    return diffraction_intensity(out)

p = np.array([1.9, 0.10])
raw = simulate(p)
print('Total intensity:', raw.sum())
# sigma=1 arbitrary toy units; no experimental precision should be inferred.
info = finite_difference_sensitivity(simulate, p, [0.02,0.02],
                                     np.ones_like(raw), np.ones_like(raw, dtype=int))
print('Sensitivity singular values:', info.singular_values)
print('Whitened-Jacobian rank:', info.numerical_rank)
print('Notice: AFM order amplitude and cant angle can be correlated.')
