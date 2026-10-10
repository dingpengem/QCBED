# Magnetic QCBED — research module integration notes

**Stage:** R0 research API only. No production magnetic refinement, GUI tab,
YAML/Resume schema or packaged app feature has been enabled.

## Repository layout

- `src/qcbed/magnetism/research.py` — SI beam relations; periodic M→A,B FFT;
  pedagogical NiO Type-II cell / Gaussian form factor; standalone scalar
  split-step and numerical Fisher information.
- `src/qcbed/magnetism/forward.py` — magnetic transmission grating in webers,
  composable with existing QCBED complex gratings without altering its
  existing propagation/FFT/masks/absorptive logic.
- `tests/test_magnetic_qcbed_research.py` and
  `tests/test_magnetic_qcbed_forward.py` — focused toy / algebra / zero-field tests.
- `examples/niO_magnetic_sensitivity.py` — explanatory toy simulation only.
- `docs/magnetic_qcbed_theory.md` — derivation, assumptions and acceptance gates.

## Use with existing QCBED engine (research mode only)

```python
from qcbed.magnetism import add_magnetic_gratings, project_az_midpoint
from qcbed.multislice import simulate_slice_sequence

# Existing baseline arrays must have been built by the established code:
# existing_gratings: tuple[complex128[ny,nx], ...] for stored slice types
# existing_propagators: same number of stored slice types
# az_at_centers_t_m: real[s,ny,nx], periodic Coulomb-gauge A_z, T*m
# stored_dz_m: positive meters per stored slice
projected_Az = project_az_midpoint(az_at_centers_t_m, stored_dz_m)
magnetic_gratings = add_magnetic_gratings(existing_gratings, projected_Az)
result = simulate_slice_sequence(
    magnetic_gratings, existing_propagators, slice_sequence_zero_based,
    beam_indices_hk=selected_beams_in_magnetic_cell,
    reciprocal_ordering='unshifted', propagator_ordering='unshifted',
)
```

**Important:** A magnetic cell can be larger than the nuclear cell; the
`selected_beams_in_magnetic_cell` indices, potential grid, cutoff and beam
catalogue must all be regenerated in the *same* magnetic simulation cell.
The snippet assumes that work was completed; it is not a one-click NiO run.

No change is made to the nonmagnetic engine's conventional calibration,
`E_g^(s)` mode, site-group target layout, stored-slice semantics, dashboard,
YAML, or checkpoint identity. The code deliberately cannot claim that a field
series has measured `m_ind` until finite-foil and `q=0` physics exist.

## Development gates

- **R0, available:** computational building blocks, unit tests, SI field/
  phase conventions, explicit q=0 exclusion and baseline-composition adapter.
- **R1, next:** integrate magnetic cell construction into `crystallography`;
  import calibrated Ni/O magnetic density; benchmark against published NiO
  diffraction and magnetic CBED; verify `A_perp`/Pauli approximation and
  `z`-boundary effects, plus sampling and probe convergence.
- **R2:** field-in-foil and vacuum solution, calibrated lens/holder field,
  field-dependent joint dataset, domain mixtures, reversible changes.
- **R3:** constrained magnetic optimizer targets with uncertainty and
  identifiability gating; then GUI/YAML/Resume/Results, full regression and
  macOS app packaging at milestone.

See [theory and validity boundaries](magnetic_qcbed_theory.md).