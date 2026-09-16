# QCBED Python

Python rewrite of a legacy Fortran multislice QCBED simulation program.

## Current status

The first legacy input-parser stage is implemented. It supports the formatted
`FOR007.DAT` title and run header, beam voltage and zone axis, mesh and slice
geometry, species and atom records, and the output and experimental-intensity
filenames. Later refinement and structure-factor sections remain deferred;
initial scan geometry has a separate parser described below. These parsers
do not calculate wavelength or physics.

Standalone electron-wavelength and grid utilities are also implemented.
`electron_wavelength_angstrom(voltage_kv)` accepts positive finite voltage in
kV and returns angstroms using float64-compatible arithmetic; convert volts
explicitly. Exact executable equivalence remains unestablished because
default-kind Fortran literal rounding may depend on the original compiler flags.

`elastic_atomic_scattering_factor(Z, stl_inverse_angstrom)` implements the
source-scaled elastic electron factor for neutral atoms Z=1–98 using the
source's tabulated cubic interpolation and Mott branches. See
[scattering conventions and validation](docs/scattering.md) for units, bounds,
and the remaining executable-equivalence limitation.

Immutable 2D lattice and reciprocal grids provide separate dimensionless H/K
indices, physical vectors and STL in inverse angstroms, using the source's
convention without a 2*pi factor. See [reciprocal geometry](docs/reciprocal_geometry.md).

`elastic_structure_factors` evaluates the raw complex128 elastic atomic sum
for fixed S atoms, with explicit applied thermal factors and vacancy fraction.
It preserves the positive atomic phase and dilated cell-area scaling; slice
shifts belong to propagation. This is not a complete SFMAP implementation,
and executable-Fortran equivalence remains unestablished. See
[structure-factor conventions and validation](docs/structure_factors.md).

Absorptive atomic factors and raw complex structure-factor sums are also
implemented, with explicit phenomenological coefficients, positive applied
Debye-Waller factors, voltage in kV, and literal legacy absorption-flag
semantics. Elastic/absorptive combination follows the source sign convention.
These absorption helpers do not compute propagation. See [absorption](docs/absorption.md)
for numerical domains, saved-output comparisons and equivalence limitations.

A serial direct-Fourier reference now reconstructs projected potential and
applies the source PGMAP exponential, with separate DSAP and Nyquist masks.
An explicit NumPy CPU FFT path is validated against that direct sum, including
odd/even ordering and the normalized mathematical inverse. Phase-grating
construction retains `method="direct"` by default; `method="fft"` selects
acceleration with identical masks. The missing readable 1D CDFT still prevents
claiming legacy FFT equivalence. See [phase gratings](docs/phase_grating.md)
and [FFT backend](docs/fft_backend.md) for validation and timing results.

The standalone PROPG kernel now provides source-backed propagation and
fractional slice-shift phases, explicit reciprocal ordering and separate
Nyquist suppression. Internal tilts are reciprocal-basis coefficients, not
angles. No RMT normalization or wave propagation is included. See
[propagator conventions](docs/propagator.md) for validation and limitations.

`propagate_one_step` connects a real-space wave, phase grating and reciprocal
propagator using the explicit Python FFT pair, without an extra RMT factor.
It returns immutable real-space and post-propagator reciprocal-wave snapshots.
Complex intermediate comparisons validate this single step against direct DFTs;
historical CDFT sign/normalization equivalence remains unresolved. See
[wave propagation](docs/wave_propagation.md).

`simulate_slice_sequence` executes precomputed slice gratings and propagators
for one tilt, using validated zero-based slice IDs and explicit reciprocal
ordering. It returns an immutable final reciprocal wave, requested complex
beam amplitudes and raw intensities, preserving beam order and duplicates.
It omits unused final transforms and uses unit plane-wave illumination by
default. Fitted output scaling remains deferred. See
[multislice sequences](docs/multislice_sequence.md).

`prepare_multislice_sample` now builds reusable immutable gratings from parsed
fixed S atoms and explicitly supplied applied parameters. `simulate_one_tilt`
builds per-slice propagators and returns raw beam amplitudes/intensities, using
the common undilated reciprocal-cell tilt coefficients. Preparation supports
direct or FFT reconstruction, a common even mesh and matching undilated lateral
cells, with unity reflection-refinement multipliers explicitly selected.
The real Al input is exercised at its recorded initial parameters; historical
CDFT/executable equivalence remains unestablished. See
[single-tilt simulation](docs/single_tilt_simulation.md) for the required
configuration, supported scope and validation.

`simulate_tilt_series_serial` and `simulate_tilt_series_batched` execute
explicit float64 `(n_tilts,2)` arrays of undilated reciprocal-cell H/K
coefficients. The chunked NumPy path transforms only spatial axes and returns
immutable `(n_tilts,n_beams)` amplitudes/intensities in request order.
Serial and independent-DFT tests validate the Python convention; historical
CDFT/executable equivalence remains unresolved. See [tilt series](docs/tilt_series.md).

`simulate_tilt_series_parallel` distributes contiguous tilt chunks to spawned
CPU processes, each using the validated batched engine with a prepared sample
loaded once per worker. Process chunk size and internal FFT batch size are
explicit; result ordering is deterministic. Windows callers need a guarded
entry point and should set numerical-library thread limits before launch.
See [CPU parallelism](docs/cpu_parallelism.md) for validation, memory costs and
benchmarks; process startup can make one-process batching faster.

## Coordinate and atom API

`parse_experimental_intensities` reads the experimental header and separate
intensity, variance and weight blocks into immutable float64 `(tilt, beam)`
arrays. Beam H/K order is supplied explicitly; optional scan validation checks
both dimensions. File values are preserved without scaling or weighting.
See [experimental intensity parsing](docs/experimental_intensity.md) for the
source mapping and unresolved preprocessing/GRADIENT semantics.

`parse_for007_scan` reads the initial scan basis/corners and only the
NUMP/NUML header of the referenced experimental file. `generate_scan_tilts`
produces float64 undilated H/K pairs in legacy column-fast order for all three
tilt-series APIs. Singleton dimensions are preserved by parsing but rejected
by generation because the source divides by zero. See [legacy scan](docs/legacy_scan.md).

- `SliceGeometry.x_shift_fractional` and `y_shift_fractional` are fractional
  direct-lattice translations, used by legacy PROPG as `H*X0 + K*Y0`.
  Cell lengths and thickness remain in angstroms; lattice gamma is in degrees.
- Exported `AtomSite` records expose `species_index_one_based`, `x_fractional`,
  `y_fractional`, `starting_occupancy`, and `record_type` (`S` or `G`). Positions
  and occupancy are input values, before gas randomization or vacancy effects.
  Use `atom.species_index_zero_based(len(parsed.species))` before Python species
  lookup; the helper and parser reject out-of-range species identifiers.
- Grid arrays use `(mesh_y, mesh_x)` with x along columns and y along rows.
  Generic PMAP conversions support odd and even dimensions and preserve dtype.
  Generic reciprocal indices for size 5 are `[-2, -1, 0, 1, 2]` centered and
  `[0, 1, 2, -2, -1]` unshifted. Legacy `oshift` requires both dimensions even.

The former angstrom shift names and ambiguous atom field names have been
removed; no compatibility aliases are provided. `SliceGeometry` and `AtomSite`
remain available from the package root. See [source provenance](docs/source_provenance.md).

## Legacy parser syntax

`parse_for007` follows each READ through the experimental-intensity filename.
It supports record continuation, comma/whitespace separators, E/D exponents,
explicit repetitions, and quoted list-directed S/G values. Numeric inputs are
finite; integer fields are parsed strictly, with field-specific bounds and
mask/flag domains. Errors include filename, section, record, and offending input.

Raw title, absorption and filename records are retained in `*_raw` fields;
`Species.label_raw` retains the full record while `label` follows FORMAT(A2).
Quotes in these formatted fields are literal. `AtomSite.record_type_raw`
preserves the list-directed token; `record_type` is its case-sensitive S/G value.
See [parser contract](docs/legacy_parser.md) for unsupported forms and the
limits of compatibility with the Fortran runtime.

## Directories

- `private predecessor sources`: original Fortran implementation
- `reference_data/`: validation inputs and outputs
- `src/qcbed/`: Python package
- `tests/`: automated tests
- `benchmarks/`: performance benchmarks
- `docs/`: technical documentation

## Current target

- NumPy/SciPy serial reference implementation
- CPU process-based parallelism over independent incident tilts
- future optional GPU backend
