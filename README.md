# QCBED Python

Python rewrite of a legacy Fortran multislice QCBED simulation program.

## Current status

An optional [macOS arm64 app build](packaging/macos/README.md) packages the existing
GUI as `QCBED.app` using PyInstaller; Python development and CLI workflows remain
available. Build artifacts are generated under `dist/` and are not tracked.

An optional [Qt configuration editor](docs/configuration_editor.md) supports Open,
validated structural transactions, Undo/Redo, compiled summaries, atomic Save/Save As,
isolated Run/Resume, cooperative Cancel, and a separate Compatibility / Legacy menu
for explicit migration and optional exports. Install with
`python -m pip install -e '.[gui]'`, then launch
`python -m qcbed.gui examples/qcbed_config.yml`. Runs use immutable canonical snapshots in a separate process;
PySide6 is not required for the core package or CLI.

[Frontend-independent APIs](docs/frontend_apis.md) provide structured validation,
physical target discovery, immutable editing, atomic YAML saving, and canonical
run snapshots. No GUI toolkit is required.

Canonical typed Python configuration is available through `QcbedConfig` and
`compile_config`, alongside the unchanged legacy loaders. Compilation performs
no file I/O and preserves existing refinement ordering and numerical mapping.
See [configuration Stages 1–5](docs/configuration.md) and the directly constructed
[reference configuration](examples/reference_configuration.py). Explicit experimental
loading supports headerless bodies or strict legacy-header validation against
configuration. The optional [strict YAML frontend](docs/yaml_configuration.md)
loads the same typed objects through `load_yaml_config`; see
[the canonical YAML example](examples/qcbed_config.yml). Optional
[legacy import/export tooling](docs/legacy_migration.md) supports migration and
compatibility bundles; canonical runs do not require FOR007.

For production runs, use `python -m qcbed.run --config examples/qcbed_config.yml`.
See the [YAML run and resume workflow](docs/yaml_runs.md); the legacy
`examples/run_full_refinement.py` entry point remains supported.

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
acceleration with identical masks. Unverified historical build/linkage still prevents
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

`compute_legacy_intensity_scaling` implements the three source YNORM branches
as one global float64 factor. `apply_legacy_intensity_scaling` returns scaled,
immutable calculated intensities without changing inputs. Unguarded singular
arithmetic propagates NaN/Inf.
See [intensity scaling](docs/intensity_scaling.md) for formulas and runtime limitations.

`compute_legacy_objective` reuses YNORM to compute the source's normalized
weighted root error in beam-major accumulation order, returning immutable
float64 scalars. Refinement and GRADIENT preprocessing remain unimplemented;
historical executable equivalence is unverified. See
[objective evaluation](docs/objective_function.md).

The [optimizer architecture audit](docs/optimizer_architecture.md) documents
legacy sequencing, parameter flow and missing source. Source-derived
`run_legacy_anneal` now provides seeded float64 adaptive Gaussian search through
the existing trial evaluator. It accepts improvements and ties, with no uphill
acceptance. Historical executable/RNG equivalence is unverified; see
[ANNEAL optimizer](docs/anneal_optimizer.md).
`build_legacy_simplex` and `run_legacy_amoeba` implement the recovered modified
Nelder-Mead algorithm through legacy-simplex trial evaluation, with explicit
final-simplex and evaluated-best diagnostics. See [AMOEBA optimizer](docs/amoeba_optimizer.md)
for source-specific stopping rules and unresolved historical stride/equivalence
limitations. `run_legacy_paramin` adds source-derived coordinate parabolic
search through physical trial evaluation, with adaptive relative ranges and
fixed or seeded random traversal. See [PARAMIN optimizer](docs/paramin_optimizer.md)
for its fallback, stopping and historical-equivalence limitations.
`run_legacy_simpara` adds simultaneous stochastic parabolic proposals with
source direction and temperature adaptation. Its startup direction/counter
are explicitly zero-initialized by agreement; historical executable equivalence
is unverified. See [SIMPARA optimizer](docs/simpara_optimizer.md).
`parse_legacy_refinement` and `run_legacy_refinement` now support the ten
source-selected standalone/hybrid modes and full-sequence LOOPS repetition.
The driver preserves retained-state handoff and the original AMOEBA scale
anchor, with immutable stage history. See [refinement driver](docs/refinement_driver.md)
for control mapping, supported scope and executable-equivalence limitations.

Refinement parameter layouts, packing/unpacking and explicit physical/simplex
trial evaluation are implemented. Trials reuse the serial numerical kernels
for vacancy-free S-atom samples without surface corrections; unresolved states
are rejected. See [refinement parameters](docs/refinement_parameters.md) for
the 18-entry reference layout, supported behavior and equivalence limitations.

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

To run a deliberately short local optimization with explicit non-legacy masked
objective and pattern checkpoints:

```powershell
Set-Location C:\data\QCBED
& "C:\Users\peng863\AppData\Local\miniconda3\envs\qcbed-python\python.exe" examples/run_local_optimization.py
```

The script prints a new output directory with checkpoint PNGs, progress JSONL
and a final HTML report. See [local optimization](docs/local_optimization.md)
for editable settings, evaluation-versus-iteration semantics and limitations.

CBED validation now includes source-derived `out.dat` I/O, immutable comparison
maps, per-beam residual diagnostics and self-contained HTML reports. Install
the optional `visualization` extra for Matplotlib. The example
`examples/validate_reference_output.py` evaluates the full reference dimensions
and writes to `tmp/end_to_end_validation/`. Its objective is undefined because
the supplied data contain zero variances; no successful reference refinement
or historical executable equivalence is claimed. See
[end-to-end validation](docs/end_to_end_validation.md) for APIs and evidence.

- NumPy/SciPy serial reference implementation
- CPU process-based parallelism over independent incident tilts
- future optional GPU backend

Wave propagation now follows the restored readable FFT signs: negative
exponent divided by the mesh size into reciprocal space, positive exponent
without normalization on return. Potential construction is unchanged.
See [FFT correction validation](docs/fft_sign_correction.md), including the
historical-output comparison and unresolved scattering-source provenance.

## Batched refinement execution

Refinement can explicitly reuse the validated NumPy batched tilt engine:
set `LegacyRefinementConfig.execution` to
`ForwardExecutionConfig("batched", batch_size=64)` (for example with
`dataclasses.replace`). Serial remains the default reference path and objective
defaults are unchanged. See [local optimization](docs/local_optimization.md#forward-execution)
for the API and [performance measurements](docs/cpu_performance_benchmark.md#integrated-refinement-trial-backend).

## Durable full-refinement runner

The editable full-refinement example now selects depth and scan coordinates
with explicit user-selected physical hard bounds, checked after source mapping
and before simulation. Package APIs default to no hard bounds. See
[refinement bounds](docs/refinement_bounds.md) for settings and resume rules.
Specific-target, category absolute and category relative rules expand over the
active layout automatically; `--dry-run` displays each resolved bound.

Production runs explicitly reject numerical-invalid proposals with an infinite
score, preserving the retained state and logging the reason. Low-level APIs
default to `invalid_trial_policy="raise"`. This is a Python safety policy, not
proven historical Fortran behavior. Optional category selection can restrict a
first run to depth and scan coordinates; see [invalid-trial policy](docs/invalid_trial_policy.md).

Run `examples/run_full_refinement.py --quick` for a short real-data pipeline
check, or `--dry-run` to inspect configuration without optimization. Edit the
settings at the top before running without flags. Use
`examples/run_full_refinement.py --resume <run_directory>` to restart from a
saved retained stage boundary. Defaults use masked objective and batched NumPy
(64 tilts); checkpoints do not claim exact mid-optimizer continuation.
See [full refinement](docs/full_refinement.md) for exact interpreter commands,
budgets, previews, output files and restart limitations.
