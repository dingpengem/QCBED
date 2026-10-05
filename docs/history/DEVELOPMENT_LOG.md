# Historical milestone record

This is an archived development record, not the current installation or GUI guide.
See [documentation navigation](../README.md) for current instructions. Measurements
and limitations below belong to their recorded milestone/configuration.

# QCBED Python

## Reconstruction milestone closure

The current reconstruction milestone is accepted. See the [release summary](../../docs/reconstruction_milestone.md)
for supported workflows, accepted test counts, limitations and recorded build provenance.
Release closure adds documentation only; no numerical or GUI workflow changes.


Python rewrite of a legacy Fortran multislice QCBED simulation program.

## Start a new project

Open `QCBED.app` (or run `python -m qcbed.gui`). The primary workflow is:

1. **A. Data Prep → New Data Preparation:** load the raw CBED pattern. Supply an
   optional aperture for calibration/threshold geometry, or draw and confirm a disk
   diameter. Select a local centrosymmetric feature and explicitly Accept/Snap its ROI.
2. Review automatically derived extraction/crop sizes and Advanced binning (default 2);
   run live preparation. Inspect beam
   boxes/labels, the zone-axis marker, and linked Pattern/Sigma/Weight. Confirm the
   origin, exclude bad pixels with weight=0, then **Validate A**.
3. **B. Atomic Model:** supply structure, zone axis, ordered physical reciprocal
   basis and site groups. Inspect Front/Top/Left projections, then **Validate B**.
4. **C. Main QCBED Parameters:** select matching signal mode, supported depth/B
   targets and bounds, workflow, batched execution and output paths. **Validate C**
   enables Run. Save the project before starting a production run.

Prepared descriptors are generated replay/cache artifacts. **Load Previous
Preparation** restores them; a saved preparation is not required to start.
Without an aperture, uncertainty is the explicitly labelled Poisson/intensity
approximation, not a calibrated detector noise model. See the [complete user
workflow and limitations](../../docs/gui_abc.md#live-experimental-preparation).

## Current status

Patch A core backend, live raw-data preparation and A/B/C initialization are complete
with documented limitations. The raw FCC Al acceptance uses live geometry and
explicit physical HK mapping; no-aperture uncertainty remains an approximation.
Historical validated-variance refinement results are documented separately below.
Patch A includes A0 specification and A1
immutable backend models in `qcbed.detector` and `qcbed.data_preparation`.
A2 adds `qcbed.detector.calibrate_psf_from_aperture()` for decoded 2-D aperture
arrays, returning PSF/OTF/MTF calibration with quality/provenance.
See the [model API and legacy workflow](../../docs/data_preparation_algorithm.md).
A3 adds `qcbed.detector.correct_cbed_with_otf()` for whole-frame direct complex-OTF
correction, separate from calibration. A4 adds `qcbed.data_preparation.fit_noise_model()` and
`variance_map_from_noise_model()` for legacy noise characterization after correction.
A5 adds `estimate_disk_geometry()` and `locate_central_disk()` for explicit profile
geometry and central-disk localization. A6 adds `detect_reciprocal_lattice_seed()` using known-radius boundary evidence and lattice
consistency, without pre-specified `inputparas` beam geometry. A7 adds global pixel-lattice fitting and a lattice-constrained beam catalogue,
keeping predicted and observed centers separate. symmetry-center adds explicit local centrosymmetry origin measurement. A9 adds beam-supported signal/variance/weight extraction. A10 supplies scan corners; `assemble_prepared_dataset()` and `adapt_prepared_dataset()` connect the canonical dataset directly to refinement. Physical HKL mapping must be supplied explicitly; applying PSF correction is not assumed or mandatory. Validated propagation/refinement is unchanged.


An optional [macOS arm64 app build](../../packaging/macos/README.md) packages the existing
GUI as `QCBED.app` using PyInstaller; Python development and CLI workflows remain
available. Build artifacts are generated under `dist/` and are not tracked.

An optional [Qt configuration editor](../../docs/configuration_editor.md) supports Open,
validated structural transactions, Undo/Redo, compiled summaries, atomic Save/Save As,
isolated Run/Resume, cooperative Cancel, and a separate Compatibility / Legacy menu
for explicit migration and optional exports. Install with
`python -m pip install -e '.[gui]'`, then launch
`python -m qcbed.gui examples/qcbed_config.yml`. Runs use immutable canonical snapshots in a separate process;
PySide6 is not required for the core package or CLI.

[Frontend-independent APIs](../../docs/frontend_apis.md) provide structured validation,
physical target discovery, immutable editing, atomic YAML saving, and canonical
run snapshots. No GUI toolkit is required.

Canonical typed Python configuration is available through `QcbedConfig` and
`compile_config`, alongside the unchanged legacy loaders. Compilation performs
no file I/O and preserves existing refinement ordering and numerical mapping.
See [configuration Stages 1–5](../../docs/configuration.md) and the directly constructed
[reference configuration](../../examples/reference_configuration.py). Explicit experimental
loading supports headerless bodies or strict legacy-header validation against
configuration. The optional [strict YAML frontend](../../docs/yaml_configuration.md)
loads the same typed objects through `load_yaml_config`; see
[the canonical YAML example](../../examples/qcbed_config.yml). Optional
[legacy import/export tooling](../../docs/legacy_migration.md) supports migration and
compatibility bundles; canonical runs do not require FOR007.

For production runs, use `python -m qcbed.run --config examples/qcbed_config.yml`.
See the [YAML run and resume workflow](../../docs/yaml_runs.md); the legacy
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
[scattering conventions and validation](../../docs/scattering.md) for units, bounds,
and the remaining executable-equivalence limitation.

Immutable 2D lattice and reciprocal grids provide separate dimensionless H/K
indices, physical vectors and STL in inverse angstroms, using the source's
convention without a 2*pi factor. See [reciprocal geometry](../../docs/reciprocal_geometry.md).

`elastic_structure_factors` evaluates the raw complex128 elastic atomic sum
for fixed S atoms, with explicit applied thermal factors and vacancy fraction.
It preserves the positive atomic phase and dilated cell-area scaling; slice
shifts belong to propagation. This is not a complete SFMAP implementation,
and executable-Fortran equivalence remains unestablished. See
[structure-factor conventions and validation](../../docs/structure_factors.md).

Absorptive atomic factors and raw complex structure-factor sums are also
implemented, with explicit phenomenological coefficients, positive applied
Debye-Waller factors, voltage in kV, and literal legacy absorption-flag
semantics. Elastic/absorptive combination follows the source sign convention.
These absorption helpers do not compute propagation. See [absorption](../../docs/absorption.md)
for numerical domains, saved-output comparisons and equivalence limitations.

A serial direct-Fourier reference now reconstructs projected potential and
applies the source PGMAP exponential, with separate DSAP and Nyquist masks.
An explicit NumPy CPU FFT path is validated against that direct sum, including
odd/even ordering and the normalized mathematical inverse. Phase-grating
construction retains `method="direct"` by default; `method="fft"` selects
acceleration with identical masks. Unverified historical build/linkage still prevents
claiming legacy FFT equivalence. See [phase gratings](../../docs/phase_grating.md)
and [FFT backend](../../docs/fft_backend.md) for validation and timing results.

The standalone PROPG kernel now provides source-backed propagation and
fractional slice-shift phases, explicit reciprocal ordering and separate
Nyquist suppression. Internal tilts are reciprocal-basis coefficients, not
angles. No RMT normalization or wave propagation is included. See
[propagator conventions](../../docs/propagator.md) for validation and limitations.

`propagate_one_step` connects a real-space wave, phase grating and reciprocal
propagator using the explicit Python FFT pair, without an extra RMT factor.
It returns immutable real-space and post-propagator reciprocal-wave snapshots.
Complex intermediate comparisons validate this single step against direct DFTs;
historical CDFT sign/normalization equivalence remains unresolved. See
[wave propagation](../../docs/wave_propagation.md).

`simulate_slice_sequence` executes precomputed slice gratings and propagators
for one tilt, using validated zero-based slice IDs and explicit reciprocal
ordering. It returns an immutable final reciprocal wave, requested complex
beam amplitudes and raw intensities, preserving beam order and duplicates.
It omits unused final transforms and uses unit plane-wave illumination by
default. Fitted output scaling remains deferred. See
[multislice sequences](../../docs/multislice_sequence.md).

`prepare_multislice_sample` now builds reusable immutable gratings from parsed
fixed S atoms and explicitly supplied applied parameters. `simulate_one_tilt`
builds per-slice propagators and returns raw beam amplitudes/intensities, using
the common undilated reciprocal-cell tilt coefficients. Preparation supports
direct or FFT reconstruction, a common even mesh and matching undilated lateral
cells, with unity reflection-refinement multipliers explicitly selected.
The real Al input is exercised at its recorded initial parameters; historical
CDFT/executable equivalence remains unestablished. See
[single-tilt simulation](../../docs/single_tilt_simulation.md) for the required
configuration, supported scope and validation.

`simulate_tilt_series_serial` and `simulate_tilt_series_batched` execute
explicit float64 `(n_tilts,2)` arrays of undilated reciprocal-cell H/K
coefficients. The chunked NumPy path transforms only spatial axes and returns
immutable `(n_tilts,n_beams)` amplitudes/intensities in request order.
Serial and independent-DFT tests validate the Python convention; historical
CDFT/executable equivalence remains unresolved. See [tilt series](../../docs/tilt_series.md).

`simulate_tilt_series_parallel` distributes contiguous tilt chunks to spawned
CPU processes, each using the validated batched engine with a prepared sample
loaded once per worker. Process chunk size and internal FFT batch size are
explicit; result ordering is deterministic. Windows callers need a guarded
entry point and should set numerical-library thread limits before launch.
See [CPU parallelism](../../docs/cpu_parallelism.md) for validation, memory costs and
benchmarks; process startup can make one-process batching faster.

## Coordinate and atom API

`compute_legacy_intensity_scaling` implements the three source YNORM branches
as one global float64 factor. `apply_legacy_intensity_scaling` returns scaled,
immutable calculated intensities without changing inputs. Unguarded singular
arithmetic propagates NaN/Inf.
See [intensity scaling](../../docs/intensity_scaling.md) for formulas and runtime limitations.

`compute_legacy_objective` reuses YNORM to compute the source's normalized
weighted root error in beam-major accumulation order, returning immutable
float64 scalars. Refinement and GRADIENT preprocessing remain unimplemented;
historical executable equivalence is unverified. See
[objective evaluation](../../docs/objective_function.md).

The [optimizer architecture audit](../../docs/optimizer_architecture.md) documents
legacy sequencing, parameter flow and missing source. Source-derived
`run_legacy_anneal` now provides seeded float64 adaptive Gaussian search through
the existing trial evaluator. It accepts improvements and ties, with no uphill
acceptance. Historical executable/RNG equivalence is unverified; see
[ANNEAL optimizer](../../docs/anneal_optimizer.md).
`build_legacy_simplex` and `run_legacy_amoeba` implement the recovered modified
Nelder-Mead algorithm through legacy-simplex trial evaluation, with explicit
final-simplex and evaluated-best diagnostics. See [AMOEBA optimizer](../../docs/amoeba_optimizer.md)
for source-specific stopping rules and unresolved historical stride/equivalence
limitations. `run_legacy_paramin` adds source-derived coordinate parabolic
search through physical trial evaluation, with adaptive relative ranges and
fixed or seeded random traversal. See [PARAMIN optimizer](../../docs/paramin_optimizer.md)
for its fallback, stopping and historical-equivalence limitations.
`run_legacy_simpara` adds simultaneous stochastic parabolic proposals with
source direction and temperature adaptation. Its startup direction/counter
are explicitly zero-initialized by agreement; historical executable equivalence
is unverified. See [SIMPARA optimizer](../../docs/simpara_optimizer.md).
`parse_legacy_refinement` and `run_legacy_refinement` now support the ten
source-selected standalone/hybrid modes and full-sequence LOOPS repetition.
The driver preserves retained-state handoff and the original AMOEBA scale
anchor, with immutable stage history. See [refinement driver](../../docs/refinement_driver.md)
for control mapping, supported scope and executable-equivalence limitations.

Refinement parameter layouts, packing/unpacking and explicit physical/simplex
trial evaluation are implemented. Trials reuse the serial numerical kernels
for vacancy-free S-atom samples without surface corrections; unresolved states
are rejected. See [refinement parameters](../../docs/refinement_parameters.md) for
the 18-entry reference layout, supported behavior and equivalence limitations.

`parse_experimental_intensities` reads the experimental header and separate
intensity, variance and weight blocks into immutable float64 `(tilt, beam)`
arrays. Beam H/K order is supplied explicitly; optional scan validation checks
both dimensions. File values are preserved without scaling or weighting.
See [experimental intensity parsing](../../docs/experimental_intensity.md) for the
source mapping and unresolved preprocessing/GRADIENT semantics.

`parse_for007_scan` reads the initial scan basis/corners and only the
NUMP/NUML header of the referenced experimental file. `generate_scan_tilts`
produces float64 undilated H/K pairs in legacy column-fast order for all three
tilt-series APIs. Singleton dimensions are preserved by parsing but rejected
by generation because the source divides by zero. See [legacy scan](../../docs/legacy_scan.md).

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
remain available from the package root. See [source provenance](../../docs/source_provenance.md).

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
See [parser contract](../../docs/legacy_parser.md) for unsupported forms and the
limits of compatibility with the Fortran runtime.

## Directories

- `private predecessor sources`: original Fortran implementation
- `reference_data/`: validation inputs and outputs
- `src/qcbed/`: Python package
- `src/qcbed/crystallography/`: B1 immutable 3-D structure/lattices, ASE frontend and diagnostic visualization
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
and a final HTML report. See [local optimization](../../docs/local_optimization.md)
for editable settings, evaluation-versus-iteration semantics and limitations.

CBED validation now includes source-derived `out.dat` I/O, immutable comparison
maps, per-beam residual diagnostics and self-contained HTML reports. Install
the optional `visualization` extra for Matplotlib. The example
`examples/validate_reference_output.py` evaluates the full reference dimensions
and writes to `tmp/end_to_end_validation/`. Its objective is undefined because
the supplied data contain zero variances; no successful reference refinement
or historical executable equivalence is claimed. See
[end-to-end validation](../../docs/end_to_end_validation.md) for APIs and evidence.

- NumPy/SciPy serial reference implementation
- CPU process-based parallelism over independent incident tilts
- future optional GPU backend

Wave propagation now follows the restored readable FFT signs: negative
exponent divided by the mesh size into reciprocal space, positive exponent
without normalization on return. Potential construction is unchanged.
See [FFT correction validation](../../docs/fft_sign_correction.md), including the
historical-output comparison and unresolved scattering-source provenance.

## Batched refinement execution

Refinement can explicitly reuse the validated NumPy batched tilt engine:
set `LegacyRefinementConfig.execution` to
`ForwardExecutionConfig("batched", batch_size=64)` (for example with
`dataclasses.replace`). Serial remains the default reference path and objective
defaults are unchanged. See [local optimization](../../docs/local_optimization.md#forward-execution)
for the API and [performance measurements](../../docs/cpu_performance_benchmark.md#integrated-refinement-trial-backend).

## Durable full-refinement runner

The editable full-refinement example now selects depth and scan coordinates
with explicit user-selected physical hard bounds, checked after source mapping
and before simulation. Package APIs default to no hard bounds. See
[refinement bounds](../../docs/refinement_bounds.md) for settings and resume rules.
Specific-target, category absolute and category relative rules expand over the
active layout automatically; `--dry-run` displays each resolved bound.

Production runs explicitly reject numerical-invalid proposals with an infinite
score, preserving the retained state and logging the reason. Low-level APIs
default to `invalid_trial_policy="raise"`. This is a Python safety policy, not
proven historical Fortran behavior. Optional category selection can restrict a
first run to depth and scan coordinates; see [invalid-trial policy](../../docs/invalid_trial_policy.md).

Run `examples/run_full_refinement.py --quick` for a short real-data pipeline
check, or `--dry-run` to inspect configuration without optimization. Edit the
settings at the top before running without flags. Use
`examples/run_full_refinement.py --resume <run_directory>` to restart from a
saved retained stage boundary. Defaults use masked objective and batched NumPy
(64 tilts); checkpoints do not claim exact mid-optimizer continuation.
See [full refinement](../../docs/full_refinement.md) for exact interpreter commands,
budgets, previews, output files and restart limitations.


A6 now uses partial-edge/arc normal-direction center voting followed by the existing
reciprocal-lattice consistency test. The old intensity-response candidate detector
and detection-scale controls are removed; A7's own response helper is preserved.
Physical radii remain approximate priors, and complete visible disk edges are not
required. Preliminary arc centers are not final beam geometry.

Canonical central-beam localization now uses known-radius partial disk-boundary
geometry rather than integrated disk intensity, because dynamical diffraction can
bias the internal intensity distribution. `locate_central_disk()` uses an image-center
or explicit neighborhood prior and shared `arc_geometry.py` boundary validation;
legacy intensity localization is isolated in `legacy_disk_localization.py` for
compatibility only, never selected automatically.

Real 256×256 `pattern.dat`, with manually estimated 38/38 px radii (matched profile
unavailable), gives central center **(127,126)** instead of (125,127). The explicit
review rerun revalidates the two prior defining centers and derives canonical
**g1=(81,9), g2=(-10,81)**. A7 confirms one new lower disk at (117,208) within its
unchanged 4 px search; left/diagonal sites remain unconfirmed. Three independent
observations are still insufficient for the four-observation global fit. The
live automatic path now selects a validated pair from three mutually consistent
arcs, recovering (81,9), (-9,82) and nine A9 beams with full ROI frame coverage.
This differs from the explicitly reviewed pair above; no archived geometry is used.
Disk radii remain explicit upstream inputs, never optimized from arcs.

A7 now enumerates every footprint-intersecting predicted site, including higher
orders and partial disks. Confirmation requires configurable usable circumference
(default 25%) plus the unchanged radius/normal/fragment criteria and 4 px search.
The expanded real search records **17 sites**: 3 observed disks, 6 missing, 7 severely
edge-clipped, and the central reference. It finds no fourth observation; global
refinement and seed-to-fit changes remain unavailable. A7 real-data completion is
not claimed. Synthetic clipped fourth-observation tests exercise the unchanged fitter.

The arc-evidence catalogue preserves recorded observations, sources, candidate IDs,
initial/final labels, confidence and residuals. Partial runs are explicitly labeled;
no intensity-only or fabricated observations fill missing sites. Run
`examples/central_geometric_validation.py` with the explicit reviewed centers as
shown in `docs/data_preparation_algorithm.md`. Updated bootstrap/beam CSVs, overviews,
17 crops per set and validation summaries are under ignored
`tmp/a7_boundary_catalogue_r38_final/`. A5/A6 science is unchanged, no inputparas
is introduced, and symmetry-center has been manually validated by the user.


symmetry-center is implemented in `qcbed.data_preparation.origin`: selected beam + selected
beam-local centrosymmetric feature → local 180° symmetry search → legacy origin
formula → zone-axis origin. The symmetry-center symmetry reference beam is selected according
to the relevant zone-axis/diffraction condition. The legacy DM convention of storing
that beam second in `inputparas` is not a physical requirement and is removed from
the canonical Python workflow.

Immutable measurement APIs preserve geometric/feature/symmetry centers separately,
use exact row/column reversal and sum-absolute mismatch, and require explicit
pairwise tolerance before optional equal-weight combination. Fractional DM subarea
rounding is unverified and rejected, never guessed. See
`docs/data_preparation_algorithm.md` for crop conventions and source-equation tests.
`python examples/a8_origin_validation.py --help` exposes explicit beam ID/indices,
feature offset and search/crop sizes. **symmetry-center has been manually validated by the user; A9 does not redetermine its origin.** Existing
A5/A6/A7 numerical behavior is unchanged. No inputparas dependency; A9 extraction is documented below.

symmetry-center now supports mouse selection with `examples/a8_origin_interactive.py`: click a
suitable diffraction-condition beam, drag its local feature ROI, inspect the raw and snapped ROI preview,
then Accept to confirm and run (or Snap to confirm without running). The final symmetry center
comes from the unchanged symmetry-center numerical backend, not the click/rectangle midpoint.
The complete beam-local selection, input hashes and result are saved as
`a8_selection.json`; `--selection FILE --output NEW_DIR` replays without a GUI.
See the interactive symmetry-center section in `docs/data_preparation_algorithm.md` for the exact
real-pattern launch command and controls. Real beam/feature selection remains manual;
A9 is now implemented; see the A9 extraction section below. A10 coordinate conversion is now implemented; see the A10 section below.

symmetry-center symmetry references may be **direct/central or diffracted beams**. The user
chooses a suitable local centrosymmetric feature for the current diffraction
condition; catalogue order is irrelevant. The central reference is eligible by
default, with the same beam-local ROI validation. Provenance records `direct` or
`diffracted`. The unchanged generalized origin formula reduces to the symmetry
center when the reference and central beam centers coincide. A9 is now implemented; see the A9 extraction section below. A10 coordinate conversion is now implemented; see the A10 section below.

symmetry-center fractional mouse ROIs now automatically preview the nearest integer-grid
rectangle alongside the raw selection. Accept explicitly confirms that proposal
and runs symmetry-center; Snap confirms it without running. Raw and accepted bounds remain in
provenance; integer crop semantics are unchanged.

A9 now provides beam-supported extraction through
`qcbed.data_preparation.prepare_experimental_qcbed`. Geometry is a binary union
(beam=1, outside=0) of explicit selected A7 ellipses at 0.9 times physical radii.
`beam_mask` requires no A4 products; `noise_characterization` requires explicit
uncertainty, final A4 mask and provenance. Signal/variance/support are separate:
gradient is eight-neighbor mean minus center, variance is center plus neighbor
sum/64, and coverage is nine-mask sum/9. Local legacy support uses strict
coverage >0.875; conventional selections subtract their minimum before binning.
Weights are support, not inverse variance. The existing objective divides by
variance separately. See [A9 algorithm and limitations](../../docs/data_preparation_algorithm.md#a9--beam-supported-signal-uncertainty-and-weight-extraction).

A reproducible beam-mask diagnostic command (choose a new output directory):

```bash
python examples/a9_extraction_validation.py \
  --pattern reference_data/pattern.dat \
  --catalogue tmp/a7_boundary_catalogue_r38_final/beam_catalogue.csv \
  --beam-ids lattice_0_0 lattice_1_0 lattice_0_-1 lattice_0_1 \
  --radii 38 38 --selection 64 64 --binning 1 --expansion 1 \
  --signal-mode conventional --weight-mode beam_mask \
  --output tmp/a9_beam_mask_review
```

Repeat with `--signal-mode gradient` and a different output path for gradient
review. The catalogue is a generated prior A7 artifact, not a packaged fixture.
Current validation extracted 4/4 beams with full ellipse coverage, 64×64 selections
and 3,207 nonzero weights per beam. These are explicit diagnostic settings.
Artifacts are under `tmp/a9_beam_mask_r38_final/{conventional,gradient}`. Real A4/noise
validation remains blocked; `sigma.dat` is not used implicitly. No uncertainty is
fabricated. A9 preserves predicted centers, stable IDs and catalogue order;
beam-major row/column blocks are available without inventing scan geometry.
Integer selection origins and odd expansion factors are supported; unresolved DM
fractional interpolation/rounding fail clearly. No inputparas or A10 is introduced.


A9 coverage follow-up: a beam need not be fully inside the detector. The default
`minimum_beam_frame_coverage=0.90` retains a selected lattice site when **at least
90%** of its 0.9R ellipse pixels are measured. Coverage uses the full ellipse, not
the enclosing rectangle. Beam=1, outside beam=0; outside detector always has zero
final weight. Fixed-size output arrays use signal0 and, when variance is supplied,
a safe variance1 placeholder only for missing/incomplete detector samples. This
avoids zero-weight 0/0 without changing the objective. Conventional minimum
subtraction uses only measured valid local beam support; gradient stencil validity
excludes incomplete edge samples individually.

Explicitly selected unconfirmed (`missing`) and edge-clipped A7 lattice sites are
eligible without changing their A7 labels. `--all-lattice-sites` audits the full
non-rejected predicted catalogue. The previous four-beam diagnostic selected only
reference/observed sites. The full audit now considers17 sites and retains9, all
with100% ellipse coverage; the remaining8 have0–10.0244% coverage. Retained indices:
(0,0), (-1,0), (1,0), (0,-1), (0,1), (-1,-1), (1,1), (-1,1), (1,-1).
Thus the increase is due to explicit prediction selection, not threshold tuning.
Current outputs: `tmp/a9_coverage90_final/{conventional,gradient}/`; the complete table
and source details are in `docs/data_preparation_algorithm.md`. Earlier four-beam
results are historical. Real A4/noise provenance remains unresolved, sigma.dat is
unused, A5–symmetry-center geometry is unchanged, and A10 coordinate conversion is now implemented; see the A10 section below.

A10 coordinate conversion now uses current symmetry-center/A7 geometry, actual A9 selection
64×64 and user-confirmed **binning2**, following the original DM point equations
and explicit basis ordering. Current corners are:

```yaml
corner_a: [-0.4241467448503984, 0.33558863328822736]
corner_b: [0.33092768004811307, 0.41948579161028415]
corner_c: [-0.33092768004811307, -0.41948579161028415]
```

Historical DM values are comparison references, **not matched-input fixtures**.
Maximum absolute difference is0.0886702083897159: significant span/basis disagreement,
not rounding. No inputs were back-solved or tuned. Full intermediates, component
comparisons and the reproducible command are in
[the A10 section](../../docs/data_preparation_algorithm.md#a10--current-geometry-compared-with-historical-dm-coordinates).
Diagnostics: `tmp/a10_current_binning2/`. APIs: `convert_scan_geometry`,
`scan_from_catalogue`, `assemble_prepared_dataset`. Current binning2 geometry is32×32;
saved A9 arrays remain binning1/64×64 and lack variance, so they are not silently
turned into a ready dataset. A2–A9 science and YAML are unchanged; no inputparas
runtime dependency or post-A10 patch is introduced.

Current A9/A10 physical-selection correction supersedes the earlier64×64 diagnostic
comparison. In the legacy DM workflow, `binning` first supplies the spatial offset
for the eight-neighbor gradient on the original grid, and later supplies final
semi-integration/downsampling. Confirmed **binning2** now controls both; conventional
mode skips the gradient but still bins2×2 after its existing minimum subtraction.

`selx/sely` mean physical QCBED selection, not a Python diagnostic crop. Explicit
`ExtractionConfig.from_disk_geometry(DiskGeometry(38,38), final_binning=2)` derives
**76×76 →38×38** without rounding. `physical_selection_shape_px` is separate from
presentation-only `diagnostic_window_shape_px`; A10 rejects diagnostic-only configs
and mismatched A9 binning. The physical half-span is(76−2)/2=37px.

Both real modes retain the same nine beams with unchanged0.9-radius/90%-coverage
policy. Corrected corners are(-0.50623966320854,0.4005412719891745),
(0.39497819876710266,0.5006765899864682),
(-0.39497819876710266,-0.5006765899864682). Maximum historical-DM difference is
0.012662801232897314 versus the old0.0886702083897159: an85.7192% reduction.
No inputs were tuned/back-solved. Full tables and commands are in
`docs/data_preparation_algorithm.md`. Outputs: `tmp/a9_physical76_binning2/` and
`tmp/a10_physical76_binning2/`. Uncertainty remains unavailable; geometry/sampling
now agrees without inventing variance. A2–symmetry-center science, A10 equations, production YAML,
refinement and optimizers are unchanged. No post-A10 work is introduced.

Opt-in A5/A6 subpixel geometry is now available through `SubpixelArcConfig` and
`refine_arc_center`/`refine_validated_candidate`. Known-radius arc detection first
establishes a robust integer/local center hypothesis. Deterministic bounded
subpixel geometric refinement then improves center precision without intensity
centroids or downstream scan-coordinate targets. Frozen validated supports, fixed
radii, ±1px component bounds and unchanged post-fit arc validation protect the fit;
failed/non-improving fits retain the original center. Original/refined coordinates
and all fit diagnostics are preserved.

The single real-data experiment retained radius38, physical selection76 and binning2.
All three defining centers improved their geometric objectives; the upper beam
reached its x bound. Unchanged A6/A7 processing preserves the nine coverage-eligible
sites and re-confirms the lower site; global refinement remains underdetermined.
Geometry-only A10 maximum historical difference changes0.01266280123→0.01153937494
(8.87% smaller), but component RMS becomes slightly worse. This is modest/mixed
agreement, not a reason to tune centers further. Historical values were used only
for downstream comparison.

**Subpixel A9 extraction remains blocked:** its existing integer-only crop guard
rejects fractional origins. No unvalidated interpolation or snapping was added.
The independently measured direct-beam symmetry-center origin remains unchanged. Results and
three center overlays are in `tmp/subpixel_arc_final/`; reproduce with
`python examples/subpixel_geometry_validation.py --output tmp/subpixel_review`.
Full methodology, shifts, basis changes and A/B/C comparisons are documented in
`docs/data_preparation_algorithm.md`. A7 policies, symmetry-center formula, A9/A10 equations,
physical radius/selection, binning and refinement algorithms are unchanged.


## Patch A closure: prepared dataset to refinement

The canonical center policy is integer/grid geometry; subpixel refinement remains
opt-in experimental diagnostics. Physical radii remain 38/38 px, physical selection
76×76, binning 2 and final beam shape 38×38. No A2–A10 equations were retuned.

`assemble_prepared_dataset()` now retains optional variance and attaches ordered
lattice coordinates, center sources, original extraction-block indices, physical
DiskGeometry, extraction settings, weight mode and provenance to the existing
immutable A1 model. `adapt_prepared_dataset()` returns the engine's existing
experimental-data and scan objects, requiring explicit stable-ID-to-HK assignments,
physical scan basis and objective mode. It does not require YAML or FOR007.
The canonical workflow no longer requires `inputparas`.

The flow is raw/aperture data → optional PSF/OTF calibration and CBED correction →
noise-characterization capability → DiskGeometry → arc central localization →
arc reciprocal seed → beam catalogue → interactive centrosymmetry origin →
signal/support extraction → scan conversion → PreparedExperimentalDataset →
existing refinement engine. Calibration availability is not evidence that a real
input pattern underwent those steps.

For the real reference, conventional and gradient each have signal/support arrays
of shape `(1444,9)` (12,996 values), with variance absent. Per-beam nonzero weight
counts are 886 and 809 respectively. The same ordered nine beams and integer-center
A10 corners are retained. Maximum historical DM corner difference remains
0.012662801232897314; exact DM equivalence is not claimed.

**Real objective execution remains blocked by missing validated variance and
physical HK/scan-basis mapping.** This was explicitly accepted as the closure
boundary. Synthetic datasets with declared variance/physical mappings validate
both modes using the unchanged Al reference configuration: conventional objective
1.4134542758264776, YNORM 33263.751315752874; gradient objective
1.8355485479505533, YNORM 498.8165729557085. These are interface smoke results, not
experimental truth; no optimization was run. Direct adaptation and legacy loading
agree exactly in arrays, amplitudes, calculated/objective intensities, YNORM and
objective in the integration tests.

Reproduce the report into a new directory, using the saved reviewed A7/symmetry-center products:

```bash
python -m examples.patch_a_closure --output tmp/patch_a_closure
```

Artifacts: `prepared_dataset_summary.md`, `prepared_dataset_provenance.json`, and
`conventional/` / `gradient/` under that directory. NPY arrays retain float64 values;
JSON records input hashes, software versions, A5–A10 provenance, beam mappings and
readiness/smoke results. Output directories must be new.

Known limitations: A7 has only three independent arc-confirmed observations and
no global fit; A4 legacy mask/gain provenance is unresolved; `sigma.dat` is a
256×256 uncertainty-related map with unresolved sigma-versus-variance meaning;
noise-characterization mode is implemented but not used for real closure;
historical geometry/calibration differences remain; subpixel is opt-in; fractional
extraction origins remain unsupported. Neither Patch B nor new GUI work is included.


## Patch B1: crystal structure and lattice geometry

`qcbed.crystallography` adds immutable `DirectLattice`, `ReciprocalLattice`,
`AtomicSite` and `AtomicStructureModel`, an ASE input adapter, and an optional
Matplotlib 3-D unit-cell/supercell viewer. Both lattices store vectors as **rows**:
direct units Å, reciprocal units Å⁻¹, **no 2π**. Fractional atom positions remain
unwrapped and ordered. Left-handed/singular cells are rejected without repair.

Install the frontend with `pip install -e '.[crystallography]'` (ASE is also in
`environment.yml`). Run `python -m examples.b1_crystal_structure` to generate a
2×2×2 FCC Al figure and print geometry. See [crystallography conventions and
limitations](../../docs/crystallography.md).

B1 is additive: Patch A `(n1,n2)` is not physical HK, and existing experimental
geometry, multislice, refinement and configuration semantics are unchanged.
B2 now adds arbitrary integer zone axes and Cartesian projection (below).
B3 now supplies exact recast-cell transforms (below); HKL enumeration and beam matching remain future work.


## Patch B2: zone-axis orientation and projected geometry

`qcbed.crystallography` now provides `ZoneAxis`, `CrystalOrientation`,
`zone_axis_cartesian()`, `zone_law_value()`, `project_atomic_structure()`,
`project_hkl()` and `plot_projected_atomic_structure()`, with immutable projected
site/vector records. Zone axes reduce by gcd and make the first nonzero component
positive. Real-space direction is `[u,v,w] @ A`, not the integer triple treated as
Cartesian. Beam-basis rows form right-handed R; row vectors transform as
`v_beam = v_cart @ R.T`, with inverse `v_cart = v_beam @ R`.

Projection preserves original site index, source site, replication offset and
depth; no columns are merged. B1 and B2 share direct-cell replication. The beam
frame is an orthonormal Cartesian view, **not** a physical integer indexing basis;
Patch A `(n1,n2)` is not mapped to HK. No reflection enumeration is included.

Run `python -m examples.b2_zone_axis_projection` for FCC Al [001]/[110]/[111]
figures and numerical diagnostics. The [crystallography documentation](../../docs/crystallography.md)
specifies tie handling, matrix/units conventions and limitations. No new dependencies,
Patch A science, multislice or refinement changes are introduced. B3 is documented below; B4 is not started.


## Patch B3: exact zone-cell recasting and original-hkl identity

`qcbed.crystallography` adds `ZoneAxisCell`, deterministic extended-Euclid
unimodular completion, exact original-hkl ↔ recast-HKL transforms, explicit internal
`ZoneReflectionHK`, original-reflection CBED basis algebra, and full structure
recasting without wrapping, replication or physical displacement.

**Users input original-cell (h,k,l), not internal recast (H,K).** The four distinct
concepts are original Miller hkl, two user original-hkl basis reflections,
experimental (n1,n2), and internal recast (H,K,0). Future B4 will explicitly form
`hkl=n1*hkl_1+n2*hkl_2` and then convert to internal HK0; B3 only implements the
integer algebra, with no Patch A wiring or configuration changes.

For direct rows `A_zone=T@A_original`, T has third row canonical [uvw] and exact
determinant +1. Reciprocal rows satisfy `B_zone=inv(T).T@B_original`, without 2π.
Indices obey `HKL=hkl@T.T`; fractional positions obey `f_zone=f_original@inv(T)`.
L equals the original zone-law integer exactly. A and B need not be perpendicular
to C; reciprocal A* and B* span its perpendicular plane.

Run `python -m examples.b3_zone_axis_cell` for [001]/[110]/[111] examples and
`tmp/b3_zone_axis_cell/zone_cell_diagnostics.json`. See [crystallography](../../docs/crystallography.md)
for derivations, sign conventions, exact [110] mappings and measured errors.
No new dependencies, Patch A/refinement changes, reflection enumeration or B4 work.


## Patch B4: explicit experimental-beam indexing

`IndexedBeamReflection` preserves stable ID, experimental (n1,n2), user-facing
original hkl and derived internal `ZoneReflectionHK`. The correspondence is explicit:
experimental g1↔original basis_hkl_1 and g2↔basis_hkl_2, giving
`hkl=n1*basis_hkl_1+n2*basis_hkl_2`, then B3's internal HK0. No sign/swap search,
automatic indexing or detector calibration is performed. Pixel g vectors and
Å⁻¹ crystal g vectors have different units.

`index_prepared_dataset()` adds a separate immutable annotation wrapper; the A1
dataset and scientific arrays remain unchanged. `adapt_indexed_prepared_dataset()`
passes internal HK to the existing adapter and derives scan HK coefficients using
the explicitly supplied A10 basis transform, preserving corners. It does not build
or validate a matching recast multislice sample, and no objective is executed.

B4a backend/adapter integration is complete. B4b production configuration integration
is explicitly deferred: current sample config contains 2-D slices, not the original
3-D lattice, and scan.basis_1/basis_2 already mean legacy HK coefficients. No YAML
fields, checkpoint identities or compiler semantics are changed. User basis input
is the validated typed `OriginalReflectionBasis`; users do not calculate internal HK.

Run `python -m examples.b4_beam_indexing` for the ordered nine-beam [110] example,
CSV/JSON records and cubic/orthorhombic/hexagonal/triclinic checks. This is a supplied
basis demonstration, not a confirmed indexing of pattern.dat. See
[crystallography](../../docs/crystallography.md) for the complete table and boundary audit.

The B4 instruction confirms sigma.dat as legacy per-pixel sigma (future variance=
sigma²), superseding the earlier uncertainty-kind question. Loading it and checking
its preparation alignment remain separate; B4 neither uses it nor invents an A4
noise mask. Real objective still needs the confirmed experimental crystal/basis,
variance integration and a compatible recast simulation sample. Four-index
Miller–Bravais input is deferred. No optimization or further stage was started.

### B5 strict simulation-cell geometry

B5 uses strict transverse-periodicity checks before explicit depth slicing. B3
`ZoneAxisCell` is an indexing/recast basis, distinct from
`MultisliceSimulationCell`. The restricted constructor supports exact integer C
shears; unsupported geometries are rejected, without full-period projection.
`build_depth_resolved_sample()` retains original site provenance and requires
explicit depth boundaries. `build_zone_axis_simulation_sample()` prepares the
supported subset through the existing propagation backend with explicit scattering
settings. Production configuration/refinement integration remains deferred. See [the B5 audit](../../docs/b5_simulation_sample_audit.md)
for conservative compatibility rules and remaining limitations.

### B4b crystallography configuration

Original structure + original zone axis + ordered original-cell `g1_hkl/g2_hkl`
can now compile through B1–B5 into the existing production sample and indexing
path. Users do not enter recast HK. See `examples/crystallography_config.yml` and
[the schema/workflow](../../docs/crystallography.md#b4b-production-configuration-integration).
Strict B5 rejection is enforced; no projection fallback exists. The legacy
pre-sliced path remains supported. GUI fields are deferred to B4c. This patch
stops at compilation/snapshot integration, without real variance or objective work.

### Real A9 variance closure

The current reference `sigma.dat` is replayed as saved per-pixel sigma, then squared
through the same uncertainty interface used by A4. Future production remains
corrected CBED → A4 noise characterization → sigma → variance → A9; saved sigma is
not a required production input. Conventional and gradient variance arrays are
both `(1444,9)`, positive and finite on active support. Signal/support and A10
corners remain unchanged. Run the no-objective diagnostic with
`python examples/a9_real_variance_closure.py --output tmp/a9_real_variance_closure`
using a new output directory. See [formulas, provenance and measured statistics](../../docs/data_preparation_algorithm.md#real-experimental-variance-closure--a4-interface-and-saved-sigma-replay).
No reference A4 refit or real objective was run. The subsequent B6a closure below uses confirmed physical indexing.

### B6a input preflight

A read-only preflight now checks B4b sample/indexing provenance, reciprocal and
zone-law invariants, Patch A array/variance/beam-order regression, A10 corners and
B4a adapter compatibility before propagation. Use
`python examples/b6a_real_input_preflight.py --synthetic --output tmp/b6a_synthetic_input_preflight`
for an explicitly synthetic demonstration. **Real B6a now PASSes** with confirmed
conventional FCC Al [001], original-cell (200)/(020), and the separate
`examples/fcc_al_real_preflight.yml`. Run `python examples/b6a_fcc_al_preflight.py
--output tmp/b6a_real_input_preflight` with a new output directory. The existing
development example remains unchanged. B6b is deferred. See [B6a gate and commands](../../docs/b6a_input_preflight.md).
No real propagation, YNORM, objective, optimization, GUI changes or app rebuild.

### B6b — fixed-input real objective closure

The confirmed conventional FCC Al [001], (200)/(020) case now executes one
primary evaluation per signal mode plus an exact deterministic repeat. Conventional
YNORM=14870390.0181119, objective=0.939675432747781; gradient
YNORM=13447206.106343817, objective=1.1896911365294378. Support remains
7974/7281 points. No inputs were tuned and no optimizer was invoked.

Run `python examples/b6b_real_objective.py --output tmp/b6b_real_objective`
with a new output directory. See [B6b diagnostics](../../docs/b6b_real_objective.md)
for provenance, residuals and limitations. B7 refinement, GUI and app rebuild
remain deferred; the current active parameter layout is empty.

### B7a — physical site-group refinement binding

The real FCC Al case now has one explicit original-cell group, `Al_fcc`, owning
four fully occupied sites with shared B=0.86 Å² and absorption(0,1,1). The existing
dynamic layout resolves one physical B target with project bounds[0.5,1.5] Å².
Baseline B6b and identical-vector reapplication match exactly. A transient
B=0.8602 perturbation changes scattering and gives objective0.9396971244307921.
No optimizer ran. Site-group Run/Resume remains explicitly disabled pending B7b.

`python examples/b7a_refinement_layout.py --output tmp/b7a_refinement_layout`
writes diagnostics and a YAML snapshot. See [B7 bring-up](../../docs/b7_refinement_bringup.md)
for the target-family audit, provenance, bounds and checkpoint boundary.

### B7b — controlled B-only refinement

Dedicated backend PARAMIN Run/Resume now validates only `Al_fcc` Debye–Waller B
in conventional mode, bounds[0.5,1.5] Å². Nine proposals/two source cycles give
best B=0.7761499999999999 Å², YNORM13982329.708364274 and objective
0.9304373947054256; this pilot is **not converged**. Real interrupted/resumed
history matches uninterrupted history exactly. Source B remains0.86.

Both conventional/gradient masks are binary inclusion gates with active weights
exactly1; variance remains the statistical denominator. Fractional support is
thresholded before objective evaluation, never used as continuous confidence.
Only conventional mode is optimized.

Use `python examples/b7b_single_target_refinement.py --output tmp/b7b_single_target_refinement`
with a new directory. See [B7 bring-up](../../docs/b7_refinement_bringup.md) for budget,
checkpoint replay, cooperative cancellation and limitations. General grouped GUI/CLI
runner integration, extra targets and application rebuilding remain deferred.

### B7c — controlled requested-depth/B refinement

Requested depth is now audited against the strict whole-slice expansion. At the
unchanged baseline (840 Å, B = 0.86 Å²), two stored slices expand to 415
occurrences and 840.1674999999964 Å effective depth. Requested depth is a physical
input; effective depth is a discrete derived quantity, with objective plateaus.
The B7c backend enumerates three adjacent sequence states and runs unchanged
one-dimensional PARAMIN for B at each (15 proposals total). It does not fit a
smooth parabola in depth or claim convergence. Source ordering is depth then B;
no additional targets or GUI controls are enabled.

```bash
python examples/b7c_two_target_refinement.py --output tmp/b7c_two_target_refinement
```

See [B7c semantics, bounds, checkpoint and validation](../../docs/b7_refinement_bringup.md#b7c--requested-depth-audit-and-controlled-bdepth-profile).
The command also checks a separate interrupted/resumed run and writes diagnostic
snapshots without changing the canonical baseline. Lattice dilation remains deferred.

### B7d — QCBED ANNEAL followed by local refinement

The backend diagnostic uses the existing QCBED `ANNEAL` (HA/global stage), then
existing `AMOEBA` for the stepped requested-depth response. Only requested depth
and shared Al_fcc B vary. The informed FCC Al window is1000–1100 Å; general legal
depth bounds remain600–1200 Å. Canonical depth840/B0.86 remains unchanged;
ANNEAL starts separately at1050/B0.86 with three explicit seeds. Requested depth
and effective whole-slice depth are always reported separately.

```bash
python examples/b7d_anneal_local_refinement.py --output tmp/b7d_anneal_local_refinement
```

The command checks baseline/masks, records global/local histories and diagnostic
snapshots, tests replay, and stops local refinement if the global result hits the
window boundary. See [settings, source semantics and limitations](../../docs/b7_refinement_bringup.md#b7d--existing-qcbed-anneal--local-amoeba).
This is multimodal thickness/B workflow validation, not a convergence claim.

Measured B7d workflow validation passed: best ANNEAL objective0.8588311455,
local objective0.8541807284. Real replay was exact for both stages;81 focused and
2816 full-suite tests passed. Seed dependence, lower-window proximity and a small
local B simplex prevent a physical-convergence claim. Much of the improvement
comes from the direct beam; four other beams worsen versus baseline. See the
B7d report for the full per-beam comparison and requested/effective depths.

The B7e lattice-dilation proposal concluded **audit-only; lattice dilation deferred**. The legacy
DILATION is vacancy-coupled, not an independent original-cell lattice scale, and
no validated bounds for the latter were recovered. The gradient baseline remains
7281 active points with binary weights. See [B7e audit](../../docs/b7_refinement_bringup.md#b7e-gradient-mode-dilation-audit--deferred)
and `python -m examples.b7e_gradient_audit --help`. No three-target optimization
is enabled by this audit.

The B7e continuation runs **gradient depth + Al_fcc B only** with the existing
ANNEAL → AMOEBA engines:
`python -m examples.b7e_gradient_depth_b_refinement --output tmp/b7e_gradient_depth_b_refinement`.
Use `--resume` for compatible stage journals. Gradient support is binary inclusion
(7281 active weights1), with variance as denominator. Independent lattice scale
remains deferred. See [refinement details](../../docs/b7_refinement_bringup.md#b7e-continuation-gradient-depthb-refinement).

Measured B7e gradient depth/B workflow: PASS; local objective **0.9121954171202451**, requested/effective depth **922.084336/923.172 Å**, B **0.5000000942 Å²**. B is lower-bound-adjacent, so optimizer stopping is not physical convergence. Both real stage replays match exactly; all9 beam contributions improve from the gradient baseline. Full suite: **2829 passed**. Detailed results and limitations are in `docs/b7_refinement_bringup.md`.


## Modern crystallography GUI and Results integration

The source GUI now opens crystallography-derived canonical configurations, shows
atomic/projection views and physical site groups, and edits requested depth,
shared B, original hkl basis, signal mode, bounds and existing optimizer controls.
Modern production accepts any nonempty supported subset of depth, site-group
Debye–Waller B / absorption, scan coordinates, and voltage without absorption.
ANNEAL, AMOEBA, PARAMIN and SIMPARA use the dynamic layout with batched, binary-support
masked objectives. Unsupported physical mappings remain explicit; legacy DILATION
is vacancy coupling, never independent lattice scaling.

Canonical snapshot Run/Resume/Cancel preserves ownership-aware journals. Results
remain available after Cancel: latest, best and retained final are distinct;
recorded objective/YNORM/parameter histories and masked experimental/calculated/
residual patterns require no numerical rerun. Legacy workflows remain separate.
See [workflow, capability audit and execution limitations](../../docs/modern_grouped_gui.md).
See packaging/macos/VALIDATION.md for the current packaged validation.

### Apple Silicon standalone application

The existing `packaging/macos/` build now includes modern crystallography (ASE
readers and resources), grouped Run/Resume/Cancel, and recorded Results after
Cancel. Paths remain relative to the external configuration document, never the
app bundle. Source/CLI/editable workflows are unchanged. See the
[build and signing procedure](../../packaging/macos/README.md) and
[measured acceptance record](../../packaging/macos/VALIDATION.md) for the local
`dist/QCBED.app` build, test counts and distribution limitations.

### A/B/C source configuration interface

The editor now groups controls into **A. Data Prep**, **B. Atomic Model**, and
**C. Main QCBED Parameters**. Data Prep links experimental pattern/sigma/weight
views, supports non-destructive manual weight exclusions and recorded Patch A
full-pattern overlays. Atomic Model provides Front/Top/Left views and physical
site groups. Fresh Run has an explicit unused-path action; Resume protection is
unchanged. See [workflow and mask persistence](../../docs/gui_abc.md).
The rebuilt standalone app includes this interface; see [packaging validation](../../packaging/macos/VALIDATION.md).

The GUI uses three sequential tabs: **A. Data Prep** (select/confirm zone axis,
inspect and mask experimental inputs), **B. Atomic Model** (validate structure),
and **C. Main QCBED Parameters** (validate run settings). Validate A → B → C
before Run; C validates simulation, targets and optimizer settings independently of
output selection. D owns output validation, Advanced execution settings and Run/Resume
readiness. Earlier scientific changes invalidate downstream confirmations. Fresh Run can
reuse any writable selected output directory: validate it separately in D, then
Fresh Run creates a unique QCBED run subfolder and preserves all existing files.
Resume separately locates a previous checkpoint and checks compatibility. See
[GUI workflow](../../docs/gui_abc.md#sequential-validation-workflow).

Modern YAML now follows `data_prep` → `atomic_model` → `qcbed`. See the compact
[FCC Al example](../../examples/fcc_al.yml) and [schema guide](../../docs/modern_yaml.md).
**User YAML = inputs; compiled snapshot = derived runtime state; checkpoint/results
= execution state.** Modern Save/Save As migrates verbose inputs automatically;
legacy inputs remain supported. Keep referenced preparation assets with the YAML.
The [rebuilt macOS app](../../packaging/macos/VALIDATION.md) supports this schema;
packaging validation uses initialization only, without rerunning refinement.

A. Data Prep now supports **New Data Preparation** from a full CBED and optional
aperture, with threshold/manual disk geometry, a confirmed symmetry ROI, live
Patch A processing in a separate process, inspection and weight-only masking.
**Load Previous Preparation** remains available. See the
[live preparation workflow](../../docs/gui_abc.md#live-experimental-preparation) for
calibrated versus Poisson fallback provenance and current detection limitations.
No refinement campaign is part of this frontend change.

Modern refinement supports dynamic subsets of implemented physical targets; see
[the catalogue and optimizer semantics](../../docs/modern_refinement_catalogue.md) and
[the user guide](../../docs/user_manual.md). Unsupported modern physical mappings remain explicit.

The modern GUI separates [Basic, Advanced and Expert controls](../../docs/gui_parameter_classification.md).
See the [user guide](../../docs/user_manual.md) for physical target selection with
recommended optimizer metadata and explicit numerical-profile validation.

**D. Run & Results** provides the visible execution dashboard, guarded Run/Resume,
Cancel and access to recorded results. Fresh projects start at requested thickness
1000 Å; loaded values and independently edited basis/thickness fields are retained.
Configured hard bounds appear before target selection; missing canonical bounds
are identified explicitly. See the [hard-bound audit](../../docs/modern_refinement_catalogue.md#hard-bound-audit-and-immediate-target-display).

Simulation Setup always shows required **Slice boundaries**. `0, 1` now selects
one stored slice per proven simulation-cell period; finer partitions remain
supported. **Transverse simulation mesh X/Y** stay Advanced and sample the
zone-axis-recast P/Q cell, with compiled dimensions and lattice-direction step
lengths shown when available. See the [backend and mesh audit](../../docs/one_slice_simulation.md).
