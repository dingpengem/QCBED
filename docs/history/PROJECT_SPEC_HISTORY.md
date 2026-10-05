# Historical milestone record

This is an archived development record, not the current installation or GUI guide.
See [documentation navigation](../README.md) for current instructions. Measurements
and limitations below belong to their recorded milestone/configuration.

# Multislice QCBED Python Rewrite

## Reconstruction milestone closure

The current reconstruction milestone is accepted. See the [release summary](../../docs/reconstruction_milestone.md)
for supported workflows, accepted test counts, limitations and recorded build provenance.
Release closure adds documentation only; no numerical or GUI workflow changes.


## Objective

Rewrite the existing Fortran multislice QCBED program in Python.

Current priorities:

1. Reproduce the Fortran numerical results.
2. Create a readable and modular implementation.
3. Parallelize independent incident tilts on CPU.
4. Optimize FFT-based propagation.
5. Preserve a path for future GPU acceleration.

## Experimental-data preparation boundary (Patch A)

Patch A core backend closure is complete: canonical dataset assembly and synthetic
refinement integration pass. Validated variance and confirmed conventional FCC Al
indexing pass real B6a; B6b fixed-input objectives now pass for both modes.
Refinement remains deferred to B7.
A0 specification and A1 typed backend models are complete. `qcbed.detector`
represents optional detector spatial-response calibration; `qcbed.data_preparation`
represents empirical noise coefficients, ordered beam/ellipse ROI records,
independent central-disk and zone-axis coordinates, scan coefficients, and prepared
intensity/variance/weight arrays. These models do not use FOR007 as an intermediate.

A2 implements aperture-only PSF/OTF/MTF calibration in `qcbed.detector.psf`,
with NumPy array input, explicit centering, candidate quality selection and safeguards.
A3 implements whole-pattern complex-OTF correction in `qcbed.detector.correction`.
A4 implements legacy noise characterization and variance maps in
`qcbed.data_preparation.noise`. A5 implements explicit profile geometry and central-disk localization in
`qcbed.data_preparation.disk_detection` and `central_arc`, using shared `arc_geometry`. A6 implements inputparas-free initial pixel-space lattice seeds in
`qcbed.data_preparation.lattice_detection`. A7 adds global lattice refinement and complete predicted beam catalogues in
`qcbed.data_preparation.lattice_refinement`. symmetry-center implements selected local centrosymmetry origin measurement. A9 extraction, A10 scan conversion, canonical dataset assembly and direct refinement
adaptation are implemented. Physical HKL mapping remains explicit caller input;
legacy compatibility is preserved. A preparation GUI is outside this closure. PSF correction is not mandatory; unresolved
legacy details remain explicit in [the specification](../../docs/data_preparation_algorithm.md).
This layer does not change validated propagation, objectives, or refinement.

## Known Fortran entry point

- Main source file: `private predecessor implementation`
- Main program name: `REPLACE_WITH_PROGRAM_NAME`
- Primary input file: `REPLACE_WITH_INPUT_FILENAME`
- Primary output file: `REPLACE_WITH_OUTPUT_FILENAME`

## Current execution procedure

Describe the existing procedure here:

1. ...
2. ...
3. ...

## Known physical workflow

The expected high-level workflow is:

1. Read microscope and specimen parameters.
2. Read atomic structure or potential parameters.
3. Construct the projected potential.
4. Construct transmission functions.
5. Construct reciprocal-space propagators.
6. Generate the incident wave for each tilt.
7. Propagate through successive slices.
8. Calculate the exit-wave diffraction intensity.
9. Integrate or save QCBED intensities.

Codex must verify this workflow against the source code.

## Parallelization target

Each incident tilt is assumed to be independent.

The intended CPU design is:

- one process handles one or several tilts;
- each process runs sequential multislice propagation;
- large tilt-independent arrays should be read-only and shared or loaded
  once per worker;
- internal FFT/BLAS thread count should be limited to one.

Codex must verify that no mutable state couples different tilts before
parallelizing them.

## Numerical validation

Reference cases are stored under `reference_data/`.

Important comparisons include:

- projected potential;
- transmission functions;
- wavefunction after selected slices;
- exit wave;
- diffraction intensity;
- integrated disk intensities.

## Known conventions

Fill in everything currently known:

- Electron energy:
- Wavelength convention:
- Coordinate convention:
- Tilt convention:
- FFT sign:
- FFT normalization:
- Array order:
- Boundary conditions:
- Absorptive potential:
- Units:

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
unrestricted automatic A6 rerun returns three candidates and remains blocked by
existing seed rules; it is not silently converted into the reviewed two-center run.

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

The symmetry-center interaction layer is now implemented separately in
`data_preparation/origin_selection.py` and `examples/a8_origin_interactive.py`.
Canonical workflow: user selects a diffraction-condition beam → boxes a local
centrosymmetric feature → confirms beam-local ROI → unchanged symmetry-center refinement →
zone-axis origin. Nearest eligible catalogue association preserves stable beam
identity and geometric center. Explicit grid conversion, validation and preview
precede acceptance. Versioned JSON records retain full selection/input provenance
and support headless replay with numerical equivalence tests. No real beam/ROI is
chosen automatically; real-data symmetry-center has been manually validated by the user. A9 is now implemented; A10 coordinate conversion is now implemented; see the A10 section below.

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

A9 beam-supported experimental extraction is implemented in
`data_preparation/extraction.py`, with optional diagnostic I/O in
`extraction_diagnostics.py`. Explicit stable-ID selection preserves A7 catalogue
ordering and predicted centers. The binary union mask means beam=1/outside=0;
default ellipse radii are 0.9× physical semi-axes. Signal modes are conventional
(with legacy local minimum subtraction) and gradient (eight-neighbor mean minus
center). Gradient variance is Vcenter+sum(Vneighbors)/64, floating support is
sum(nine masks)/9, and final local support threshold is strictly >0.875 after
source expansion/semi-integration. Explicit dimensions/binning/expansion replace
inputparas. Integer origins and odd expansion factors are supported; unresolved
DM interpolation/half-tie rounding is rejected, not guessed.

`WeightMode.BEAM_MASK` uses geometric support without requiring A4 provenance.
`NOISE_CHARACTERIZATION` requires explicit A4 uncertainty and final mask;
weights remain support because the engine already divides by variance.
Absent variance remains absent. A1 beam/ROI models are reused; intermediate
products preserve beam-major row/column ordering without a fabricated
PreparedExperimentalDataset scan. The real 256×256 pattern validation extracted
four selected beams in both modes, full ellipse coverage, 64×64 diagnostic windows,
3,207 nonzero weights each. Outputs: `tmp/a9_beam_mask_r38_final/`. Real noise mode remains
blocked by missing A4 provenance; sigma.dat is not implicitly consumed. Details,
source limitations and numerical results are in `docs/data_preparation_algorithm.md`.
symmetry-center is manually validated. Beam geometry, refinement and numerical kernels are
unchanged. No inputparas dependency or A10 scan/corner conversion is introduced.


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

A10 (`data_preparation/scan_conversion.py`) converts current symmetry-center origin, explicitly
ordered A7 reciprocal geometry and actual A9 selection dimensions using the DM
half-width `(selx-binning)/2` and2×2 solve. Point1/2/3 map to corner_a/b/c, protected
against the production scan generator. Explicit integer basis transforms preserve
sign/order; singular and ill-conditioned bases are rejected.

Current selection64×64 with confirmed binning2 yields corners
(-0.4241467448503984,0.33558863328822736),
(0.33092768004811307,0.41948579161028415),
(-0.33092768004811307,-0.41948579161028415).
Historical DM values are comparison references, not a matched-input regression
fixture. Maximum absolute difference0.0886702083897159 is significant; geometry
was not back-solved or tuned. Diagnostics and all component comparisons are in
`tmp/a10_current_binning2/` and `docs/data_preparation_algorithm.md`.
The confirmed binning2 is recorded separately from saved A9 binning1; arrays remain
unchanged. Dataset assembly requires matching geometry and preserves absent variance;
refinement readiness separately requires validated variance at evaluated points.
A2–A9 scientific behavior is unchanged. There is no inputparas runtime dependency,
physical HKL requirement, or post-A10 implementation.

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

### B4b — production crystallography configuration

Implemented mutually exclusive original-structure and pre-sliced input modes,
optional ASE file resolution, pure canonical compilation, strict B5 validation,
original-hkl CBED basis validation, ordered beam/HK derivation, B4a dataset adapter
integration and YAML Save/Save As. Compiled execution uses the existing sample,
layout and runner; no numerical or checkpoint algorithm changes. See
[the B4b contract](../../docs/crystallography.md#b4b-production-configuration-integration).
The current GUI rejects this mode explicitly pending B4c. Crystallographic legacy
export is deferred; legacy imports do not invent original structure/basis data.
The synthetic example is not real-pattern indexing. Real variance and objective
work remain separate and have not been started.

### Real experimental variance closure

Completed the shared immutable sigma→sigma² interface, exposed existing A4 sigma
map evaluation, and added generic saved-text-sigma replay. Reference sigma.dat
is treated as legacy sigma (historical save statement unrecovered; source/numerical
support), not as the future production input. No A4 gain-threshold change/refit.
Existing A9 conventional/gradient propagation was reused: 2×2 mean variance
sum(V)/16; gradient variance Vcenter+sum(Vneighbors)/64 with ±2 stencil.
Both modes produce `(1444,9)` variance, all active entries positive/finite, with
unchanged 7974/7281 support totals, beam order, signal and A10 corners.
This supersedes earlier reference variance-unavailable status; archived products
remain historical. B6 real objective/optimization is not started. Confirmed
physical zone/basis/sample inputs remain required. Details and replay command:
[variance closure](../../docs/data_preparation_algorithm.md#real-experimental-variance-closure--a4-interface-and-saved-sigma-replay).

### B6a — input closure/preflight

Implemented the read-only production-compiled sample + prepared-dataset gate,
exact basis/beam provenance checks, integer zone law and hkl basis algebra,
float reciprocal/direction checks, baseline array/corner regression and existing
refinement-adapter acceptance. JSON/Markdown PASS/FAIL reports include nine-beam
and per-mode/per-beam diagnostics. Synthetic and confirmed real preflight pass.
The real case uses conventional FCC Al, a0=4.049 Å, full four-site motif, B=0.86 Å²,
[001] and (200)/(020); see `examples/fcc_al_real_preflight.yml`.
Original indices refer to the standard conventional crystallographic cell where
practical; recast/legacy indices remain derived compatibility representations.
No objective or propagation is invoked. See [B6a](../../docs/b6a_input_preflight.md).
B6b, GUI B4c and app rebuilding remain deferred.

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

### B7c — requested depth semantics and mixed refinement bring-up

The shared strict-< cyclic slice mapper preserves sequential addition order.
Requested depth remains continuous input; effective depth and occurrence count
are derived discrete geometry. Baseline requested 840 Å expands two stored slice
types into 415 occurrences, effective 840.1674999999964 Å. Plateaus and exact
transition boundaries are audited before depth activation.

Exactly two targets are enabled in the scoped backend diagnostic: sample requested
depth and Al_fcc shared B, in that source-defined numerical order. Existing
reference requested-depth bounds are [600,1200] Å; B bounds are [0.5,1.5] Å².
Three explicit depth candidates (414/415/416 occurrences), each with five unchanged
PARAMIN B proposals, avoid smooth interpolation of the stepped depth response.
Baseline and a nine-point B/depth grid precede search. No fractional slicing,
lattice dilation, occupancy, absorption, scan, voltage or other target is enabled.

B7c checkpoints preserve mixed-target identity, exact sample geometry and replay
journals. Separate diagnostic snapshots retain the best sampled physical state;
canonical inputs stay B0.86/depth840. Binary support, variance, beam/scan indexing
and objective semantics are invariant. GUI/app work remains deferred. Full details
and measured validation are in `docs/b7_refinement_bringup.md`.

### B7d — controlled ANNEAL → AMOEBA workflow

Existing QCBED ANNEAL provides the HA/global search; existing AMOEBA provides
local refinement without a smooth parabolic assumption for requested-depth
plateaus. Source algorithms, objective, variance/binary support and propagation
are unchanged. Exactly requested depth and Al_fcc B remain free. Current informed
window1000–1100 Å is distinct from legal requested-depth bounds600–1200 Å;
B bounds remain0.5–1.5 Å². Canonical depth840/B0.86 is preserved; explicit midpoint
initialization1050/B0.86 is used for seeds2026–2028. Relative temperatures.05/.0001,
cooling.9 and MAXITS59 are documented diagnostic overrides; no new optimizer.

The lowest valid recorded ANNEAL objective selects handoff. Boundary-touching
sample geometry blocks the local stage pending window review. Requested/effective
depths, occurrence identity, global best/last/retained state, per-beam contributions
and exact deterministic checkpoint replay are recorded. AMOEBA MXF40 is a soft
source budget; physical convergence and basin reproducibility remain separate
from its stopping criterion. See `docs/b7_refinement_bringup.md` for full settings.
No additional target families, GUI or packaging work is enabled.

B7d measured integration status: PASS. All three ANNEAL seeds completed; the
interior lowest-objective state handed off to AMOEBA. Best global/local objectives
are0.8588311454653902 /0.8541807283629408. Both real interrupted/resumed histories
and final optimizer metadata matched exactly.81 focused and2816 full-suite tests
passed. Physical convergence is not established: seeds retain different B states,
two runs are boundary-adjacent, local B motion is small, and aggregate improvement
is dominated by the direct beam while four beam contributions worsen.

### B7e gradient dilation gate

B7e uses gradient mode: 7281 active points, active weights `[1.0]`, canonical
YNORM 13447206.106343817 and objective 1.1896911365294378. The dilation audit is
complete with an audit-only PASS: independent global lattice dilation remains
DEFERRED. Legacy DILATION controls vacancy-induced volume contraction; its 0–1
coefficient domain is not a bound for an original-cell lattice scale. No validated
bounds for that new physical target were found. No three-target optimization or
new checkpoint format is introduced. See `docs/b7_refinement_bringup.md` for
source provenance, gate rationale, and the reproducible fixed-baseline command.

### B7e continuation: gradient depth and group B only

The controlled case workflow now selects gradient mode explicitly through the
existing two-target physical adapter and ANNEAL → AMOEBA engines. No independent
lattice-scale target is introduced; legacy DILATION remains vacancy coupling.
Gradient support >0.875 is binary inclusion (7281 active weights1), with variance
the statistical denominator. Baseline, data/ordering invariants, explicit local
simplex steps, checkpoint replay, and masked diagnostics are checked. Historical
conventional checkpoint identities and numerical behavior are preserved. The
case-specific window is900–1100 Å; general depth bounds600–1200 Å and B bounds
.5–1.5 Å² remain unchanged. See `docs/b7_refinement_bringup.md` and
`examples/b7e_gradient_depth_b_refinement.py`.

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

### Modern macOS packaging closure

The additive native arm64 PyInstaller specification explicitly collects ASE input
readers/spacegroup data, modern crystallography/grouped execution and Results.
Frozen acceptance covers the real gradient binary support, a short ANNEAL/AMOEBA
run, cooperative cancellation, ownership-compatible resume, masked Results,
Save/Save As and legacy compatibility. Export of modern crystallographic provenance
remains explicitly unsupported; its GUI report must not attempt unresolved pure
compilation before reporting that limitation. No physics or additional refinement
families are enabled by packaging. Build fingerprints, ad-hoc signature audit,
acceptance outcomes and remaining manual/distribution limits are recorded in
`packaging/macos/VALIDATION.md`; build artifacts remain outside source control.

### A/B/C configuration presentation

The source GUI organizes existing controls into Data Prep, Atomic Model and Main
QCBED Parameters. Manual exclusions are persisted as new standard experimental
bodies with zeroed weights; signal/variance and source files remain unchanged.
Companion provenance enables restoration and explicit Patch A overview attachment;
no GUI-only numerical mask is applied. Existing data identity protects Resume.
Beam order/duplicates, dynamic target layout and strict B5 behavior remain intact.
Orthogonal atomic views use the existing beam-frame projection. Unique fresh-run
paths preserve overwrite protection. No physics, objective, optimizer or app build
change is included. See `docs/gui_abc.md` for asset retention and unavailable-data
behavior.

### GUI ordered preparation workflow

A/B/C are actual tabs with explicit ordered validation receipts. A includes manual
centrosymmetric-feature selection through the existing symmetry-center/A10 equations, linked
signal/sigma/weight inspection, weight-only exclusions and an immutable prepared
snapshot. B reuses strict crystallography validation; C reuses execution readiness.
Run requires all three; upstream edits invalidate downstream receipts. Existing
output-only paths are accepted by the GUI with explicit backup-and-replace
confirmation. Canonical runner/CLI checks and resume identities remain unchanged.
The packaging smoke for this patch is GUI-only, with no optimizer execution.

### Modern compact input documents

The version-2 A/B/C YAML frontend constructs the unchanged immutable QcbedConfig.
Modern Save emits physical/user inputs, named dynamic targets and an immutable
prepared-input descriptor reference. Mask/origin provenance remains external and
single-source; validation flags are not scientific state. Original sites and
conventional hkl remain user-facing; internal HK, recast structure, effective depth,
layout, YNORM and results remain compiled/run/checkpoint state. Verbose modern and
legacy schemas remain readable; legacy Save behavior is unchanged. Versioned
optimizer defaults preserve exact explicitly supplied settings on migration.
See docs/modern_yaml.md. No numerical algorithm, checkpoint identity or refinement
campaign changes are included, and app rebuilding is deferred.

## Live Patch A frontend

New preparation consumes current raw CBED inputs, not a historical manifest.
An optional aperture supplies threshold geometry and detector calibration through
separate existing functions; an explicit geometry/profile override takes precedence.
Without an aperture, manual geometry and the explicitly labelled Poisson/intensity
uncertainty approximation are supported. Geometry, ROI snapping and refined origin
require confirmation. A5–A10, propagation and optimization kernels are unchanged.
The spawned worker publishes exclusive replay assets; masking modifies only weight.
A validation is independent of an incomplete Atomic Model; B/C validation still
gates execution. Raw YAML drafts are frontend documents, not fabricated compiled
configurations. See [workflow and limitations](../../docs/gui_abc.md#live-experimental-preparation).
The real pattern now yields nine usable live A9 beams via mutually consistent sparse
arc seed selection, with explicit upstream geometry and no reference catalogue
fallback. A7 missing remains nonfatal. No B6/B7 campaign is repeated.

Modern refinement supports dynamic subsets of implemented physical targets; see
[the catalogue and optimizer semantics](../../docs/modern_refinement_catalogue.md) and
[the user guide](../../docs/user_manual.md). Unsupported modern physical mappings remain explicit.
