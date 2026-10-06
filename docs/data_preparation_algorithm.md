# Experimental CBED data preparation — Patch A, A0–A10

## Current closure status

Patch A core backend is closed with documented limitations: A5–A10 products
assemble into the existing immutable `PreparedExperimentalDataset`, and the direct
adapter passes synthetic end-to-end tests for both modes. Real objective evaluation
now **passes B6b fixed-input integration** for both modes. Real variance and
confirmed FCC Al physical indexing pass B6a without invented variance or HK.
See [B6b](b6b_real_objective.md); no refinement is started. Earlier stage reports
remain historical evidence.

## Scope and evidence

A0 is complete. A1 implements immutable typed models for a new
**preprocessing/backend layer**, with a GUI to follow later. A2 now implements
aperture-derived detector PSF/OTF/MTF calibration. A3 implements whole-pattern
CBED OTF correction. A4 implements legacy noise characterization and variance
maps. A5 implements disk-profile geometry and direct-disk localization. A6 infers
an initial pixel-space reciprocal lattice without pre-specified beam coordinates.
A7 refines the global pixel lattice and generates a predicted beam catalogue. No step
changes canonical configuration semantics. Validated
multislice, propagation, objective, optimizer, refinement, and checkpoint code is
outside Patch A's scope. **FOR007 is not an intermediate representation.** No
automatic physical HKL indexing is proposed. symmetry-center implements explicitly selected local centrosymmetry origins; symmetry-center has been manually validated by the user. A9 extraction is now implemented; A10 coordinate conversion is now implemented.
Detector-response metadata is optional: **PSF correction is not established as a
mandatory CBED-processing stage**. Other A0 ambiguities below remain unresolved; A3 records the user-confirmed correction method.

This specification distinguishes the requested scientific workflow from behavior
actually visible in three supplied DigitalMicrograph (DM) source files. Two
workflows are involved: aperture calibration and CBED preparation; preparation
has both a combined script and a standalone extraction variant.

Source identifiers below refer to the inspected external files, not repository
fixtures. Paths are relative to the supplied `DMScripts` directory. Line numbers
are for these exact versions; no DM scripts were copied or modified.

| ID | Source | SHA-256 |
|---|---|---|
| P | `CURRENTADQCBEDScripts/AMMPSFDM/AMMPSF.s` (387 lines) | `a811a0480ecee5ba8f48151f88510719e0775d1e5475203dd45ceafddecc62b6` |
| N | `simplified_finaledition_withipp.s` (847 lines) | `044162c8a38d64a7c0d7d12fda12f23e204cd4feac0550624868fcdcc73e4d2d` |
| E | `CURRENTADQCBEDScripts/IGQCBEDdataprep.s` (235 lines) | `ccb0b697df44cdf63c0ccd37fd1c82778cbd7cf29a59583d11b23d2228af1e32` |

P describes aperture calibration; N combines noise characterization, disk/origin
localization and extraction; E accepts already supplied signal, sigma and mask
images. These are source observations, **not executed DM/Python equivalence**.
In particular, neither N nor E visibly applies P's calibrated detector OTF to
the experimental CBED frame. The user recalled the historical upstream step as
IFFT(FFT(raw CBED)/complex OTF), establishing the A3 operation. This is workflow
testimony, not an operation located in N/E or inferred from variable names.

## Intended end-to-end flow and stage contracts (calibration, correction, noise and central-disk localization implemented)

The correction branch below describes the intended scientific workflow, not a
required A1 state transition. A dataset may carry calibration metadata without
claiming that it was applied.

```text
Focused aperture image (arbitrary shape)
  → PSF / OTF / MTF calibration
  → Detector spatial-response calibration

Raw CBED frame + detector calibration
  → PSF correction / deconvolution
  → Corrected CBED frame
  → Valid-pixel / gain masking
  → Noise characterization → sigma(I) model → Variance map
  → CBED disk geometry estimation → Central disk localization
  → Beam-catalogue alignment
  → Zone-axis / reciprocal-origin determination
  → Reflection ROI / mask construction
  → Conventional or gradient-QCBED signal generation
  → Intensity / variance / weight extraction
  → Pixel-space scan geometry
  → Reciprocal-basis coefficient conversion
  → corner_a / corner_b / corner_c
  → Canonical prepared experimental dataset
```

| Stage | Explicit inputs | Outputs |
|---|---|---|
| Detector calibration | Focused aperture image, pixel sampling, candidate/acceptance settings | PSF, complex OTF, MTF, normalization/centering conventions, fit diagnostics |
| Detector correction | Raw CBED frame and compatible detector calibration; correction policy | Corrected frame with provenance and unchanged coordinate registration, or an explicit documented transform |
| Mask/noise | Corrected frame, optional supplied validity mask, gain threshold, optional dark sigma, sampling settings | Validity/blemish masks, empirical noise samples and fitted coefficients, sigma and variance maps |
| Disk geometry | Disk-profile image or supplied physical radii, working CBED signal, optional search prior | Physical semi-axes and geometric central center in pixels; boundary support |
| Catalogue alignment | Ordered beam records and old/new central center | The same stable beam IDs/HK assignments, with consistently translated pixel centers |
| Origin | Selected Bragg-condition region, aligned catalogue, central disk center | Independently determined zone-axis origin and symmetry diagnostics |
| ROI/signal | Ordered ellipse ROIs, valid mask, signal, variance, mode, sampling stride | Conventional or gradient signal, propagated variance, coverage information |
| Extraction | ROI centers/sizes, signal/variance/coverage, binning and expansion settings | Ordered per-beam intensity, variance, weight maps |
| Geometry | Zone origin, two reference endpoints/vectors, scan sample-center corners | Basis coefficients for A/B/C; rows/columns and coordinate metadata |
| Prepared dataset | Extracted arrays, ordered beam IDs/HK, scan geometry, mode, provenance | Immutable frontend dataset adaptable to existing canonical experimental inputs |

Images use `(row, column)` arrays; geometric points use `(x, y)` with x increasing
with column and y with row. Do not introduce an implicit y-axis flip, transpose,
rotation or origin shift. Pixel centers versus edges and fractional interpolation
need an explicit convention before implementation. Pixel lengths are in pixels;
A/B coefficients are dimensionless. Intensities retain declared detector units,
sigma has those units, variance their square, and validity/weights are dimensionless.

## Scientific behavior to preserve

### Aperture calibration is separate from CBED correction

P accepts a focused aperture of **any shape**, not necessarily circular. Its
recursive procedure (P:108–387) is:

1. Estimate a sharp aperture by thresholding subpixel-shifted images. The initial
   threshold is randomized; later iterations use half the estimated maximum.
   The sharp estimate is scaled to the measured aperture's total intensity.
2. Form a rough PSF by Fourier division of the measured aperture by the sharp
   estimate, inverse transform, normalize its sum, and handle Fourier centering.
3. Regularize through rotations/reflections and coverage normalization. Separate
   central/core and outer treatments are blended (source radii approximately
   5 and 25 pixels). The result is **not necessarily perfectly radial**.
4. Construct an OTF from the PSF; MTF is its modulus. The source also regularizes
   a modulus-derived transfer function and adjusts an OTF with a negative real
   minimum. These are calibration policies, not a generic Wiener filter.
5. Deconvolve the aperture with that OTF and feed it back into the sharp-aperture
   estimation. Reconvolve the estimate to assess its fit to the measurement.
6. Retain improving candidates and repeat independent measurements until the
   requested number passes acceptance; average the accepted PSF solutions.

For measured aperture `M`, reconstructed aperture `R`, and the source's edge-emphasis
weights `w = abs(1 - abs(R - max(R)/2)/(max(R)/2))`, P:328–341 uses

```text
GOF = sum(w * (R - M)^2) / sum(w)
Poisson mismatch = sqrt(GOF / (sum(w * M) / sum(w)))
```

Within a measurement, improving solutions require mismatch ≤50; inclusion in the
average requires ≤25. This Poisson mismatch is a **calibration diagnostic**, not
justification for replacing the downstream empirical noise model with Poisson
variance. Iterative PSFs are sum-normalized, whereas displayed/averaged PSFs are
peak-normalized; the final `avepsf` is peak-normalized. The calibration object must
identify these conventions explicitly before producing a correction OTF.

Applying the calibration to a raw CBED frame is a distinct operation. Its padding,
OTF-zero handling, regularization, flux/DC behavior, detector sampling, and negative
output policy are not established by N/E. A3 explicitly documents its Python
safeguards and boundary convention below; exact historical settings remain unverified.

### Empirical noise, validity, and gain

The scientific model is

```text
sigma(I) = a*I + b*sqrt(I) + c
V(I) = sigma(I)^2
```

N:315–404 internally fits `ay*(I + bee*sqrt(I)) + cee`, so the public coefficients
are `a=ay`, `b=ay*bee`, `c=cee`, not `b=bee`. The final sigma image evaluates this
expression at `abs(front)`. A4 explicitly preserves this user-confirmed legacy
policy for negative corrected intensities; it is not general-purpose clipping.

N:44–180 accepts a supplied mask or constructs one from positive signal and a
gain-normalization threshold; shifted masks exclude borders. It rejects point
blemishes at 10 standard deviations of the estimated high-frequency residual.
A positive supplied dark-noise value fixes `c`; a negative sentinel requests
fitting. The exact zero-dark-noise case is ambiguous in the source control flow; A4
explicitly treats supplied zero as a fixed zero constant, rather than N’s initialized one.

Noise estimation is not simply `sqrt(I)`: N constructs a normalized 3×3 smoothing
kernel (center 1, axial 0.5, diagonal 0.25), combines repeated smoothing as
`4*S(front) - 6*S²(front) + 4*S³(front) - S⁴(front)`, and scales the residual using
a Fourier-derived reference RMS. It groups masked intensity samples and their
RMS noise with population weights before fitting. The comment about the
“high frequency half” must be reconciled with the actual Fourier crop/centering.
N also warns that for a difference image A−B, its noise input `front` should be A+B.

Crucially, N's later signal is the smoothed image `in`, while sigma is constructed
from `front`; E instead consumes a supplied `front` and sigma image. The intended
corrected-versus-smoothed input to the future pipeline must be decided explicitly.
The internal smoothing kernel is **not** an application of the aperture-calibrated
PSF. The separate image named `Psf` in N is thresholded for disk dimensions, not
used there for detector deconvolution.

### Disk center, beam catalogue, and reciprocal origin

N:417–579 thresholds the supplied profile image, projects its support onto each
axis to estimate disk dimensions, and searches sum-binned images from coarse to
fine resolution. The localization score is integrated signal **inside an
elliptical mask**, not the brightest pixel. The final fine search uses 0.9 times
the estimated semiaxes. Automatic and preset/manual central positions are supported.
This is historical source behavior only. Canonical A5 now uses boundary geometry; the intensity search survives only in the compatibility helper.

`central_disk_center_px` is the center of the direct disk. In N:580–590 the measured
offset from its stored position is added to **every ordered beam position**.
Future records must have stable beam IDs, explicit role selection, and an ordered
`(h, k)` assignment that need not be unique. No indexing is inferred; reordering,
HK deduplication, and using HK as a dictionary key would lose scientific identity.

`zone_axis_origin_px` is independent. N:593–630 locally searches a selected
Bragg-condition region (legacy second beam) for minimum centrosymmetry mismatch:

```text
score(center) = sum(abs(Cut(center) - rotate_180(Cut(center))))
```

The exact crop uses `width-1-x`, `height-1-y` for the 180° comparison. Importantly,
the best symmetry center `s` is not directly assigned as the zone origin. The
source combines it with direct disk center `c` and the aligned selected Bragg
beam center `b` as `origin = s + (c-b)/2`. Basis endpoints are then translated
relative to their prior origin. Record all these points separately.

There is a visible record-layout discrepancy: N writes the new origin at
`2*numbeams+3,+4` and reference endpoints at `+5…+8`, but N's final geometry and E
read origin at `+1,+2` and endpoints at `+3…+6`. Do not reproduce positional DM
records in Python or choose one layout without the actual parameter image.

### ROIs, signal, variance, and weights

The first implementation needs **ellipses only**: circular, polygonal and arbitrary
masks are not required extensions in A0. N:632–643 sums slightly shrunken elliptical
disk masks and multiplies by the validity mask. Summation means overlapping disks
can produce coverage above one; do not silently substitute a Boolean union.
The extraction additionally applies a radial mask of radius `selx/2` centered at
`(selx/2, sely/2)` (N:681–683; E:60–62). This is distinct from the ellipse support.

At the signal-generation stage conventional mode is the working corrected
experimental intensity `I`, with variance `V` and its validity mask. Gradient mode
uses the eight neighboring shifts at **the chosen raw-pixel binning stride** `q`:

```text
G(p) = sum(I(p + d) for d in eight offsets)/8 - I(p)
V_G(p) = sum(V(p + d) for d in eight offsets)/64 + V(p)
```

The source obtains sigma as `sqrt(sum(V_neighbors))/8`, then
`sqrt(sigma_neighbors² + sigma_center²)` (N:685–730; E:64–109).
This is the legacy independent-error propagation, omitting covariance. Deconvolution
and resampling can correlate noise; do not silently add covariance or replace this
formula while claiming legacy equivalence. DM interpolation and boundary behavior
remain to be established. Do not assume they match the backend's zero padding.

Extraction variants differ scientifically:

| Operation | E | N |
|---|---|---|
| Gradient coverage before extraction | Threshold sum of nine masks at 8.9999999 | Sum of nine masks divided by 9 |
| Final coverage threshold | 0.999999 | 0.875 |
| Conventional working signal | Supplied `front` | Smoothed `in` |

Both scripts subtract the minimum of the expanded conventional ROI before
aggregation (E:135–138; N:749–752); gradient ROIs are not minimum-subtracted.
This is a real offset operation, **not** hidden normalization to introduce in the
rewrite. The requested conventional signal definition and this legacy extraction
step must be reconciled via an explicit, documented policy before implementation.

For expansion factor `e` and binning factor `q`, extraction replicates samples via
`Trunc`, aggregates q×q blocks, then takes a center-offset decimation at
`e*index + Round((e-1)/2)`. If `S_I` is the sum of expanded signal samples and
`S_V` the sum of their squared sigma samples within a block, the written operations
are:

```text
intensity = (S_I/e²)/(q²/e²) = S_I/q²
sigma_out = sqrt(S_V/e²)/(q²/e²)
variance_out = sigma_out²
coverage = sum(expanded valid/radial mask samples)/q²
weight = source-specific threshold(coverage)
```

Do not reinterpret replicated samples as independent detector observations or
replace this “semi-integrating” procedure with generic resize/mean defaults.
Equality at mask thresholds follows the literal absolute-value expressions;
rounding, fractional Warp interpolation, cropping and divisibility need fixtures.

The final blocks are **all beam intensity maps, all beam variance maps, then all
beam weight maps**, preserving beam and row order. Variance is not sigma; weights
are not `1/variance`. No additional flux normalization, histogram normalization,
variance flooring, mask repair, or objective scaling is part of the adapter.

### Scan geometry and canonical interface

For central disk center `c`, zone origin `o`, ROI width/height `(W,H)` and sampling
stride `q`, the printed A/B/C locations (E:194–235; N:810–847) are sample centers:

```text
pA = c - o + (-(W-q)/2, -(H-q)/2)
pB = c - o + ( +(W-q)/2, -(H-q)/2)
pC = c - o + (-(W-q)/2,  +(H-q)/2)
v1 = reference_endpoint_1 - o
v2 = reference_endpoint_2 - o
p = A*v1 + B*v2
```

These are not the ROI's outer pixel edges. Store A/B/C as dimensionless coefficients
in that **ordered** basis, with B traversed along columns and C along rows. The
legacy determinant formula for A and division by `v2.y` for B can be replaced by a
2×2 linear solve of the same equation. Reject degenerate/ill-conditioned geometry
explicitly; do not require `v2.y != 0` when the basis is otherwise nonsingular.
Detector-plane vectors in pixels must remain distinct from the canonical scan's
reciprocal H/K basis vectors; the association between each reference vector and
its reciprocal identity must be supplied, not inferred by indexing.

The existing contracts fit this boundary:

- [ExperimentalIntensityData](../src/qcbed/experimental.py) stores immutable
  float64 arrays of shape `(scan_rows*scan_columns, beam_count)`. Flatten as
  `tilt = row*scan_columns + column`; preserve duplicate ordered HK entries.
- [Experimental-data loading](../src/qcbed/experimental_config.py) reads body-only
  intensity/variance/weight blocks with metadata from compiled configuration.
  It performs no preprocessing or normalization; keep it that way.
- [ScanConfig](../src/qcbed/configuration.py) and
  [LegacyScanDefinition](../src/qcbed/scan.py) already express the basis, A/B/C and
  scan dimensions. Their current singleton restrictions remain unchanged.
- The existing [trial gradient](../src/qcbed/refinement_trial.py) acts on calculated
  data when gradient mode is selected. Prepared experimental gradients must not
  be differentiated again by the frontend or input reader.

The intended prepared-data boundary includes stable beam IDs, central center,
origin, reference pixel vectors, masks, noise/calibration provenance, mode and
sampling policies. A1 implements the core models below; full mask/policy payloads
remain for later steps once their scientific semantics are resolved. An adapter can provide the existing
experimental arrays and canonical scan/experiment fields without a temporary
FOR007 or a new physics input path. A1 implements the independent dataset model below, but no canonical adapter or
serialization schema. Detector metadata never asserts that correction occurred.

## Implemented A1 API and future boundaries

A1 adds only `__init__.py` and `models.py` in each of
`qcbed.detector` and `qcbed.data_preparation`. Public imports are from those
packages; nothing is added to the root `qcbed` namespace. No empty numerical
modules are created. The models use NumPy and standard-library frozen dataclasses,
with no YAML, Qt, file I/O, or FOR007 model dependency.

| Public model | Fields / semantics |
|---|---|
| `qcbed.detector.PSFCalibrationMetrics` | Optional `goodness_of_fit`, `poisson_mismatch`, `accepted_measurements`; finite nonnegative diagnostics, not acceptance rules |
| `qcbed.detector.DetectorResponseModel` | `calibration_method`; optional `psf`, `otf`, `mtf`, `artifact_reference`, `metrics`; immutable textual `provenance` |
| `qcbed.data_preparation.NoiseModel` | Explicit `a`, `b`, `c`, `intensity_units`; finite coefficients only, with no fitting/evaluation or sign-policy inference |
| `DiskGeometry` | Positive finite `radius_x_px`, `radius_y_px`: axis-aligned ellipse semi-axes |
| `BeamObservation` | Nonempty stable `id`, finite `(x,y)` `center_px`, optional integer `indices` of length 2 (HK) or 3 (HKL), strict Boolean `enabled`, optional `label` |
| `ReflectionROI` | Stable `id`, `beam_id`, independent `center_px`, `geometry: DiskGeometry`, `enabled`; ellipse only, with room for separate future ROI variants |
| `ReciprocalFrame` | Separate `central_disk_center_px`, `zone_axis_origin_px`, and nonzero `basis_vector_1_px`, `basis_vector_2_px`; basis vectors are pixel displacements, not endpoints |
| `ScanGeometry` | `corner_a`, `corner_b`, `corner_c` finite coefficient pairs, positive integer `scan_rows`, `scan_columns` |
| `SignalMode` | `CONVENTIONAL` / `GRADIENT` enum; required explicitly on a dataset |
| `PreparedExperimentalDataset` | Three arrays, ordered `beams`, ordered `rois`, `reciprocal_frame`, `scan_geometry`, `mode`, `intensity_units`; optional `detector_response`, `noise_model`, and textual `provenance` |

Example model construction (no correction or geometry calculation):

```python
from qcbed.data_preparation import BeamObservation, ReciprocalFrame, ScanGeometry

beam = BeamObservation('direct', (100.25, 120.5), indices=(0, 0))
frame = ReciprocalFrame(
    central_disk_center_px=(100.25, 120.5),
    zone_axis_origin_px=(101.0, 118.0),
    basis_vector_1_px=(50.0, 2.0), basis_vector_2_px=(-3.0, 49.0),
)
scan = ScanGeometry(
    corner_a=(-0.508156, 0.407641), corner_b=(0.407641, 0.508156),
    corner_c=(-0.407641, -0.508156), scan_rows=39, scan_columns=39,
)
```

A dataset has finite float64 arrays of shape
`(scan_rows*scan_columns, len(beams))`. A single defensive byte-backed copy makes
each array immutable even through `setflags(write=True)`; no dtype conversion,
value adjustment, normalization or inverse-variance weighting is performed.
Signed intensities and zero variance/weight are allowed; negative variance/weight
are rejected. Input lists of records/coordinates are defensively converted to
tuples. Array-bearing models use identity equality, consistent with existing
array-bearing models; compare their arrays explicitly rather than relying on
ambiguous NumPy equality.

Beam IDs and ROI IDs are unique within their respective collections; HK/HKL may
repeat. ROI references must resolve to a beam ID. Multiple ROIs per beam and empty
ROI collections are representable. Disabled beams still occupy columns: A1 does
not silently filter data, reassign IDs, or reorder either collection. A prepared
dataset requires at least one beam. Coordinates may be negative or fractional.

PSF/MTF arrays, when supplied, must be nonempty finite float64 2-D grids; OTF is
complex128. Supplied arrays must have matching shapes. This is a full-grid storage
contract, not an FFT implementation or a check that one array is the transform of
another. Compact half-spectrum/unequal-grid representations need a later explicit
extension; metadata-only or opaque artifact references avoid mandatory large
arrays. References are stored verbatim, with no existence check or CWD resolution.
Calibration normalization/sign conventions are not imposed by the model.

`ScanGeometry` follows canonical corner/dimension naming but intentionally does
not reuse the entire `ScanConfig`: that type also requires reciprocal HK basis
assignments, which cannot be inferred from pixel displacements. Singleton scans
can be represented as data; this does not relax engine readiness restrictions.
Likewise, finite nonzero collinear measured vectors are representable; singularity
and conditioning checks belong to the later conversion/readiness stage, not an
implicit geometry calculation in A1. The meaning of origin records, smoothing,
baseline subtraction, noise-fit constraints, mask thresholds and PSF application
is not decided by any constructor.

The remaining module split is a proposal, except `detector/psf.py` (A2), `detector/correction.py` (A3), and
`data_preparation/noise.py` (A4), `data_preparation/disk_detection.py` and `central_arc.py` (A5), and
`data_preparation/lattice_detection.py` (A6), and
`data_preparation/lattice_refinement.py` (A7), now implemented:

| Boundary (models and `detector/psf.py` exist) | Responsibility |
|---|---|
| `qcbed/detector/models.py` | Immutable calibration data, sampling/normalization conventions and diagnostics |
| `detector/psf.py` | Aperture-derived PSF/OTF/MTF estimation and candidate acceptance |
| `detector/correction.py` | Explicit calibration application to experimental frames |
| `data_preparation/models.py` | Frames, beam-ID catalogue, geometry, policy/provenance and prepared dataset |
| `data_preparation/noise.py` | Valid/gain/blemish masking, empirical fitting, sigma-to-variance construction |
| `data_preparation/disk_detection.py`, `central_arc.py`, `arc_geometry.py` | Profile dimensions and canonical known-radius boundary localization; shared arc evidence |
| `data_preparation/origin.py` | Bragg-region centrosymmetry search and independent origin determination |
| `data_preparation/roi.py` | Ordered ellipse masks and explicit coverage/overlap semantics |
| `data_preparation/signal.py` | Conventional/gradient signals and documented variance propagation |
| `data_preparation/extraction.py` | Sampling, binning, radial coverage and ordered experimental-array adapter |
| `data_preparation/geometry.py` | Pixel sample corners and reciprocal-basis coefficient conversion |

Use explicit inputs and returned snapshots, not shared DM image state. Detector
FFT centering/normalization needs its own specification/tests; do not repurpose
or change the validated multislice FFT conventions to satisfy calibration.

## Replaceable DM mechanisms

Front-image stack order, `HideImage`/`ShowImage`, dialogs, results-window printing,
image names, packed one-row parameter images and scratch-image allocation are
UI/storage mechanisms, not scientific requirements. Replace them with typed
arguments, IDs, diagnostics and named intermediate results. Vectorization may
replace pixel loops only after preserving sampling and arithmetic conventions
within stated tolerances. Literal real*4 storage, tiny denominator constants,
random-search control flow and rotation loops are not automatic Python design
requirements: their numerical/regularization effects still require explicit
validation before substitution. In particular, do not replace thresholding,
FFT padding, interpolation, PSF normalization or noise propagation merely because
the implementation uses an inconvenient DM primitive.

## Historical coordinate comparison reference

The user-supplied historical DM comparison values are:

```yaml
corner_a: [-0.508156,  0.407641]
corner_b: [ 0.407641,  0.508156]
corner_c: [-0.407641, -0.508156]
```

These six-decimal values are historical comparison references, not a matched-input
regression fixture. A10 reports differences from current geometry without tuning
inputs or imposing exact equality.
The corners alone cannot prove intensity, variance, weight or calibration
correctness. Future evidence must include the matching raw/calibration images,
beam/parameter records, sampling settings and prepared arrays. Check stage
contracts, known synthetic geometry, PSF/OTF normalization, variance propagation,
mask coverage, duplicate-beam order and conventional/gradient final outputs.
Do not demand pixel-identical reproduction of every temporary DM image.

## Open decisions and legacy evidence still required

1. Obtain the actual CBED correction script/settings and clarify what N's `Psf`
   input contains. Confirm whether supplied `front` is already corrected.
2. Select the authoritative N/E variant and resolve smoothed versus corrected
   signal, conventional minimum subtraction, and the two coverage thresholds.
3. Obtain the parameter-image schema; resolve the `+1/+2` versus `+3/+4` origin
   offsets, Bragg-beam role, basis-endpoint definitions and projection x/y naming.
4. Establish DM Warp/FFT centering, padding, fractional sampling, pixel centers,
   crop bounds, rounding/tie rules, ellipse overlap, and noninteger dimension
   handling with small exported fixtures. Source comments alone are insufficient.
5. Specify detector calibration transfer across image sizes/pixel samplings,
   DC/flux conventions, OTF-zero/negative handling and failure/iteration limits.
   P's acceptance loop can retry indefinitely; random seeds are not recorded.
6. Resolve noise fitting constraints/seed, zero dark sigma, absolute-value handling
   of negative intensities, difference-image inputs, the Fourier noise reference,
   and the reused `binning` variable (N prompts for it twice before sampling).
7. Decide whether legacy covariance omission is the first explicit supported
   policy; any covariance-aware alternative must be separately selected/tested.
8. Archive approved source/DM exports as future fixtures with checksums and data
   provenance; current scripts are external, and no executable DM equivalence
   or raw-image-to-golden-corners result has been established here.

A1 model contracts are tested in `tests/test_data_preparation_models.py`. A2
adds calibration, A3 adds whole-pattern correction, and A4 adds noise characterization
below. No new dependencies, GUI, automatic indexing, canonical
adapter, or propagation/refinement changes are included. A5 localization and the new A6 lattice-seed architecture are described below; A7 global fitting/catalogue generation is described below; symmetry-center is now documented below; its real-data workflow has been manually validated by the user.


## Implemented A2 aperture calibration

Public imports are from `qcbed.detector`:

```python
import numpy as np
from qcbed.detector import (
    PSFCalibrationConfig, PSFCalibrationError, calibrate_psf_from_aperture,
)

# aperture is an already decoded numerical array, not a filename.
response = calibrate_psf_from_aperture(
    aperture, config=PSFCalibrationConfig(), seed=2026,
)
# response.psf, response.otf, response.mtf, response.metrics, response.provenance
```

`qcbed/detector/psf.py` contains the algorithm. It depends only on NumPy and the
standard library. No image-format loader, HyperSpy dependency, reciprocal-space
CBED logic, or experimental-frame correction is added. The inverse filtering in
this module applies **only to the measured aperture during calibration**.

### Input, frame size and normalization

The input must be a nonempty 2-D real numerical NumPy array of shape `(ny,nx)`.
Integer and real floating dtypes are copied into float64 working storage;
complex, Boolean, object, string/path, non-finite, and non-2-D inputs are rejected.
Nonconstant contrast and positive finite summed intensity are required. Signed
pixels are not silently clipped or background-subtracted. Unusable weighted
means/intermediates fail explicitly. The original array is never mutated.
Square, rectangular, even and odd frames and arbitrary aperture shapes are
supported. The output arrays inherit the frame shape and are immutable.

A centered circle of the configured larger radius must fit within the pixel-center
bounds. With `cy=ny//2`, `cx=nx//2`, the available radius is
`min(cy, ny-1-cy, cx, nx-1-cx)`. This is a **Python safety check derived from the
spatial core radii**, not a DM image-format restriction. The legacy 24.99 px radius
requires each dimension to be at least 51; smaller configured radii permit smaller
frames. Small radius must be positive and no larger than the large radius.

The returned PSF has unit sum. As in P, each accepted candidate contributes a
peak-normalized PSF to the ensemble; the ensemble is then converted to unit sum
for a consistent PSF/OTF/MTF triplet. P instead ends with a peak-normalized display
image. This explicit final scalar conversion is not an assertion of bit-identical
DM output. OTF DC and MTF DC are approximately one. Negative PSF lobes or local
MTF values above one are not clipped; strict monotonicity is not imposed.

### FFT and resampling conventions

- Spatial array indices are `(row,column)`. The PSF center is the pixel
  `(ny//2,nx//2)`; a unit impulse there yields an all-ones OTF.
- FFT origin/DC is `(0,0)`. NumPy's default backward convention is used:
  forward FFT unscaled, inverse FFT divided by `ny*nx`.
- `OTF = fft2(ifftshift(centered_psf))`; `MTF = abs(OTF)`. Both outputs are in
  **unshifted** frequency order, complex128 and float64 respectively.
- Rough-PSF inverse transforms are origin-ordered and shifted to the PSF center.
  Magnitude transfer functions are explicitly shifted for geometric rotational
  averaging about DC, then unshifted for inverse transformation. P rotates
  `RealFFT` output without an explicit shift at that point; its exact DM frequency
  storage/display convention has not been verified. This is a documented Python
  interpretation, not evidence of exact legacy-array equivalence.
- Real image reconstruction explicitly uses the real part of the full-complex
  inverse FFT, corresponding to the real-image role of DM `RealIFFT`.
- Subpixel shifts and rotations use bilinear interpolation with zero extension,
  including fractional coverage at frame edges. Reflections are about the centered
  pixel (`2*cx-x`, `2*cy-y`), not about the half-pixel array edge. For even sizes this
  matches the source's `xrange-x`/`yrange-y` formulas. Odd frames use the declared
  integer center; DM's half-integer center behavior still requires an exported fixture.
- Detector calibration does not call or change the multislice FFT backend.

### Preserved stages and parameters

`PSFCalibrationConfig` exposes the following fields:

| Field | Default | Meaning |
|---|---:|---|
| `small_core_radius_px` | 4.99 | Less-averaged central core, from P |
| `large_core_radius_px` | 24.99 | Core/outer blend boundary, from P |
| `accepted_measurements` | 10 | Required accepted candidate measurements, not merely ten attempts |
| `max_attempts` | 100 | Python cap on the otherwise unbounded candidate retry loop |
| `max_iterations` | 100 | Python cap per candidate, including repeated subdivision levels |
| `fourier_relative_floor` | 1e-12 | Explicit relative Fourier-denominator safeguard |

The algorithm retains these source stages:

1. Subdivision dimension starts at 1 and normally advances by 2 through 13.
   At dimension 1, the measured aperture is restored and the threshold is
   `max(measured)*(1+U)/3` for uniform `U`. Otherwise the threshold is half the
   previous sharp-estimate maximum. Sum the shifted threshold masks, using
   shifts from +0.5 to −0.5 pixels for dimension >1, then match measured total intensity.
2. Fourier-divide measured by sharp aperture to estimate the rough PSF, sum-normalize
   it, and center it. Threshold equality is excluded, matching the source sign expression.
3. Apply successive angle additions `360/2**k` degrees with x/y reflections of
   each updated sum. Levels 1–2 provide the less-averaged core; levels 3–16 provide
   its strongly averaged complement. A rotated/reflected constant field provides
   coverage normalization. This is not a radial-bin average.
4. Rotationally average the rough PSF's magnitude transfer function at levels 1–16,
   inverse transform, then rotationally average that PSF. Blend it with the earlier
   core at the larger radius. P's center-pixel value of 2 at dimension 1 is retained,
   followed by sum normalization; it is not an error-recovery delta-PSF fallback.
5. For a negative minimum real OTF in the source-supported range (≥−1e30), add
   `2*abs(minimum)` to the **complex OTF**, divide by maximum real value, reconstruct
   the PSF, and repeat the same subdivision dimension. A strictly positive minimum
   advances the subdivision by two. The measured aperture is deconvolved by this
   OTF and becomes the next current estimate; measured, sharp and deconvolved
   arrays remain distinct.
6. Reconvolve the sharp aperture and use P's edge-emphasis weights and the GOF and
   Poisson mismatch equations above. Retain only GOF improvements with mismatch
   ≤50; a completed candidate is accepted into the ensemble only if its retained
   mismatch is ≤25. No Wiener, Richardson–Lucy, Gaussian fitting or replacement
   blind-deconvolution method is used.

### Safeguards, failures and provenance

Nonzero Fourier denominators below `relative_floor*max(abs(denominator))` retain
phase but have their magnitude raised to that floor (bounded below by float64's
smallest normal magnitude). Exactly zero bins have quotient zero: those
noninvertible components are explicitly omitted rather than assigned a made-up
phase. An entirely zero denominator fails. This adds a minimal, documented inverse
safeguard to effectively unguarded source division; it is not a noise-dependent
filter. Counts of guarded bins are recorded for successful candidate attempts.

Non-finite intermediate arrays, invalid sums/peaks/fit weights, and unusable OTFs
raise `PSFCalibrationError`. An exactly zero minimum real OTF or a value below
−1e30 is rejected: P's reset/control path does not define a usable positive-OTF
advance there. There is no epsilon-based change to the positive/negative branch.
Candidate failures may be retried up to `max_attempts`. A candidate which exhausts
`max_iterations` before advancing beyond dimension 13 is rejected, even if it had
a promising intermediate state. Fewer than the requested accepted measurements
causes an informative failure, never a partial-success or Gaussian/delta fallback.

Pass either an integer nonnegative `seed` or a NumPy `Generator`, not both. A
supplied generator is consumed; otherwise a local generator is created. Global
NumPy RNG state is untouched. Provenance stores initial generator state, seed,
parameters, grid/centering conventions, guard counts, OTF-adjustment counts, and
per-attempt accepted/rejected status with selected iteration and quality or failure
reason. This does not promise reproduction of DM's `UniformRandom` sequence.

`PSFCalibrationMetrics` reports the **best accepted candidate's** GOF and Poisson
mismatch plus the accepted count. It does not claim a newly computed quality score
for the averaged response. The detailed selections are JSON text in the existing
model's provenance tuple; no new giant-array serialization format is introduced.

### A2 validation status

`tests/test_detector_psf.py` checks independent source equations for threshold
mask combinations, quality weights, OTF adjustment and Fourier guards; centered
impulses and displaced-impulse FFT phase; core/outer averaging stages; shape and
finite-value rejection; radius-derived frame checks; candidate acceptance/rejection
and budget failures; RNG reproducibility; and output consistency/immutability.

Synthetic tests blur a noncircular union of rectangles with a known simple PSF.
On square, rectangular and odd frames, they require centered recovered response,
unit DC, lower high-frequency MTF, and reconvolved agreement substantially better
than the unblurred baseline (RMSE <1 intensity unit at aperture intensity 100).
The synthetic PSF is used only to generate test data, never as a calibration fit.
The default ten-accepted-measurement ensemble is also exercised.

No matching aperture image plus exported DM PSF/OTF/MTF pair was found in the
repository's reference/test data. Exact DM-versus-Python equivalence therefore
remains **unverified**, especially interpolation, Fourier display ordering and
odd-frame centering. Historical scan corners are now compared with current geometry in A10 below;
this is not an A2 regression result. Extraction, reciprocal-origin calculation and
adaptation were implemented in later stages below; preparation GUI remains outside this patch.


## Implemented A3 whole-pattern OTF correction

A0, A1, A2 and A3 are complete within their stated scope. The workflow is:

```text
Aperture image → A2 detector PSF/OTF/MTF calibration
Raw CBED + calibrated OTF → A3 whole-pattern Fourier correction
                         → corrected CBED → A4 noise characterization / variance map
```

Public API (also exported by `qcbed.detector`):

```python
corrected = correct_cbed_with_otf(
    raw_pattern, detector_response,
    fourier_relative_floor=1e-12,
    imaginary_relative_tolerance=1e-12,
    imaginary_absolute_tolerance=1e-12,
)
```

`raw_pattern` is a finite, nonempty real numerical NumPy array in `(ny, nx)`
order. It is copied to float64. Square, rectangular, odd and even frames are
supported. `detector_response` is an existing `DetectorResponseModel` with an
OTF of exactly the same shape. Missing OTFs and mismatched shapes fail; there is
no calibration, loading, interpolation or resampling here. Inputs remain unchanged.
The return is a new float64 array of the original whole-frame shape.

The operation is `ifft2(fft2(raw_pattern) / complex_otf)`, using A2's unshifted
`otf = fft2(ifftshift(centered_psf))`. No shifts are applied to either spectrum
in A3. The PSF center remains `(ny//2, nx//2)`. This is the discrete periodic
Fourier convention, with no padding. The complex OTF preserves phase; MTF is not
used as a divisor. No intensity normalization, minimum subtraction, clipping,
ROI processing, masks, variances or weights are applied. DC transfer determines
intensity scale, including for supplied non-unit-normalized OTFs.

The existing A2 Fourier-division helper is reused unchanged. Nonzero magnitudes
below `max(max(abs(OTF))*fourier_relative_floor, float64.tiny)` are raised to that
floor while retaining phase. Exactly zero bins contribute zero corrected spectrum
(null-space omission); an all-zero OTF fails. Well-conditioned bins use direct
division. This deterministic safeguard is not Wiener filtering and cannot recover
lost frequencies. The relative floor must be finite and strictly between 0 and 1.
Non-finite/overflowing intermediate results fail with informative errors.

The inverse transform is checked before taking its real part. Its maximum
absolute imaginary component must be at most `imaginary_absolute_tolerance +
imaginary_relative_tolerance * max(abs(real_part))`. Both tolerances are finite
and nonnegative; the absolute tolerance has input intensity units. Larger residuals
raise `ValueError` mentioning OTF phase/Hermitian symmetry. Thus roundoff can be
discarded explicitly without silently discarding a materially complex result.

`tests/test_detector_correction.py` checks identity and intensity scale, input
immutability, asymmetric known-blur recovery, complex phase versus magnitude-only
correction, odd/even rectangular grids, shifted-impulse centering, zero/small OTF
bins, imaginary residuals, overflow, and input/settings validation. Synthetic
recovery is checked at absolute tolerance 1e-12 and relative tolerance 1e-13.
No matching historical raw/corrected CBED fixture is available, so exact DM
correction equivalence is not claimed. The correction formula is user-confirmed
historical workflow evidence, separate from the archived preparation scripts;
original safeguard, padding and roundoff settings remain unverified.

symmetry-center implements explicitly selected local centrosymmetry origins; symmetry-center has been manually validated by the user. A9 extraction is now implemented; A10 coordinate conversion is now implemented. Validated propagation,
refinement, objectives, optimizers and canonical configuration remain untouched.


## Implemented A4 noise characterization and variance maps

A0–A5 are complete within their stated scope. A4 consumes the **A3-corrected whole
CBED pattern**, never repeats OTF correction, and never changes the supplied image.
Public exports from `qcbed.data_preparation` are:

```python
config = NoiseFitConfig(
    gain_normalization_limit=40000., blemish_sigma=10.,
    intensity_slices=100, sampling_stride=16,
)
model = fit_noise_model(corrected_pattern, valid_mask=None, dark_noise=None,
                        config=config, intensity_units="counts")
variance = variance_map_from_noise_model(corrected_pattern, model)
# Optional immutable diagnostics, including the same model:
result = characterize_noise(corrected_pattern, config=config)
```

Inputs are finite, nonempty real numerical NumPy arrays, copied to float64 in
`(ny, nx)` order. Fitting accepts rectangular frames with each dimension at least
12 and divisible by four. This explicit Python restriction ensures integer crop
coordinates and a nonempty four-pixel interior for N’s half-frame Fourier reference;
it is not a detector/file-format constraint. Unsupported dimensions fail rather
than silently crop/resample. Variance evaluation has no such dimension restriction.
Caller-supplied masks must be Boolean arrays of the same shape.

### Masks, smoothing, and Fourier reference

Without a supplied mask, retain positive input pixels at or below the named gain
normalization limit (source default 40000). With a supplied mask, N bypasses the
gain threshold; A4 does likewise. Both branches apply N’s successive diagonal
mask shifts: shifted = mask(x−4,y−4), then shifted * shifted(x+4,y+4), with explicit
zero extension. This excludes four-pixel borders and retains both original and
shifted mask holes; it is not a generic symmetric erosion. Nonpositive input
pixels are excluded. No corrected intensities are clipped or gain-normalized.

`_smooth` applies the normalized 3×3 kernel (center 1, axial neighbors 0.5,
diagonal neighbors 0.25; sum 4) four times and forms
`4*S(I)-6*S²(I)+4*S³(I)-S⁴(I)`. Separable periodic convolutions implement N’s
Fourier-convolution boundary convention. Tests compare against an independently
constructed FFT kernel. This is noise-estimation smoothing, not a detector PSF.
The unnormalized residual is `I-smooth`.

`_reference_rms` implements the literal N:119–143 crop rather than interpreting
the comment “high frequency half” as a high-pass filter: fftshift the full FFT,
take the central half-width/half-height rectangle, ifftshift that rectangle,
take the real inverse FFT divided by two, and compute RMS over mask[::2,::2].
This crop includes DC under NumPy conventions (a constant I gives reference RMS
2*I). That consequential behavior is regression-tested and exposed, not silently
replaced by high-pass filtering. DM RealFFT storage/crop interpretation remains
unverified. No claim is made that this is an unbiased noise estimator.

The residual mean and population standard deviation are computed over the initial
valid mask. Blemishes satisfy abs(residual−mean) > blemish_sigma*standard_deviation
(default 10); equality is retained. As in N, the normalization denominator still
uses the initial valid mask, before blemish rejection. The residual alone is
multiplied by reference_RMS/residual_standard_deviation. A degenerate residual,
empty mask, or non-finite intermediate fails explicitly. The final noise mask is
valid_mask & blemish_mask & (smooth > 0). Neither corrected intensity nor final
variance is normalized by this internal factor.

### Slicing and fitting

N reuses `binning`: the initial piecewise default 1 is overwritten by the later
final-pattern default 16 before noise sampling. A4 names this independent control
`sampling_stride`, default 16; callers may explicitly select 1 for dense sampling.
This is stride sampling, not pixel averaging. Masked smooth signal and normalized
residual are sampled at that stride, then divided into 100 descending equal-width
intensity intervals by default. The masked zero background participates in the
range. Intervals are (low, high]; each records mean smooth intensity, residual RMS
(not mean-subtracted slice variance), and valid sample count.

`_rebin` accumulates adjacent descending slices until population reaches the
largest original slice population. Mean intensity and squared sigma are
population-weighted. Despite N’s comment, this does not guarantee monotonic sigma.
The source stores a final remainder group at index xcoord, but allocates/copies
only xcoord groups for fitting. A4 preserves this **last-group exclusion** and
returns the excluded remainder as a diagnostic. More than three populated effective
fit groups and a full-rank design are required; otherwise fitting fails clearly.

The physical coefficients in `NoiseModel` remain a (linear), b (square root), and
c (constant), with sigma(I)=a*I+b*sqrt(I)+c. N’s internal `ay,bee,cee` correspond to
`a=ay, b=ay*bee, c=cee`. A4 uses deterministic nonnegative weighted least squares
(SciPy NNLS), minimizing sum(population*(observed_sigma−predicted_sigma)²). This
preserves the source objective and nonnegative physical-coefficient constraints,
without reproducing its random-walk trajectory or approximate coordinate updates.
Column scaling only improves solver conditioning. No RNG is used.

For a>0 the physical parameterization is exactly equivalent to N’s product form.
The solver includes the closed boundary a=0,b>0, approached in N as ay→0 and
bee→infinity; this explicit limiting case avoids an artificial division by zero.
Thus solver outputs need not match the legacy stochastic stopping point. The
reported goodness of fit retains N’s denominator sum(weights)*(1−3/n_groups),
even with fixed dark noise. Fit populations are **not QCBED objective weights**.

`dark_noise=None` fits c; a supplied finite nonnegative dark_noise fixes c exactly.
Positive values follow N. Supplied zero explicitly means fixed zero (the source
zero sentinel branch leaves c initialized to 1); negative sentinels are rejected
in favor of None. The existing immutable NoiseModel is unchanged.

Variance evaluation is exactly `[a*abs(I)+b*sqrt(abs(I))+c]²`. Sigma is standard
deviation, whereas this output is variance, including at negative corrected
intensities. Evaluation rejects negative coefficients for this legacy policy and
non-finite/overflowing results. No mask, weight generation, offset subtraction,
intensity normalization, or downstream extraction occurs.

### Diagnostics, validation, and remaining evidence

`NoiseFitResult` contains model/config, immutable validity/blemish/final masks,
smooth estimate, normalized high-frequency image, original slice points,
effective fit points and excluded remainder (columns intensity/sigma/population),
reference RMS, normalization factor, and goodness of fit. It does not alter
PreparedExperimentalDataset. Units are recorded on NoiseModel; a is dimensionless,
b has square-root intensity units, and c and sigma have intensity units.

`tests/test_noise_preparation.py` covers input/mask validation, rectangular frames,
immutability, exact kernel and fourth-order Fourier comparison, constant residual,
impulse/blemish rejection, source crop/DC scaling, threshold intervals, population
rebinning/remainder exclusion, nonnegative weighted fits, known coefficient recovery,
fixed/free/zero dark noise, absolute-intensity variance, deterministic diagnostics,
and failure/overflow cases.

The repository has 256×256 pattern.dat and sigma.dat arrays. A held-out consistency
test fits the empirical model to half their paired pixels and checks the other
half within 1 stored intensity unit (the sigma text is rounded). This checks a
numerical association only: correction state, gain/mask settings, fitted coefficients,
Fourier intermediate images, and provenance are not supplied. They therefore cannot
validate the full DM noise-characterization pipeline. Exact DM-versus-Python noise
equivalence remains unverified. DM Warp borders and FFT storage conventions,
source crop interpretation, and later smoothed-versus-corrected extraction input
remain open. Historical scan corners are comparison references for current A10 geometry below.

symmetry-center local centrosymmetry measurement is implemented; A9 extraction is now implemented; A10 geometry conversion is now implemented.
Propagation, refinement, optimizer, objective, and canonical configuration code
remain untouched. There are no new dependencies, file loaders, GUI, or adapters.


## Implemented A5 disk geometry and canonical boundary localization

Canonical central-beam localization now uses known-radius partial disk-boundary
geometry rather than integrated disk intensity, because dynamical diffraction can
bias the internal intensity distribution.

```python
from qcbed.data_preparation import DiskGeometry, CentralArcConfig, locate_central_disk
result = locate_central_disk(corrected_pattern, DiskGeometry(38,38),
    search_config=CentralArcConfig(search_radius_px=8, search_step_px=.5))
central_disk_center_px = result.central_disk_center_px
```

`estimate_disk_geometry(profile)` still uses the separate matched disk-profile
image: select above its mean, then above half the selected mean; maximum row/column
populations give x/y diameters divided by two. This is not detector PSF estimation.
For current real validation the matched profile is unavailable: radii 38/38 px
are manually estimated physical semi-axes, not an automatic A5 result.

`central_arc.locate_central_disk()` is the only canonical central localization API.
Its immutable `CentralDiskResult` stores one central center, physical geometry,
search prior/configuration, support evidence and status. The default search prior
is the pixel-frame midpoint `((nx-1)/2,(ny-1)/2)`; optional `initial_center_px` is a
neighborhood prior, not a preset accepted result. No intensity estimator is called.
The local Euclidean search uses radius 8 px and step 0.5 px by default; radius must
not exceed half the smaller semi-axis. Full physical footprints must fit the frame.
An empty support search raises an explicit error without a fallback; an optimum
within half a step of the search limit is flagged `search_boundary` for review.
No automatic search expansion is performed.

`arc_geometry.py` owns shared edge/ridge extraction, polarity-invariant normals,
normal-direction voting, nominal-radius validation, angular coverage and fragment
aggregation. Central localization, A6 and A7 bootstrap use this single implementation.
Central localization enables central edges and uses a local threshold mask extending
by physical semi-axes + search radius + tolerance + 1 px. Filtering remains on the
whole input frame (Gaussian sigma 1 px, Sobel/8, ridge thinning). This smoothing is
for edge evidence only: it does not alter the scientific intensity array.
Each local center hypothesis gathers nominal-radius and normal-compatible pixels
across disconnected fragments. Default tolerance is 2 px, normal angle 12 degrees,
minimum support 12 pixels, angular coverage .15 and arc length 8 px. The score is
shared with A6 and includes support, coverage, fragments, radius/normal errors and
center-vote RMS. Ties prefer smaller prior displacement then x/y. Physical radii
are never shrunk or freely fitted. Complete circles and continuous arcs are not required.
A6 exclusion defaults, arithmetic and A7 acceptance/search rules are unchanged.

The former intensity algorithm and its sum-binning/integrated-ellipse helpers are
isolated in `legacy_disk_localization.py`; its explicit API is
`locate_central_disk_intensity()`. They are absent from canonical package exports
and never called automatically. Source-equation regression tests remain as
compatibility coverage. `compare_central_disk_centers()` and
`examples/central_disk_center_check.py` are explicitly compatibility comparisons;
they do not select or replace the canonical center.

`tests/test_central_arc.py` covers canonical clean, asymmetric, strong-fringe,
disconnected-arc, missing-edge, wrong-radius, ellipse/polarity, noisy/background,
dynamic-frame, deterministic and immutable results. A monkeypatched legacy
initializer raises if called by canonical localization. Diagnostic tests verify
one canonical center and non-overwriting output. The previous intensity algorithm
is still tested only through its compatibility module. No DM golden localization
claim is made, and the reciprocal/zone-axis origin is still a separate quantity.

## Implemented A6: inputparas-free partial-arc reciprocal-lattice seed

Canonical A6 now uses **partial edge / arc detection → known-radius normal-direction
center voting → reciprocal-lattice consistency**. The former intensity-response
candidate detector and its scale/percentile/background configuration have been
removed; there is no fallback detector. The standalone `disk_response_map` helper
remains solely because A7 still uses it for its existing local-evidence algorithm.
A7 algorithms and all simulation/refinement kernels are unchanged.

Dynamical diffraction strongly redistributes disk-interior intensity. A response
maximum inside a disk need not be its center, whereas visible portions of its
boundary can retain useful geometric information. Full circles, connected contours
and central-disk intensity templates are not required. The physical `DiskGeometry`
semi-axes remain an approximate prior, never a detection-kernel size. An arc-voted
center is a preliminary observation, not the final beam position. A6 tests global
lattice consistency; A7 then refines global geometry and predicts the catalogue.
No `inputparas`, HKL or crystallographic reference table is used.

### Public API and coordinates

```python
votes, arrays, diagnostics = arc_center_votes(pattern, geometry, central_center,
                                             config=LatticeDetectionConfig())
candidates = detect_lattice_candidates(pattern, geometry, central_center)
seed = infer_lattice_seed(candidates, central_center, geometry)
# Or the combined path:
seed = detect_reciprocal_lattice_seed(pattern, geometry, central_center)
```

Input is a finite real `(ny,nx)` array. Coordinates are `(x,y)` = (column,row), in
pixels; y increases downward. The original intensity array is never overwritten.
`arc_center_votes` exposes gradient magnitude, selected edge mask, edge coordinates,
unit normals, initial-gradient and ridge masks, proposal peaks, eligible-center mask,
vote threshold and radius offsets. Candidates
are immutable and record ID, center, vote strength, supporting-edge count,
combined angular coverage, radius residual RMS, fragment count, center-vote RMS,
mean/std radius, mean normal residual, geometric confidence, and exact supporting pixels. The existing `response_score`
field is retained as an A7 compatibility slot: in A6 it now means **vote strength**
and has the explicit `vote_strength` property. It is not an intensity score. A7's
own local-response values still have their original meaning.
IDs use integer vote-peak coordinates, remain deterministic for identical input,
and do not promise identity across edited images. Result assignments preserve
abstract integer indices, residuals and inlier/outlier flags.

### Edge estimation and center voting

1. Gaussian smoothing (sigma 1 px by default) is used for derivatives only.
   Sobel x/y derivatives divided by 8 define Gx/Gy; magnitude is hypot(Gx,Gy).
   Nearest extension is confined to derivative computation, and a border of
   `ceil(4*sigma)+1` pixels is excluded from edge selection. No reflected/padded
   scientific intensities or wrapped votes enter the calculation.
2. Gradient strength only selects initial pixels: use the 85th percentile of eligible
   magnitudes plus a scale-relative floating-point guard. Exclude the known central
   ellipse enlarged by 1.2. Thin gradients to ridges using bilinear samples one pixel
   along both normal directions. This prevents weak gradient tails of a small,
   wrong-radius feature from masquerading as a boundary at the physical radius.
3. For every ridge normal n and sampled radius, use the exact axis-aligned ellipse
   displacement `(a² nx,b² ny)/sqrt((a nx)²+(b ny)²)`. Both signs vote. Five offsets
   uniformly cover physical semi-axes ±2 px by default; the same offset applies to
   both axes and must keep them positive. Zero tolerance is supported.
4. Equal edge/radius weights and bilinear deposition create an in-frame vote map,
   smoothed at sigma 1 px. Outside votes are discarded. Proposal maxima pass both
   0.2 effective votes (reduced from 1 after ridge thinning) and 30% of the largest
   eligible vote. This is a permissive proposal gate, not disk-boundary acceptance.
5. Treat each vote peak as a **known-radius center hypothesis**. Gather all ridge
   pixels satisfying `abs(distance(pixel, center) - R) <= radius_tolerance_px`
   for a circular disk, and normal-angle residual ≤12° modulo 180°. R is the nominal
   physical radius; it is not fitted. For an ellipse, radial error is distance to
   the nominal ellipse along the pixel's radial ray, and the expected normal is
   proportional to `(dx/a²,dy/b²)`. No additional support tolerance expands the prior.
6. Validate the combined support across **all disconnected fragments** around that
   same center. Require ≥12 supporting pixels, ≥8 px total projected arc length,
   and ≥0.15 circumference (54°) combined angular coverage. There is no minimum
   count, length, or angular span for an individual fragment and no closed-contour
   requirement. Coverage is the union of projected pixel angular intervals, with
   periodic wrap handled explicitly: missing gaps and overlapping support are not
   counted. Arc length integrates nominal ellipse tangent speed over this union.
7. Count spatially independent fragments with the existing 8-neighbor/25°
   sign-invariant normal-continuity graph. This count describes spatial fragments,
   not statistical independence; connectivity is a diagnostic, never an acceptance
   gate. Rank accepted centers using total count N, combined coverage C, fragment
   count F, radius RMS r, mean normal angle theta and center-vote RMS s:
   `N*C*(1+0.1*log(F)) / (1+s/center_vote_rms_scale_px+r+theta/maximum_normal_residual_deg)`.
   This geometric score is not a probability. Center-vote RMS is measured about the
   proposed center and penalizes diffuse support; it is no longer a hard rejection
   threshold. Suppress accepted centers within 6 px by this score, cap at 64, and
   sort by center/ID before the unchanged lattice fit.

The architecture is **known-radius center hypothesis → gather compatible edge
support → validate combined arc evidence**. It is not fragment validation followed
by center inference. The physical radius is the primary constraint; disconnected
weak/short boundary sections may jointly provide sufficient support. Strong internal
fringes are rejected when they cannot satisfy the physical radius and normal geometry
around the same center. A coherent internal feature with the same geometry remains
intrinsically ambiguous; accepted centers are still preliminary observations.

No local curvature derivative is required. Radial/normal consistency supplies the
geometric test, and center-vote spread contributes to confidence rather than forcing
one long smooth fragment. The tests include eight individually insufficient fragments
which succeed together, rejection of inconsistent radii/normals/centers, and an image
with disconnected visible boundary sections. Existing straight/parallel/wrong-radius
internal-feature rejection, polarity, partial-arc and lattice tests remain.

For elliptical priors, `mean_radius_px`/`radius_std_px` use circularized nominal radius
`sqrt(a*b)*sqrt((dx/a)^2+(dy/b)^2)`. `radius_residual_rms_px` uses radial distance to the
nominal ellipse. For circles these are ordinary Euclidean radius statistics. The
compatibility field `radius_offset_px` is now always zero: support is gathered at the
nominal radius with one explicit tolerance, not at a second fitted radius. Discrete
radius samples remain only in hypothesis voting. `support_tolerance_px` was removed;
`maximum_center_vote_rms_px` was replaced by scoring scale `center_vote_rms_scale_px`.

Controls are validated finite scalars; count fields require strict integers.

| Detection setting | Default |
|---|---:|
| edge_smoothing_sigma_px | 1 |
| edge_percentile | 85 |
| radius_tolerance_px | 2 |
| radius_samples (odd, ≥3) | 5 |
| vote_smoothing_sigma_px | 1 |
| minimum_vote_strength | 0.2 |
| relative_vote_threshold | 0.3 |
| candidate_merge_radius_px | 6 |
| minimum_support_pixels | 12 |
| minimum_arc_length_px | 8 |
| center_vote_rms_scale_px | 2.5 |
| maximum_normal_residual_deg | 12 |
| maximum_neighbor_normal_change_deg | 25 |
| minimum_angular_coverage | 0.15 |
| central_exclusion_scale | 1.2 |
| max_candidates | 64 |

### Minimal seed from two validated centers

`infer_lattice_seed()` (and therefore `detect_reciprocal_lattice_seed()`) now has two
explicit paths. Exactly two geometrically validated diffracted centers use direct
initial displacements from the known central center. A sufficiently large candidate
set still uses the unchanged displacement-voting path below. Three to five candidates
retain the existing insufficient-support error; we do not silently select a pair
from a larger inconsistent set or fall back after a failed multi-candidate fit.

Minimal mode requires both candidates' retained physical-radius support evidence,
minimum pixel count/combined coverage/arc length, positive geometric confidence,
fragment count, and normal/radius residual criteria. Raw coordinate-only candidates
are rejected. Recorded supporting pixels are rechecked against the physical radius.
Arc validation thresholds and candidate detection are not weakened.

Safeguards: initial displacement lengths and reduced basis lengths must be at least
`minimum_lattice_vector_scale * min(radius_x,radius_y)` (default 1.5×radius, 57 px for
this dataset), `abs(sin(angle)) >= min_basis_sine` (default 0.2), and initial/reduced
basis condition numbers must not exceed 1e6. The visual 85 px spacing is not an input.
The existing Gauss reduction/sign convention canonicalizes the direct basis. Both
observations must reconstruct consistently under its integer labels. Pairwise
family voting is not called, and no additional observations are invented.

Seed provenance records `seed_method="two_validated_arc_centers"` or
`"multi_candidate_displacement_voting"`, the reference center, original candidate
records (IDs, centers, all arc metrics/support), initial vectors, canonical vectors,
canonical included angle in degrees and condition number. `primitive_basis_certain`
is explicitly false. Minimal support=1 and residual≈0 describe reconstruction of
two defining observations, not independent validation or primitive-basis certainty.

**A7 handoff:** the free-origin global fit still requires four independent candidate
observations by default. Minimal seeds now enter the bootstrap stage below first;
the known central reference and the two defining disks are never counted as new
independent confirmation. A failed bootstrap retains predictions and crops without
claiming a refined lattice.

### Multi-candidate displacement voting and basis testing

All pairwise displacements among candidates plus the fixed reference center are
computed. Vectors shorter than minimum_lattice_vector_scale*min(radius_x,radius_y) are excluded.
Deterministically sorted vectors are greedily clustered by Euclidean pixel distance,
allowing sign reversal. Each family's center is the running mean of aligned vectors;
its vote count is the number of displacement pairs, not independent reflections.
Families passing minimum votes are ranked by votes, length, and components and
bounded by max_families. Every non-collinear family pair is a basis hypothesis.

For each hypothesis, solve G*n = candidate−reference, round n using NumPy's nearest
integer convention, and calculate Euclidean reconstruction residual. A candidate
is supported only within the residual tolerance and with nonzero n. At most one
candidate supports each integer site; the lowest-residual candidate wins. Other
candidates remain diagnostic outliers. Three iterations of equal-weight inlier
least squares refine **only the basis**, keeping the A5 reference fixed. Candidate
vote strength does not weight lattice support, so weak detected disks can contribute.
No physical HKL meaning is assigned to integer coordinates.

Accepted hypotheses meet minimum inlier count/fraction, length and angle constraints.
Ranking is lexicographic: most unique supported sites, shortest sum of reduced basis
lengths, highest pair-family votes, lowest RMS, then deterministic vector components.
This combines repeated geometry and global consistency; it does not blindly choose
the nearest two spots or the two shortest raw displacements. Singular/ill-conditioned
hypotheses are rejected. Insufficient, collinear or inconsistent data raises an
informative error; no basis or legacy fallback is fabricated.

Basis vectors are 2-D Gauss-reduced by integer subtraction/swapping, first vector
shortest, first vector's x positive (or y positive for exactly zero x), and positive
determinant. Remaining exact ties follow sorted family enumeration/vector ranking.
This deterministic representative is invariant to input candidate ordering, but
near symmetry, changing image noise can select another equivalent basis. The API
does not promise invariant orientation under physical rotation or perturbed data.

Unchanged lattice controls: `minimum_lattice_vector_scale=1.5` (renamed from the
old shared candidate-separation field, preserving its inference value),
`displacement_tolerance_px=1.5`, `lattice_residual_tolerance_px=2`,
`min_family_votes=3`, `max_families=32`, `min_inliers=6`,
`min_inlier_fraction=0.6`, `min_basis_sine=0.2`.

### Synthetic validation and limits

Tests explicitly reject stronger straight/parallel internal fringes and wrong-radius
curved features while preserving 20–40% true outer arcs. Tiny arcs are rejected;
polarity reversal preserves centers and supporting pixels. Diagnostics tests compare
CSV metrics and exact pixel membership against immutable candidate records.
Tests also cover full arcs, weak boundary sectors, unequal elliptical
semi-axes, bright displaced interior features, moderate noise/background gradients,
missing lattice disks, unrelated straight edges and a false non-lattice disk,
radius error, positive intensity scaling, deterministic IDs/basis, and 38 px radii
with 85 px separation. A7's direct response-sum regression remains. Existing
inference outlier/collinearity/ordering tests remain. Obsolete response-detector,
background-annulus and detection-scale tests and their example study were removed.

The center-vote map can still contain false peaks from internal edges. Sparse arcs,
very small disks, inaccurate radius priors, heavy overlap, excessive false candidates
and candidate/family budgets remain limitations. Support/coverage are diagnostics,
not calibrated probabilities. A failed fit must not be presented as final geometry.
The physical radius is never adjusted to force a lattice. No symmetry-center work is included.

## A7 bootstrap: visible boundary support and independent observations

```text
A6: known-radius arc geometry → minimal reciprocal seed
A7: predicted lattice sites → local arc confirmation
    → sufficient independent observations → existing global lattice refinement
    → complete arc-evidence beam catalogue
```

`bootstrap_lattice(pattern, seed, geometry, config=BootstrapConfig())` preserves the
seed and uses only predicted-site searches. A5/A6 scientific algorithms and shared
`arc_geometry.py` criteria are unchanged. No whole-image redetection is performed.
The central reference is never counted as an independent observation.

`boundary_sites.enumerate_boundary_sites()` enumerates all ellipses with nonzero
in-frame circumference, including higher-order sites and centers outside the frame.
It solves detector bounds expanded by physical semi-axes into integer coordinates;
it neither truncates the enumeration nor imposes a fixed beam count. Ordering is
central first, squared displacement, then integer pair, with the existing enumeration
budget. `enumerate_seed_sites()` and `enumerate_lattice_sites()` remain available
with their original fully-contained semantics for the older response catalogue.

Visible boundary fraction is arc-length-weighted midpoint quadrature over 4096
uniform ellipse phases. This is a numerical approximation, not image interpolation.
The pixel-center frame is [0,nx-1]×[0,ny-1]. Records distinguish physical in-frame
fraction from usable fraction excluding the shared derivative margin
`ceil(4*sigma)+1` (5 px by default). No invented boundary pixels, image padding,
wrapping or reflected data enter confirmation. Existing shared filtering is unchanged;
only its stencil-safe interior edge pixels can support an observation.

`BootstrapConfig.minimum_boundary_fraction` defaults to **0.25**, is explicit and
must lie in (0,1]. A local hypothesis must meet it inside the stencil-safe frame.
The nominal physical ellipse is never shrunk. Local searches remain **4 px radius,
1 px step**, strictly less than a quarter of seed spacing. They use the same
radius tolerance (2 px), normal consistency, support count, combined angular coverage,
fragment aggregation and center-vote ranking as before. Disconnected fragments may
jointly confirm an ellipse; local intensity maxima do not confirm beams.

A site is `edge_clipped` when no hypothesis has sufficient usable circumference;
otherwise failed arc validation is `missing`. Thus a predicted fraction below .25
can still be searched if an allowed nearby hypothesis reaches .25. Both fractions
in the table refer to the prediction, not a substituted observed center.
A strongest solution within half a step of the search limit is explicitly
`boundary_limited=True`, `weak_evidence`; it is never silently accepted. Shared
support pixels cannot be reused by another independent confirmation. Duplicate
seed IDs, positions or integer sites are rejected. Defining observations have
source `a6_defining`; accepted new observations have source
`a7_bootstrap_confirmed`. Prediction, local hypothesis and accepted observation
remain separate, including their shifts. The central reference is not fabricated
as an independent observed disk.

The existing `refine_reciprocal_lattice()` runs only with enough independent
observations (default **four**). Rank, trimming, reassignment, convergence and
canonicalization checks are unmodified. Insufficient support retains a partial
catalogue, with no invented origin, vectors or fit residuals.

`generate_arc_catalogue()` uses the existing fit's assignments when available and
regenerates all footprint-intersecting predictions. It preserves candidate IDs,
initial/final integer labels, residuals, rejection status and exact recorded arc
observations. Unobserved newly entering sites remain unobserved. It does not call
the old intensity-response confirmation API. The older response-based catalogue
API remains unchanged for compatibility. Without a fit the new catalogue is
explicitly labeled **partial seed catalogue**, not final refined geometry.

Artifacts include `bootstrap_overview.png`, `bootstrap_sites.csv`, individual
`bootstrap_beams/` crops, plus `beam_catalogue_overview.png`, `beam_catalogue.csv`
and `beams/` crops. Every enumerated site is retained, even missing/clipped ones.
Rows include physical footprint bounds, visible/usable fractions, bounded-search
status, source, candidate ID, confidence, support metrics, prediction/observation,
local-hypothesis shift, integer mapping and available fit residual. Images show
physical ellipses and supporting arcs. `refinement_comparison()` reports actual
seed-to-fit changes and before/after RMS; these are null without a successful fit.

Focused tests cover higher-order sites, clipped fourth-observation fitting,
insufficient circumference, no double counting/fabrication, boundary-limited
solutions, ellipse fractions, enumeration limits, full post-fit arc catalogue,
deterministic IDs/order and recorded-coordinate preservation. All prior A5/A6/A7
regressions are retained; old full-footprint bootstrap count assertions now include
the newly requested edge-intersecting sites and explicit provenance names.

## Implemented A7: global pixel lattice and predicted beam catalogue

A7 consumes an A6 `ReciprocalLatticeSeed`; neither fitting nor catalogue generation
reads inputparas. The fitter and older response catalogue are implemented and synthetically tested; real-data A7 completion is still blocked by independent support. The older API flow is:

```text
A6 seed → fit origin + g1 + g2 over candidate observations
        → geometric rejection/reassignment → refined global lattice
        → all footprint-contained integer sites → constrained local response evidence
        → predicted BeamObservation catalogue with observed/weak/missing/rejected status
```

This is a **global lattice origin**, not a physical or centrosymmetry-refined zone
axis. The original A5 central-disk evidence remains in `lattice.initial_seed` and
is not constrained to equal the fitted origin. No `ReciprocalFrame` is produced.

### Public APIs and model boundary

```python
config = LatticeRefinementConfig()
lattice = refine_reciprocal_lattice(seed, config=config)
catalogue = generate_beam_catalogue(corrected_pattern, lattice, disk_geometry,
                                   valid_mask=None)
# Combined convenience API:
catalogue = refine_lattice_and_generate_catalogue(
    corrected_pattern, seed, disk_geometry, config=config, valid_mask=None)
# Pure geometric enumeration, independent of image intensity:
sites = enumerate_lattice_sites(lattice, corrected_pattern.shape, disk_geometry)
```

New frozen models are `LatticeRefinementConfig`, `RefinedCandidateAssignment`,
`RefinedReciprocalLattice`, `PredictedBeam` and `LatticeBeamCatalogue`, with enum
`BeamStatus`. The existing `BeamObservation` is unchanged. Each `PredictedBeam`
wraps a BeamObservation whose center **always equals predicted_center_px** and whose
HK/HKL `indices` is None. Its separate `lattice_indices` holds abstract `(n1,n2)`.
`observed_center_px` is optional independent evidence, never substituted for the
prediction. Stable IDs are `lattice_n1_n2` within the canonical integer frame;
these are not persistent crystallographic identities across arbitrary reindexing.
The catalogue exposes `.observations` and a fresh `.status_counts` dictionary.
Weak/missing entries are retained and enabled; this stage makes no extraction or
objective-weight decision. Arrays use `(ny,nx)` and positions/vectors use `(x,y)`.

### Deterministic robust affine fit

For fixed candidate integer coordinates, float64 least squares solves the design
matrix `[1,n1,n2]` against both measured coordinates simultaneously. All distinct
A6 inlier sites initialize the fit. Duplicate sites resolve by smallest residual,
then greatest response, then candidate ID, retaining all records in diagnostics.
At least four inliers are required by default (configurable minimum >=3), with
full affine design rank. Singular/ill-conditioned lattices fail clearly.

After each fit, fractional coordinates are solved and rounded with NumPy's nearest
integer rule. Reassignment occurs only if its reconstruction residual improves by
more than 0.1 pixel by default. Candidates within the residual threshold (default
1.5 pixels) become inliers, including previously rejected candidates if now
consistent; duplicate-site filtering is reapplied. Fit and classification repeat
until integer labels and membership are unchanged. A repeated nonconverged discrete
state or exhaustion of 20 iterations raises an error rather than returning an
oscillating/incomplete result. This is transparent threshold trimming, not a
high-breakdown estimator; a badly contaminated seed can still bias or defeat it.

A6's deterministic Gauss-reduction/sign/positive-determinant convention is reused.
Integer labels are transformed into the resulting basis; predictions do not change.
The `reassigned` flag records geometric reassignment, not mere basis canonicalization.
Every assignment retains ID/center, initial and final labels, residual and inlier
state. Initial RMS uses initial unique A6 inliers; final RMS uses final inliers, so
the populations can differ. The initial seed, final geometry, iteration count and
configuration are retained for comparison. The initial center is not pinned in the fit.

### Complete footprint-contained enumeration

The image rectangle is reduced by the x/y disk semi-axes. Its corners are transformed
into lattice coordinates to bound integer enumeration. Every site is then tested
against that rectangle. A comparison tolerance `64*float64.eps*max(ny,nx)` handles
roundoff at exact boundaries; positions are never snapped, clipped or normalized.
The central `(0,0)` footprint must fit, otherwise enumeration fails. Exactly one
central site is returned. No visible peak is needed to predict any site.

Order is central first, then increasing squared reciprocal displacement from the
origin, then integer-coordinate tuple. Counts are dynamic and duplicate sites are
impossible by construction. The configurable enumeration budget (100000 bounding-box
sites by default) raises on excessive allocation/work; it never truncates a catalogue.
The valid mask does not remove geometric predictions: unsupported observations remain
as rejected entries for review.

### Local evidence and statuses

A7 local detection retains the tapered elliptical response minus its background annulus;
this helper is no longer used by canonical A6.
The A3 corrected pattern is not altered or renormalized. Searches are confined to
radius `min(local_radius_scale*min(rx,ry), 0.2*shortest_reduced_basis_length)`, default
scale 0.5. A circular grid includes the prediction and uses spacing no greater than
`local_step_px` (default 0.25). Bilinear response sampling, not intensity interpolation
for downstream science, provides local scores. Ties prefer smallest displacement,
then x/y offset. There is no new global disk search.

Defaults `observed_response_threshold=1` and `weak_response_threshold=0.1` are
**absolute response thresholds in input intensity units**, not universal detector
calibrations or probabilities. They must be chosen for the experiment's units.
Below the weak threshold means missing; between thresholds means weak. Scores above
the observed threshold yield an observed center only if the best point is strictly
inside the local search boundary. A boundary maximum is rejected as insufficiently
constrained. Weak/missing/rejected beams retain predictions and no observed center.

The complete response annulus plus search margin must fit the frame; otherwise the
site is retained as rejected. A Boolean same-shape valid mask optionally requires
at least `minimum_mask_coverage=0.95` over that bounding support rectangle. Coverage
includes background support, not just the disk. This conservative criterion avoids
claiming padded edge responses as observations; a nonzero tolerated fraction can
still influence a response and must be selected appropriately. Mask and image are
never mutated. Score, coverage, search radius, classification, predicted position
and optional observed position are separately inspectable.

### A7 tests and limits

`tests/test_lattice_refinement.py` tests exact affine recovery, noisy centers,
unique geometric outliers, duplicate observations, corrected provisional assignments,
missing-observation enumeration, brute-force oblique enumeration, rectangular frame
bounds, weak/missing status retention, displaced strong intensity maxima, reliable
small local offsets, mask rejection, central identity, unchanged physical index
fields, candidate-order determinism, immutable results, convergence/rank/budget
failures, and a complete synthetic A6-detection-to-A7-catalogue run.

This new geometric inference is not an executed legacy DM equivalence claim. True
primitive indexing may remain ambiguous from missing reflections. Biased intensity
centers and contaminated seeds can bias the global origin; symmetry/noise can select
equivalent canonical bases. Repeated identical inputs are deterministic, but IDs and
ordering are not guaranteed invariant under changed data or a different basis label.
A later physical mapping and centrosymmetry analysis must resolve these meanings.
symmetry-center measurement is documented below, with real-data validation confirmed by the user. No physical HKL,
reflection ROI extraction, signal/variance/weight extraction, scan corners, GUI,
canonical adapter or simulation/refinement changes. Legacy compatibility is unchanged.


## A7 real-data validation extension (no numerical redesign)

The user identifies `reference_data/pattern.dat` as already legacy PSF-corrected
and noise-prepared. The diagnostic command uses its intensities directly and does
**not run A2, A3 or A4**. The command uses geometric A5 localization and the unchanged A6/A7 rules.

The file inspected on 2026-09-28 is **256×256**, with 65,536 whitespace-delimited
values, not the requested 2048×2048 (4,194,304 values). `numpy.loadtxt(dtype=float64)`
reads physical text rows as y and columns as x, without flattening, transposing,
reshaping or resampling. No dedicated pattern.dat image loader exists in the current
preparation package. Existing logical-record experimental-data readers serve a
different data format and are not reused for this rectangular image text.

The initial attempt had no matched A5 disk-profile fixture or independently
validated radius pair in these pixel units. The subsequent manual-radius run below
uses the user-supplied geometry instead. The aperture calibration file `apertureimage.dm3` is
not established as N's disk-profile `Psf`; neither it nor pattern.dat is silently
substituted. That initial attempt was **blocked before A5**.
A5 centers, A6/A7 bases, beam counts and crop figures are unavailable, not zero.
No thresholds were tuned and no real geometry accuracy conclusion is made.

Reproduce the input check from the repository root (choose a fresh output directory):

```bash
python -m qcbed.data_preparation.diagnostics \
  --pattern reference_data/pattern.dat --output tmp/a7_pattern_validation
```

The command writes `validation_summary.json` and `validation_summary.md`, recording
shape, value count, float64 interpretation, SHA-256, configuration and blockers,
and exits nonzero when blocked/failed. Generated artifacts stay under ignored `tmp/`.
Existing output directories are refused. Reference inputs are never modified.

Once the matching image and profile are supplied, pass `--pattern PATH --profile
PROFILE_TEXT_PATH --geometry-provenance "description of the matched reference"`.
Alternatively use `--radii RX RY --geometry-provenance "independent validation and
pixel units"`. Profiles use the same rectangular numeric-text loader. No DM file
loader or fabricated replacement is introduced. Shape is now dynamic by default, with optional `--expected-shape NY NX` as an
explicit assertion. The confirmed current dataset is 256×256. Default A5 geometric search, A6 configuration, and A7 configuration are recorded
and used unchanged. Stage failures retain completed stage information. Non-default
profile/radii, expected shape and crop scale are recorded; there is no auto-tuning.

On a completed geometry run, outputs are:

- `beam_catalogue.csv`: ID, n1/n2, predicted x/y, observed x/y, shift magnitude,
  status, response score and matching inlier-candidate lattice residual. Missing
  observed coordinates/residuals are empty CSV cells, never substituted predictions.
  The lattice residual is absent if no fitted inlier corresponds to that site; it
  is not a renamed local shift or an invented zero.
- `beam_catalogue_overview.png`: original prepared intensity, A5 center, A6 candidates,
  A7 predicted and reliable observed positions, with integer labels for <=150 sites.
- `beams/beam_NNN_nN1_N2.png`: deterministic names (`m` denotes negative indices),
  raw local intensity and separate predicted/observed markers joined by a line,
  status, response and **expected footprint, not detected edge** ellipse.
- JSON/Markdown summaries: initial/final geometry, assignments/outliers, RMS,
  counts, shift statistics, largest shifts, low-response observed sites and missing
  sites. These are review diagnostics, not automatic scientific conclusions.

Crop half-width/height are `crop_scale*radius_x/y`, default 1.25. Bounds round
outward (floor lower, ceil upper plus one for inclusive pixel centers), are clipped
without padding/reflection/wrapping, and are recorded with a clipped flag. An empty
intersection is skipped with a reason, retaining its numerical beam record. Ellipses
have width/height 2*radius_x/y and use full-frame x/y coordinates; `imshow` extents
align pixel centers correctly. Plotting uses raw intensities; it never modifies fit
results. Matplotlib is imported only in diagnostic plotting, not numerical APIs.

`tests/test_a7_diagnostics.py` checks text row/column interpretation, footprint-derived
crop bounds, overlay coordinates, deterministic filenames, absent observations,
CSV precision/semantics, clipping, PNG generation, image immutability and blocked
input reporting. The expensive real geometry run is a reproducible command, not a
unit test with unreviewed hard-coded beam counts. A7 scientific behavior remains
unchanged. symmetry-center local symmetry APIs are documented below; A9 is now implemented; see the A9 extraction section below. A10 coordinate conversion is now implemented; see the A10 section below.


### Current real-data validation: canonical geometric center, manual 38 px radii

The matched A5 profile remains unavailable; physical `DiskGeometry(38,38)` is a
manual visual estimate. The actual pattern is dynamically loaded as 256×256.
Canonical geometric localization from image midpoint (127.5,127.5) gives
**(127,126)**, versus historical intensity (125,127): delta **(+2,-1)** px,
magnitude **2.2360679775** px. Boundary evidence: 235 pixels, 10 fragments,
92.6703% angular coverage, radius RMS .551999 px, mean normal residual 5.759964°,
center-vote RMS 4.465548 px. Status `estimated`, not a search-limit optimum.
This is a user-selected scientific method change, not a claim of Fortran equivalence.

The unrestricted diagnostic rerun is recorded in `tmp/a7_geometric_central_r38/`.
The new exclusion mask changes the adaptive threshold and automatic proposals:
(118,208), (136,44), (208,135). Existing A6 rules reject this three-candidate set
(insufficient multi-candidate support; minimal seeding requires exactly two).
No candidates are silently dropped and no A6 rule is relaxed to force a seed.

The separately requested review rerun explicitly revalidates the previously accepted
(137,45) and (208,135) centers using the NEW center's edge mask and shared boundary
criteria. Both pass. It is reproducible with:

```bash
python examples/central_geometric_validation.py \
  --output tmp/a7_geometric_reviewed_r38 --radii 38 38 \
  --geometry-provenance "Manually estimated physical 38 px radius; matched A5 profile unavailable" \
  --defining-center 137 45 --defining-center 208 135
```

Direct source displacements are **(10,-81)** and **(81,9)**, lengths
81.61494961 and 81.49846624 px; angle 89.30225098°. Existing deterministic basis
canonicalization swaps/signs these to **g1=(81,9), g2=(-10,81)**, angle
90.69774902°, condition 1.01233727. This indexing change is expected, not a new
indexing algorithm. Reconstruction RMS is zero by construction, not fit validation.

### Expanded A7 real-data search (current result)

The starting geometry is unchanged: center (127,126), g1=(81,9), g2=(-10,81),
manual physical semi-axes 38/38 px. No old intensity center or inputparas is used.
The current reproducible command is:

```bash
python examples/central_geometric_validation.py \
  --output tmp/a7_boundary_catalogue_r38_final --radii 38 38 \
  --geometry-provenance "Manual R=38 px; matched profile unavailable" \
  --defining-center 137 45 --defining-center 208 135
```

All 17 footprint-intersecting sites are recorded below. Sites with no sufficiently
usable local hypothesis are classified without applying arc validation to inadequate
circumference. The rest are checked with the unchanged local radius/geometry criteria.

| Indices | Prediction (x,y) | Status | Usable boundary fraction | Reason |
|---|---|---|---:|---|
| (0, 0) | (127.0, 126.0) | reference | 1.0000 | Known A5 center; not independent validation |
| (-1, 0) | (46.0, 117.0) | missing | 1.0000 | No validated local arc support |
| (1, 0) | (208.0, 135.0) | defining | 1.0000 | Defines seed; not new validation |
| (0, -1) | (137.0, 45.0) | defining | 1.0000 | Defines seed; not new validation |
| (0, 1) | (117.0, 207.0) | confirmed | 1.0000 |  |
| (-1, -1) | (56.0, 36.0) | missing | 0.8037 | No validated local arc support |
| (1, 1) | (198.0, 216.0) | missing | 0.8525 | No validated local arc support |
| (-1, 1) | (36.0, 198.0) | missing | 0.8037 | No validated local arc support |
| (1, -1) | (218.0, 54.0) | missing | 0.8188 | No validated local arc support |
| (-2, 0) | (-35.0, 108.0) | edge_clipped | 0.0000 | Insufficient usable circumference |
| (2, 0) | (289.0, 144.0) | edge_clipped | 0.0000 | Insufficient usable circumference |
| (0, -2) | (147.0, -36.0) | edge_clipped | 0.0000 | Insufficient usable circumference |
| (0, 2) | (107.0, 288.0) | edge_clipped | 0.0000 | Insufficient usable circumference |
| (-2, -1) | (-25.0, 27.0) | edge_clipped | 0.2036 | Insufficient usable circumference |
| (2, 1) | (279.0, 225.0) | missing | 0.2236 | No validated local arc support |
| (-1, 2) | (26.0, 279.0) | edge_clipped | 0.2051 | Insufficient usable circumference |
| (1, -2) | (228.0, -27.0) | edge_clipped | 0.1812 | Insufficient usable circumference |

The lower observation is still (117,208), 1 px from prediction (117,207).
There is **one additional independent confirmation**, plus the two defining disks:
**three observations, below the unchanged four-observation minimum**. Global
refinement was not invoked. Refined origin/g1/g2, RMS before/after and seed-to-fit
changes are **unavailable**, not zero. No final refined geometry is claimed.

The partial catalogue contains 3 observed, 6 missing, 7 edge-clipped, 0 weak,
0 rejected, and 1 central reference (17 records). Seven sites are limited by field
of view/stencil-safe circumference. Six searched sites fail boundary evidence;
this does not prove zero diffraction intensity. The diagonal (198,216) is now
searched despite partial clipping, but still fails arc validation. Higher-order
(279,225) has sufficient circumference at some local hypotheses, but no validated
arc. No thresholds or local windows were loosened to obtain a fourth observation.

Outputs are in `tmp/a7_boundary_catalogue_r38_final/`:
`validation_summary.md/json`, `bootstrap_overview.png`, `bootstrap_sites.csv`,
`bootstrap_beams/`, `beam_catalogue_overview.png`, `beam_catalogue.csv`, and `beams/`.
Both crop sets contain all 17 sites. JSON keeps exact settings, support pixels,
input hash, source provenance and absent refinement quantities. A5/A6 science is
unchanged; inputparas remains absent; symmetry-center is now documented below; its real-data workflow has been manually validated by the user. **A7 real-data completion
is not claimed.**


## symmetry-center: explicitly selected local centrosymmetry and zone-axis origin

The symmetry-center symmetry reference beam is selected according to the relevant
zone-axis/diffraction condition. The legacy DM convention of storing that beam
second in `inputparas` is not a physical requirement and is removed from the
canonical Python workflow.

```text
selected beam + selected beam-local centrosymmetric feature
→ local 180° symmetry search → legacy origin formula → zone_axis_origin_px
```

`qcbed.data_preparation.origin` provides immutable `SymmetryReference`,
`SymmetrySearchConfig`, `SymmetryMeasurement` and `OriginConsistency` models.
`select_symmetry_reference(arc_catalogue, beam_id=..., feature_offset_px=(dx,dy))`
resolves a stable ID; `lattice_indices=(n1,n2)` is an alternative. Supplying both
requires them to identify the same beam. The central/direct `(0,0)` reference and observed diffracted beams are
accepted by this catalogue adapter. A directly constructed `SymmetryReference`
is an explicit externally validated beam/feature declaration, not automatic beam
quality assessment. No reference or feature is inferred from catalogue ordering.
An observed disk is not automatically suitable for centrosymmetry.

The reference stores its geometric center and mandatory beam-local feature offset.
The feature seed equals beam center + offset; its seed must lie inside the physical
ellipse. Re-resolving the same selection against updated beam geometry keeps the
feature attached through this offset. The result separately retains the beam
center, feature seed, best symmetry trial center, crop pixel midpoint and derived
zone-axis origin. The selected A7 geometric center is the lattice prediction, not
an automatically substituted local arc-observation displacement. Abstract lattice
indices are not physical HKL.

### Source-defined search and explicit sampling limits

Inspected source: `simplified_finaledition_withipp.s`, lines 594–623 (N).
The source uses `in` as its input image; symmetry-center accepts the caller's chosen working
signal without performing smoothing, correction, normalization or interpolation.
With disk width/height = twice the physical semi-axes:

- default search half-size = (width/8, height/8);
- default crop half-size = (width/4, height/4);
- y-major outer loop, x inner loop; one-pixel increments;
- lower-inclusive, upper-exclusive search intervals;
- crop `[y-half_y:y+half_y, x-half_x:x+half_x]`, half-open;
- strict smaller mismatch updates the result, preserving the first minimum on ties.

An explicit `SymmetrySearchConfig(search_half_size_px=(sx,sy),
crop_half_size_px=(hx,hy))` overrides sizes. Every crop must fit the original frame;
no clipping, padding, wrapping or interpolation is allowed. Search is capped at
100000 trials. Fractional DM subarea conversion is not established by the available
source: **fractional crop bounds raise an error rather than being rounded**.
Half-integer trial centers with half-integer crop sizes are supported when bounds
are integral. With R=38 the default search half-size is 9.5 and crop half-size is
19; an integer feature seed therefore needs explicitly supplied compatible sizes.
This limitation is not a claim about what DM does with fractional bounds.

`rotate180(Cut)` is exactly `Cut[::-1, ::-1]`: columns and rows both reverse.
The mismatch is the unnormalized float64 sum `sum(abs(Cut-rotate180(Cut)))`.
No radial or whole-disk intensity score replaces this local comparison. For an
even crop, the actual pixel midpoint is (trial_x-.5, trial_y-.5). The legacy formula
uses the **trial coordinate**, not an undocumented half-pixel correction:

```python
zone_axis_origin_px = symmetry_center_px + 0.5 * (
    central_disk_center_px - reference_beam_center_px)
```

For odd crops, the sampled midpoint is also recorded explicitly. Tests implement
the integer-bound legacy loops independently and compare every mismatch, as well
as even/odd rectangular reversals and the origin equation. No DM executable/golden
output is available: source-equation tests do not establish bitwise DM equivalence
or its float32/fractional-coordinate behavior.

The result retains every trial score, best crop bounds, tie count, crop contrast
and whether the best is on the sampled search boundary. Flat crops can give zero
mismatch, so zero mismatch alone is not proof of a useful diffraction-condition
feature. No arbitrary acceptance threshold is invented; diagnostics require review.

### Multiple references and real validation interface

Call `measure_symmetry_origin()` separately for each explicitly supplied beam and
feature. `combine_origin_measurements(measurements, tolerance_px=..., combine=False)`
reports maximum Euclidean pairwise origin separation and consistency. Combination
requires `combine=True`; only consistent measurements produce an equal-weight
arithmetic mean. Inconsistent origins are retained individually and never averaged.
Consistency itself does not establish feature quality or statistical independence.

`examples/a8_origin_validation.py` accepts an existing A7 arc `beam_catalogue.csv`,
`--beam-id` and/or `--indices`, required `--feature-offset`, `--central-center`,
`--radii`, optional `--search-half-size`/`--crop-half-size`, and a new output directory.
For the real frame, supply the known central center 127 126 and radii 38 38;
the actual reference selector, feature offset and compatible ROI/search sizes
remain **user inputs still required**. No example chooses those scientific values.
Use `python examples/a8_origin_validation.py --help` for the complete command.

`write_symmetry_diagnostics()` writes `symmetry_diagnostic.png` and
`measurement.json` into a fresh directory. Panels show the detector and selected
full disk, beam geometry, feature/search seed region, best trial, crop rectangle,
original crop, exact rotated crop, difference and derived origin. The CLI also
records input/catalogue paths and hashes, shape and preprocessing provenance.

`tests/test_origin.py` covers non-positional selection, suitable versus asymmetric
local regions, off-center features in an asymmetric disk, geometry-relative offsets,
legacy ranges and crop sizes, all trial scores, exact origin relation, even/odd
rotation, ambiguous selections, invalid bounds, deterministic ties, immutable
results, multiple-origin disagreement and diagnostic overwrite protection.

**Real pattern.dat symmetry-center validation has been confirmed by the user.** A9 does not
redetermine its symmetry beam, feature region or origin.
A5/A6/A7 algorithms and canonical configuration remain unchanged. symmetry-center uses no
inputparas file. A9 is now implemented; see the A9 extraction section below. A10 coordinate conversion is now implemented; see the A10 section below.

### symmetry-center interactive real-data selection and replay

The canonical manual workflow is now:

```text
user selects suitable diffraction beam
→ user boxes approximate centrosymmetric feature
→ accepted selection stored in beam-local coordinates
→ existing symmetry-center legacy centrosymmetry refinement
→ zone-axis origin and saved diagnostics
```

The beam is chosen according to the relevant diffraction condition, not catalogue
position. The useful feature may occupy only part of the disk. Mouse coordinates
provide coarse initialization only; the final symmetry center is determined by
`origin.py`, whose numerical implementation is unchanged.

From the repository root, in the existing GUI-capable Python environment:

```bash
MPLBACKEND=QtAgg python examples/a8_origin_interactive.py \
  --pattern reference_data/pattern.dat \
  --catalogue tmp/a7_boundary_catalogue_r38_final/beam_catalogue.csv \
  --central-center 127 126 --radii 38 38 --selection-radius 38 \
  --output tmp/a8_manual_selection
```

This opens the image and waits for the user's beam click and feature rectangle.
It does not preselect either. QtAgg uses the optional PySide6 GUI environment;
a native interactive matplotlib backend may be used instead. No new dependencies
are added. Headless rendering backends cannot perform mouse selection and produce
an explicit launch error. The current A7 CSV is a partial catalogue with three
observed diffracted disks; observed status does not imply suitability for symmetry-center.

1. Click the intended diffraction disk. Association uses Euclidean distance to
   catalogue geometric centers, with configurable `--selection-radius`; ties use
   stable beam ID. The click never replaces the beam center. Default eligible
   statuses are `observed` and the central `(0,0)` `reference`. An explicit
   `--allow-unreliable` flag permits weak/missing/edge-clipped predictions and is
   recorded; it does not certify their scientific reliability.
2. Release, then drag the coarse feature rectangle. Preview shows stable ID,
   lattice indices, geometric center, click distance/status, detector/local ROI
   bounds, size, center and feature offset. Percentile contrast is display-only.
3. A mouse rectangle normally has fractional coordinates. Its nearest integer
   pixel-center ROI (ties to even) is automatically previewed with a green dashed
   outline alongside the orange raw rectangle. This proposal is not applied or run.
   **Accept** explicitly confirms the displayed proposal and runs symmetry-center; **Snap ROI to
   grid** confirms it without running. Raw mouse bounds remain in provenance.
   Collapsed or otherwise invalid previews require reselection.
   No clipping, silent rounding or implicit change to physical beam geometry occurs.
4. **Accept** validates and runs symmetry-center. **Reselect ROI**, **Reselect beam**, and
   **Cancel / Close** allow correction or abandonment. No analysis output is written
   until acceptance. Output directories must be new.

ROI bounds use `(x0,y0,x1,y1)`, never row/column interchange. Beam-local bounds
subtract `(beam_x,beam_y)` from both corners; the inverse adds the current beam
center. This lets a local feature definition follow changed beam geometry when
creating a new selection. Replay instead preserves and verifies the originally
accepted geometry; it does not silently update an old measurement.

Default search half-size is half the accepted ROI width/height; its center is the
ROI center, with the backend's upper-exclusive unit-step search. The default crop
half-size remains physical disk width/height divided by four. A snapped integer ROI
therefore supplies integer trial centers with the current R=38 crop half-size 19.
Advanced CLI options `--search-half-size SX SY` and `--crop-half-size HX HY` are
optional. Advanced search bounds must remain inside the coarse ROI.

Validation rejects zero-size ROIs, detector-exterior or ellipse-exterior ROI corners,
nonintegral crop bounds, and insufficient room for every trial crop. The coarse
ROI bounds the search centers, not every crop: surrounding crop pixels may extend
outside the ROI, but must remain in the detector. Invalid previews show the reason
and permit reselection; no numerical thresholds or crop semantics are changed.
Inputs changed on disk since display require restarting selection.

`origin_selection.py` implements pure coordinate/selection helpers plus provenance
I/O; all mouse handling is in the example script. It calls existing
`measure_symmetry_origin()` with an explicit `SymmetryReference` and configuration.
The new helper does not implement a second symmetry scorer or origin formula.

Every accepted run saves:

- `a8_selection.json` (schema version 1): absolute pattern/catalogue paths and
  SHA-256 hashes, selected stable ID/indices/status/geometric center, click/distance,
  eligibility/radius policy, raw and accepted detector ROI, beam-local ROI, physical
  geometry, central center, search/crop settings and full numerical measurement;
- `selection_context.png`: accepted coarse ROI, beam footprint/center, feature seed,
  final symmetry center, central center and derived zone-axis origin;
- `symmetry_diagnostic.png` and `measurement.json`: existing backend crop, rotated
  crop, absolute difference and numerical diagnostics;
- `interactive_preview.png`: accepted interactive view, including displayed values.

Non-interactive replay needs no window or mouse input:

```bash
python examples/a8_origin_interactive.py \
  --selection tmp/a8_manual_selection/a8_selection.json \
  --output tmp/a8_manual_replay
```

Replay selects Agg internally, verifies input hashes and catalogue identity/geometry,
checks local/detector ROI consistency, calls the unchanged backend, and writes to a
fresh directory. Absolute input paths make replay independent of process CWD.
Repeated launches can create independent beam/feature records; their measurements
can be passed to existing `combine_origin_measurements()` with explicit tolerance.
No automatic multi-reference combination is performed by this UI.

`tests/test_origin_selection.py` covers click association and rejection, explicit
status overrides, local/detector round trips and geometry shifts, invalid ROI/crop
bounds, provenance and tampering checks, direct/backend replay equality, headless
CLI replay, and construction of the UI without auto-selecting a real beam.
The user has confirmed real symmetry-center manual validation; A9 leaves this result unchanged. A9 is now implemented; see the A9 extraction section below. A10 coordinate conversion is now implemented; see the A10 section below.


### Direct or diffracted symmetry-center symmetry references

symmetry-center permits either the direct/central beam or a diffracted beam. The only physical
selection criterion is that the beam contains a local feature suitable for 180°
centrosymmetry under the current zone-axis/diffraction condition. Catalogue order
has no role. The central `(0,0)` entry with `reference` status is eligible without
an unreliable-beam override; existing diffracted eligibility safeguards remain.
The beam-local ROI and crop validation are unchanged.

Measurements and saved `a8_selection.json` records label `reference_kind` as
`direct` or `diffracted`, retaining stable ID, lattice indices, geometric center,
local ROI, best symmetry center and derived origin. Direct references are identified
by `(0,0)` or equality of the explicitly supplied reference and central centers. No center is silently replaced. The generalized formula
remains `origin = symmetry + 0.5*(central - reference)`; when the reference center
is exactly the central center, it reduces exactly to `origin = symmetry`. No
separate scientific branch is used. A9 is now implemented; see the A9 extraction section below. A10 coordinate conversion is now implemented; see the A10 section below.

## A9 — beam-supported signal, uncertainty and weight extraction

A9 consumes canonical A5–symmetry-center geometry; it never determines beam positions. symmetry-center has
been manually validated by the user. A9 does not need the symmetry-center origin to extract
pixel-grid selections: origin/scan conversion is implemented in A10 below.

Public APIs in `qcbed.data_preparation.extraction` are `ExtractionConfig`,
`WeightMode`, `build_beam_geometry_mask`, `construct_signal_fields`,
`extract_beam_data`, and `prepare_experimental_qcbed`. Immutable `SignalFields`,
`BeamExtraction`, and `ExtractionResult` retain full-frame intermediates and
beam-local products. They reuse A1 `BeamObservation` and `ReflectionROI` rather
than inventing scan geometry for `PreparedExperimentalDataset`. Abstract `(n1,n2)`
are stored separately from physical HK/HKL. No `inputparas` is read.

### Selection, geometry and ordering

An explicit set of unique stable beam IDs is required. Results follow the A7
catalogue's order, irrespective of selection-list order or brightness. Reference,
A6 defining and A7 observed sites are distinguishable by retained status/source.
Weak/edge-clipped and unconfirmed (`missing`) lattice predictions require explicit
selection and pass the same coverage gate; rejected/unknown statuses are refused. Predicted centers are always
used; observed centers are preserved separately. No rounding or center switching
occurs. Duplicate lattice labels do not merge beam IDs.

The ellipse uses strict `((x-cx)/rx)**2 + ((y-cy)/ry)**2 < 1`, with radii
`beam_radius_factor * physical_radius`, default factor **0.9** (38 becomes 34.2,
not 17.1). The full-frame mask is a **binary union**: beam=1, outside=0. It is not
the background complement. This intentionally differs from N's additive mask at
overlaps; repeated/overlapping ellipses never create weights greater than one.

Coverage counts all detector-grid pixels in the ellipse, including those outside
the image for the expected count. In-frame/expected coverage must meet
`minimum_beam_frame_coverage` (default0.90), inclusively (`>=`). A rejected beam remains in the ordered results
with status `rejected_coverage` and zero final weights, and does not contribute to
the retained full-frame union. Noise-valid/in-frame coverage is reported separately.
The explicit rectangular selection is independent of this ellipse coverage count;
its bounds and cropped ellipse mask are also retained.

### Signal, variance and weight are separate

Conventional full-frame signal is a copy of the input. The local expanded
selection has its minimum subtracted **before** semi-integration. The minimum
is calculated only over measured, valid local beam support (ellipse and local
radial mask), excluding placeholders and nonparticipating background.
`conventional_minimum_subtraction=True` is the compatibility default; disabling
it is an explicit alternative recorded in provenance. No smoothing, gain
correction, intensity normalization or detector correction is implicitly applied.

Gradient offsets are `(-q,0),(q,0),(0,-q),(0,q),(-q,-q),(q,q),(-q,q),(q,-q)`,
where q is `final_binning`. At detector pixel p:

```
G = sum(eight shifted intensities)/8 - I(p)
V_G = V(p) + sum(eight shifted variances)/64
coverage = (base_mask(p) + sum(eight shifted base masks))/9
```

There is no gradient minimum subtraction. The coverage remains floating until
local masking/rebinning. Integer shifts use zero extension, never wrapping or
reflection; incomplete detector/noise stencils receive zero final weight. This
hard validity gate is explicit and stricter than permissive fractional legacy
mask coverage: invalid detector/noise values cannot influence a weighted sample.
The independent-error variance formula does not model covariance introduced by
prior processing; exact DM border equivalence remains unverified.

`WeightMode.BEAM_MASK` needs no uncertainty or A4 mask. Optional explicitly supplied
Boolean detector validity can restrict support. `NOISE_CHARACTERIZATION` requires
explicit A4 variance (or declared sigma), final Boolean A4 noise mask and nonempty
provenance. Callers may use `variance_map_from_noise_model(pattern, fit.model)` and
`fit.noise_mask`; A9 does not refit noise or guess thresholds. The mode combines
noise validity with geometry before shifted support construction, and excludes
invalid stencil terms. A noise mask cannot silently enter beam-mask mode.

The validated objective is `sqrt(sum(w*(s*c-x)**2/V)/sum(w*x*x/V))`.
Thus w is support, **not** inverse sigma or inverse variance. Noise uncertainty
belongs in V. `support_weight` includes the coverage threshold and hard validity;
`noise_weight` is one (no additional uncertainty factor); `final_weight` is their
product, or zero for a coverage-rejected beam. Zero variance remains zero; A9 does
not repair the engine's singular zero-variance semantics or fabricate variance.

Explicit `sigma_map` is squared once; explicit `variance_map` is copied. Supplying
both is rejected. Provenance records this distinction. Without uncertainty,
variance and sigma remain absent, and these products are not yet a complete
refinement input.

### Extraction and source normalization

`selection_width`, `selection_height`, `final_binning` and `expansion_factor`
replace the legacy `selx/sely`, `binning`, and expansion settings. Selection origin
is `(cx-width/2,cy-height/2)`. Only integer origins are currently supported;
fractional origins fail rather than silently snapping or choosing an unverified
DM interpolation. Dimensions must divide by q. Odd expansion factors are
supported; even factors fail until the DM `Round((e-1)/2)` tie convention is
validated.

The local radial mask preserves N's **inclusive** circle of radius width/2,
centered at `(width/2,height/2)`. Rectangular selections explicitly retain this
legacy x-radius rule; no ellipse generalization is silently introduced.
Expansion repeats each pixel e times per axis. q-by-q sums are accumulated in
source row/column loop order, then signal/support divided by q². Variance sums are
scaled by e²/q⁴, matching N's semi-integration normalization. Support is thresholded
**strictly greater than `minimum_weight_coverage=0.875`**, then sampled at
e-spaced coordinates with offset `(e-1)/2`, as are signal and variance. This is not
a replacement with a detector-grid block mean. The source's treatment of repeated
samples and omitted correlations remains a compatibility convention, not a new
uncertainty theory. N's threshold is supported; E's different prethresholding is
not silently mixed into this mode.

Both signal modes yield `(height/q,width/q)` per beam; effective support may
differ. `blocks()` returns beam-major, then row, then column arrays, with stable
IDs retained alongside them. `tilt_beam_arrays()` provides the same values as
`(local row-major pixel,beam)` for eventual A1 integration. A10 must still supply
real scan geometry and determine any physical HKL mapping. No corner values or
canonical refinement adapter are constructed here.

### Initial four-beam validation (historical, superseded by full-catalogue audit below)

`examples/a9_extraction_validation.py` requires explicit beam IDs, selection size,
binning, expansion and physical radii. It writes new directories only; uncertainty
files require an explicit `--uncertainty-kind sigma|variance` plus provenance.
The current run used pattern.dat as supplied (256×256), manual physical radii38/38,
factor0.9, canonical bootstrap catalogue, diagnostic selection64×64, q=1, e=1.
These selection settings are diagnostic choices, not recovered legacy metadata.
Both conventional and gradient modes extracted all four selected sites without
coverage rejection, in catalogue order:

| Beam ID | Predicted center | Ellipse pixels / in-frame | Final nonzero weights | Conventional mean | Gradient mean |
|---|---|---|---|---|---|
| lattice_0_0 | (127,126) | 3681 / 3681 | 3207 | 692695.3601074219 | -637.2065734863281 |
| lattice_1_0 | (208,135) | 3681 / 3681 | 3207 | 319421.93994140625 | -304.9551696777344 |
| lattice_0_-1 | (137,45) | 3681 / 3681 | 3207 | 311369.90283203125 | -348.34735107421875 |
| lattice_0_1 | (117,207) | 3681 / 3681 | 3207 | 302773.99169921875 | -313.557373046875 |

Means describe all output selection pixels, not only positive-weight pixels.
The lower observed center (117,208) was not substituted. Figures, exact NPY arrays,
CSV tables and JSON provenance are under `tmp/a9_beam_mask_r38_final/conventional/` and
`tmp/a9_beam_mask_r38_final/gradient/`. They include whole-frame union, overlay, signal,
support fraction and individual signal/mask/weight crops. No variance images are
fabricated when uncertainty is unavailable.

Real noise-characterization validation remains blocked: sigma.dat is ASCII,
256×256, 65,536 finite values and numerically associated with the sigma model,
but its sigma-versus-variance provenance is unconfirmed. It was **not used**.
Default A4 gain exclusion40000 retains only884 initial/882 final-mask pixels and
fails fitting for too few effective fit points. No gain limit or Noicemask was
invented. Synthetic tests exercise genuine A4-produced variance and final masks.
Source-equation tests protect masks, gradient sign/stride/variance, strict threshold,
minimum subtraction, normalization, ordering, clipping, missing uncertainty,
immutable arrays and diagnostics. No claim of complete DM end-to-end equivalence
is made without fractional Warp and historical noise metadata fixtures.


### A9 follow-up: partial detector coverage and full lattice-site audit

A9 does not require a beam ROI to be fully contained within the detector. A beam
is retained when at least the configured fraction of its physical ROI is measured;
the default is 90%. Any out-of-frame pixels receive zero experimental weight and
therefore do not contribute to refinement. Specifically, the ROI is the existing
0.9× physical-radius ellipse (34.2px semi-axes for the current radius38px input).
`minimum_beam_frame_coverage=0.90` is inclusive. The denominator is the integer
pixel count of the complete ellipse, including out-of-frame positions; the numerator
counts its measured detector pixels. `window_frame_coverage` describes the separate
rectangular crop and never determines beam retention.

Beam region=1, outside beam=0, outside detector=0 final weight. All beam-local
outputs keep the configured fixed shape. Missing detector samples are not wrapped,
reflected, interpolated or extrapolated. Missing or incomplete detector-stencil/bin
outputs have signal0 and weight0. If uncertainty was supplied, those positions use
variance1 as an explicitly nonphysical safe placeholder, avoiding 0/0 in the
unchanged legacy objective. If uncertainty was absent it remains absent. Measured
variance values, including their existing zero-variance semantics, are not repaired.
A4 validity/provenance requirements remain unchanged.

Conventional minimum subtraction now uses measured valid pixels within the beam
ellipse AND local radial mask, not all rectangular pixels. Placeholder zeros and
background cannot change this minimum. This intentionally changes some conventional
means reported in the initial four-beam diagnostic; it is not intensity normalization.
Gradient mask contributions from missing neighbors are zero; the hard detector
stencil validity gate remains necessary because eight of nine contributions (8/9)
would still exceed0.875. It zero-weights individual incomplete samples rather than
rejecting a sufficiently covered beam. The gradient formula/sign is unchanged.

The explicit CLI `--all-lattice-sites` policy considers all non-rejected A7
predictions, including sites labelled `missing` or `edge_clipped`. Alternatively,
`--beam-ids` selects a subset. These policies never relabel A7 observations or
interpret `missing` as scientifically absent. IDs and catalogue order are preserved.
The initial four-beam result selected only reference/observed sites and refused
`missing` sites; it was NOT caused by boundary loss of the five additional nearby
predictions. All nine nearby ellipses are fully contained.

The current rerun uses unchanged 256×256 pattern.dat, canonical predicted centers,
manual radii38/38, factor0.9, threshold0.90, beam-mask mode, selection64×64, binning1,
expansion1. Both conventional and gradient modes considered17 sites, retained9,
and rejected8 for insufficient ellipse coverage. Each retained beam has3207 final
nonzero weights. No retained real beam is partially clipped in this catalogue;
95.38%, exactly90% and88.89% clipping scenarios are tested synthetically.

| (n1,n2) | Predicted center | A7 status | In-frame / complete ellipse pixels | Coverage | A9 |
|---|---|---|---|---|---|
| (0,0) | (127,126) | reference | 3681 / 3681 | 100.0000% | retained |
| (-1,0) | (46,117) | missing | 3681 / 3681 | 100.0000% | retained |
| (1,0) | (208,135) | observed | 3681 / 3681 | 100.0000% | retained |
| (0,-1) | (137,45) | observed | 3681 / 3681 | 100.0000% | retained |
| (0,1) | (117,207) | observed | 3681 / 3681 | 100.0000% | retained |
| (-1,-1) | (56,36) | missing | 3681 / 3681 | 100.0000% | retained |
| (1,1) | (198,216) | missing | 3681 / 3681 | 100.0000% | retained |
| (-1,1) | (36,198) | missing | 3681 / 3681 | 100.0000% | retained |
| (1,-1) | (218,54) | missing | 3681 / 3681 | 100.0000% | retained |
| (-2,0) | (-35,108) | edge_clipped | 0 / 3681 | 0.0000% | rejected: <90% |
| (2,0) | (289,144) | edge_clipped | 7 / 3681 | 0.1902% | rejected: <90% |
| (0,-2) | (147,-36) | edge_clipped | 0 / 3681 | 0.0000% | rejected: <90% |
| (0,2) | (107,288) | edge_clipped | 24 / 3681 | 0.6520% | rejected: <90% |
| (-2,-1) | (-25,27) | edge_clipped | 320 / 3681 | 8.6933% | rejected: <90% |
| (2,1) | (279,225) | missing | 369 / 3681 | 10.0244% | rejected: <90% |
| (-1,2) | (26,279) | edge_clipped | 369 / 3681 | 10.0244% | rejected: <90% |
| (1,-2) | (228,-27) | edge_clipped | 228 / 3681 | 6.1940% | rejected: <90% |

Current diagnostic outputs are `tmp/a9_coverage90_final/conventional/` and
`tmp/a9_coverage90_final/gradient/`. CSV/JSON audit every considered site, status, counts,
coverage, decision/reason and final nonzero weight count. Overviews show all sites
(cyan retained/red rejected) and the detector boundary. Per-beam arrays retain the
full local ellipse and separate detector-in-frame mask. Retained clipped windows
also receive an edge-support figure distinguishing measured/missing beam portions;
the synthetic diagnostic test exercises this since the current real retained set
has no edge beams. Noise-mode real validation remains blocked; sigma.dat is unused.
A5–symmetry-center geometry and refinement kernels are unchanged. A10 coordinate conversion is now implemented; see the A10 section below.

## A10 — current geometry compared with historical DM coordinates

A10 is coordinate conversion only: current symmetry-center origin + explicitly ordered
legacy-equivalent A7 vectors + actual A9 selection dimensions + user-confirmed
**binning=2** → original DM point1/2/3 equations → 2×2 solve → existing
`ScanGeometry.corner_a/b/c`. A2–A9 science and saved data are unchanged. No
historical input geometry is reconstructed or back-solved from output coordinates.
The historical values below are comparison references, not a matched-input
regression fixture, pass/fail target, or universal constants.

### Source mapping

N=`simplified_finaledition_withipp.s`:413–414 specifies central beam first in the
historical record. N:810–841 (also E:198–231) uses that point as point4 and subtracts
the zone origin. Python resolves the central beam by explicit stable ID and (0,0)
label, not row order. symmetry-center `zone_axis_origin_px` and the geometric center remain
separate inputs even when their measured values coincide.

The source uses ordered reference endpoints minus the zone origin as vec1/vec2.
N:625–630 explicitly expresses endpoint reanchoring (endpoint−old origin+new
origin), preserving displacements. The canonical helper likewise reanchors ordered
predicted displacement vectors at the measured symmetry-center origin, rather than subtracting
a shifted origin from unshifted raw beam centers. The low-level API also accepts
already resolved absolute endpoints. No observed center is silently substituted.

Current canonical ordered basis IDs `lattice_1_0`, `lattice_0_1` give columns
g1=(81,9), g2=(-10,81). The historical reference ordering is rightward then upward;
the canonical ordering is rightward then downward. The explicitly supplied integer
transform `T=[[1,0],[0,-1]]` gives `scan_basis=canonical_basis@T`, vec1=(81,9),
vec2=(10,-81). Reference endpoints are origin+vec1 and origin+vec2. The transform
is never inferred from catalogue order and must be integer/unimodular. Basis swaps
and sign changes are tested. Determinant sign is reported in detector coordinates
(x right, y down), without implying physical 3D handedness.

### Equations and implementation

`convert_scan_geometry` builds point4=reference−origin and
half_x=(selx−binning)/2, half_y=(sely−binning)/2. Point1 adds(-half_x,-half_y),
point2 adds(+half_x,-half_y), point3 adds(-half_x,+half_y). Matrix columns are
vec1 and vec2; `numpy.linalg.solve` yields coefficients(A,B) without swapping axes.
Tests independently compare against literal DM determinant/A/B algebra. The matrix
solve handles vec2y=0 without the removable division singularity of that expression.

DM output order and the existing engine agree: point1→corner_a, point2→corner_b,
point3→corner_c; columns advance b−a, rows advance c−a. A test exercises the actual
production `generate_scan_tilts` mapping. Scan rows=sely/binning and columns=selx/binning.
Expansion is recorded but does not enter the source corner equations. Singular,
nonfinite or condition>1e10 bases are rejected; the maximum condition is an explicit
API argument. Positive dimensions must divide by binning.

APIs: `convert_scan_geometry`, `scan_from_catalogue`, `assemble_prepared_dataset`
in `scan_conversion.py`. Immutable `ScanConversionResult` reuses `ScanGeometry`
and records only geometry/diagnostics, not duplicate image arrays. Dataset assembly
preserves retained A9 beam order and values, requires matching sampling dimensions,
and uses A1's existing owned immutable arrays. Missing variance is an explicit
readiness error, not an invented noise model. No physical HKL is required.

### Current measured inputs and results

Inputs are `tmp/a8_manual_selection/measurement.json`, the existing canonical A7
catalogue, and `tmp/a9_coverage90_final/conventional/validation_summary.json`.
No A5–A9 algorithm is rerun. Values used:

- symmetry-center origin=(127,126); central/reference point=(127,126); point4=(0,0).
- Ordered endpoints=(208,135), (137,45); vec1=(81,9), vec2=(10,-81).
- Actual A9 selection64×64; confirmed coordinate-conversion binning2; expansion1.
- Point1=(-31,-31), point2=(31,-31), point3=(-31,31).
- Determinant≈-6651; condition≈1.0123372719941535; negative determinant.

```
corner_a: [-0.4241467448503984,   0.33558863328822736]
corner_b: [ 0.33092768004811307,  0.41948579161028415]
corner_c: [-0.33092768004811307, -0.41948579161028415]
```

Historical DM comparison values only:

```
corner_a: [-0.508156,  0.407641]
corner_b: [ 0.407641,  0.508156]
corner_c: [-0.407641, -0.508156]
```

| Component | Current Python | Historical DM | Signed difference | Absolute difference |
|---|---|---|---|---|
| corner_a.A | -0.42414674485039838 | -0.508156 | 0.0840092551496 | 0.0840092551496 |
| corner_a.B | 0.33558863328822736 | 0.407641 | -0.0720523667118 | 0.0720523667118 |
| corner_b.A | 0.33092768004811307 | 0.407641 | -0.0767133199519 | 0.0767133199519 |
| corner_b.B | 0.41948579161028415 | 0.508156 | -0.0886702083897 | 0.0886702083897 |
| corner_c.A | -0.33092768004811307 | -0.407641 | 0.0767133199519 | 0.0767133199519 |
| corner_c.B | -0.41948579161028415 | -0.508156 | 0.0886702083897 | 0.0886702083897 |

Maximum absolute difference=**0.0886702083897159**. This is significant (roughly
17–19% per-component magnitude difference), not six-decimal rounding. Origin
translation is zero in the current geometry; corner signs/order agree. The current
scan span relative to its reciprocal basis is smaller. Selection dimensions and
basis geometry set that span; historical output values alone cannot uniquely
separate their contributions. Inputs were not tuned to force agreement. No
matched-input golden fixture or historical reconstruction is used.

The saved A9 diagnostic has binning1 and64×64 arrays. The newly confirmed binning2
is explicit in A10 provenance and yields32×32 scan geometry; the saved A9 arrays
are not silently rebinned or relabelled. They cannot yet be assembled with these
corners as a ready dataset. Missing real variance is a separate readiness issue.

Reproduce into a NEW directory:

```bash
python examples/a10_scan_geometry_validation.py \
  --a8-measurement tmp/a8_manual_selection/measurement.json \
  --a9-summary tmp/a9_coverage90_final/conventional/validation_summary.json \
  --catalogue tmp/a7_boundary_catalogue_r38_final/beam_catalogue.csv \
  --central-beam-id lattice_0_0 \
  --canonical-basis-ids lattice_1_0 lattice_0_1 \
  --basis-transform 1 0 0 -1 --binning 2 --output tmp/a10_current_review
```

Current diagnostics: `tmp/a10_current_binning2/scan_geometry.json` and
`scan_geometry.md`. They include all points/vectors, transform, dimensions,
determinant/condition, source hashes, historical comparison and every signed/absolute
difference. There is intentionally no matched-input pass/fail flag. No inputparas
runtime dependency or post-A10 work is introduced.

## Corrected physical A9/A10 selection and shared binning (current workflow)

This section supersedes the earlier64×64 diagnostic coordinate comparison.
In the legacy DM workflow, the final pattern-matching `binning` value is used
first as the spatial offset for the eight-neighbor gradient calculation and later
again as the final semi-integration/downsampling factor. The confirmed current
setting is **binning2** in BOTH roles. The prior real A9 run used binning1; the
core already bound both operations to the same value, but that validation setting
was wrong for this workflow. The image must not be downsampled before the gradient.

`selx/sely` represent the physical QCBED selection size, not the Python diagnostic
crop window. Following the user's explicit physical-diameter selection policy,
`ExtractionConfig.from_disk_geometry(DiskGeometry(38,38), final_binning=2)` derives
76×76 exactly. This is an explicit current policy, not a claim that every historical
DM parameter record was automatically generated from diameter. Noninteger diameters
are rejected rather than rounded. Explicit physical sizes may alternatively be
provided in configuration; no coordinates are back-solved from historical targets.

Configuration distinctions (all shapes use rows,columns):

- `physical_selection_shape_px`: declared scientific selection. Current76×76.
- `selection_width/height`: actual numerical local extraction; when physical size
  is declared, these must match it exactly.
- `diagnostic_window_shape_px`: optional presentation metadata only. The old64×64
  choice may be recorded here, but does not size numerical extraction or scan geometry.
- `final_pattern_matching_shape_px`: derived selection/binning. Current38×38.
- `gradient_step_px`: explicit read-only legacy binding to `final_binning` (2).

A10's catalogue adapter refuses configurations without a declared physical selection.
Its CLI rejects a binning value different from the saved A9 extraction. It no longer
silently substitutes diagnostic crop dimensions or overrides sampling. Low-level
DM point construction and 2×2 solve are unchanged.

Current dimensional flow:

```
256×256 prepared input, unchanged
→ full-grid gradient with offsets (-2,0),(2,0),(0,-2),(0,2),
   (-2,-2),(2,2),(-2,2),(2,-2) [gradient mode only]
→ physical/local extraction76×76
→ expansion1 and2×2 semi-integration
→ final pattern-matching arrays38×38
→ A10 half-span (76−2)/2 =37 detector pixels
```

Variance uses the identical ±2 stencil, Vcenter+sum(Vneighbors)/64; support uses
sum(center+eight masks)/9. Existing normalization, boundary validity and strict
coverage>0.875 are unchanged. Conventional mode bypasses the gradient, keeps its
valid-support minimum subtraction, and still undergoes2×2 semi-integration. The
0.9-radius support ellipse remains34.2×34.2 semi-axes, independent of the76×76
selection. Nine-beam inclusion and minimum frame coverage0.90 remain unchanged.

Both real reruns in beam-mask mode considered17 sites and retained the same9.
Every retained beam has38×38 arrays and100% ellipse coverage. Conventional output
has886 nonzero final weights per retained beam; gradient has809. No sigma.dat/A4
mask was inferred and variance remains unavailable. A9/A10 sampling now agrees,
although variance readiness is still unresolved.

Origin/reference=(127,126), point4=(0,0), vec1=(81,9), vec2=(10,-81), and explicit
transform diag(1,-1) are unchanged. Points are(-37,-37), (37,-37), (-37,37);
determinant≈-6651 and condition≈1.01233727. Corrected corners:

```
corner_a: [-0.50623966320854,    0.4005412719891745]
corner_b: [ 0.39497819876710266, 0.5006765899864682]
corner_c: [-0.39497819876710266,-0.5006765899864682]
```

| Component | Previous64 selection | Corrected76 selection | Historical DM | Signed difference | Absolute difference |
|---|---|---|---|---|---|
| corner_a.A | -0.42414674485 | -0.506239663209 | -0.508156 | 0.00191633679146 | 0.00191633679146 |
| corner_a.B | 0.335588633288 | 0.400541271989 | 0.407641 | -0.00709972801083 | 0.00709972801083 |
| corner_b.A | 0.330927680048 | 0.394978198767 | 0.407641 | -0.0126628012329 | 0.0126628012329 |
| corner_b.B | 0.41948579161 | 0.500676589986 | 0.508156 | -0.00747941001353 | 0.00747941001353 |
| corner_c.A | -0.330927680048 | -0.394978198767 | -0.407641 | 0.0126628012329 | 0.0126628012329 |
| corner_c.B | -0.41948579161 | -0.500676589986 | -0.508156 | 0.00747941001353 | 0.00747941001353 |

Maximum absolute difference decreased from0.0886702083897159 to
**0.012662801232897314**: absolute reduction0.07600740715681858, **85.7192%**
improvement (about7× smaller). Residual disagreement remains; no exact equivalence
is claimed and no origin/basis/size was tuned to the historical targets.

Reproduce into NEW directories:

```bash
python examples/a9_extraction_validation.py \
  --pattern reference_data/pattern.dat \
  --catalogue tmp/a7_boundary_catalogue_r38_final/beam_catalogue.csv \
  --all-lattice-sites --radii 38 38 --physical-selection-from-disk \
  --diagnostic-window 64 64 --binning 2 --expansion 1 \
  --minimum-beam-frame-coverage 0.90 --signal-mode gradient \
  --weight-mode beam_mask --output tmp/a9_physical_review_gradient
```

Repeat for conventional mode with a different output directory. Feed its summary
to `examples/a10_scan_geometry_validation.py`, with the same existing symmetry-center measurement,
A7 catalogue, central/basis IDs, `--basis-transform 1 0 0 -1 --binning 2`.
Actual saved outputs are `tmp/a9_physical76_binning2/{conventional,gradient}/` and
`tmp/a10_physical76_binning2/`. The latter JSON/Markdown contains all component
differences; this section records the old/new comparison.

Regression tests independently construct ±2 gradients and variances on the original
grid and then integrate2×2; they distinguish ±1 and the erroneous bin-first path.
Other tests protect conventional minimum-before-integration, exact diameter
selection without rounding, final dimensions, explicit physical-selection readiness,
and invariance of A10 under changes to diagnostic-window metadata. A2–symmetry-center science,
A7 catalogue, symmetry-center origin, A10 equations and refinement/optimizer code are unchanged.
No post-A10 work is introduced.

## Opt-in fixed-radius subpixel center refinement

Known-radius arc detection first establishes a robust integer/local center
hypothesis. A deterministic bounded subpixel geometric refinement then improves
center precision without using intensity centroids or downstream scan-coordinate
targets. Existing detection and acceptance thresholds are unchanged. A10 comparison
is downstream validation only and is not part of the center-refinement objective.
This path is opt-in because A9 still explicitly rejects fractional crop origins.

`subpixel_arc.py` exposes `SubpixelArcConfig`, immutable `SubpixelArcResult`,
`refine_arc_center` and `refine_validated_candidate`. `locate_central_disk`,
`detect_lattice_candidates` and `detect_reciprocal_lattice_seed` accept optional
`subpixel_config`. With None, existing grid localization remains unchanged. Results
retain the original center, refined center, shift, objective values, configuration,
frozen support coordinates, acceptance/fallback reason and active-bound flag.

The two-variable float64 trust-region bounded least-squares fit has default bounds
±1px per coordinate, fixed axes, deterministic tolerances1e-12 and at most200
function evaluations. For displacement d to a support point and axes(a,b), its
signed radius residual is `norm(d)-norm(d)/norm(d/(a,b))`, exactly the existing
radial-ray ellipse convention. Expected normals are normalized `d/(a²,b²)`.
Additional residuals are the sign-independent normal cross product divided by
sin12degrees and center-vote components divided by2.5px. These scales match the
existing default validator's angle and vote scales; they are explicit fit settings,
not relaxed validation thresholds. Normal polarity for center votes is frozen at
the initial center. Reported objective is summed squared residuals/support count.
There is no free radius, intensity centroid, brightest-pixel, scan-corner or
historical-DM term.

Only already validated support enters the optimizer; membership never changes
inside a fit. Once afterward, the unchanged arc validator runs on that frozen
support (it may discard now-incompatible points, never add new ones). A failed
validator, unsuccessful/ill-conditioned solve, out-of-bounds proposal or increased
objective retains the original center. Accepted bound-limited solutions are flagged,
not silently given larger bounds. Exact objective equivalence is allowed.

### Real-data result: one predeclared ±1px run

The real seed previously used explicitly reviewed/revalidated A6 candidates,
not the entire automatic proposal list. The experiment consumes that saved
selection/support from `tmp/a7_boundary_catalogue_r38_final/validation_summary.json`
and checks the pattern hash. A5 reruns its unchanged robust grid search first.
The physical radius remains38/38. No parameters were changed after the downstream
historical comparison.

| Center | Initial | Refined | dx,dy | Shift magnitude | Objective before → after |
|---|---|---|---|---|---|
| Central | (127,126) | (126.7480175314,125.9368156961) | (-0.2519824686,-0.0631843039) | 0.2597834112 | 3.8055152514 →3.7604634183 |
| Upper defining | (137,45) | (136.0000000000,44.2018423363) | (-1.0000000000,-0.7981576637) | 1.2794747579 | 5.4027823240 →4.1773368242 |
| Right defining | (208,135) | (207.8935050874,135.0988322037) | (-0.1064949126,+0.0988322037) | 0.1452892663 | 4.0128022087 →3.9996324694 |

The upper center reaches its lower x bound. It is not an unconstrained optimum;
no bound expansion or tuning was performed. Euclidean shift can exceed1px while
each component remains within the configured ±1px box.

Existing A6 canonicalization naturally changes g1=(81,9) to
(81.1454875560,9.1620165075), and g2=(-10,81) to
(-9.2519824686,81.7349733598). Component changes are
(+0.1454875560,+0.1620165075) and (+0.7480175314,+0.7349733598).
Lengths increase0.1626181065px and0.6419956868px; included angle changes−0.6815472088°.

Unchanged A7 bootstrap/catalogue logic re-confirms the lower site at
(117.4960350628,207.6717890559); the previous observation was(117,208).
The catalogue still lacks sufficient independent observations for a global fit.
The same nine sites pass the unchanged A9 coverage rule, each at100% coverage.
No new confirmations or beam counts were forced.

**A9 extraction is blocked for the refined geometry.** Both conventional and
gradient attempts raise the existing fractional-selection-origin error. A9's
integer crop semantics do not establish a DM fractional interpolation rule.
No rounding, interpolation, resampling or intensity/variance/weight changes were
introduced to bypass it. Therefore this experiment does not provide new subpixel
A9 intensity arrays or claim an end-to-end ready dataset.

Geometry-only A10 still runs unchanged with76×76 selection, binning2 and the
existing explicit scan-basis transform. The independently validated symmetry-center direct-beam
symmetry origin remains(127,126); it is not replaced by the geometric center. The
symmetry-center direct-beam relation reduces to its symmetry center regardless of the equal
central/reference geometric-center shift. Point4 is now(-0.2519824686,-0.0631843039).
No symmetry-center search or formula is changed.

| Stage | corner_a | corner_b | corner_c | Maximum historical difference |
|---|---|---|---|---|
| A: diagnostic64 | (-0.4241467449,0.3355886333) | (0.3309276800,0.4194857916) | (-0.3309276800,-0.4194857916) | 0.08867020839 |
| B: physical76 integer | (-0.5062396632,0.4005412720) | (0.3949781988,0.5006765900) | (-0.3949781988,-0.5006765900) | 0.01266280123 |
| C: physical76 subpixel | (-0.5043324934,0.3969228880) | (0.3961016251,0.4978563308) | (-0.4024079398,-0.4970171545) | 0.01153937494 |

Signed C-minus-historical differences, in(A,B) order, are:
corner_a=(+0.003823506573,−0.010718112034),
corner_b=(−0.011539374938,−0.010299669244),
corner_c=(+0.005233060236,+0.011138845493).
Maximum error decreases8.87186%, but component RMS increases from0.0090060600 to
0.0093112502. Thus the change is modest/mixed, not evidence that localization
precision explains the remaining discrepancy. Other geometry/calibration
differences remain plausible; this task makes no further adjustments.

Reproduce: `python examples/subpixel_geometry_validation.py --output tmp/subpixel_review`.
The saved run is `tmp/subpixel_arc_final/`: JSON/Markdown three-stage comparison,
per-center overlays (integer/refined centers, fixed38px boundary and original
support), unchanged A7 catalogue figures/CSV, A9 coverage and explicit extraction
blockers. Tests cover exact subpixel circular/elliptical evidence, disconnected
arcs, strong internal features, bounded shifts, solver/validator fallback,
determinism and opt-in A5 behavior. No A7 thresholds, symmetry-center formula, A9/A10 scientific
equations, physical radii, physical selection, binning, refinement or optimizer
algorithms were changed.


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

### Adapter objective audit and ordering contract

The engine consumes `x=prepared intensities`, `v=scientific variance`, and
`w=support`. The objective remains
`sqrt(sum(w*(s*c-x)^2/v) / sum(w*x*x/v))`, accumulated beam outermost, tilt
innermost. YNORM remains fixed for IREF>0, `sum(w*c*x/v)/sum(w*c*c/v)` for IREF=0,
and `sum(w*abs(x))/sum(w*abs(c))` for IREF<0. Even the last branch does not remove
variance from the objective. `beam_mask` therefore does not eliminate the engine's
variance requirement. The gradient adapter sets GRADIENTOPT=1 so the existing
engine transforms calculated intensities once; already prepared experimental
signals are not differentiated again.

Masked mode selects `w != 0` before division and requires positive variance there.
Legacy mode does not short-circuit zero weight: zero variance can still produce
0/0. The adapter preserves all values and rejects missing variance or invalid
variance at evaluated points. It uses **no placeholder**, including at zero-weight
points. Existing zero-weight/zero-variance entries remain valid for masked mode.
Legacy mode readiness instead requires positive variance everywhere. Objective
and YNORM implementations were not changed.

Assembly retains extracted beams in A9 catalogue order, omits rejected diagnostic
slots explicitly, and records both original extraction-block and retained column
indices. Within each beam, `tilt=row*scan_columns+column`; beam-major flattening is
`array.T.ravel(order="C")`. All source arrays are copied into owned read-only
storage. Physical HK duplicates are allowed; stable IDs must remain unique.
Geometric `(n1,n2)` are stored separately from optional physical indices.

| Refinement column | Stable ID | (n1,n2) |
|---|---|---|
| 0 | lattice_0_0 | (0,0) |
| 1 | lattice_-1_0 | (-1,0) |
| 2 | lattice_1_0 | (1,0) |
| 3 | lattice_0_-1 | (0,-1) |
| 4 | lattice_0_1 | (0,1) |
| 5 | lattice_-1_-1 | (-1,-1) |
| 6 | lattice_1_1 | (1,1) |
| 7 | lattice_-1_1 | (-1,1) |
| 8 | lattice_1_-1 | (1,-1) |

Each real beam has 5,776 raw selection pixels and 1,444 final samples.
Conventional support totals 7,974 nonzero points; gradient totals 7,281.
Zero-weight entries remain in their original positions. The symmetry-center zone origin is
(127,126), and scan corners remain:

- A = (-0.50623966320854, 0.4005412719891745)
- B = (0.39497819876710266, 0.5006765899864682)
- C = (-0.39497819876710266, -0.5006765899864682)

Physical radius factor .9, beam frame coverage >=.90, strict final support >.875,
±2-pixel gradient stencil, variance propagation, minimum subtraction and 2×2
semi-integration all remain the existing A9 implementation. Input image shape
remains dynamic. Historical corners never enter assembly or adaptation.

## Real experimental variance closure — A4 interface and saved-sigma replay

This supersedes the earlier closure's *variance unavailable* readiness status for
the current reference dataset. Historical variance=None outputs remain untouched.
No Patch A geometry, crystallography, propagation, objective or optimizer equation
has changed.

### Canonical future production source

```text
raw CBED -> OTF correction -> corrected CBED -> A4 characterize_noise
         -> sigma_map_from_noise_model -> variance_map_from_sigma_map -> A9
```

A4 remains the canonical production source. `characterize_noise()` returns the
fitted `NoiseModel` and its actual final `noise_mask` (detector/gain support,
10-sigma blemish rejection and positive smoothed intensity). In noise-characterization
weight mode, pass this mask explicitly; it is not reconstructed from variance.
`variance_map_from_noise_model()` retains its public signature and now delegates
through the exposed sigma stage. No fitting thresholds or fit equations changed.

The existing evaluated model is `a*abs(I)+b*sqrt(abs(I))+c`; on the positive-intensity
fit domain this is `a*I+b*sqrt(abs(I))+c`. Its established negative-intensity
extension is deliberately preserved, as are nonnegative coefficient requirements.
This is not a replacement by `sqrt(I)`. The current default gain threshold 40000
is inappropriate for the reference pattern; no reference refit was attempted.

`variance_map_from_sigma_map()` accepts a finite, real, nonempty 2-D NumPy sigma
map, rejects negative values (including tiny negatives), and returns same-shape
owned immutable float64 `sigma**2`. No abs correction, epsilon floor or nonfinite
cleanup is used. Squaring overflow is rejected. Zero sigma remains zero, requiring
separate positivity validation before using that pixel on active support.

### Legacy/reference source

`load_legacy_sigma_map(path, expected_shape=...)` reads whitespace-delimited saved
sigma images, retaining row/column orientation and dynamic shape. It does not
identify an arbitrary file's physical units by filename and does not square its
output. `tests/fixtures/legacy_reference/sigma.dat` is accepted as a saved legacy per-pixel **sigma**
map: 256×256, 65,536 finite positive values, range 46.5183–122206.0. Its variance
is the square, not the unsquared file and not sigma squared twice.

The exact historical save statement was not recovered. Source semantics and
numerical validation support this interpretation; it is not claimed as recovered
save-statement provenance. Saved sigma is a compatibility/replay input only;
future production does not require a sigma.dat file.

Both A4 evaluation and saved-file replay converge through the same sigma-to-variance
API. A9's existing `sigma_map=` entry also uses that API. The alternative
`variance_map=` entry accepts variance directly, with explicit provenance.

### Exact A9 propagation reused

For the current `expansion_factor=1`, `binning=2`, signal semi-integration is an
average of four original pixels, not a sum:

```text
S_out = (S00 + S01 + S10 + S11)/4
V_out = (V00 + V01 + V10 + V11)/16
```

Conventional minimum subtraction affects signal only. Variance is cropped on the
same physical 76×76 window and averaged with squared coefficients; it is never
minimum-subtracted. This retains the established treatment of the subtracted
minimum as a signal baseline, without adding a new baseline-estimation covariance
model.

Gradient fields are constructed before cropping, with the eight neighbors at
(±2,0), (0,±2), and (±2,±2):

```text
G_p = sum(I_neighbors)/8 - I_p
V_G,p = V_p + sum(V_neighbors)/64
V_gradient_out = sum_{p in 2×2 block}(V_p + sum(V_neighbors(p))/64)/16
```

For independent detector-pixel uncertainties, the four ±2 stencils in a 2×2
block occupy distinct parity classes, so they have no shared input pixels within
that block. This validates the marginal output-variance formula. Covariance
between neighboring output bins is not added to the existing objective model.
The more general expansion-factor path keeps its existing legacy `e²/q⁴` variance
convention; this patch does not establish new interpolation statistics for it.

Support remains independent: radius 38 px, beam radius factor 0.9, 76×76 physical
selection, final 38×38, frame coverage >=0.90, strict support >0.875. No variance
inversion enters the support map. Out-of-frame/incomplete outputs retain signal=0,
variance placeholder=1 and weight=0. No mask or variance floor was fabricated.

### Reference diagnostic and measured results

```bash
python examples/a9_real_variance_closure.py --output tmp/a9_real_variance_closure
```

Use a new output directory. The example reads the original pattern/sigma and the
immutable archived A7/symmetry-center/A10 geometry in `tests/data/real_variance`; hashes check
input association. It does not rerun geometry detection. It writes ordinary NPY
signal/variance/weight arrays plus JSON and Markdown provenance/statistics.

| Quantity | Conventional | Gradient |
|---|---:|---:|
| Shape | (1444,9) | (1444,9) |
| Finite values | 12996 | 12996 |
| Active values | 7974 | 7281 |
| Support per beam | 886 | 809 |
| Active minimum variance | 4911091.31189375 | 6251955.962016113 |
| Active maximum variance | 2527281206.4593754 | 2410438153.199834 |
| Active median variance | 49783785.29718749 | 58697909.21734375 |
| Active nonpositive/NaN/Inf | 0 / 0 / 0 | 0 / 0 / 0 |

Beam order stays `(0,0),(-1,0),(1,0),(0,-1),(0,1),(-1,-1),(1,1),(-1,1),(1,-1)`.
Both modes match the old variance=None signal/support arrays exactly, including
archived signal and weight hashes. Corners remain exactly:

```text
A = (-0.50623966320854,  0.4005412719891745)
B = ( 0.39497819876710266, 0.5006765899864682)
C = (-0.39497819876710266,-0.5006765899864682)
```

The real dataset now contains variance. Real adapter execution is not claimed:
physical zone axis/original hkl basis remains unconfirmed. A known synthetic
adapter test verifies variance acceptance and safe weight-zero input preparation,
without objective evaluation. Real propagation, YNORM, chi-square and optimization
were not run. B6 remains pending confirmed physical indexing and a strict-compatible
physical sample; historical A4 calibration/mask provenance remains separate.

Validation for this closure: **89 focused tests passed**, **2694 full-suite tests
passed**. Focused checks include sigma validation/immutability, unchanged A4
values, hand-calculated crop/2×2 and ±2 gradient variance, placeholder handling,
real archived signal/weight/corner regression, and synthetic adapter preparation.

## B6a downstream preflight

The completed variance closure now feeds a read-only B6a input gate. Preparation
manifests provide exact signal/variance/weight hashes, beam metadata and A10
geometry for regression checks. The real CLI reuses the existing closure rather
than rerunning or changing A7–A10. Physical indexing must be explicitly supplied
and confirmed before a real preflight can pass; the confirmed FCC Al assignment
below now satisfies that requirement. Synthetic and real preflight pass. No real propagation or
objective is evaluated. See [B6a input preflight](b6a_input_preflight.md).

### Confirmed real B6a closure

The later confirmed conventional FCC Al [001], (200)/(020) assignment closes the
previously missing physical-indexing input. Real B6a now passes with the existing
sigma-derived variance and unchanged A7–A10 arrays/geometry. See
[b6a_input_preflight.md](b6a_input_preflight.md). Earlier blocked milestones above
are historical, not the current input status. B6b objective remains unexecuted.

## Live experimental preparation

A. Data Prep opens with **New Data Preparation**. Import a full, single-channel
CBED image (numeric text, NPY, TIFF or PNG). An aperture image is optional.
Geometry precedence is explicit/manual override → optional disk-profile override
→ threshold geometry from the supplied aperture. With no image-derived geometry,
draw a diameter on the CBED and press Accept. A separate disk-profile image is
not required. No radius comes from array dimensions or the historical 38 px case.

With an aperture, existing A2/A3/A4 calibrate the detector PSF/OTF, correct the
whole pattern, and fit noise. Calibration and geometry remain separate backend
operations; their use of the same supplied aperture is explicit in provenance.
Without an aperture there is no correction or calibrated noise fit: the selected
Poisson/intensity approximation is sigma=sqrt(max(I,0)), variance=max(I,0), with
positive-intensity validity. Calibration failures never silently select this fallback.

Draw a rough centrosymmetric ROI, inspect its integer-grid preview, and press
Accept/Snap. Enter explicit physical extraction dimensions/binning and symmetry-center crop
half sizes. Run Data Preparation invokes the existing A5–A10 routines in a spawned
process. Inspect the resulting beam boxes and refined Zone axis marker, explicitly
confirm the origin, inspect Pattern/Sigma/Weight, then Exclude/Restore bad pixels
and Validate A. Masking changes only weights. A may validate before B is complete;
an empty project starts with a clearly incomplete Atomic Model draft, never the
FCC Al reference model. Validate B and C remain required before Run.

**Load Previous Preparation** is the secondary replay path. Successful preparation
writes an exclusive `.qcbed-data/preparation-…/` cache containing ordinary
experimental blocks, full corrected image, uncertainty arrays, a Patch A manifest,
selection/noise provenance and a compact prepared descriptor. Manual edits create
new immutable snapshots. Existing compact/verbose/legacy documents remain supported.

Raw-input drafts can be saved and reopened without prepared data:

```yaml
schema_version: 2
data_prep:
  pattern: full_cbed.dat
  aperture: aperture.tif  # optional
```

Existing atomic_model/qcbed sections may accompany the draft. After preparation,
normal modern Save retains raw references plus `data_prep.preparation` for replay;
selections belong to the cache, not an expanded YAML schema. An unprepared draft
cannot be compiled or run by the CLI. Paths are document-relative on input and
resolved on Save As.

Limitations are explicit: the current OTF backend requires matching array grids;
A4 retains its frame/fit requirements; A9 requires integer crop origins and does
not silently round refined beam centers. Fractional geometry or insufficient
lattice evidence produces an error. Live sparse A6 evidence (three to five arcs)
now reuses the existing two-validated-arc seed, requiring every detected candidate
to agree within the existing lattice residual tolerance. Deterministic selection
uses geometric confidence; detector thresholds and upstream radii are unchanged.
For `pattern.dat`, explicitly supplied 38 px radii give center (127,126), arc centers
(118,208), (136,44), (208,135), and seed (81,9), (-9,82). Sixteen predicted sites
produce nine A9 beams with full physical ROI coverage (76×76 selection, binning 2).
A7 insufficient/missing evidence remains nonfatal; seven sites fail A9 coverage.
No archived centers, catalogue, or origin are imported. The one-pixel seed difference
from the reviewed historical pair is retained, not tuned away. This is geometry/
extraction validation, not validation of the no-aperture Poisson noise approximation.


### Live A → B → C handoff

Validate A publishes an immutable, content-addressed `*.prepared.json` descriptor
in the document's `.qcbed-data` directory, using the same descriptor publisher as
Save. It links the confirmed symmetry-center origin, mask snapshot, beam identities/order,
Pattern/variance/Weight, ROIs, live lattice provenance and current A10 corners.
Sigma remains the square root of the recorded variance; full-frame uncertainty
products are retained by live preparation. B/C consume the validated A fields;
changes to beam order/labels, scan transform, scan geometry or source artifacts
invalidate the A receipt. Live A validation rejects mismatched A10 corners.
The descriptor is a preparation output/cache and can be opened with Load Previous
Preparation; it is not required to begin with raw CBED.

The real 38 px-radius live smoke uses g1=(81,9), g2=(-9,82), central/symmetry-center origin
(127,126), physical selection 76×76 and binning 2. Its identity-basis corners are
A=(-0.5008180871634687,-0.3962516733601071),
B=(0.40175516882344187,-0.4953145917001339),
C=(-0.40175516882344187,0.4953145917001339). These use the current live basis,
not the historical DM sign/basis convention or archived corners. The explicit
[001], (200)/(020) B controls map every live (n1,n2) to (2n1,2n2,0).
The smoke stops at stage dispatch with zero optimizer evaluations: its checkpoint
is marked failed by the intentional test exception, not a completed refinement.
No-aperture Poisson uncertainty remains labelled as an approximation.
