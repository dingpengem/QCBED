# CBED output and end-to-end validation

> Clean-import provenance update: the [canonical source inventory](legacy_source_inventory.md)
> supersedes source-location, binary-only, duplicate and missing-provider claims below.
> Driver: private predecessor implementation; ONEIT: private predecessor implementation; ONEITPARA: private predecessor implementation; optimizers now
> occupy their matching .f files; RNG: private predecessor implementation; ATOMICSCATFACT:
> private predecessor implementation (all tables confirmed). The four recently corrected core
> file bodies are unchanged. Historical test results below remain historical.


> FFT convention update (2026-09-16): historical measurements below used the
> earlier wave signs. The corrected full reference run and before/after results
> are recorded in [FFT sign correction](fft_sign_correction.md). Numerical
> output is comparable for this Al case, but that symmetry does not establish
> correct signs for general asymmetric models. Objective formulas are unchanged.

An explicitly selected non-legacy `objective_mode="masked"` now permits local
optimization of the reference dataset without changing legacy arithmetic.
See [masked objective](masked_objective.md) and [local runner](local_optimization.md).
The NaN/blocker discussion below continues to apply to the default legacy mode.
Reports prominently name the selected mode; masked reports compute per-beam
contributions using only the same nonzero-weight set as YNORM and the objective.

## Executable source evidence

The readable executable statements are authoritative. No conclusion below
depends on the related-version `private predecessor implementation` declarations.

* `private predecessor implementation` reads `BMNAME`, opens formatted unit 6 with that
  filename, and reads the requested beam list. The supplied FOR007 names
  `out.dat`; that name is not hard-coded. Beam H/K indices follow input order.
* `private predecessor implementation` reads the experimental header and separate intensity,
  variance and weight blocks. These experimental values are not copied into
  unit 6.
* `private predecessor implementation` extracts squared complex beam amplitudes, calls
  `GRADIENT` when requested, then copies `BEAM` into beam-major `THEORINTS`.
* `private predecessor implementation` calls `CHISQCAL`, and on `CSQ <= CHISMIN` rewinds
  unit 6 and writes `YNORM*THEORINTS(M+K+N)`, with beam outermost, scan row
  next, and scan column innermost. A blank WRITE follows each scan row.
  `private predecessor implementation` has the corresponding callback output path.
* `private predecessor implementation` defines the three global YNORM branches; lines 772-782
  accumulate the normalized weighted root error. No display normalization
  occurs in these writes.
* `private predecessor implementation` opens/rewinds `firstout.dat` on unit
  19 and writes pre-stage snapshots. The same write pattern recurs throughout
  the optimizer dispatch. The file can contain multiple concatenated snapshots,
  rather than one universally identifiable initial state.
* `private predecessor implementation` applies `PARMIN` and writes the final parameter
  record to unit 9 (`PARS.DAT`), using FORMAT 543, `(I6,X,200(E13.5))`.
  That formatting does not apply to CBED unit 6. `FOR0013.DAT` contains
  logging; `OCCPLOT` and `VGVSH` are structural diagnostics, not CBED arrays.

Unit 6 has **no header**, beam identifiers, scan coordinates, experimental
values or residuals. Its numbers are global-YNORM-scaled, optionally
GRADIENT-processed simulated values. They can be negative and are not generally
raw intensities. `WRITE(6,*)` is list-directed: physical wrapping, decimal
precision and spacing depend on the Fortran runtime, with no fixed field width.
An overwrite records a callback best/tie in the current bookkeeping context;
it is not guaranteed to be an explicit final-state reevaluation or an all-stage
global minimum. Initial unit-19 snapshots use the same numeric ordering, but
can be concatenated. The reader deliberately accepts only one dataset.

## Python representation and output

`SimulatedCbedData` owns immutable float64 snapshots with `(tilt, beam)` shape:
`raw_intensity` is squared amplitude, `objective_intensity` is after optional
GRADIENT and before scaling, and `scaled_intensity` is exactly
`ynorm*objective_intensity`. Tilts are undilated fractional reciprocal-cell
`(H,K)` coefficients. All map values and YNORM must be finite. An undefined
objective remains NaN/Inf as a diagnostic; it is never silently repaired.

`write_legacy_out_dat(path, data)` writes beam/row/column values with blank
row separators and portable `.17E` precision. This is numeric compatibility,
not byte-for-byte reproduction of a particular Fortran runtime.
`parse_legacy_out_dat(path, scan=..., beam_indices_hk=...)` requires external
scan and beam metadata because the file contains neither. It accepts physical
line wrapping and D exponents, rejects empty/truncated/extra datasets and
nonfinite values, and returns `LegacySimulatedOutput`. It cannot reconstruct
raw intensities, YNORM, GRADIENT flags or beam identities from file contents.

## Orientation and comparison

`reshape_scan_map(values, scan)` maps `values[row*scan_columns+column]` to
`map[row,column]`. Explicit loop-oracle tests establish column-fast ordering;
there is no transpose. This agrees with the source output address
`(beam-1)*NUML*NUMP + (row-1)*NUMP + column`.

`CbedComparison` contains experimental, initial simulated and final simulated
data and defines residual as **final scaled simulation minus experimental file
values**. Experimental input is already in the objective's comparison space;
GRADIENT is not applied again to it. Raw simulation remains available separately.
Beam order, scan dimensions and GRADIENT option must agree. Fitted scan
coordinates can differ from initial coordinates: acquisition indices remain
aligned and the report states this, without interpolation.

`plot_cbed_comparison` produces four panels for a selected beam, lower origin,
with one shared visible-point min/max for experimental/initial/final and a symmetric diverging
residual scale. Constant panels get an explicit nonzero color interval solely
for rendering. No independent normalization hides poor fits. Pixel axes are
scan column/row with their fractional H/K step vectors: a general oblique scan
must not be mislabeled as horizontal H and vertical K. A report includes one
readable figure for every beam rather than an overcrowded overview.

Per-beam diagnostics accumulate, in source tilt order,
`numerator_b = sum(w*residual*residual/variance)` and
`denominator_b = sum(w*experimental*experimental/variance)`.
The objective is `sqrt(sum(numerator_b)/sum(denominator_b))`, subject to
floating-point differences from regrouping source sequential accumulation.
Unweighted RMS is a separately labeled convenience statistic, not the objective;
its default uses only nonzero-weight display points. All-point and masked-only
RMS are also available, with NaN for an empty subset.
Zero weight does not mask zero variance in the executable equations.

## Driver validation and reports

`validate_refinement(state, experimental, config)` uses the existing driver for
positive LOOPS, independently reevaluates initial/final vectors, and verifies
objective and YNORM consistency. Equal vectors reuse the same initial map.
Diagnostic reevaluations are not counted as optimizer history.

For LOOPS <= 0 it directly evaluates one trial without running the guarded
optimizer driver. This explicit diagnostic path preserves undefined objectives
while allowing finite CBED maps to be inspected. It is not a successful driver
run, and its final map is the unchanged initial map. No refinement improvement
is claimed. Neither path reduces scan size, mesh or specimen thickness.

`plot_refinement_convergence` and `plot_parameter_trajectories` show the existing
stage endpoints, optimizer names, loop numbers and physical parameter metadata.
No dense evaluation history is invented. Undefined objectives are annotated.
`generate_refinement_validation_report` creates self-contained HTML with embedded
PNG figures, configuration/input summary, initial/final values, parameter table,
stage history, all-beam comparisons, contributions and limitations. Matplotlib
is optional and imported only when plotting; install the `visualization` extra.

## Reference case and unresolved validation

Run `examples/validate_reference_output.py` using the qcbed-python interpreter.
It writes only ignored `tmp/end_to_end_validation/`: `report.html`, `out.dat`,
`comparison.npz` and `summary.json`. It uses the real FOR007/in.dat pair with
18 active parameters, 39 by 39 = 1521 tilts and nine beams, original wave mesh
and thickness, and seed 2026. LOOPS is explicitly set to zero.

The full driver attempt failed its existing CHISMIN finite-objective guard.
The supplied experimental data contain **117 zero variances**. The unchanged
source objective divides by variance even at zero weight, producing NaN.
Evaluation-only diagnostics preserve this failure instead of altering weights,
variances or the objective. A meaningful full refinement therefore remains
blocked pending source/data evidence resolving this discrepancy. Finite YNORM,
raw/processed/scaled maps, map reconstruction, unchanged baseline and numeric
output round-trip can still be checked. Small synthetic integration tests also
exercise actual ANNEAL/AMOEBA driver stages with finite objective data.

`tests/fixtures/legacy_reference/out.dat` is empty, so historical final-output equivalence cannot
be tested. `firstout.dat` is a separate pre-stage artifact: the script attempts
a strict single-dataset parse and records maximum absolute and relative L2
differences when possible, without treating it as trusted final output.
Generated Python output is source-derived. No historical executable, original
compiler/RNG equivalence or production refinement reproduction is claimed.

Measured full-size evaluation-only results (2026-09-16, qcbed-python):

| Check | Result |
| --- | --- |
| Initial/final YNORM | 4723684.377058765 (same unchanged state) |
| Initial/final objective | NaN, from zero-variance source arithmetic |
| Baseline unchanged | yes |
| Generated output round-trip maximum absolute error | 0 |
| Initial scaled map vs single `firstout.dat` maximum absolute difference | 0.33890555339166895 |
| Same comparison, relative L2 difference | 4.120887163251928e-7 |
| Evaluation plus output/report wall time | 162.03 seconds; not a scaling benchmark |

The firstout differences have no established executable-equivalence tolerance.
The run produced all nine beam comparisons and an unchanged final diagnostic
map, not a refined solution. The focused tests passed (14); the complete suite
passed (1525) using
`C:\Users\peng863\AppData\Local\miniconda3\envs\qcbed-python\python.exe`.

## Weight audit and explicit display masking

The existing model fields are `intensities`, `variances`, and `weights` (plural),
immutable float64 arrays in `(tilt, beam)` order. The names in the input block
do not imply an objective skip policy.

All readable uses of experimental `WEIGHTS` were audited:

* `private predecessor implementation` reads variance then weight rows, without transforming
  either. It sums `TOTWEIGHTS` and prints "The total number of experimental
  data points being fitted is". This supports a point-inclusion interpretation
  for binary weights, but imposes no binary restriction on general weights.
* `private predecessor implementation` multiplies weight into both least-squares YNORM sums before
  division by variance. Lines 758-770 multiply weight into equal-count sums
  without using variance. Lines 772-782 multiply weight into objective terms
  before division by variance. **There is no weight-zero test or skip branch.**
* The readable `GRADIENT` definition is in the misleadingly named
  `private predecessor implementation`, called by `private predecessor implementation` and
  `private predecessor implementation`. Its executable assignments modify temporary maps
  and calculated `BEAM` only. It does not modify weights or variance, or repair
  zero variance. Other weight mentions in `private predecessor implementation` are declarations only.
* `private predecessor implementation` and `private predecessor implementation` write every row/column of each
  beam, without consulting weights. No readable legacy display-mask routine
  was found. No source rule requires zero-weight points to have zero variance.

Thus weight is effectively a binary point mask **in the supplied dataset**,
supported by the fitted-point count and actual 0/1 values. General source
arithmetic supports multiplicative weights; it does not establish guarded
exclusion of masked points or prescribe legacy masked visualization.

### Measured reference statistics

| Quantity | Count |
| --- | ---: |
| Total points (1521 times 9) | 13689 |
| Zero variance | 117 |
| Zero weight | 6140 |
| Zero variance AND zero weight | 117 |
| Zero variance AND nonzero weight | 0 |
| Nonzero variance AND zero weight | 6023 |
| Nonzero weight / displayed points | 7549 |

Unique weights are exactly `{0.0, 1.0}`. There are no fractional, negative,
NaN or infinite weights. Masks are not identical across beams.

| Beam H,K | Masked | Valid display | Zero variance | Overlap |
| --- | ---: | ---: | ---: | ---: |
| 0,0 | 682 | 839 | 0 | 0 |
| 1,1 | 685 | 836 | 0 | 0 |
| 1,-1 | 681 | 840 | 0 | 0 |
| -1,-1 | 683 | 838 | 0 | 0 |
| -1,1 | 684 | 837 | 0 | 0 |
| 2,0 | 682 | 839 | 0 | 0 |
| 0,-2 | 681 | 840 | 78 | 78 |
| -2,0 | 681 | 840 | 39 | 39 |
| 0,2 | 681 | 840 | 0 | 0 |

In the reference dataset, zero variance coincides with zero display weight/mask.
More precisely, the zero-variance set is a strict subset of the masked set:
beam `(0,-2)` has zeros at rows 0 and 1 across all 39 columns; beam `(-2,0)`
has zeros at column 0 across all 39 rows (zero-based acquisition indices).
No nonzero-weight zero-variance locations exist. This spatial observation does
not prove how the experimental preprocessing generated the zeros.

### Python display and diagnostics policy

`display_mask.valid_display_mask(experimental)` returns immutable booleans
`weights != 0`. Fractional and negative nonzero weights remain visible; they
are not thresholded, rounded, or multiplied into displayed intensities.
`experimental_mask_statistics` provides totals, per-beam counts, unique weights
and detailed locations of any visible zero-variance points.

`plot_cbed_comparison(..., mask_zero_weight=True)` applies the same per-beam
mask to experimental, initial, final and residual panels. Grey masked pixels
are Matplotlib masked-array entries, not zero intensities. Shared color limits
use only visible values; no independent panel normalization is applied. An
all-masked beam remains grey with fallback plotting bounds. Use
`mask_zero_weight=False` for an explicitly labeled unmasked debugging view.
Stored raw/scaled/residual arrays and output files are unchanged.

`CbedComparison.rms_by_beam` now excludes weight-zero pixels by default;
`rms_all_by_beam` and `rms_masked_by_beam` expose the other subsets. These are
unweighted convenience summaries only. Source objective contributions still
use every point, so zero-variance contributions remain nonfinite. An empty RMS
subset is NaN rather than a fabricated zero.

`plot_weight_diagnostics` displays unmodified weight, zero-variance indicator,
and combined zero-weight/zero-variance indicator using exact column-fast map
ordering and lower origin. The HTML report includes these three spatial maps
for every beam, total/per-beam mask counts, subset RMS, and explicit separation
of display policy from objective arithmetic. It lists beam H/K, tilt index,
row/column and intensity/variance/weight if visible zero-variance points exist.

Updated report: `tmp/end_to_end_validation/report.html`, regenerated by the
existing full-size evaluation-only example. The legacy objective remains NaN;
neither display masking nor GRADIENT resolves the arithmetic blocker. The next
step is to recover/audit the experimental preprocessing or matching historical
executable's handling of masked zero variances. Any proposed skip or variance
repair would require separate evidence and an explicitly approved numerical
policy; it must not be inferred solely from binary weights or these pictures.

Mask-update validation: 20 focused tests and 1531 complete-suite tests passed
using the exact qcbed-python interpreter above. The full-size report was
regenerated with 22 embedded figures; YNORM, source objective status, unchanged
baseline, and zero output round-trip error agree with the preceding run.
