# B6a — input consistency before real propagation

B6a stops immediately before propagation. It does not compute calculated
intensities, YNORM, chi-square, objective or optimization. Existing B4b compilation
may prepare phase gratings; no incident wave is propagated. GUI work and app
rebuilding remain deferred.

## Current real-input status

**Confirmed real preflight PASS.** The user confirmed FCC Al in its conventional
cubic cell, zone [001], ordered original-cell g1=(200), g2=(020). The separate
`examples/fcc_al_real_preflight.yml` and `examples/fcc_al_conventional.extxyz`
implement this assignment; the user-editable canonical example is unchanged.

Physical provenance:

- `reference_data/FOR007.DAT` record 4 gives layer spacing 2.0245 Å. For the
  confirmed FCC [001] half-period planes, a0=2*2.0245=4.049 Å. This is a
  crystallographic derivation from the reference spacing, not an independently
  measured lattice constant. Its a0/sqrt(2)=2.863075359... Å agrees with the
  separately rounded legacy lateral length 2.8630754 Å. Do not use that lateral
  number as cubic a0 or adjust the two values to force bitwise equality.
- The four standard FCC sites have full chemical occupancy. Their two sites per
  half-period slab preserve FCC density. The legacy one-slice projected motif
  has two 0.5-weight records in half the conventional lateral area; these
  projected weights must not be copied as 50% chemical occupancy into the four
  conventional sites. This is not a claim of identical finite-mesh propagation.
- The same reference provides Al B=0.86 Å², 200 kV, 16×16 mesh, depth840 Å,
  DSAP2.5 Å^-1 and independent absorption (0,1,1). No new thermal fit is used.
- Two explicit slabs [0,.5), [.5,1) each have thickness2.0245 Å. Strict B5
  accepts the transverse-periodic cell. The full-period projection shortcut
  remains unsupported for production. Four original atoms remain four atoms.

B6b remains deferred. This confirms input consistency, not experimental calibration
or numerical equivalence to a different legacy projected cell.

## API and invariants

`qcbed.input_preflight.preflight_inputs(compiled, datasets,
baseline_manifests=..., input_kind=...)` accepts the existing B4b `CompiledConfig`
and a mapping of mode names to immutable `PreparedExperimentalDataset` objects.
The manifest mapping comes from the preceding validated preparation stage.
`require_physical_config(config)` rejects configurations without explicit
crystallography. `PreflightError` provides a machine-readable error `code`.

The gate verifies:

- The B5 sample, B4 indexing and B2 orientation share the same original lattice
  and exact `ZoneAxisCell` provenance; an equivalent but distinct zone cell fails.
- Strict integer-C-shear B5 compatibility, model species/coordinates/occupancy,
  slice geometries, depth and B metadata agree with the compiled execution state.
- All original hkl follow the user-defined ordered basis and satisfy the integer
  zone law exactly, including the central zero reflection.
- Independently evaluated original and internal reciprocal vectors agree; the
  experimental-basis linear combination agrees; all vectors lie in the beam plane.
- B2, B3 and B5 positive propagation directions agree.
- Beam IDs, coordinates, HK and array ordering match; arrays and A10 geometry
  match preparation manifests exactly; active variance is positive and finite.
- B4b and direct B4a adapter scans agree with the compiled scan exactly. Arrays
  remain unchanged, and adapter input preparation safely handles zero weights.

Float geometry uses `atol=rtol=1e-14`, as in B4b noncubic reciprocal regressions;
actual maximum vector/component errors are reported. Integer relations, beam
order, corners, configuration metadata and array hashes use exact comparisons.
Reciprocal errors are Å^-1; direction errors are dimensionless. The legacy
engine's printed-pi metric approximation is unchanged; this gate checks the
attached crystallographic bases rather than replacing kernel trigonometry.

A coherent user-selected basis swap produces a different table and scan basis;
it is not automatically rejected as scientifically wrong or corrected. Mixing
one basis with another sample/index table fails. Algebraic consistency cannot
independently confirm an experimental indexing assignment or detector calibration.

Both prepared signal modes can be inspected without overriding configuration.
The report records the configured gradient option and each adapter's gradient
option. A later B6b run must select its intended mode explicitly; a B6a report is
not a run snapshot or an authorization to execute an objective.

## Commands and reports

Explicit synthetic demonstration (no reference pattern or real indexing):

```bash
python examples/b6a_real_input_preflight.py --synthetic \
  --output tmp/b6a_synthetic_input_preflight
```

Real input, only after the assignment is confirmed:

```bash
python examples/b6a_real_input_preflight.py --config /absolute/path/real.yml \
  --physical-input-confirmation "Reference/source confirming this structure, zone and ordered hkl basis" \
  --output tmp/b6a_real_input_preflight
```

The confirmation records user provenance; it is not evidence generated by the
program. The CLI never chooses a development crystal as a real default.
Missing input or incompatibility writes FAIL diagnostics and exits nonzero.
Outputs must use a new directory; existing artifacts are not overwritten.
JSON and Markdown contain sample cells/T/slices/depth, complete beam tables,
per-mode/per-beam support and variance statistics, and maximum geometry errors.
No full experimental arrays are embedded in Markdown.

The real branch reuses `a9_real_variance_closure.prepare_reference_variance()` and
its archived A7/symmetry-center/A10 geometry, pattern/sigma hashes and signal/weight regression
checks. It does not rerun A7–A10 or invent a mask. The real sample must have the
same beam IDs/order and scan geometry; the gate never edits config to make it fit.

## Measured synthetic demonstration

This is explicitly synthetic: a two-site cubic Si demonstration, zone `[001]`,
original basis `(1,0,0)/(0,1,0)`, identity T, two 2 Å slices, requested depth 5 Å,
effective whole-slice depth 6 Å. It is not a real-pattern assignment.

Nine ordered experimental coordinates map to original `(n1,n2,0)` and internal
`(n1,n2)`. Maximum basis-algebra, original/internal reciprocal, beam-direction and
beam-parallel errors were all zero for this fixture. Both modes have shape
`(1444,9)` and unchanged corners. Conventional support is 886 per beam / 7974;
gradient support is 809 per beam / 7281. Synthetic active variance min=max=median
is 1.0 conventionally and 1.125 for gradient, with no nonpositive/NaN/Inf entries.

The earlier blocked real report is superseded by the confirmed-real result below.
Synthetic PASS must never be substituted for real PASS. B6b has not started.

Validation: **42 focused tests passed; 2704 full-suite tests passed**. Negative
cases cover missing confirmation/configuration, invalid zone basis, strict-B5
rejection, coherent versus mixed basis swaps, zone-cell provenance mismatch,
reciprocal corruption, original-site replacement, invalid variance, changed
signal, beam order and A10 corners. Execution sentinels reject propagation or
objective calls during the synthetic compilation/preflight test.

## Confirmed real command and measured result

```bash
python examples/b6a_fcc_al_preflight.py --output tmp/b6a_real_input_preflight
```

Use a new directory; the command refuses overwrites. The previous blocked report
was archived before creating the requested updated report. This case wrapper
calls the existing source compiler, unchanged real-variance preparation and B6a
gate. It adds FCC/legacy diagnostics only. It writes `preflight.json`,
`preflight.md`, `legacy_for007_mapping.csv`, and `conventional_body.dat` (17-digit,
headerless beam-major signal/variance/weight blocks). The YAML points at this
last artifact at the standard output location. For a different output directory,
update that explicit path before any future separately authorized run. No runner
or objective is invoked to serialize the already prepared arrays.

The nine conventional hkl / internal HK pairs, in retained Patch A order:

| (n1,n2) | Original hkl | Derived HK | FCC diagnostic |
|---|---|---|---|
| (0,0) | (0,0,0) | (0,0) | direct |
| (-1,0) | (-2,0,0) | (-2,0) | allowed |
| (1,0) | (2,0,0) | (2,0) | allowed |
| (0,-1) | (0,-2,0) | (0,-2) | allowed |
| (0,1) | (0,2,0) | (0,2) | allowed |
| (-1,-1) | (-2,-2,0) | (-2,-2) | allowed |
| (1,1) | (2,2,0) | (2,2) | allowed |
| (-1,1) | (-2,2,0) | (-2,2) | allowed |
| (1,-1) | (2,-2,0) | (2,-2) | allowed |

All zone-law values are exactly zero. Maximum basis-algebra error, original versus
internal reciprocal error, beam-direction mismatch and beam-parallel component
are all **0** in the measured preflight. Reciprocal units are Å^-1. All beam IDs,
order, signal, variance/support and A10 checks pass exactly. Both modes retain
(1444,9), with 7974 conventional / 7281 gradient active samples, all finite and
positive active variance. The geometry remains center(127,126), pixel
(81,9)/(-10,81), selection76×76, binning2, final38×38. A10 corners remain:
A=(-0.50623966320854,0.4005412719891745),
B=(0.39497819876710266,0.5006765899864682), C=-B.

No input blocker remains for a separately authorized B6b investigation. That next
stage still needs an explicit signal-mode choice (this YAML selects conventional)
and its own numerical validation. This patch executes no real propagation,
calculated intensities, YNORM, objective or optimization.

Confirmed-real validation: **258 focused tests passed; 2727 full-suite tests
passed in 129.09 seconds**. New tests cover the conventional nine-beam assignment,
FCC diagnostics, strict integer/parity legacy round trips, no-propagation/objective
sentinels, and exact serialized experimental-body round trip. `git diff --check`
is clean. The normal suite exercises existing numerical routines as regression
tests; no real B6b calculation was executed.

## Subsequent B6b status

The separately authorized [B6b fixed-input evaluation](b6b_real_objective.md) now
passes for both signal modes, reusing this preflight on the exact attached sample.
The B6a command itself remains read-only with respect to propagation/objective.
Earlier statements that B6b was deferred describe the B6a milestone. B7 remains
deferred; no optimization has run.
