# B6b — one deterministic fixed-input real objective evaluation

B6b closes the forward/objective interface for the confirmed real case. It is
**not refinement or optimization**. Each signal mode has one primary evaluation
and one identical repeat for determinism. Tests separately exercise the same path.
No physical parameters, indexing, scan geometry, variance or support are tuned.

## Inputs and execution path

`examples/fcc_al_real_preflight.yml` and `fcc_al_conventional.extxyz` provide the
confirmed conventional FCC Al cell: a0=4.049 Å, four fully occupied sites, B=0.86 Å²,
[001], ordered original g1=(200), g2=(020). Voltage200 kV, mesh16×16,
absorption(0,1,1), DSAP2.5 Å^-1, requested region depth840 Å, applied dilation1,
zero vacancy, and no explicit reflection multiplier records (unity) are retained.
The conventional cell contains (0,0,0), (0,.5,.5), (.5,0,.5), (.5,.5,0).
The active parameter layout is empty; all inputs are fixed. No optimizer is called.

The command loads the production YAML frontend and source compiler, regenerates
the established prepared datasets from `pattern.dat` and saved `sigma.dat`
(variance=sigma²), and runs B6a on those exact objects. No noise refit or geometry
redetection occurs. The attached B5 prepared sample is passed **by identity** to
`simulate_tilt_series_batched`; it is not manually reconstructed. Its direct
phase-grating method and unshifted reciprocal ordering remain unchanged.

`CompiledCrystallography.adapt_dataset` supplies native experimental arrays and
scan. `generate_scan_tilts`, the existing batched propagation, the existing
`refinement_trial._gradient`, and `compute_selected_objective` perform the work.
The general refinement trial path normally rebuilds gratings after parameter
mapping. B6b deliberately uses the already-approved attached sample instead;
it does not change that general path or claim a new optimizer integration.
The diagnostic currently requires configured batched execution and masked
objective policy; unsupported policies are rejected, never silently replaced.

Strict B5 retains two stored slices of2.0245 Å each, 415 occurrences, effective
depth840.1674999999964 Å. This whole-slice overshoot is existing behavior.
The detector geometry remains center(127,126), pixel basis(81,9)/(-10,81),
selection76×76, binning2, final38×38, and unchanged A10 corners. B3-derived
HK=(h,k) for this identity recast. Nine-beam order stays:
000, -200, 200, 0-20, 020, -2-20, 220, -220, 2-20.
The direct beam participates in both forward propagation and the objective.
The B6a FCC diagnostic remains reporting-only, with no channel filtering.

## Objective and decomposition

Both modes use their own unchanged observed signal, variance and weights.
Calculated raw intensity is nonnegative. Gradient mode applies the existing
zero-padded eight-neighbor-minus-center calculated-signal transform; its signed
values are valid and are saved separately from raw calculated intensities.

The selected policy is **masked**, IREF=-1. Existing global YNORM is
sum(w*abs(observed))/sum(w*abs(calculated objective signal)) over active points.
The unchanged engine returns sqrt(E/N), where
E=sum(w*(YNORM*calculated-observed)^2/variance), N=sum(w*observed^2/variance).
Accumulation remains beam-major. No independent per-beam scales are fitted.

Diagnostic point contributions select nonzero weights *before* division. Inactive
contributions are exactly zero, including placeholder zero variances. Per-beam
E contributions sum to E; E_beam/N sums to objective². For readable additive
root-error allocation, (E_beam/E)*objective sums to the root objective (within
floating-point summation tolerance). These allocations are not separate fitted
beam objectives. Reported ratios of signed sums may cancel in gradient mode;
a zero denominator is reported as null, never repaired.

## Measured real results

| Mode | Active points | Global YNORM | Objective |
|---|---:|---:|---:|
| Conventional | 7974 | 14870390.0181119 | 0.939675432747781 |
| Gradient | 7281 | 13447206.106343817 | 1.1896911365294378 |

Both raw calculated intensity arrays are(1444,9), min1.4349943290287628e-6,
max0.2559831234260208, median0.0035248862840709053; NaN/Inf/negative counts are0.
Complex amplitudes are(1444,9), with12996 finite real and12996 finite imaginary
components, maximum magnitude0.5059477477230439. All determinism repeat errors
(amplitudes, raw intensities, objective signals, YNORM, objective) are **exactly0**.

| Beam HK | Conventional E contribution | Gradient E contribution |
|---|---:|---:|
| (0,0) | 813831.8994988562 | 4518.866397694295 |
| (-2,0) | 618711.0369727366 | 2708.062600686007 |
| (2,0) | 108670.98823923779 | 2820.268697685904 |
| (0,-2) | 262145.96457844594 | 2969.499524109343 |
| (0,2) | 252213.91480925283 | 2855.4842761014825 |
| (-2,-2) | 86296.05363928458 | 1416.5267383530904 |
| (2,2) | 381510.55456545285 | 1207.8345535502035 |
| (-2,2) | 73024.68379721984 | 1164.3528622744184 |
| (2,-2) | 76653.51434294181 | 1135.9686185472635 |

These substantial residuals are reported without tuning. PASS concerns finite,
consistent integration, not goodness of fit. Historical `firstout.dat` is not a
matched regression: cell/motif, beam/scan definitions and experimental objective
inputs have not been shown equivalent. No legacy numerical equivalence is claimed.

## API, command and artifacts

`qcbed.objective_diagnostics.evaluate_preflight_objectives(compiled, datasets,
baseline_manifests=..., input_kind=..., physical_input_confirmation=...)` returns
report metadata and arrays; it performs B6a before any propagation.
`residual_diagnostics` provides the masked pointwise/beam decomposition without
running another objective or fitting a scale.

```bash
python examples/b6b_real_objective.py --output tmp/b6b_real_objective
```

The output directory must be new. It contains `objective_summary.json/.md`,
`beam_residuals.csv`, and separate conventional/gradient directories with full
complex amplitudes, raw intensities, objective signals, globally scaled signals,
observations, residuals, weighted squared residuals, variance, weights and tilts.
Each beam also has38×38 NumPy images. Figures use common linear observed/scaled
calculated limits per beam, and symmetric residual limits. Display limits do not
change arrays or support. Outside support unweighted residuals can be nonzero;
weighted contributions stay zero. JSON retains exact floating-point provenance,
source/resolved configurations, slicing and parameter layout.

## Next boundary

B7 is not started. Before refinement, choose/review explicit physical targets,
bounds and signal mode: the current configuration has **zero active parameters**.
Any integration of future parameter updates with this crystallographic path needs
its own tests. GUI and app rebuilding are deferred. B6b introduced no kernel,
optimizer, checkpoint or configuration-semantics changes.

## Validation results

58 focused tests passed. The full suite passed **2731 tests in148.95 seconds**.
`git diff --check` is clean. Tests cover both modes, exact sample identity, no
optimizer invocation, finite outputs, repeat determinism, beam ordering,
zero-weight exclusion, residual decomposition, and saved-array round trips.
