# Mapped-physical refinement bounds

**Hard bounds are a user-selected production refinement policy, not legacy
Fortran behavior.** They prevent expensive, physically unwanted proposals from
entering simulation. They do not imply a guaranteed physical or numerical domain.

Three separate operations remain distinct:

1. Existing source coordinate conversion and physical mapping/sign folding.
2. Optional hard-bound rejection.
3. Optional numerical-invalid rejection during evaluation.

Readable `private predecessor implementation` defines UPDPARS and assigns
`DEPTH(I)=DABS(PT(COUNT+1))` at line 817. Python keeps that existing mapping.
D1/D2 search ranges are not reinterpreted as hard limits. No declarations from
the related-version private predecessor implementation are used to invent constraints. Optimizer
implementations, mapping, FFT, multislice, GRADIENT, YNORM and objective formulas
are unchanged by this addition.

## API and identification

```python
from qcbed import RefinementBound, RefinementBounds, evaluate_refinement_trial

bounds = RefinementBounds((
    RefinementBound("depth", (0,), lower=600.0, upper=1200.0),
    RefinementBound("scan_coordinate", (0, 0), lower=-0.65, upper=-0.35),
))
trial = evaluate_refinement_trial(
    state, vector, experimental,
    bounds_policy="reject", refinement_bounds=bounds,
    invalid_trial_policy="raise",
)
```

These frozen dataclasses use the existing layout's category and zero-based
physical target tuple, not transient active-vector positions. Bounds support
lower-only, upper-only, or both. Unspecified active parameters remain unbounded.
Both endpoints are **inclusive**. Bounds must be finite real values, not bool;
lower must not exceed upper. Empty-sided entries, duplicate targets, unknown
categories and nonexistent physical targets are errors. Known but inactive
targets remain metadata and are not enforced, including after
`select_refinement_categories(...)`.

The check reads mapped active values using the existing `_entries` semantics,
immediately after `apply_refinement_parameters` and before `_prepare_trial`,
`prepare_multislice_sample`, structure-factor preparation or forward simulation.
For example, raw depth -27925 maps to physical +27925 before comparison.
AMOEBA simplex coordinates are converted before the same check. No values are
clipped, reflected, repaired or replaced.

`bounds_policy="none"` disables enforcement. `"reject"` returns a
`BoundRejectedTrial` containing the mapped state, supplied evaluation vector,
coordinate mode, violation records and numeric `objective=+inf`. It is a
specialization of the existing rejecting-result protocol so unchanged optimizer
wrappers can process it; its reason is `bound_rejection`, not a numerical
exception. The low-level default (`bounds_policy=None`) inherits an explicitly
scoped driver policy if one exists and otherwise means `"none"`. An explicit
`"none"` disables bounds for that evaluation.

`run_legacy_refinement` accepts `bounds_policy="none"` and `refinement_bounds`.
The driver scopes the policy with a ContextVar and restores it on exit, including
exceptions. It enables the wrappers' existing infinite-score comparison path
without modifying any optimizer implementation. The evaluator's numerical
policy remains independently selected: bounds can reject while genuine
numerical errors still raise. Existing comparison rules cannot promote an
infinite score to retained state after a valid starting anchor.

For direct optimizer-wrapper use, the same public scope is available:

```python
from qcbed import refinement_bound_policy

with refinement_bound_policy(state, bounds_policy="reject",
        refinement_bounds=bounds, invalid_trial_policy="raise"):
    # Existing wrapper needs reject to permit +inf comparison; the scope keeps
    # actual numerical failures on the separately selected raise policy.
    result = run_legacy_anneal(..., invalid_trial_policy="reject")
```

Prefer the driver API for normal use. Bounds are checked in the parent before
any forward backend scheduling; no optimizer or worker globals are patched.

## Production example

`FullRefinementConfig` adds `bounds_policy` (package default `"none"`) and
`refinement_bounds` (default empty). The editable example explicitly chooses:

```python
REFINEMENT_CATEGORIES = ("depth", "scan_coordinate")
BOUNDS_POLICY = "reject"
DEPTH_MIN_ANGSTROM = 600.0
DEPTH_MAX_ANGSTROM = 1200.0
SCAN_RELATIVE_BOUND = 0.20
SCAN_ABSOLUTE_FALLBACK = 0.05
```

Depth bounds are in angstrom. Scan coordinates and fallback are dimensionless,
undilated fractional reciprocal-cell H/K coefficients, not angles.
`scan_relative_bounds(state, relative_fraction=.20, absolute_fallback=.05)`
creates bounds for active scan coordinates centered on their starting values,
with half-width `max(abs(start)*.20, .05)`. Negative starts are handled by
absolute magnitude. At zero, width is .05 rather than zero; the absolute floor
also dominates when abs(start) < .25 with these settings. Both settings must
be finite and positive. This formula and all example numbers are explicit user
choices, not inferred scientific limits. The helper does not adapt during a run.

The current edited input starts at depth 1000 A; the example retains that
value and selects seven active parameters. Neither selection nor bounds edit
FOR007.DAT. Run `examples/run_full_refinement.py --dry-run` with the documented
qcbed-python interpreter to inspect the full bound table and active layout.

## Rules and metadata-driven expansion

`RefinementBoundRules` accepts dictionaries and stores immutable copies. Resolve
them once against the initial active state, not anew against each trial:

```python
from qcbed import RefinementBoundRules, resolve_refinement_bounds

rules = RefinementBoundRules(
    specific_bounds={("depth", (0,)): (600.0, 1200.0)},
    category_bounds={
        "absorption": (0.5, 1.5),
        "debye_waller": (0.5, 1.5),
        "reflection_multiplier": (0.5, 1.5),
    },
    category_relative_bounds={"scan_coordinate": 0.20},
    relative_bound_absolute_fallback={"scan_coordinate": 0.05},
)
resolved = resolve_refinement_bounds(state, rules)
controls = FullRefinementConfig(driver_config, bounds_policy="reject",
    refinement_bounds=resolved, bound_rules=rules)
```

The example exposes these dictionaries as `SPECIFIC_BOUNDS`, `CATEGORY_BOUNDS`,
`CATEGORY_RELATIVE_BOUNDS`, and `CATEGORY_RELATIVE_FALLBACK`. They are explicit
user-selected example constraints, not authoritative physical or legacy limits.

Precedence is **specific target > category absolute > category relative >
unbounded**. An override replaces the whole interval; it is not an intersection
with lower-priority rules. Category rules expand to every matching active
parameter in source layout order. Eight active reflection multipliers receive
eight bounds from one category rule. `REFINEMENT_CATEGORIES=None` expands these
rules over all eighteen current parsed parameters without per-target edits.
The controlled seven-parameter subset resolves only those seven active entries.
Inactive specific rules remain in the original rule metadata for persistence.

Relative rules use half-width `max(abs(start)*fraction, absolute_fallback)`.
The fallback dominates exactly when `abs(start)*fraction < absolute_fallback`;
at equality both widths are identical. Positive fractions and fallbacks must be
explicit, finite, and supplied for the same categories. There is no implicit
near-zero threshold, inferred fallback, or adaptive update. Both positive and
negative starts use their magnitude for width. Each parameter has its own
center. Unconfigured categories have no bound.

`resolve_refinement_bounds(layout, rules)` also accepts a
`RefinementParameterLayout` whose initial values already represent mapped
physical starts. Prefer supplying `RefinementState`: the resolver then uses
the existing pack/map functions to derive physical starts and validates specific
targets against the complete model, including inactive targets. Do not supply
simplex coordinates or assume raw input signs are already physical. No
18-parameter layout or per-reflection index list is hard-coded.

Dry-run prints `bounds=[lower, upper]` or `bounds=unbounded` next to every active
parameter. Run config/checkpoints/HTML save both the original normalized rule
configuration and the fully resolved per-target table. The runner verifies that
re-resolution against the unchanged mapped initial state exactly matches the
supplied table. Resume restores the saved rules/table, ignoring edits to the
example's current dictionaries; changed rules or resolved limits are rejected.
Direct per-target configuration remains supported with `bound_rules=None`.

## Progress, reporting and restart

Bound rejections consume attempted evaluations but not successful evaluations
or retained-state updates. They produce `progress.jsonl` rows with:

- `classification` and `reason`: `bound_rejection`;
- global evaluation index, stage, loop and stage-local index;
- `active_vector` as supplied to the evaluator, with `coordinate_mode`;
- mapped physical parameters and all violations: category, target, value,
  lower and upper;
- `objective="inf"` in strict JSON, `ynorm=null`, `accepted=false`, and
  `retained_state_unchanged=true`.

For AMOEBA/SIMPARA, `raw_trial_vector` is also available from the callback.
ANNEAL/PARAMIN project their vector before the evaluator, so the pre-projection
raw vector is unavailable there and this field is null. The supplied projected
vector and mapped state are still recorded; no raw values are fabricated.

Numerical-invalid events retain their separate `invalid_trial` classification
and counters. Summaries and HTML include bound policy, a full physical-bound
table with active/inactive indicators, bound rejection counts by stage, and
separate numerical-invalid counts/reasons. Rejected events have +inf markers
at the top of the trial plot, not finite replacement objective values.

Configuration and checkpoints serialize the bound table. Resume reconstructs
it with `RefinementBounds.from_dict(...)` and requires exactly the same policy
and bound configuration. Policy changes, reordered/different bound tables or
changed bounds are rejected; unchanged bounds can accompany increased external
evaluation/time budgets. Inactive bounds are preserved. Old checkpoints lacking
these fields mean no bounds; they cannot silently resume with bounds enabled.
The original stage-boundary restart limitations still apply.

Initial and retained-state verifications must satisfy bounds before simulation.
An out-of-bound starting/retained model is a configuration failure, not a viable
optimization seed. Source mapping errors, unsupported physics, bad experimental
variance and proposal-generation failures outside evaluation are not converted
to bound violations. A proposal inside bounds can still fail numerically; its
separate `invalid_trial_policy` controls that outcome.

## Reproduction and timing

`benchmarks/validate_refinement_bounds.py` reads the current reference inputs,
restores historical depth 840 A **in memory only**, and reconstructs all eighteen
proposal values using seed 2026 and source driver stage seed
11945904881327827335. The current file's depth 1000 A is preserved. It checks the
same 27925.729123683803 A mapped proposal diagnosed previously.

Recorded output:
`tmp/refinement_bounds/20260917_074026_900583/summary.json`.

| Evaluation | Time | Result |
|---|---:|---|
| Inclusive depth 600..1200 A bound | 0.0008605 s | Bound rejection; zero preparation/forward calls, no warnings |
| No hard bounds; numerical rejection enabled | 0.8998835 s | Propagation multiply overflow |

The single-call measured difference was 0.8990230 s. This is a diagnostic
comparison, not a repeated throughput benchmark. Plot/report generation is
excluded. The same real eighteen-parameter ANNEAL initialization plus two
proposals took 12.3604531 s: three evaluations, one forward preparation, two
bound rejections. The second proposal was reached normally. Retained vector
was exactly unchanged and objective remained 1.4154022826663093. No long fit or
scientific improvement is claimed. The script saves both raw and mapped vectors.

## Validation recorded 2026-09-17

All validation used exactly
`C:\Users\peng863\AppData\Local\miniconda3\envs\qcbed-python\python.exe`.

- New hard-bound/rule tests plus the checkpoint retry: **36 passed in 5.28 s**
  (35 new tests, one existing checkpoint test).
- Extended focused policy/driver/four-optimizer suite: **204 passed, 1 failed**
  in 23.03 s. The failure was intermittent Windows `WinError 5` replacing
  `checkpoint_latest.json` in `test_snapshot_cadence`; that test passed in the
  subsequent focused retry. An earlier run before adding rule expansion had
  198 focused passes. Sandboxed temporary-directory access also required running
  tests with the existing validation permissions; no filesystem failures were
  hidden as bound or numerical rejections.
- Final complete suite: **1728 passed, 5 failed in 103.18 s**. All five failures
  are existing assertions expecting the old 840 A depth while the user-edited
  FOR007.DAT now contains 1000 A:
  `test_reference_exact_layout[FOR007.DAT]`,
  `test_reference_problem_short_trial_supported_path`,
  `test_reference_input_preparation_and_zero_tilt`,
  `test_reference_thickness_only_layout_and_fixed_physics`, and
  `test_real_al_subset`. Neither the reference file nor those tests was changed
  to force a green suite. The current reference-file SHA256 remained
  `5bd48ba770f8390bdda108296658717877b71ec1360e3631d2eb680a1091d3fb`
  throughout this task.
- Both controlled-subset and all-parameter dry-run behavior are tested, including
  unbounded display, all eight reflection rules, rule precedence, sign-folded
  physical starts, zero fallback, immutable configuration, retained-state safety,
  all four optimizers/hybrids, and saved-rule/resolved-table resume consistency.
- `git diff --check` and separate new-file whitespace checks passed. Nothing
  was staged or committed. No protected source/reference file was modified.
