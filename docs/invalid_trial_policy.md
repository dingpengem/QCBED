# Numerical-invalid refinement proposals

**Python production-run numerical safety policy; not proven historical Fortran behavior.**

`evaluate_refinement_trial`, the four optimizer wrappers, and
`run_legacy_refinement` accept `invalid_trial_policy="raise"` (default) or
`"reject"`. Raise preserves the original low-level exceptions and arithmetic.
`FullRefinementConfig` and the editable production example explicitly use
`"reject"` by default. Neither objective formula changes.

Reject runs the evaluation with NumPy overflow, invalid-operation and division
errors enabled. Numerical exceptions and existing nonfinite-array validation
errors produce an immutable `InvalidRefinementTrial` with objective `+inf`,
exception type, reason, original trial vector/coordinate mode, and mapped
physical state when available. Nonfinite YNORM/objective results are rejected
too. There are no fabricated amplitudes or YNORM. No parameter clipping, epsilon,
variance repair, array replacement, or weakened FFT checks are introduced.

Structural input errors, unsupported physics, bad retained experimental
variances, unexpected exceptions, callbacks and budget interruptions still
propagate. Classification of existing ValueErrors uses their finite-validation
messages, so new kernel validation messages must be covered when introduced.
Initial/retained-state verification must succeed: the policy rejects proposals,
not an unevaluable starting model. Arithmetic errors in proposal generation
outside the evaluator are not converted to trial scores.

ANNEAL, AMOEBA, PARAMIN and SIMPARA permit a positive-infinite proposal score
only after a valid anchor. Such a score cannot update the retained vector or
objective. Existing comparison rules and parabolic fallbacks remain unchanged.
SIMPARA skips callback parameter writebacks for invalid evaluations. Required
retained-state reevaluations fail explicitly if they become invalid.

## Progress, reports and subset selection

An invalid proposal consumes an attempted evaluation, but not a successful
evaluation or retained-state update. Its `progress.jsonl` classification is
`invalid_trial`, with global/stage indices, vector, mapped physical parameters,
exception type/reason, `accepted=false`, `retained_state_unchanged=true`, and
`ynorm=null`. Strict JSON records the objective as the string `"inf"`; the
optimizer receives numeric positive infinity. No traceback is stored per row.
A concise console warning identifies the stage, evaluation and error.
Checkpoints, summary JSON and HTML report include invalid counts, reasons and
counts by stage. Invalid trials use separate markers at the top of the progress
plot, not invented finite objective values. Resume restores rejection accounting.

`select_refinement_categories(state, config, ("depth", "scan_coordinate"))`
returns immutable state/config copies: excluded masks are disabled, physical
starting values are preserved, and D1/D2 ranges and PLAMDA are selected in the
original source parameter order. This is a user-selected refinement subset,
not a change to input-file semantics. Set `REFINEMENT_CATEGORIES` in
`examples/run_full_refinement.py`; `None` keeps all parsed active parameters.
The reference subset contains seven rather than eighteen parameters. Resume
reconstructs the same subset from the saved category list.

## Reproduction of the seed-2026 failure

`benchmarks/diagnose_invalid_trial.py` reconstructs the first ANNEAL proposal
using the existing source-derived driver seed and RNG, then temporarily
instruments the existing forward path. It does not replace physics kernels.
The stage seed is `11945904881327827335`. Controls are masked objective,
NumPy batched/64, all eighteen parameters, maximum temperature 100, minimum 1,
cooling .9, and the original reference data. Full raw diagnostics are saved in
`tmp/invalid_trial_audit/diagnosis.json`.

The raw proposal is below. The existing UPDPARS/packing mapping folds the first
four entries to positive values before evaluation; all other entries retain
their signs. Those folds predate this policy.

| Parameter (source order) | Starting value | Raw proposal |
|---|---:|---:|
| Depth, region 0 | 840 | -27925.729123683803 |
| Absorption (0,1) | 1 | -17.54921103197698 |
| Absorption (0,2) | 1 | -9.32128452070585 |
| Debye-Waller (0,0) | .86 | -2.2999633002671773 |
| Scan (0,0) | -.508156 | 27.851284058182173 |
| Scan (0,1) | .407641 | -19.76376396906165 |
| Scan (1,0) | .407641 | -55.86352057579832 |
| Scan (1,1) | .508156 | 15.569703843908052 |
| Scan (2,0) | -.407641 | -18.507550624072188 |
| Scan (2,1) | -.508156 | -30.616442883553447 |
| Elastic reflection 0 | 1 | -77.95416700889349 |
| Elastic reflection 1 | 1 | 117.9261054735951 |
| Elastic reflection 2 | 1 | -48.31972193033225 |
| Elastic reflection 3 | 1 | -8.70507505824593 |
| Elastic reflection 4 | 1 | -190.43547276670452 |
| Elastic reflection 5 | 1 | -42.89339826093701 |
| Elastic reflection 6 | 1 | -114.88858466286261 |
| Elastic reflection 7 | 1 | 54.919273665329044 |

Mapped parameters and both structure-factor arrays were finite. Their maximum
magnitudes were 34.39632671380606 (elastic, before reflection scaling) and
11.234432168737746 (absorptive). The base grating magnitude ranged from
.8995501944072085 to 1.0000012952754551; the reflection-adjusted grating remained
finite but ranged from .18121905668627927 to **5.51555028243405**.

The first nonfinite array was the wave **after grating multiplication at the
738th slice occurrence**, zero-based occurrence 737, region 0, slice type 0.
The preceding finite wave magnitude ranged from 8.594253637231415e69 to
6.170327186333842e307. Two complex entries became nonfinite in the first
64-tilt block. Prior FFT outputs were finite; the next FFT input check raised
`real waves must be finite`. Calculated intensities and YNORM were not reached.
Reject mode instead catches the arithmetic overflow at the multiplication.

The mapped depth expands the sequence from 415 to 13,794 occurrences. Restoring
one category at a time, retaining the others at the proposed values, gave:

| Restored category | Existing forward path, first 64 real tilts |
|---|---|
| Reflection multipliers | All 13,794 occurrences finite; maximum intensity 2.9782855995816603e-30 |
| Absorption | Nonfinite reciprocal wave |
| Debye-Waller | Nonfinite reciprocal wave |
| Depth | Wave/intensity path reaches nonfinite intensity within the 415-occurrence calculation |
| Scan coordinates | Nonfinite grating product at occurrence 458 |

This implicates the large signed reflection multipliers and the resulting
amplifying grating, rather than absorption magnitude alone. It does not establish
universal causality or validate the restored trial over all 1,521 tilts. No
objective was computed for these isolation cases.

## Bounded production validation

The original seed/configuration with rejection enabled and an external cap of
five total evaluations produced
`tmp/invalid_trial_audit/continuation_20260917_004442/report.html`.
It rejected the reproduced fourth evaluation, then completed a fresh fifth
evaluation of the retained valid model. There were five attempts, four successful
evaluations, one invalid ANNEAL trial, and no completed optimizer stage.
Initial/final objective was exactly **1.4154022826663093**, YNORM exactly
**4723684.377058742**, and all eighteen retained parameter changes were zero.
Total runner time, including diagnostics, was **80.5886 seconds**.

An earlier exploratory cap-seven run was stopped while the subsequent costly
proposal was in flight; its unfinished directory is
`tmp/invalid_trial_audit/continuation_20260917_004001`. Use the completed bounded
report above. The real validation proves rejection and valid finalization;
focused tests prove subsequent proposal and hybrid-stage continuation. No long
fit or scientific convergence is claimed. Budgets remain cooperative: an
expensive finite proposal can exceed the wall-time limit before returning.

Validation used only
`C:\Users\peng863\AppData\Local\miniconda3\envs\qcbed-python\python.exe`:

- Focused evaluator/runner/four-optimizer tests: **143 passed in 13.38 s**.
- Complete suite: **1698 passed in 67.39 s**.
- Example `--dry-run` passed with the user's existing 25/150/1800 settings.
- Final bounded-run vector/objective/YNORM and logged proposal were checked
  against saved diagnostics. The invalid-trial plot was visually inspected;
  its x limits explicitly include rejected evaluations.
- A focused attempt encountered Windows `WinError 5` during checkpoint atomic
  replacement; the identical rerun passed. This intermittent filesystem error
  is not converted into a numerical rejection or silently ignored.
- `git diff --check` and new-file whitespace checks passed. No protected source
  or reference file was changed; nothing was staged or committed.
