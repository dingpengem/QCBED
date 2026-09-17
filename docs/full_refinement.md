# Durable iterative refinement

`examples/run_full_refinement.py` is an editable local runner for the existing
source-backed driver. It adds persistence, budgets and diagnostics, not a new
optimizer or scientific convergence criterion. Use exactly:

```powershell
Set-Location C:\data\QCBED
& "C:\Users\peng863\AppData\Local\miniconda3\envs\qcbed-python\python.exe" examples/run_full_refinement.py --quick
& "C:\Users\peng863\AppData\Local\miniconda3\envs\qcbed-python\python.exe" examples/run_full_refinement.py --dry-run
# Edit the settings near the top of the script before a full run:
& "C:\Users\peng863\AppData\Local\miniconda3\envs\qcbed-python\python.exe" examples/run_full_refinement.py
& "C:\Users\peng863\AppData\Local\miniconda3\envs\qcbed-python\python.exe" examples/run_full_refinement.py --resume "C:\data\QCBED\tmp\full_refinement\YOUR_RUN"
```

Outputs use unique timestamped directories beneath `tmp/full_refinement`.
Dry-run parses the full problem and prints settings and semantic parameter
labels without evaluating an objective or creating a run directory. Matplotlib
is required for this runner's mandatory initial/final previews and report.

## Configuration and defaults

The example reads the physical/refinement problem from `FOR007.DAT`, and
experimental intensity, variance and weight arrays from `in.dat`. The latter
does not define the whole refinement problem. Input paths, file SHA256 hashes,
the parsed problem fingerprint and layout are saved and checked on resume.

Default example settings:

- Objective: `masked`; execution: NumPy `batched`, maximum batch size 64.
- Mode: `ANSI`, the existing ANNEAL then AMOEBA sequence, one loop, seed 2026.
- ANNEAL/AMOEBA/SIMPARA controls: 20; PARAMIN cycle control: 3.
- Snapshot request every 25 objective evaluations, checkpoint every 25.
- Cumulative evaluation limit 150; cumulative wall-time limit 1800 seconds.
- Invalid proposed numerical states: explicit `invalid_trial_policy="reject"`.
- The example selects `REFINEMENT_CATEGORIES=("depth", "scan_coordinate")`;
  `None` preserves all parsed active categories.
- The example enables user-selected hard bounds: depth 600..1200 angstrom and
  scan half-width `max(abs(start)*.20, .05)`. Category rules automatically cover
  absorption, Debye-Waller and reflection multipliers when activated; specific
  targets override category rules. Dry-run shows every resolved bound. Package configuration defaults to
  no bounds. See [refinement bounds](refinement_bounds.md) for mapped-value
  checks, inactive targets, separate rejection counts and exact resume policy.

`ANNE`, `SIMP`, `PARA`, `SPAR`, `PAAN`, `ANPA`, `SIPA`, `PASI`, `SIAN`, `ANSI`
are the existing driver modes. `SPAR` preserves the input's SIMPARA strategy as
an explicit option. The example exposes AMOEBA and SIMPARA limit names for
clarity, but they map to the **same legacy MAXFSIMP** field as ANNEAL. Settings
for selected hybrid stages must agree; conflicting limits raise. No new
stage-specific budget or optimizer sequence is invented. Source soft/inclusive
limits remain soft/inclusive, not strict universal evaluation caps.

Quick mode overrides the example to ANNEAL, one loop, MAXFSIMP=0 (one source
inclusive proposal), temperatures .001/.0001, cooling .5, seed 2026, masked
objective, batch 64, no intermediate snapshots, checkpoint every evaluation,
20 total evaluations and 600 seconds. With no category subset it retains all 18 active parameters,
1521 tilts, nine beams and original reference physics. It is pipeline validation,
not a scientifically converged fit. A normal/full run is never launched as part
of pytest or validation automation.

Masked objective is an explicit non-legacy numerical policy in which zero-weight
points are excluded before variance division and before YNORM computation.
The same retained points enter YNORM and the weighted objective. Variances and
legacy arithmetic are not changed. `legacy` remains selectable, but this real
dataset's zero-weight/zero-variance entries yield an undefined legacy objective;
the runner fails explicitly rather than replacing it.

## Reusable APIs

```python
from dataclasses import replace
from qcbed import ForwardExecutionConfig, FullRefinementConfig, run_full_refinement

driver_config = replace(driver_config,
    execution=ForwardExecutionConfig("batched", batch_size=64))
controls = FullRefinementConfig(
    driver_config, objective_mode="masked", snapshot_every=100,
    checkpoint_every=25, max_total_evaluations=1000, max_wall_time=3600.,
)
run = run_full_refinement(
    state, experimental, controls,
    output_directory="tmp/a_new_run",
    input_paths=("reference_data/FOR007.DAT", "reference_data/in.dat"),
)
```

The immutable full-run config honors the supplied driver's execution selection;
the example explicitly recommends batching. `run_full_refinement(..., resume=True)`
uses the same baseline problem/config and directory. `load_refinement_checkpoint`
loads versioned JSON with basic validation; the runner additionally validates
input, configuration and layout identity. Existing APIs keep their previous
defaults. New optional driver hooks expose pre-stage retained boundaries and a
before-evaluation budget check. Optimizer wrappers only forward the latter;
their algorithms are unchanged.

## Retention, progress and snapshots

Before the first optimizer proposal, the runner performs and saves an independent
initial evaluation, active vector, physical state, raw intensity, objective,
YNORM, preview, configuration and checkpoint. The driver's existing pre-stage
evaluation still runs; it is not silently replaced by the initial diagnostic.

`progress.jsonl` is append-only, flushed and fsynced per record. Rows contain
global evaluation index, loop, optimizer stage, stage-local evaluation index,
classification, objective, YNORM, source-ordered active values, elapsed seconds
and restart count. Trial vectors are **mapped physical active values**, even
when AMOEBA internally evaluates simplex coordinates. Semantic labels/targets
are stored once in `run_config.json` and checkpoints. Stage-local indices count
pre-stage evaluations and reset when restarting a stage.

Classifications distinguish `initial_retained`, `trial`, `invalid_trial`,
`retained_boundary`, `stage_retained`, `snapshot_verification`, and final
`verification`. Invalid trials use `accepted=false` and retain their original
trial coordinate mode/vector alongside mapped physical values where available.
For finite trials `accepted` is null: callbacks expose evaluations, not source-backed
proposal acceptance. A trial is never relabeled as retained just because its
objective is lower. Retained plots show only confirmed driver states; rejected
or acceptance-unknown proposals are not drawn as accepted trajectories.

No full CBED arrays are kept at every trial. Small vector/scalar YNORM lookup
metadata are kept for the current stage only. Initial and final snapshots always
exist after successful evaluation/finalization. A snapshot cadence hit requests
the **next completed retained stage boundary**, where a fresh diagnostic
evaluation creates the snapshot. Multiple cadence hits in a long stage coalesce;
there is no invented mid-stage retained state. Snapshot metadata records global
evaluation, retained boundary index, stage/loop, vector, objective and YNORM.

## Checkpoints and restart

Checkpoints are `qcbed-retained-stage`, version 1, JSON only (no pickle).
Numbered files include evaluation, completed-stage count, restart count and
status; `checkpoint_latest.json` is published last by atomic same-directory
replacement. Temporary files are flushed/fsynced before replacement. NPZ
snapshots are similarly replaced atomically. Checkpoints include:

- Retained raw active vector, objective/YNORM, confirmed-boundary evaluation.
- Separate next-stage packed vector, loop and stage-position cursor.
- Completed stage summaries, evaluation counters, cumulative elapsed time.
- Original seed, optimizer/execution/objective settings and active layout.
- Input paths/hashes and parsed-problem digest, restart count and policy.
- `rng_state=null`, because exact internal continuation is not available.

Only the initial, pre-stage and completed-stage retained states are restartable.
After interruption within a stage, restart begins that stage from its last
confirmed boundary using the original derived stage seed. In-stage accepted
improvements may be lost. This is **not exact stochastic or simplex continuation**.
At completed-stage boundaries, the next stage starts with the source-required
UPDPARS/packing transition. The original baseline remains unchanged: rebuilding
the baseline from the retained model would incorrectly change AMOEBA scaling.
Original loop/stage indices are preserved for seed derivation.
Resume with the same code and numerical environment; checkpoints are not a
cross-version optimizer-state interchange format.

Budgets are cumulative across sessions. Increase them explicitly when resuming
a limited run, for example:

```powershell
& "C:\Users\peng863\AppData\Local\miniconda3\envs\qcbed-python\python.exe" examples/run_full_refinement.py --resume "C:\data\QCBED\tmp\full_refinement\YOUR_RUN" --max-total-evaluations 2000 --max-wall-time 7200
```

One evaluation slot is reserved for mandatory fresh final verification. Counts
include initial, driver, periodic verification and final evaluations; attempts
interrupted in flight also consume a slot. Diagnostic evaluations do not count
as optimizer decisions. Wall time is cooperative, checked before evaluations;
an in-flight numerical call and safe finalization/report generation can exceed
it. It is not a hard process watchdog.

Ctrl+C checkpoints the latest confirmed boundary and attempts final verification
and output. A second Ctrl+C during finalization saves a checkpoint and exits;
resume regenerates outputs. Interruption before a successful initial evaluation
cannot produce a valid CBED snapshot; the checkpoint explicitly retains an
unevaluated initial vector. Other errors save a failed checkpoint/summary and
propagate, except numerical-invalid proposals under the explicit rejection
policy: these receive a rejecting positive-infinite score and do not alter the
retained state. Initial and retained-state verification must remain valid.
See [invalid-trial policy](invalid_trial_policy.md) for classification, logging,
subset selection and the reproduced seed-2026 failure. Physics is never repaired.

A `.run.lock` prevents concurrent writers. After a hard process kill, verify the
old process is gone before manually removing a stale lock. Complete progress
records newer than a periodic checkpoint advance counters on resume while its
safe retained boundary is preserved. A torn/corrupt JSONL record is rejected for
explicit recovery; it is never silently truncated. Exact accounting of an
unlogged in-flight call after a hard kill cannot be guaranteed. These are local
filesystem durability measures, not a distributed transaction system.

## Outputs and interpretation

```text
run_config.json                 initial settings, labels, input hashes
summary.json                    status, objectives, improvements, counts, timing
progress.jsonl                  append-only scalar/vector events
checkpoint_latest.json
checkpoint_<evaluation>_<stages>_<restart>_<status>.json
initial_parameters.json        mapped physical model
final_parameters.json
initial_active_vector.npy
final_active_vector.npy
initial_simulated_intensity.npy raw |amplitude|² before GRADIENT
final_simulated_intensity.npy
final_scaled_intensity.npy      post-GRADIENT, YNORM-scaled
initial_snapshot.npz/.json
final_snapshot.npz/.json
retained_*_snapshot.npz/.json    optional boundary diagnostics
initial_preview/index.html     one image per beam
final_preview/index.html
retained_*_preview/             optional boundary galleries
final_out.dat
convergence.png
trial_progress.png              trials (acceptance unknown) vs retained markers
parameter_trajectories.png
report.html
```

The existing plotter uses a four-panel layout: experiment, initial simulation,
current retained simulation, signed residual. Every panel uses the same weight
mask, measured beam H/K labels, and aligned acquisition indices. Data panels
share intensity limits; residual limits are symmetric around zero. No independent
normalization hides mismatch. Grey means weight-zero, not zero intensity. The
real reference's GRADIENT=1 means the fit view contains processed objective-space
values; raw intensities are saved separately. Each gallery shows all nine beams
as separate readable images. Previews and plotting occur outside numerical trial
calls; total runner wall time does include diagnostics and report generation.
Keep the output directory together: the main report embeds its standard CBED
plots and also links the preview galleries and trial-progress PNG.

`final_out.dat` reuses `write_legacy_out_dat`: headerless, beam-major, scaled
calculated values after optional GRADIENT, portable `.17E` output. This is numeric
format compatibility, not compiler-byte equivalence. Successful finalization
freshly evaluates the retained raw vector with the selected backend and demands
exact equality of saved/recomputed objective and YNORM before publishing final
outputs. Reports include per-beam diagnostics, initial/final parameters and
changes, stage/loop boundaries, improvements and restart limitations. An objective
increase or no improvement is reported honestly; completing configured stages
does not imply scientific convergence.

## Source evidence and scope

The current readable `private predecessor implementation` implements ANSI, including ANNEAL at
1441 and AMOEBA at 1510–1511, restoration of PARMIN and ARRANGEPARS between
stages. `private predecessor implementation` supplies ANNE/SPAR paths; `private predecessor implementation` places
GRADIENT before objective/output and writes scaled beam-major values. These
confirm reuse of the existing driver and output APIs. Checkpoint files, budgets,
atomic writes and restart policy are explicit Python orchestration features,
not reconstructed historical behavior. private predecessor implementation is only a related-version
declaration reference. No new executable-Fortran equivalence is claimed.

Parallel refinement, new FFT backends, GPU, new search algorithms, spherical
waves, covariance estimation and automatic scientific-convergence claims are
outside this runner. Existing unsupported trial-state errors remain in force.


## Validation recorded 2026-09-17

Using only `C:\Users\peng863\AppData\Local\miniconda3\envs\qcbed-python\python.exe`:

- Focused runner/driver/backend suite: **111 passed in 16.86 s**.
- Complete pytest suite: **1669 passed in 55.07 s**.
- Quick and normal `--dry-run` both parsed the real 18-parameter problem without
  evaluating or creating a run directory.
- The final quick-run artifacts are in
  `tmp/full_refinement/20260917_001808_806073/`; open `report.html` there.
- Quick controls: ANNE, LOOPS=1, MAXFSIMP=0, seed 2026, temperatures .001/.0001,
  cooling .5, masked objective, batched NumPy 64. All reference physics and scan
  points retained. No tuning was performed to force improvement.
- **5 total completed evaluations**, comprising one initial diagnostic, three
  driver calls and one fresh final verification. Total runtime **57.25211 s**,
  including persistence, previews and report generation.
- Initial and final objective: **1.4154022826663093**; absolute/relative
  improvement: **0 / 0**. Initial and final YNORM: **4723684.377058742**.
- All 18 retained parameter changes were zero. The one perturbed proposal had
  objective 1.4157744956478286; the final driver-retained state was unchanged.
  No scientific convergence is claimed.
- Final objective/YNORM verification was exact. `final_out.dat` parsed back
  exactly equal to `final_scaled_intensity.npy`. Both initial/final galleries
  contain nine labeled beam PNGs. The trial-versus-retained plot was inspected.
- Tests cover exact driver agreement for all four standalone optimizers,
  deterministic repeats, hybrid/loop-boundary restart, mid-stage and pre-initial
  Ctrl+C, cumulative evaluation limits, corrupted checkpoint versions, input
  fingerprint mismatch, checkpoint cadence, snapshots, final verification
  failures, output round-trip, report generation and mapped initial-scan geometry
  on resume.
- `git diff --check` and separate new-file whitespace checks passed. Protected
  source/reference files are unchanged. Nothing was staged or committed.

The earlier smoke run in `tmp/full_refinement/20260917_000922_224086/` used
identical numerical settings and returned the same objective/YNORM and five
calls (60.72 s); it preceded final preview metadata/progress-report corrections.
Use the newer directory above. Only short quick runs were performed, not a full
production fit.
