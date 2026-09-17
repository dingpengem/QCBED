# Running local optimization and inspecting patterns

The runner calls the existing `run_legacy_refinement` and optimizer algorithms.
It adds explicit objective selection and observation, not a second optimizer.
This is not claimed to reproduce a historical Fortran executable.

## Run the short reference example

```powershell
Set-Location C:\data\QCBED
$env:MPLBACKEND = "Agg"
& "C:\Users\peng863\AppData\Local\miniconda3\envs\qcbed-python\python.exe" examples/run_local_optimization.py
```

The script prints its unique output directory beneath `tmp/local_optimization/`.
Open `report.html` there when finished. While it runs, inspect `progress.jsonl`
and the checkpoint `*_beam_*.png` files. Initial, selected trial and final maps
use identical weight masks across experimental/simulated/residual panels.
PNG files update at checkpoints; the complete HTML report is written at the end.
Matplotlib is provided by the package's optional `visualization` extra.

Edit the settings in the example for subsequent runs. It initially uses the
unaltered 18-parameter, 1521-tilt, nine-beam reference problem, explicit masked
objective, seed 2026, ANNEAL only, LOOPS=1, MAXFSIMP=0, maximum temperature .001,
minimum temperature .0001 and cooling .5. MAXFSIMP=0 permits **one proposal**
because the existing ANNEAL source loop is inclusive. This is deliberately short
validation, not an adequate production fit. Other optimizer limits also retain
their source soft/inclusive behavior; limits are not universal evaluation caps.
Do not infer improvement from successful execution: inspect the measured scores.

## API

```python
from qcbed import run_local_optimization

run = run_local_optimization(
    state, experimental, config,
    objective_mode="masked",       # default is "legacy"
    snapshot_every=10,              # None disables intermediate snapshots
    record_progress=True,
    output_directory="tmp/my_new_run",  # must not exist
    reference_input_name="my FOR007.DAT",
)
```

Pass parsed state/experimental data and `LegacyRefinementConfig` as in the example.
Omit `output_directory` for an in-memory run without Matplotlib/file output.
Protected source/reference directories are rejected as output destinations.
Directories must be new so previous runs are not silently overwritten. Failures
leave already-written checkpoint files for inspection; no resume mechanism is
provided. Invalid trial states still raise; no objective penalties are invented.

Masked mode is **an explicit non-legacy numerical policy for datasets whose
zero-weight masked points contain zero variance**. It excludes weight==0 before
variance division in both data-dependent YNORM and weighted root error. It keeps
fractional weights, rejects invalid retained variances, and never adds epsilon.
The reference case includes 7549 points, masks 6140, and all 117 zero variances
are masked. Legacy mode retains its NaN result for this case and cannot refine.

## What progress means

`LocalOptimizationRunResult` contains validation data, immutable progress and
snapshots, objective mode and measured elapsed time. `final_physical_parameters`
exposes the driver-returned immutable physical state.

Each `RefinementProgressRecord` contains evaluation index, stage index/name,
loop index, objective, YNORM and a read-only physical parameter vector:

* `kind="evaluation"` is the **evaluated trial**, possibly rejected. These are
  actual objective evaluations, including pre-stage, initialization and reset
  calls. They are not universally equivalent to optimizer iterations.
* `kind="retained"` is the exact driver-retained vector at a stage boundary.
  YNORM is matched to the corresponding evaluated physical state, without an
  extra objective call. No inference of acceptance for each trial is made.

The convergence plot distinguishes trial scores from retained stage endpoints
and labels optimizer/loop boundaries. Parameter trajectories use retained stage
states, not an invented dense sequence of accepted iterates. Full evaluation
vectors are available in JSONL. Stage selection retains source behavior; it need
not be monotonically improving across every source reset or stage.

Observers run after an existing trial computation and after stage retention.
They draw no random numbers and make no objective calls. No optimizer branch,
mapping equation, RNG order, or forward physics changes. Tests compare complete
retained results with recording enabled/disabled for every optimizer and a
multi-stage sequence. A final diagnostic evaluation runs **after** the driver
on both recorded and unrecorded paths and checks the retained score. It is
excluded from driver evaluation counts.

## Snapshot and file semantics

For positive integer N, `snapshot_every=N` selects global driver evaluation
indices divisible by N. Evaluation 1 is always the initial snapshot and is not
duplicated when N=1. The final retained state is always stored independently,
even if it equals an earlier state. Intermediate snapshots are explicitly trial
states, not necessarily improvements. None disables only intermediate snapshots.

No full intensity array is retained for an unselected evaluation. Small
scalar/vector progress records are O(evaluations times active parameters);
selected snapshots retain the raw/objective/scaled arrays needed by the existing
validation representation. `record_progress=False` disables detailed persisted
history, not mandatory initial/final snapshots or the existing driver history.

Output files:

* `progress.jsonl`: appended evaluated/retained scalar and parameter records.
* `initial.npz`, `evaluation_000NNN.npz`, `final.npz`: selected map arrays,
  tilts, beam indices, parameters, scores, YNORM and objective mode.
* Matching per-beam PNGs: experimental, initial, checkpoint/final, and residual.
* `report.html`: self-contained masked/legacy-labeled report, counts, final
  comparisons, initial residuals, optional trial gallery, convergence, retained
  parameter trajectories, and configuration/parameter tables.
* `summary.json`: initial/final scores, YNORM, parameter changes, seed, limits,
  evaluation count and measured wall time including final report generation.
* `final_physical_parameters.json`: complete driver-returned physical state.

Individual figures share a scale across their data panels. Color limits can
change between checkpoint figures; compare the displayed numbers/colorbars,
not apparent brightness alone. Initial and final data are compared on a common
scale within the final report figure. No display normalization changes values.

## Measured reference validation (2026-09-16)

One run was performed with the short configuration above, without retuning:

| Quantity | Result |
| --- | --- |
| Mode / stages / seed | masked / ANNEAL / 2026 |
| LOOPS / MAXFSIMP / MAXFPARA | 1 / 0 / 1 (MAXFPARA unused) |
| Initial objective | 1.4154022826662773 |
| Proposed trial objective | 1.4157557004972587 (worse; not retained) |
| Final objective | 1.4154022826662773 |
| Objective improvement | 0 |
| Initial/final YNORM | 4723684.377058765 |
| Proposed trial YNORM | 4726252.953286001 |
| Retained parameter changes | all 18 exactly zero |
| Driver evaluations | 3: pre-stage, ANNEAL initialization, one proposal |
| Additional diagnostic evaluations | 1 final-state verification |
| Snapshot cadence / labels | 3 / initial, evaluation_000003, final |
| Wall time including original report generation | 695.7369672 seconds |
| Baseline mutation check | passed |

Artifacts are in `tmp/local_optimization/20260916_202242_824447/`; open
`report.html`. The rejected proposal is included in the gallery rather than
hidden. Checkpoint objectives and YNORM were independently recomputed exactly
from saved objective-space arrays. Passing those same arrays through the legacy
objective still produces NaN. Rendering labels were refreshed from saved
checkpoints after final review; no numerical evaluations or optimizer reruns
were used for that refresh, and its rendering time is not included above.

Validation: 195 focused tests passed; the full suite passed 1558 tests with
the exact qcbed-python interpreter. Legacy objective/scaling files and protected
source/reference files were unchanged. This short run validates execution and
recording, not optimizer efficacy or a historical production refinement.
