# Optional configuration editor and runner (Stages 6B, 6C and 6D)

> Historical milestone/audit record. For current controls and the A/B/C/D workflow,
> use the [User Manual](user_manual.md); dated counts and replay inputs below retain
> their original provenance.

Install the optional frontend into your development environment and launch it:

```sh
python -m pip install -e '.[gui]'
python -m qcbed.gui examples/qcbed_config.yml
```

Launching without a filename opens an empty window with an Open action. The
numerical package, CLI, and editor document model do not import PySide6. Run and Resume execute canonical snapshots in a separate worker process.

## Workflow

Open loads YAML through `load_yaml_config`. The left navigation selects editors
for project title, microscope energy (kV) and integer zone axis, independent
region depths (angstrom), scan dimensions/bases/corners, experimental path,
refinement bindings, execution backend/batch size, and output paths. All other
canonical fields are visible under **All fields (read-only)**. Structural records and bounds are edited through **Structure / Bounds**. Optimizer
settings and category filters remain read-only.

Numeric fields use text inputs, with comma- or whitespace-separated vector
components. Values display with round-trip precision. Typing changes only a
draft. **Apply edits** parses and validates the complete draft, then replaces the
immutable `QcbedConfig` snapshot using the public Stage 6A editing APIs. Invalid
field or compile edits remain visible and do not change the committed snapshot.
**Revert drafts** restores that snapshot, not the last saved file. No values are
wrapped, clipped, or rounded by the editor. Batched execution requires an explicit
positive batch size; switching to serial requires clearing the batch size.

The persistent summary compiles the committed snapshot and shows regions, stored
slices, tilts, beams, active parameter count, execution preflight status, and resolved
hard bounds. **Validate** also checks experimental files and output readiness.
Configuration validity and execution/file readiness are separate: valid
configurations with missing data can be edited and saved. Preflight does not run
numerical trials. The issue panel shows structured validation phases and locations;
selecting an issue navigates to its section and field where a matching editor exists.

## Refinement bindings

Targets come from `available_refinement_targets`, including unsupported targets
and their explanations. Region/species/reflection stable IDs distinguish physical
instances, including duplicate reflections. Sorting the table changes presentation
only. The editor never reorders physical records.

**Binding** selects a target binding. **Active** additionally reflects the existing
category filter, which is shown above the table and is read-only in this milestone.
For example, selecting an elastic multiplier in a depth-plus-scan configuration
does not silently remove that filter. New bindings require explicit D1, D2, and
PLAMDA values. These retain legacy range-factor/step semantics and signs; they are
not hard bounds. Linked or unsupported targets remain identified by the catalogue.

## Saving and paths

Save atomically replaces the current YAML; Save As atomically creates a new file,
with confirmation before replacing an existing destination. Both delegate to
`save_yaml_config`. Unsaved edits prompt before closing/opening another document.
Invalid drafts prevent saving. YAML comments and formatting are not preserved.

Relative path edits resolve against the current document directory, never the
incidental process working directory. Save As preserves resolved file identities;
subsequent relative edits use the new document directory. Stable IDs, full numeric
precision, scientific ordering, and all unedited canonical fields are preserved.
When title or experimental/output filename changes, only its associated legacy raw
spelling is cleared because it describes the old value; unrelated export metadata
is retained. Canonical saves do not overwrite protected input/reference files.

## Implementation and tests

`qcbed.gui.document.EditorDocument` manages drafts/committed snapshots without Qt.
`qcbed.gui.window.ConfigWindow` presents those APIs with Qt Widgets. Neither owns
scientific validation, target enumeration, compilation, or YAML serialization rules.
The frontend delegates execution and checkpoint semantics to the canonical runner.

```sh
python -m pytest tests/test_gui_document.py
QT_QPA_PLATFORM=offscreen python -m pytest tests/test_gui_qt.py
python -m pytest --ignore=tests/test_gui_qt.py
```

Qt tests skip when PySide6 is absent. They use QApplication directly, without a
pytest-qt dependency. Offscreen tests cover widget editing, issue navigation,
sorting, compile summaries, and canonical Save As/reload equivalence. Native desktop
file-dialog styling and window-manager behaviour still warrant manual platform checks.

The issues panel and summary distinguish **Configuration validity**, **Compile
validity**, **Execution readiness**, and **External file/output readiness**.
Existing run output paths are displayed as execution-readiness **Run blockers**,
not configuration errors: the snapshot remains structurally valid, compilable,
and saveable. Resume is suggested only when canonical compatibility validation
succeeds. Directory readiness distinguishes:

- Missing output directory: fresh Run eligible; Resume disabled.
- Existing directory with a compatible checkpoint: Run blocked; Resume eligible.
- Existing directory without a checkpoint: both blocked; choose a new directory
  or remove/rename the existing directory.
- Incompatible or corrupt checkpoint: both blocked; the compatibility or parsing
  reason is shown.

Other configuration/data readiness checks must also pass. Missing checkpoints
are checked before attempting resume validation, without creating any artifacts.
Underlying validation severity, overwrite protection, and checkpoint checks are unchanged.


## Run and Resume (minimal Stage 6D)

Apply edits and select **Validate** to enable eligible Run/Resume actions. Their
availability is checked separately: an existing output directory blocks a fresh
run, while Resume requires a compatible checkpoint and experimental data. Reasons
are shown in the run log and action tooltips. Any edit invalidates readiness until
Validate is used again. Launch repeats readiness checks; the canonical runner
remains authoritative about locks, output protection and checkpoint compatibility.

The GUI captures configuration and experimental data using `prepare_run_snapshot`
and launches `run_snapshot` in a spawned process, without a temporary YAML file.
The editor, Open/Save actions, and window closing are disabled while a worker is
active. Cancel requests safe completion without unlocking the editor early. This intentionally conservative first version permits no edits during a
run; later edits cannot mutate the captured snapshot. Only cooperative cancellation is provided (see Stage 6E below). Preflight/data loading is synchronous; numerical refinement runs
outside the Qt thread. Native/library process failures are isolated from Qt.

The persistent run panel shows Idle, Running, Completed, Failed, Cancelled, or Interrupted
(budget-limited runs also display Interrupted, with the original status/reason in
the log). It polls progress events and worker stdout/stderr, preserving the
runner's stage, loop, evaluation, objective, rejection and snapshot records as
reported. Elapsed time continues updating during long evaluations. Final status
includes objective, evaluation count, output directory and checkpoint/final-state
presence. Worker exceptions include tracebacks; abnormal exits show the exit code.
The bounded display retains the latest 3,000 log blocks. Canonical progress logs,
checkpoints, reports and final outputs remain the durable production artifacts.
The GUI does not add optimizer callbacks or reinterpret checkpoint identities.

Process tests compare actual spawned execution with direct canonical output and
cover snapshot independence, compatible/incompatible resume, worker failures,
progress delivery, and Qt timer responsiveness. Run these with
`python -m pytest tests/test_gui_run.py tests/test_gui_qt.py` (Qt tests are optional).

## Structural transactions and history (Stage 6C)

Choose **Structure / Bounds** after applying or reverting ordinary field drafts.
The structural window edits canonical objects directly. Select a record in its
navigation tree, or double-click a row in a collection table. Tables/tree sorting
only changes display order. **Scientific up/down** explicitly reorders an ID-bearing
collection; the stable IDs and associated refinement targets are retained.

The editor provides:

- Species labels and atomic numbers; add, duplicate, delete, scientific reorder.
- Regions with independent depths and symmetry; ordered stored slices, reflections,
  and species-associated Debye–Waller entries.
- Stored-slice thickness, fractional shifts, cell lengths (angstrom), gamma (degrees),
  record type, and ordered atoms with species, fractional coordinates and starting
  occupancy. Values are not folded, clipped or normalized.
- Ordered reflections with h/k and elastic/absorptive multipliers. Duplicates are
  retained as distinct stable-ID records.
- Species-associated absorption coefficients with independent/linked-C mode.
- Specific-target, category-absolute and category-relative bounds, and policy.
  The main compiled summary displays resolved values and the responsible rule,
  including unbounded targets and whether rejection policy is enabled.

**Add** opens a typed form with starting values for review, creates a fresh ID for
ID-bearing records, and appends in scientific order. New region defaults include
a displayed unit-cell slice and zero thermal entries for current species; review
and edit these before committing. **Duplicate** uses the public stable-ID cloning
helper, copying nested physical records with fresh IDs, but not copying bindings
or bounds automatically. Adding/duplicating a species requires explicitly adding
its thermal entries in every region and absorption entry if enabled. Missing
associations remain draft validation issues until supplied; they are never guessed
from another species. Blank optional bound endpoints mean no endpoint.

**Stage fields** updates only the transaction draft. The tree can hold incomplete
cross-object relationships while you make coordinated edits. For example, enabling
absorption requires adding its coefficient records before Commit. **Commit
transaction** uses the shared structured field/compile validator; only a valid
complete snapshot enters the document. Execution-readiness issues remain separate.
Cancel discards the draft. The available-target catalogue refreshes after valid
structural changes; existing bindings remain attached to their stable identities.
New bindings are enabled in the main Refinement panel after structural commit.

Deletion previews the removed subtree and the transitive external references:
atoms, thermal/absorption entries, refinement bindings, and specific bounds. **No**
cancels; **Yes** explicitly authorizes removing exactly the reported dependencies
in the draft. Alternatively cancel and reassign/remove references individually.
Deletion never silently retargets a reference. Linking absorption C to B asks for
explicit permission to disable an existing independent-C binding; declining leaves
the draft unchanged. Stored C remains intact and is read-only while linked.

**Undo/Redo** restores exact immutable snapshots: values, IDs, scientific ordering,
bindings, and bounds. One valid structural Commit is one history entry. Common
Apply edits also participates. Incomplete widget text and intermediate structural
staging do not enter history. Apply or revert ordinary drafts before Undo/Redo;
a new committed edit clears the redo branch. Save does not erase history; Open
starts a new document/history. All structural/history actions remain locked during
Run/Resume. The numerical worker's captured configuration/data cannot be mutated.

## Cooperative Cancel (Stage 6E)

**Cancel** is enabled only during an active Run/Resume. The first request sets a
process-shared event, disables repeated requests, and displays cancellation pending.
The Qt timer continues receiving logs while the editor remains locked.

The worker checks the event before canonical execution and after completed,
durably logged progress records (initial retention, completed trial/rejection,
retained stage boundaries, and verification). It never injects an exception into
propagation or an atomic checkpoint write. At the first observed request it invokes
the runner's existing `KeyboardInterrupt` handling once. That handling preserves
the confirmed retained boundary and may perform final verification and produce
reports. Cancellation therefore can take as long as the current calculation plus
safe finalization; there is no force-kill fallback. Optimizer algorithms, checkpoint
identity, evaluation semantics and CLI handling remain unchanged.

The GUI reports **Cancelled** only when the worker actually observed a request at
a safe boundary. It logs both the request and the effective boundary/evaluation.
Canonical on-disk artifacts retain their existing interrupted/limited status and
resume policy; Cancelled is the frontend's additional explanation, not a new
checkpoint format. Failures, including finalization failures after cancellation,
remain **Failed** with tracebacks. Unrequested interruptions/budgets remain
**Interrupted**. A request arriving after the last observable safe boundary cannot
undo a completed run: completion is reported and the late request is explained.

Cancellation before execution creates no output directory or checkpoint. On an
early-cancelled Resume, existing checkpoints remain untouched. Subsequent Resume
uses the same canonical compatibility checks; in-stage work resumes from the last
confirmed boundary according to the existing restart policy, not an invented
mid-optimizer state. Run controls and editor/history operations are restored only
after the worker exits. Immutable configuration/data snapshot isolation is retained.

## Compatibility / Legacy menu (Stage 6F)

These optional tools are separate from normal Open/Save/Run. Neither Run nor
Resume creates or consumes an exported FOR007 file.

**Import Legacy Project…** selects FOR007 and its header-bearing experimental
file, then requires explicit migration choices: seed, serial/batched execution,
batch size if batched, objective mode, invalid-trial policy, no imported hard
bounds, canonical output paths, snapshot/checkpoint intervals and evaluation/time
budgets. Zero explicitly disables snapshots or budgets where indicated; no choice
is inferred from absent legacy fields. Add new bounds later through Structure /
Bounds. Relative paths resolve against the selected FOR007 directory. The existing
Stage 5 migration API owns parsing, ambiguity rejection and stable-ID generation.
A failed import leaves the current document unchanged. A successful import opens
an unsaved canonical document; Save invokes Save As YAML. Input files are never
rewritten by import.

**Export FOR007…** asks for a new filename and explicitly verified target-build
EXPINT/BMNAME filename widths. It uses the current resolved canonical filenames;
standalone export requires existing header-validation data mode. For short bundled
names use Export Legacy Bundle; explicit relocation aliases remain available in
the compatibility CLI. No legacy executable capacities are guessed.

**Export Legacy Bundle…** asks for a new directory and delegates to the Stage 5
bundle exporter. It writes `FOR007.dat`, a separate header-bearing `in.dat` and
`compatibility.json`. Canonical experimental data are read, never overwritten.
Both export modes are exclusive: existing destinations are refused, even if an
operating-system file picker would otherwise offer replacement. Choose a new
path. Protected source/reference directories remain protected by the exporter.

Both export actions first show the read-only compatibility report and require
**Export acknowledging omissions**. The report separates representable fields,
FOR007 omissions/compatibility caveats with their actual canonical values, and
structured unsupported/unexportable diagnostics. Errors disable export. Python
execution, seed, hard bounds, objective/rejection policies, output/checkpoint/
snapshot/budget settings, stable IDs and inactive-category metadata are explicitly
reported according to `legacy_export_report`. Representability is not a promise
that an arbitrary legacy executable accepts its dimensions or filenames. Missing
DSMAX is an error; the GUI does not invent it. Data-read/write failures are reported
without changing the document.

**View Compatibility Report…** reviews the committed configuration in bundle mode
without writing files or changing history. Pending field drafts are not committed
by report/export actions; apply or revert them before export. Compatibility actions
are disabled during Run/Resume, including while cancellation is pending. Normal
YAML saving, history, snapshot execution and checkpoint compatibility are unchanged.

## Read-only Results (Stage 6G)

Results is a persistent workspace tab, not a completion popup. The macOS
**Results** menu provides **Show Results**, **Open Existing Run…**, **Open HTML
Report**, and **Open Output Folder**. Run/Resume updates the workspace without
changing the selected page. An open workspace refreshes in place about once per
second. **Show Current / Latest Run** returns from an externally opened directory
to the application's active/latest run; opening an existing directory never
changes the editor document.

The workspace contains Overview, Refined Parameters, χ² History, Pattern Matching,
and Report / Files. Overview identifies the displayed evaluation and the latest
valid completed evaluation separately. “Latest valid” means a finite, completed,
newline-terminated durable progress record, not a currently computing trial.
It is not necessarily an accepted optimizer state or a resumable checkpoint.
Checkpoint status remains explicitly separate.

The parameter table uses the physical active vector recorded for that evaluation,
with initial values and resolved run bounds. The χ² plot preserves recorded
evaluation numbers, stage/loop/resume boundaries, and rejected-bound markers.
Rejected events have no invented objective value. Click a completed plot point or
use the evaluation selector to inspect historical parameters; selecting
**Latest valid (follow live)** resumes following new events. A cancelled or
interrupted run remains partial even if finalization produced a verified snapshot.

Pattern Matching displays stored experimental, scaled objective-space calculated
(including GRADIENT when enabled), and difference maps. Scan columns vary fastest;
new snapshots retain the actual scan basis and corners, used only for display
coordinates (undilated reciprocal H/K). Older snapshots without physical geometry
use clearly labeled scan row/column axes if dimensions exist, or numeric values
only otherwise. Common or independent intensity color scaling affects rendering
only; residual colors are centered on zero. Beam choices retain ordered duplicate
H/K entries. No propagation, objective evaluation, or normalization is performed
by the viewer.

Output additions are confined to reporting: completed progress events include
mapped physical active values; snapshots include scan geometry; one atomically
replaced **latest_evaluation.npz** stores the most recent completed trial's
patterns and evaluation identity. This is published before its normal durable
progress event. The viewer requires a matching event before displaying it, so an
unlogged/in-flight pattern is never considered valid. Only one current pattern
file is kept, avoiding a full pattern archive at every evaluation. Historical
patterns are available only at existing snapshot evaluations; other historical
selections explicitly say unavailable rather than recomputing them. The additional
snapshot write adds output I/O per evaluation, with no extra evaluations,
callbacks, or checkpoint changes.

Old runs remain inspectable, with missing names/geometry/patterns clearly marked.
Some older events store only optimizer vectors; the viewer explicitly reports
that mapped physical values were not separately recorded. Corrupt or truncated
history preserves the valid prefix with a warning. File actions use native
opening and never modify run artifacts.

During the current GUI session, Running, Cancel pending, Cancelled, Interrupted,
and Completed are shown explicitly. Offline artifacts may not distinguish
cooperative cancellation from another interruption; their recorded canonical
status is shown as partial. Before the first completed event, Results shows no
valid evaluation rather than presenting starting values as refined results.
