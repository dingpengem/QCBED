# Modern crystallography GUI and grouped refinement

Launch `python -m qcbed.gui` in the development environment. PySide6 remains
optional. The current QCBED.app includes this integration; see packaging/macos/VALIDATION.md.

## Workflow

1. Open `examples/fcc_al.yml` for the raw-input teaching project and complete
   [the FCC Al study](EXAMPLE_FCC_AL.md). The verbose `fcc_al_grouped_gui.yml`
   remains a historical 840 Å/two-slice replay input requiring its generated body.
   Legacy/pre-sliced documents retain their compatibility editor and run path.
2. In **B. Atomic Model**, choose the structure file, zone axis and original/conventional
   `g1_hkl`/`g2_hkl`. These are not internal HK indices. Inspect the atomic-cell
   view, zone-axis projection and strict B5 compatibility summary.
3. Inspect **Site groups**: stable group identity, element, original members,
   chemical occupancy and absorption. Edit shared B; membership is explicit in
   canonical YAML and is never inferred from simulation occurrences.
4. Select **Conventional** or **Gradient** and its matching prepared experimental
   file. Changing this selector does not transform experimental data. Their
   objective values are not directly comparable. Gradient support is binary
   inclusion, active weight one, with variance in the denominator.
5. In **C / Refinement Targets**, select physical targets and edit hard bounds.
   D1/D2/PLAMDA are optional Expert metadata.
   The catalogue is dynamic; current modern production execution supports
   any nonempty subset of supported depth, site-group thermal/absorption, voltage
   (without absorption), and scan targets; see [the audit](modern_refinement_catalogue.md).
6. Choose ANNEAL, AMOEBA, PARAMIN, SIMPARA or an available stage sequence in **C / Optimizer**. No optimizer mathematics changed. Explicit loops, seed, temperature,
   cooling, tolerance and soft evaluation budget are used without example overrides.
7. **Validate C** validates simulation, targets and optimizer independently of output
   selection. In **D. Run & Results**, select and **Validate Output Directory**.
   **Advanced execution settings** holds execution mode, batch size and run limits;
   modern grouped runs retain the canonical batched backend. Save/Save As YAML,
   then **Run**, **Resume**, or cooperative **Cancel**. Paths
   retain the explicit document base. Runs use immutable snapshots in a separate
   process. Unsupported geometry remains an error, never a full-period projection.
8. Open **Results**, including after Cancel. Inspect the latest completed state,
   best-so-far and retained final separately. Select recorded evaluations, view
   objective, requested depth/B and YNORM histories, and masked patterns.
9. Compare each beam's experimental/calculated/residual images using its recorded
   global normalization. Excluded pixels are hidden, not zero measurements.
   Details show stable beam identity, original hkl/internal HK and residual
   contribution. Viewing results never recomputes propagation or the objective.

The reference group is `Al_fcc`, Al, original sites `(0,1,2,3)`, initial
B = 0.86 Å² and chemical occupancy = 1.0. Diagnostic optimized values do not
replace this baseline. Requested depth and whole-slice effective depth are
separately displayed. A stopped optimizer is not evidence of physical convergence.

## Parameter-family capability audit

These statuses describe the **modern** path; supported legacy bindings remain intact.

| Family | Status / scope / units | Initial and bounds source | Binding and GUI policy |
|---|---|---|---|
| Requested depth | READY; sample; Å | Explicit crystallography depth; explicit target bounds (reference 600–1200 Å) | Original-cell rebuild and strict cyclic whole-slice expansion. Effective depth is derived; never fitted independently. |
| Shared Debye–Waller B | READY; physical site group; Å² | Explicit group B; explicit target bounds (reference 0.5–1.5 Å²) | Original/recast/slice occurrences share a physical owner. Editable; catalogue target uses group ID. |
| Occupancy/vacancy | READ-ONLY; site group; fraction; refinement DEFERRED | Structure/group chemical occupancy; no approved modern refinement bounds | Chemical occupancy is not the projected boundary half-weight. Legacy VAC is a vacancy fraction, not a lossless alias for every modern site/group model. Nonzero legacy vacancy trials remain unsupported. |
| Absorption a/b/c | READY; physical site group; dimensionless weights | Explicit group coefficients and nonnegative rejection bounds | Independent coefficients map to the existing absorptive grating; modern linked-C bridge remains unsupported. |
| Accelerating voltage | READY without absorption; kV | Committed microscope input and positive bounds | Existing energy-dependent sample rebuild; absorption-enabled voltage refinement remains unavailable. |
| Scan A/B/C | READY; six signed reciprocal-basis coefficients | Committed scan corners and explicit bounds | Explicit opt-in fitting; A10 conventions and experimental arrays remain unchanged. |
| Reflection multipliers | NOT APPROPRIATE for current modern path; reflection; dimensionless | Unity in strict B5, legacy ordered records otherwise | No validated original-hkl multiplier target bridge. Legacy duplicate/order semantics remain unchanged. |
| Atom x/y | NOT APPROPRIATE; original sites; fractional coordinates | Structure file; no symmetry-mode bounds | Independent occurrence coordinates are not symmetry-aware physical variables. Display structure only. |
| Legacy DILATION | DEFERRED; vacancy coupling; dimensionless | Legacy coefficient; domain 0–1, not lattice bounds | Legacy `(1 - VAC*DILATION)^(1/3)` coupling in its supported domain; zero vacancy gives unit applied scale. Never label as strain or independent lattice constant. |

Implementation evidence: `refinement_mapping.py` / `refinement_trial.py` retain
legacy vacancy, absorption and voltage restrictions; `crystallography/site_groups.py`
and `crystallography_refinement.py` resolve original physical ownership;
`depth_refinement.py` preserves strict whole-slice depth semantics;
`crystallography/simulation_sample.py` supplies the strict bridge. Existing B7 binding
and propagation equivalence tests remain the scientific authority.
The GUI **Capabilities** table exposes these distinctions without offering fake controls.

## Execution and checkpoint contract

`prepare_run_snapshot` → `run_snapshot` is shared by YAML, canonical Python and
GUI callers. `grouped_run.py` dispatches the existing ANNEAL/AMOEBA/PARAMIN/SIMPARA engines through a dynamic journal adapter with explicit configuration-derived bounds and experimental arrays.
It does not call example scripts or generate temporary YAML. The current gate is
batched execution, masked objective, reject bounds, a nonempty supported target
subset and binary support. Depth need not be active. See [the complete audit](modern_refinement_catalogue.md).

Grouped canonical checkpoints have their own explicit `qcbed-grouped-run` format.
Identity includes the existing ownership-aware canonical identity, compiled
physical/data digest, optimizer/bounds/execution settings and occurrence provenance.
Diagnostic B7 journals and legacy production checkpoints are not silently
reinterpreted. Compatible interrupted runs replay recorded optimizer prefixes;
previous completed evaluations are not recomputed. Cancellation is checked before
new evaluations, and completed child journals survive a parent update interruption.
Results read those committed child journals as well as the parent checkpoint.

Every completed physical evaluation retains its comparison NPZ before journal
commit. This intentionally uses disk space for reliable latest-state viewing after
cancellation. Current modern execution checkpoints each evaluation, independent of
legacy snapshot intervals; interval-based retention is not yet supported. Wall-time
and total-evaluation safety budgets apply per invocation (including Resume).
A bound rejection is recorded separately and does not create a fabricated pattern.
Final output files are written only on completion. The GUI accepts any writable
output directory, creating it during explicit validation if needed. Each Fresh Run
uses a unique managed subfolder, preserving unrelated files and previous results.
Resume separately checks the newest managed checkpoint (or a legacy direct
checkpoint) and requires matching identity. A stale run lock requires explicit
investigation, never automatic deletion.

Modern structure errors, zone/basis/B5 errors, unavailable target selections,
bounds and checkpoint mismatches appear as structured diagnostics. Configuration
validity, compilation, execution capability and external/output readiness remain
separate. Save/Save As preserves stable IDs, ordering and numerical precision.
Legacy import does not invent modern structure/site groups; modern FOR007 export
continues to be blocked where provenance cannot be represented safely.

## Validation and remaining limits

Small synthetic conventional/gradient jobs exercise fresh execution, cancellation,
resume, immutable state, output loading and uninterrupted equivalence. A real
baseline evaluation reproduces gradient active count 7281, YNORM
13447206.106343817 and objective 1.1896911365294378 exactly. No new fitting campaign
is part of this integration.

Before distributing a rebuilt app, run the separate packaging smoke checklist
against this new source: structure-reader dependencies, spawned worker, cancellation,
resume and retained NPZ Results. Signing/notarization remain separate. The current dynamic target and optimizer capability limits are documented in
[the catalogue](modern_refinement_catalogue.md); mutable group membership and the
unavailable physical families remain outside this GUI workflow.

## Standalone macOS application

The Apple Silicon packaging path now explicitly includes ASE readers/spacegroup
resources and modern grouped worker/Results modules. External YAML, structure and
experimental files retain their document-based paths. The app does not require
conda activation or a source checkout as an import path. Use the normal Open,
Run/Resume/Cancel and Results controls described above. READY/READ-ONLY/DEFERRED
capabilities are identical to source execution; packaging adds no scientific modes.

See [build instructions](../packaging/macos/README.md) and the
[measured local acceptance record](../packaging/macos/VALIDATION.md). Ad-hoc signing
is for local use, not Developer ID signing or notarized public distribution.

## A/B/C source workspace

Controls are now grouped under **A. Data Prep**, **B. Atomic Model**, and
**C. Main QCBED Parameters**. Structure/site-group controls remain in B; requested
depth, refinement and microscope/signal mode are in C; execution/output controls
are exclusively in D. A
provides linked pattern/sigma/weight inspection, non-destructive mask snapshots
and explicitly attached Patch A full-pattern geometry. See
[selection, persistence and overview details](gui_abc.md). Front/Top/Left replace
the previous two structure panels. The rebuilt standalone application includes these views; see the latest
[packaging acceptance record](../packaging/macos/VALIDATION.md).

## Compact project documents

Use [examples/fcc_al.yml](../examples/fcc_al.yml) for the raw-input teaching project expressed
as A/B/C/D user inputs. Modern Save/Save As writes this schema automatically;
[modern_yaml.md](modern_yaml.md) documents descriptor assets, defaults and migration.
The existing verbose example remains a compatibility input. No internal HK,
compiled lattice, effective depth, YNORM, objective or optimizer history is dumped
into the editable document. Compiled/run/checkpoint provenance is unchanged.
