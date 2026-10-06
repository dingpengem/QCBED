# QCBED User Manual

[Documentation index](README.md) · [Installation](INSTALLATION.md) ·
[Quick Start](QUICK_START.md) · [FCC Al study](EXAMPLE_FCC_AL.md)

## 1. Purpose and scientific scope

QCBED integrates automated CBED preparation with user-reviewed geometry and physical
inputs, plane-wave multislice simulation, visual refinement-target selection, and
interactive results inspection. Its Python workflow runs on CPU: modern refinement
uses batched evaluation, while simulation APIs also support process-based parallelism
across independent tilts. Normal use requires neither a GPU nor a Fortran compiler.

Select supported physical parameters from the [target catalogue](modern_refinement_catalogue.md) and inspect their starting values, hard bounds
and activation status before choosing ANNEAL, AMOEBA, PARAMIN, SIMPARA or a supported
multi-stage workflow. Checkpointing and compatible-run Resume support continued
work; recorded parameter/objective histories and per-beam pattern comparisons let
you inspect available latest, best and retained states, including after cancellation.
Basic / Advanced / Expert controls keep routine work accessible while preserving
access to numerical, optimizer and execution settings.

The A → B → C → D workflow separates measurement preparation, the atomic model,
simulation/refinement settings, and execution. It preserves established numerical
conventions but does not decide whether a model, mask, bound interval, mesh or
convergence budget is scientifically appropriate for your specimen. Structure and
site-group choices still require user review; the teaching example is instructional.

The maintained application is an Apple Silicon macOS package. Source installation
also provides numerical APIs and a CLI. See the installation guide for the tested
platform and ad-hoc-signing limitations. Normal app use needs no Fortran compiler.

## 2. The project and its files

**One project YAML** is the user-facing configuration entry point. It records input
references and explicit choices; it is not a second optimizer history or a dump of
all derived matrices. The normal actions are **Load Configuration…**, **Save**, and
**Save As**. There are two meaningful input states:

- A raw project points to a CBED pattern (and optional calibration/profile inputs).
  Loading it opens interactive A preparation; it cannot run refinement yet.
- A prepared project points to an immutable preparation descriptor, prepared arrays
  and mask/origin provenance. Loading it restores the recorded inputs without
  reconstructing measurements from a previous fit.

The generated `.qcbed-data` folder, prepared arrays, source pattern and structure
are dependencies of the YAML. **Save As does not copy every dependency** into the
new folder. Keep referenced paths available when moving/sharing a project; this
release has no self-contained modern project-export command. A loaded project
must be revalidated; GUI validation receipts are not serialized as scientific data.

Use **Save As** for your own working copy of a shipped example. Saving a modern
prepared project requires A/B/C validation and committed edits. A raw-input project
can be saved through its raw-project save path before preparation. During an active
worker, editing and project replacement are locked.

## 3. Starting and loading

A blank launch opens **A. Data Prep**, in **New Data Preparation** mode. That is the
new-project entry; there is no additional "New Project" button. **Load Configuration…**
opens either a raw template or a saved prepared project. **Load Previous Preparation…**
replays existing preparation/provenance. Starting a new preparation clears prior A
measurements; retained B/C settings still require validation against the new data.

The four workflow tabs are:

| Tab | Owns | Validation |
|---|---|---|
| A. Data Prep | Experimental signals, geometry, origin and masks | Validate A |
| B. Atomic Model | Structure, zone/basis, site groups and physical factors | Apply Site-Group Parameters; Validate B |
| C. Main QCBED Parameters | Simulation Setup, Refinement Targets, Optimizer | Validate Simulation Setup; Validate Refinement Targets; Validate Optimizer |
| D. Run & Results | Output directory, execution settings, Run/Resume readiness and actions | Validate Output Directory; Validate Execution / Run Readiness |

Basic controls stay visible. **Advanced** exposes numerical/execution choices;
**Expert / Legacy** exposes internal search metadata and compatibility details.
Opening a section never resets a loaded value to an example value.

## 4. A. Data Prep

### Input and optional detector calibration

Browse **Full CBED** to a finite, single-channel two-dimensional image. Supported
readers include numeric text (`.dat`), NumPy `.npy`, and image formats handled by
Pillow (`.tif`, `.tiff`, `.png`). Colour images are not silently converted.
**Aperture (optional)** supplies detector-response calibration/geometry when
available. **Disk-profile override (optional)** is under Advanced settings.
No aperture file is needed for the manual-diameter FCC Al study.

Detector/aperture correction, geometric disk measurement and signal/uncertainty
preparation are different operations. Do not use an unrelated aperture or profile
merely to make validation pass. Inputs and hashes are recorded; changed source
files invalidate replay integrity checks.

### Disk geometry and centrosymmetry

Select **Draw Disk Diameter**, drag across a diffraction disk and click
**Accept geometry / Snap ROI preview**. The diameter defines its physical radius.
Then select **Select Centrosymmetric Feature**, draw a rectangle around a local
symmetry feature, and accept the snapped ROI. The backend searches the local
centrosymmetry origin. A disk centre and the zone-axis/symmetry origin need not be
the same position.

The displayed extraction size and symmetry crop are derived from accepted geometry,
ROI and binning. Advanced binning changes require re-preparation. Array extents may
be expanded for valid bin multiples; that does not redefine the measured radius.
Positions use detector **x = column**, **y = row**; inspect the displayed axes rather
than treating a screenshot's top-left corner as an origin convention.

### Prepare and inspect

Choose **conventional** or **gradient**, then **Run Data Preparation**. In a blank
project, choose a project folder for its cache. Preparation runs in a separate
process; wait for completion. Inspect the full-pattern overview, beam catalogue and
selected reference beam. Use **Confirm recorded zone-axis position** once correct,
or select the symmetry centre on the full pattern and inspect the revised mapping.
Click **Validate Data Prep** only after these measurements and masks are accepted.

Changing the signal-mode selector does not transform a saved experimental file.
A gradient simulation requires its matching gradient-prepared signal and variance;
conventional and gradient objective values are not directly comparable.

### Pattern, Sigma, Variance and Weight

| Quantity | Meaning |
|---|---|
| Pattern | The configured experimental signal for the selected mode/beam |
| Sigma | Standard uncertainty, displayed as the square root of stored variance |
| Variance | Squared uncertainty used by the objective, not raw intensity squared |
| Final Weight | Binary inclusion/exclusion for the modern masked objective |

Preparation estimates/propagates uncertainty using the implemented noise and mode
rules. Do not substitute the raw `sigma.dat` reference for a newly prepared gradient
variance, or interpret display normalization as a change to objective inputs.
For details and units see [preparation algorithms](data_preparation_algorithm.md).

### Mask layers

**Final Weight**, **Geometric Mask**, **Processing Validity** and **Manual Mask** are
separate inspection layers. Geometry reflects beam-radius/extraction support;
processing validity excludes numerically invalid points within that support;
manual exclusions record the researcher's decisions. Final Weight combines them.
Even conventional preparation checks positive/valid input intensities.

Draw a region or use **Advanced settings / Exact pixel selection**, then **Exclude**,
**Restore Selected**, or **Reset Manual Mask**. Stops in the exact pixel box are
exclusive. Restore cannot invent validity outside original geometric/processing
support. Masking changes weight only; it does not repaint Pattern or Sigma.
The immutable mask update is recorded in project provenance and invalidates later
validation. Legacy nonbinary weights can be inspected but not edited by this tool.

## 5. B. Atomic Model

### CIF, zone axis and reciprocal basis

Load **Structure file** and remain on Structure. Inspect lattice/atoms and projections, confirm orientation and g1/g2, then use **Next: Configure Site Groups**. CIF and extxyz are covered by packaged acceptance; other
ASE formats are available through ASE but not all are certified by release tests.
The source cell and original site order are the physical reference. Inspect the
atomic-cell/projection views and the strict transverse-periodicity diagnostics.

Set **Zone axis [u,v,w]**, **g1 hkl**, and **g2 hkl** in the original/conventional
cell's indexing. These define the meaning of experimental lattice coordinates;
they are not the derived internal Fourier indices. They must be compatible with
the chosen zone and periodic simulation cell. Validate B preserves your committed
basis; it must not be inferred from the apparent pixel directions alone.

### Site groups, thermal factors, occupancy and absorption

A site group ties explicit original sites to a physical parameter. The current
bridge supports one physical group per element. Source membership, occupancy,
thermal factors and enabled absorption are validated before compilation.

**Debye–Waller B (Å²)** is an explicit isotropic group input, not U and not a length.
The present GUI does not infer it from absent/unsupported CIF thermal metadata.
When it is missing, enter a justified value; B = 0.86 Å² in the FCC Al study is an
example/reference value, not a universal default. Group ownership applies to all
recast occurrences without inventing separate B values for boundary copies.

Chemical **occupancy** is a fixed input in the modern GUI. It is not equivalent to
legacy vacancy and is not a supported modern refinement target. Absorption is
optional and uses explicit supported coefficients; do not enable it without an
appropriate model. Modern voltage refinement with absorption, linked-C absorption,
and multiple groups of the same element remain restricted.

Click **Apply Site-Group Parameters**, then **Validate Atomic Model**. A pending table edit is
not yet the model used for refinement. Source-site indices and derived mappings
are available in Expert/diagnostic views.

## 6. C. Main QCBED Parameters

C has exactly **Simulation Setup**, **Refinement Targets**, and **Optimizer**.
It has no output-directory or execution-mode editor. C can be Valid while D still
says **Output directory not selected**.

### Simulation Setup

Check **Energy (kV)** (accelerating voltage), **Signal mode**, **Requested depth (Å)** and
**Slice boundaries**. A fresh modern draft starts at **1000 Å**; opening a saved
project preserves its saved depth. Requested depth is a physical starting value,
not an instruction to achieve an exactly continuous propagated thickness: the
existing cyclic whole-slice rule determines the effective simulated depth.

**Slice boundaries** are visible basic inputs, fractional along the recast
simulation-cell propagation direction. Start at 0, end at 1, use finite strictly
increasing values:

- `0, 1` means **one stored slice** per cell period.
- `0, 0.5, 1` means **two stored slices** per cell period.

These are stored slices, not the number of slice propagations through the full
sample. One-slice support retains strict transverse-periodicity checks and does
not certify slicing convergence.

**Use recommended numerical settings** retains the current configured mesh/cutoff;
it is not a specimen-specific convergence test. In **Advanced numerical settings**,
**Transverse simulation mesh X (points)** and **Transverse simulation mesh Y (points)**
sample the projected-potential/wavefield simulation cell perpendicular to the zone
axis. X/columns follow P and Y/rows follow Q; the directions can be oblique.
They are **not detector pixels** or the experimental scan dimensions. Committed
cell lengths, transverse angle and sampling steps are displayed after validation.
The **Reciprocal cutoff (Å⁻¹)** controls reciprocal-space sampling support.
Keep initial settings for the tutorial, but establish mesh/cutoff/slice convergence
before interpreting a scientific refinement.

Click **Validate Simulation Setup** to compile the physical model and expose its
supported target catalogue. Geometry errors remain errors; the GUI does not
silently substitute a non-periodic or full-period projection.

### Projected Elastic Coefficients

**Projected elastic coefficient $E_g^{(s)}$**: the slice-resolved transverse
Fourier coefficient of the elastic potential in the zone-axis-oriented
multislice simulation cell. Reflections are displayed using their original
crystallographic hkl identities.

QCBED currently refines the slice-resolved projected elastic Fourier coefficients
used directly by the multislice calculation. The dedicated C → Refinement Targets
section identifies each coefficient by original hkl and physical slice interval.
Choose Real only, Complex (Re + Im), or Magnitude only; supply physical component
bounds and use Reset to IAM to restore the current atomic model. Elastic Friedel
pairs are linked; absorption remains separate. Explicit overrides become stale
after upstream model changes. These are not conventional 3D unit-cell structure
factors. See the [complete E_g workflow and Results guide](projected_elastic_coefficients.md).

### Refinement Targets and hard bounds

The table shows **Physical target**, **Refine**, **Committed starting value**, **Hard lower
bound**, **Hard upper bound**, and availability/status. Starting values come from
committed physical controls and applied site groups. Changing a starting value
means editing its owning physical control, not a second copy of the parameter.

Choose any nonempty supported subset: requested depth, group B, supported enabled
absorption coefficients, voltage without absorption, or six scan-coordinate
components, and slice-resolved **Projected elastic coefficient E_g** targets. Depth
is optional. Unavailable quantities are not secretly fitted.
The [refinement catalogue](modern_refinement_catalogue.md) lists exact domains and
limitations. Chemical occupancy, original-site positions, reflection-specific
elastic/absorptive multipliers (`reflection_multiplier`), and physical lattice/strain
refinement are unavailable in the modern path. Reflection multipliers are supported
by the pre-sliced compatibility backend/editor, not the modern physical-target GUI;
they scale computed structure factors for fixed reflection records, rather than
refining absolute structure factors.

Every selected target needs finite ordered **physical hard bounds** containing the
starting value. Domain rules still apply (positive depth, nonnegative thermal
factors, and so on). Proposals outside bounds are rejected, not clipped to the
boundary. Bounds are scientific choices, not universal software defaults.
Click **Validate Refinement Targets**. Incomplete drafts remain visible but do not
become active bindings.

### Optimizer and recommended settings

Select a supported workflow: **AMOEBA**, **ANNEAL**, **ANNEAL → AMOEBA**, or the
supported Expert/Legacy alternatives **PARAMIN** and **SIMPARA**. The teaching study
uses AMOEBA. **Use recommended optimizer settings** supplies existing per-target
search metadata; it does not choose scientifically optimal budgets or bounds.
**Advanced optimizer settings** holds loops, temperatures, tolerances, seed and
stage/search budgets. Irrelevant controls are disabled for the selected workflow.
Loaded explicit overrides are preserved until you request recommended values.
Click **Validate Optimizer** after completing setup, targets and optimizer configuration.

In **Expert / Legacy settings — per-target search metadata**, D1/D2 are dimensionless
search/scaling factors and PLAMDA is optimizer metadata, **not hard bounds**.
For AMOEBA the physical simplex displacement is
`PLAMDA × (D2 − D1) × starting value`; a zero starting value or equal factors can
make scaling invalid. PARAMIN uses multiplicative coordinate factors; SIMPARA uses
an absolute uncertainty width. Signed/zero values follow the existing engine rules.
Consult the catalogue before changing these values. Most users need only physical
starting values, hard bounds and recommended metadata.

## 7. D. Run & Results

### Output directory and readiness

Enter **Output directory** or use **Browse…**, then **Validate Output Directory**.
Validation checks only selection, directory type, creatability and writability.
Missing directories are created. Existing/nonempty folders, including old QCBED
outputs, are acceptable. There is no user-editable output filename.

D shows independent statuses: **Output directory: Valid**, **Run: Ready**, and, when
appropriate, **Resume: Not available — no previous run found**. A/B/C validation
alone does not enable Run. D must validate output, and actual runtime inputs must
still pass checks. Editing the output path invalidates its D receipt without
invalidating scientific A/B/C configuration. **Validate Execution / Run Readiness**
applies execution edits and rechecks runtime readiness.

Each fresh modern GUI run owns a unique `qcbed-run-<timestamp>-<id>` subfolder.
Existing files in the selected directory are not replaced, renamed or deleted.
**Open Results** follows the actual run subfolder. The CLI uses an exact output
path instead; see [installation/CLI instructions](INSTALLATION.md).

### Advanced execution settings

The D panel scrolls; expand **Advanced execution settings** to see:

| Control | Behavior |
|---|---|
| Execution mode | Canonical `serial` or `batched`; modern grouped refinement currently requires batched |
| Batch size (blank for serial) | Number of independent tilts scheduled together; affects memory/execution, not physical inputs |
| Maximum total evaluations | Safety budget per invocation; blank means unlimited |
| Maximum wall time (seconds) | Safety limit checked at safe boundaries; blank means unlimited |
| Checkpoint / result cadence | Fixed after every completed evaluation for modern grouped runs |

Fresh modern GUI projects use **Batched / 64**. Closing or never opening Advanced
preserves that default or the loaded explicit setting. The low-level generic
execution API retains its own serial default; it is not a new GUI default.
Worker/process count is not exposed because the canonical run configuration has
no such field; standalone tilt APIs support separate process execution. Legacy
`checkpoint_every` and `snapshot_every` metadata round-trip in YAML but do not
control modern grouped cadence. The UI labels this fixed behavior explicitly.

### Run, Cancel and Resume

**Run QCBED** uses an immutable validated snapshot in a spawned worker. Editing is
locked while it is active. **Cancel** is cooperative: allow the current safe unit
of work to finish; do not repeatedly force-quit or remove a lock file. A very early
stop can leave no evaluated pattern, and a tiny run may finish before Cancel arrives.
A budget limit yields an interrupted partial run rather than a converged result.

**Resume QCBED** is independent of output writability. It requires a checkpoint,
matching configuration/data identity and no active/stale writer lock. The GUI
looks for the newest managed child with a checkpoint, falling back to an older
checkpoint directly in the selected directory. It does not silently choose an
older compatible run when the newest is incompatible; select that older run folder
explicitly. Checkpoints are not interchangeable across edited data, masks,
parameters, settings or output identity. Investigate a stale `.run.lock` before
removing it; never delete a lock belonging to an active writer.

## 8. Results viewer

**Open Results** opens recorded outputs after completion, cancellation or a safety
interruption. It does not launch new simulations. The tabs are **Overview**,
**Refined Parameters**, **χ² History**, **Physical / YNORM history**, **Beam details**,
**Pattern Matching**, and **Report / Files**.

- **Initial:** committed starting vector for this run.
- **Latest:** most recently evaluated valid physical vector.
- **Best:** lowest recorded valid objective; it may differ from latest.
- **Retained:** optimizer-retained/final state when available; do not substitute
  the best or latest state when a completed retained state is absent.

The objective/χ² measures agreement for the recorded signal, variance, inclusion
mask and normalization. YNORM is the recorded intensity normalization. A smaller
objective can reflect a different mask/noise model; compare runs only with their
provenance. A flat history or bound hit is not sufficient evidence of convergence.

Inspect experimental/calculated/residual comparisons beam by beam. Excluded pixels
are not fitted and are displayed as excluded/NaN where applicable. Recorded
patterns and values may be partial after Cancel. Before the first evaluation,
parameter/history/compare views can legitimately contain no data. Check status,
budget, active points and source identity before interpreting a plot scientifically.

## 9. Saving, loading and troubleshooting

**Save** updates the associated YAML; the first save of an unassociated prepared
project uses **Save As**. **Load Configuration…** replaces the working project only
after unsaved/preparation checks. Revalidate loaded A/B/C and D output before Run.
Changing A clears downstream confirmations; B changes invalidate physical setup
and targets; setup changes preserve compatible committed bindings but require
revalidation. A visible draft has not necessarily been committed.

| Message/symptom | Action |
|---|---|
| Select/confirm zone-axis position | Inspect the full-pattern reference beam and confirm or select its symmetry position in A |
| Region too close to image boundary | Choose a slightly smaller local symmetry ROI; do not alter source pixels to bypass it |
| Debye–Waller B required | Enter a justified B for each required group, Apply Site-Group Parameters, Validate B |
| Strict transverse-periodicity/basis error | Check the supplied cell, zone and g1/g2; unsupported geometry is not repaired automatically |
| Incomplete refinement binding / bounds required | Select a supported target, enter ordered physical bounds containing its start, validate targets |
| Zero starting value / invalid search scaling | Inspect the chosen optimizer's domain and metadata; do not reinterpret D1/D2 as bounds |
| Output directory not selected / not checked | Select it in D and click Validate Output Directory |
| Output path is a file / cannot write/create | Choose a directory in a writable location; inspect the reported filesystem reason |
| Resume: no previous run found | Normal for a fresh folder; use Run QCBED when ready |
| Resume identity mismatch | Restore the exact original inputs/settings or start a new fresh run |
| Missing prepared descriptor / hash mismatch | Restore referenced cache/source files or rerun A with the actual data |
| No pattern in Results | Check whether at least one evaluation completed; zero-evaluation initialization has none |
| Controls below the visible area | Scroll D or the expanded panel; Advanced/Expert sections start collapsed |

Known limits include Apple Silicon-only packaged acceptance, ad-hoc rather than
notarized distribution, strict periodic geometry, fixed modern grouped checkpoint
cadence, and the target restrictions above. The app does not prove mesh/slice
convergence, experimental calibration quality, parameter identifiability or fit
uncertainty. See the [release report](RELEASE_REPORT.md) for current validation and
[developer architecture](ARCHITECTURE.md) for implementation boundaries.

## Equality groups and termination reporting

In Projected Elastic Coefficients, select owner rows and use **Share One Refinement Parameter** / **Make Independent** for optional user equality groups. Real-only shares Re while keeping each Im fixed; Complex shares Re and Im; Magnitude-only shares magnitude while retaining each IAM phase. Groups require compatible starts, intersecting bounds and the same slice. Results and D report **Iterations**, **Objective evaluations** and **Stop reason** separately; see the [termination audit](25_iteration_audit.md).


### Resolved refinement parameter layout

After **Validate Refinement Targets**, the **Resolved Refinement Parameter Layout**
shows each shared group, its original hkl/slice members, mode, optimizer coordinates,
starting values and intersected hard bounds. Independent targets and mandatory
Friedel relationships are listed separately. The summary and D dashboard distinguish
**Active physical targets** from **Optimizer parameters** using the compiled layout.
Depth plus two shared Real-only coefficients therefore has three physical targets
and two optimizer coordinates. Multiple groups remain separate; Complex groups
contribute two coordinates each. Draft edits require revalidation.

Coefficient-only selections are supported, including category-filtered configurations:
the model preview may contain zero ordinary coordinates before elastic coordinates
are appended. Sharing remains enforced in the optimizer layout and writeback, not
by copying independently optimized results afterward. YAML and Resume use the
existing physical identities and group representation.


### Sharing through the validation buttons

**Validate Refinement Targets** commits the visible coefficient checkboxes, modes,
starting values and bounds before validating the complete layout. No separate
**Apply Projected Elastic Coefficients** step is required. **Share One Refinement
Parameter** and **Make Independent** immediately commit the sharing change, clear
refinement/C validation and disable Run/Resume readiness. Their feedback reports
the actual compiled vector length before and after the change and asks for
revalidation. For depth plus two Real-only coefficients this is **3 → 2 → 3**.


### Navigation, drafts and validation state

The A/B/C/D workflow strip uses icons, text and color for Not started, Ready to
validate, In progress, Valid, Needs revalidation and Blocked. The active stage has
an emphasized border. The persistent **Next step** message guides the next action.

Site Groups explicitly belongs to the loaded structure. It shows Applied/Modified
state, keeps the user in place after Apply, and provides Back: Structure and
Validate Atomic Model. Applying thermal/occupancy/absorption inputs preserves
orientation/basis drafts and safe projection-view limits.

C is sequential: Simulation Setup → Refinement Targets → Optimizer. Local Validate
buttons report in place. Next buttons become enabled after validation and navigate
only. Sharing, Make Independent, Reset to IAM and target validation retain the
Refinement Targets tab. Compatible target drafts are restored by physical identity
when widgets regenerate; missing hkl/slice identities are flagged for review.

The compact shared-parameter summary shows membership, mode and counts. Expand
Shared parameter groups for starts/bounds and component semantics; expand Resolved
Refinement Parameter Layout for the validated coordinate list. Draft counts are explicitly estimates; validated counts and the
resolved layout come from compilation. Changing parameterization retains compatible
starting components; new/incompatible component bounds must be supplied explicitly.
D exposes readiness states, physical/optimizer counts and return links to each
incomplete configuration section. See the [quick-start screenshots](QUICK_START.md#workflow-guidance-and-shared-parameters).

Compatible target drafts also survive the temporary absence of the catalogue during B → C revalidation. Explicit **Revert drafts** and **Clear saved coefficient targets** intentionally discard the affected drafts; ordinary navigation does not.


### Independent versus shared coefficient parameters

For independent Real-only `(200)` and `(020)` refinement:

1. Check **Refine** for both coefficient targets and enter valid starting values/bounds.
2. Leave their sharing state **Independent**, then **Validate Refinement Targets**.
3. Each reflection has its own real optimizer coordinate: two coordinates, or three
   with sample depth. Choose this when you want the measured reflections to vary independently.

To impose a user-selected equality constraint:

1. Check **Refine**, then separately select both rows using Ctrl/Cmd or Shift.
2. Click **Share One Refinement Parameter**, name the group and review confirmation.
3. Confirm the group members/mode and dimension change beside the buttons.
4. **Validate Refinement Targets** using the persistent bottom action area.

Both Real-only targets now use one shared real coordinate: one parameter, or two
with sample depth. Choose sharing when your refinement model intentionally assumes
those real components should agree; it reduces dimensionality, but is not automatic
space-group symmetry enforcement. **Make Independent** restores separate coordinates.
Real-only members retain their own fixed imaginary parts; Complex shares Re and Im
(two coordinates); Magnitude-only shares magnitude with each member's phase fixed.
The nearby **Example** expands these distinctions and the 3 → 2 depth example.
**Friedel partner — automatic physical constraint** remains separate and mandatory.

C keeps long explanations, IAM/component details, numerical information and the
resolved layout collapsed by default. Select a coefficient and expand **Selected
coefficient details** to inspect IAM, phase, Friedel identity or reset it to IAM.
Simulation Setup keeps its physical inputs visible; Optimizer keeps workflow,
budgets/convergence and the recommended-settings toggle visible. Advanced settings
remain available. Bottom action rows stack when space is insufficient; tables
scroll horizontally for remaining columns without moving the primary actions.
