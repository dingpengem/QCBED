# QCBED release and validation records

## Optional coefficient equality groups and termination audit — accepted

**25-ITERATION AUDIT: PASS**. The two matching real runs stopped at AMOEBA's
soft search-evaluation budget: 12 iterations, 25 objective evaluations (4 initial
vertices + 21 search), `maxfsimp=20`, no convergence. The existing objective
journal was replayed without new forward evaluations. Global-limit values were
not readable in those historical records and are not inferred; new checkpoints
record them explicitly. See the [complete audit](25_iteration_audit.md).

**OPTIONAL COEFFICIENT EQUALITY GROUPS: PASS**. Independent remains the default.
User groups share real, complex or magnitude coordinates within the same physical
slice, separately from mandatory Friedel conjugacy. GUI linking/unlinking,
explicit unequal-start resolution, bound intersections, compact YAML, Resume,
shared histories and individual beam diagnostics are integrated. An optional
exact per-stage iteration cap is separate from unchanged legacy soft budgets.

| Final check | Result |
|---|---|
| Focused milestone and related tests | **169 passed**, 32 warnings |
| Complete suite | **3,007 passed**, zero failures/skips, 225 warnings, 420.77 seconds |
| Source and frozen GUI smoke | **PASS**, zero optimizer evaluations; linked Complex pair, A/B/C/D, YAML, Results counts and cancellation reason, Resume |
| Native app | Rebuilt; **329 arm64 Mach-O binaries** |
| Signature | Ad-hoc; **deep/strict verification PASS** |
| ZIP | Refreshed; **integrity PASS** |

Two historical hash assertions were updated to use the original six-field RunConfig
schema, excluding the newly optional unset cap. The 35 compatibility tests passed,
and an existing pre-milestone frozen checkpoint remained identity-compatible.
The complete suite was then rerun to the final green result above.

The complete suite used an isolated copy of the current source/tests with private
maintainer fixtures for five existing provenance checks. No private source was
restored to the shared tree or packaged. No new refinement campaign ran. Run
artifacts used by the diagnosis were read-only. Local logs are in
`tmp/equality_milestone/` (ignored); artifact hashes are recorded in
[release provenance](../packaging/macos/RELEASE_PROVENANCE.json).

Limitations: same-slice user equality only; no automatic space-group constraints;
no independent magnitude-plus-phase mode or conventional 3D F_hkl mapping. Legacy
soft budgets can overshoot; the optional exact iteration cap does not silently
reinterpret them. Changed cap/group identity requires a new compatible run.
Historical missing counters/limits remain unavailable. Apple Silicon-only app,
ad-hoc signed and **not notarized**; FCC Al remains instructional.

## Projected elastic coefficient milestone — accepted

**PROJECTED ELASTIC COEFFICIENT REFINEMENT: PASS**, using the agreed **Projected elastic
coefficient E_g** terminology. QCBED currently refines the slice-resolved projected
elastic Fourier coefficients used directly by the multislice calculation. This is
not conventional 3D unit-cell structure-factor refinement.

| Check | Final accepted result |
|---|---|
| Focused E_g tests | **19 passed**; zero/very-low optimizer evaluations |
| Complete regression suite | **2,984 passed, zero failures/skips**, 225 existing Qt/Matplotlib warnings, **578.49 s** |
| Source GUI and frozen GUI smoke | **PASS**, zero optimizer evaluations |
| Single-project YAML | **PASS**: hkl/slice identity, mode, complex starting value, physical bounds and existing A/B/C/D state restored |
| QCBED.app build | **PASS**, native arm64; all **329 Mach-O files** verified arm64 |
| Signing | **Ad-hoc signed**, deep/strict verification passed after frozen smoke |
| QCBED-arm64.zip | Rebuilt; ZIP integrity **PASS** |

Final checklist closure includes the exact definition below, selected-coefficient
IAM/starting Re, Im, magnitude and phase details, and individually labeled real
and imaginary bounds in Complex mode. The focused GUI regression covers these
controls; the full suite and app were rerun/rebuilt after these final changes.

**Projected elastic coefficient E_g^(s): the slice-resolved transverse Fourier
coefficient of the elastic potential in the zone-axis-oriented multislice
simulation cell. Reflections are displayed using their original crystallographic
hkl identities.**

The full suite ran in a temporary validation copy with byte-identical current
Python source/tests and the private predecessor fixtures required by five existing
provenance checks. Only external predecessor evidence and the original private source inventory
were overlaid there; no private source was restored to the shared repository or
included in the app. A shared checkout without those fixtures cannot reproduce
those five checks. The unchanged assertions were executed, not skipped or weakened.
The first run found one outdated Results-tab assertion; it was updated for the new
tab. Final review also added the zero-refined-magnitude/IAM-phase regression; the
reported complete suite includes that correction.

Accepted functionality: original-hkl/physical-slice identities; Real-only,
Complex (Re + Im), and Magnitude-only modes; elastic Friedel linking; separate
absorptive A_g; explicit bounds; Reset to IAM and stale-reference detection;
all four dynamic optimizers; physical E_g Results states/history and beam links.
IAM settings reproduce baseline gratings, amplitudes and intensities exactly in
the focused cases. The equivalent legacy multiplier comparison differs by at most
**8.1e-18** in complex amplitudes (floating-point rounding), within the test's
1e-14 relative / 1e-16 absolute tolerance. This is not a claim of unrestricted
historical executable equivalence or scientific convergence.

Frozen acceptance drove the real preparation workflow and packaged E_g controls,
then saved/reopened one project YAML and initialized a cancelled checkpoint before
any optimizer evaluation. It also retained one-slice [0,1], g1/g2, explicit B,
three C tabs, D-only output validation, Batched/64, Resume and Open Results behavior.
See [E_g workflow](projected_elastic_coefficients.md),
[FCC Al instruction](EXAMPLE_FCC_AL.md), and
[build provenance](../packaging/macos/RELEASE_PROVENANCE.json).

Known limitations: Apple Silicon-only packaged app; ad-hoc signed, **not notarized**;
FCC Al is instructional; conventional 3D F_hkl mapping and automatic space-group
constraint inference remain future work. Elastic (000) is read-only. Magnitude-only
mode requires nonzero IAM phase. Modern occupancy, atom-coordinate, lattice/strain,
vacancy/DILATION and linked-C absorption refinement remain deferred; voltage
refinement still requires absorption disabled. Legacy multiplier YAML remains on
its compatibility path rather than being silently reinterpreted.

## Previous release acceptance

**Historical release acceptance: PASS (before private-source separation).**

The source-distribution policy has since changed. Private predecessor sources and
detailed inventories are no longer distributed. Five unchanged source-dependent
checks require private maintainer assets; this shared checkout alone cannot
reproduce the full-suite count below. Python behavior and scientific reference data
remain unchanged. See [source provenance](source_provenance.md).

This closure verifies the final working tree after the earlier readiness-test
correction. It adds no features or physics changes and runs no refinement campaign.
The only executable-file change in this closure extends the opt-in packaging smoke
assertions; production Python and test sources remain unchanged.

## Project overview

QCBED combines automated, user-reviewed experimental preparation with CPU
plane-wave multislice simulation, visual physical-target selection, modular
optimization, checkpoint/Resume support, and interactive per-beam results.
The Python application presents a Basic / Advanced / Expert workflow; modern
refinement uses batched execution, and simulation APIs additionally support
independent tilt-level process parallelism. No GPU backend or Fortran compiler
is needed for normal use. See the [software highlights](../README.md#highlights)
for the concise capability overview and the [modern target catalogue](modern_refinement_catalogue.md)
for supported physical parameters. Reflection-specific elastic/absorptive multipliers
are supported only in the pre-sliced compatibility workflow. The modern GUI now
adds [slice-resolved projected elastic coefficients E_g](projected_elastic_coefficients.md)
as a distinct physical parameterization, not conventional 3D structure factors.
The acceptance evidence and limitations below remain unchanged. The FCC Al study demonstrates the workflow, not scientific
convergence or automatic selection of all crystallographic site groups.

## Tests

The complete automated suite passed: **2,965 passed, 225 warnings,
422.22 seconds; zero failures/skips**. One full-suite invocation was needed
for this closure. This supersedes the previous 2,964-pass/one-stale-assertion run;
the earlier record remains in [packaging validation history](../packaging/macos/VALIDATION.md).
The warnings are the existing Qt/Matplotlib deprecated high-DPI enum notices.
No test removal or numerical-reference modification was needed.

Command: `QT_QPA_PLATFORM=offscreen python -m pytest -q` in the documented
`qcbed-python` environment, with a writable Matplotlib cache.
Local log: `tmp/release_acceptance/full.log`.

## Frozen GUI and project YAML

Both the pre-rebuild packaged app and final rebuilt/signed app passed their
zero-evaluation acceptance smoke. The final smoke started with **New Data
Preparation**, the actual GUI's new-project entry point, without a YAML. It used
`reference_data/pattern.dat` and `reference_data/Al_FCC.cif` and covered:

- A → B → C → D navigation and real data preparation/validation.
- Editable visible `0, 1` boundaries accepted as one stored slice; 1000 Å depth.
- Explicit B = 0.86 Å², site-group Apply/Validate B, preserved g1 = 200 and g2 = 020.
- Populated Refinement Targets with depth/B selections and explicit study bounds.
- Optimizer page with explicit objective tolerance `2e-6` and stage budget `23`;
  both values were restored in widgets and committed configuration after loading.
- D-only output-directory creation/validation, Run readiness and the independent
  fresh-folder message **Resume: Not available — no previous run found**.
- Visible D Advanced execution controls, Batched mode and batch size 64.
- One compact project YAML created through Save As's document-save path; actual
  Load Configuration restoration of A/B/C/D state, target selections/bounds,
  thermal factors, basis, thickness, slicing, optimizer and output/execution values.
  No duplicate YAML, derived matrices, validation receipts or optimizer histories
  were serialized; no user-facing output filename was added.
- Open Results disabled for the fresh project; enabled after opening the cancelled
  initialization checkpoint, with **zero completed optimizer evaluations**.
- Compatible Resume after initialization and preservation of unrelated folder files.

The preparation and initialization smoke does not produce or imply a converged
fit. No optimizer evaluation was performed. Local evidence and screenshots are in
`tmp/release_acceptance/final-frozen/`, including `result.json` and `project.yml`.

## Repository and documentation

`git status` was inspected without resetting/discarding the accumulated work.
`git diff --check`, conflict-marker checks, Python syntax, local documentation
links and user-guide/example absolute-path checks passed. Generated acceptance
logs, projects and build products remain in ignored `tmp/`, `build/` and `dist/`;
no unintended generated/debug files were added to source control. Existing valid
uncommitted work remains intact.

All required guides exist and match the implemented controls:

- [Installation](INSTALLATION.md)
- [Quick Start](QUICK_START.md)
- [User Manual](user_manual.md)
- [FCC Al example study](EXAMPLE_FCC_AL.md)
- This release report
- [Modern refinement catalogue](modern_refinement_catalogue.md)
- [Packaging validation](../packaging/macos/VALIDATION.md)

C has exactly Simulation Setup, Refinement Targets and Optimizer. Output selection
and Advanced execution belong only in D. There is no separate New Project button;
New Data Preparation starts the project. The example paths exist, and the FCC Al
study explicitly identifies its settings and short budget as **instructional,
not a converged scientific refinement**. Previous cleanup details remain in
[the repository audit](REPOSITORY_AUDIT.md).

## Accepted package

The final build was performed after the full suite and pre-build smoke checks
passed. The signed build then passed the extended frozen smoke.

| Check | Accepted result |
|---|---|
| `QCBED.app` rebuild | PASS |
| Native executable architecture | arm64 (Apple Silicon) |
| Signing | Ad-hoc |
| `codesign --verify --deep --strict` | PASS |
| ZIP rebuilt after signing | PASS |
| Every ZIP entry CRC / decompression | PASS |
| Archived executable matches signed app | PASS |

Artifacts: `dist/QCBED.app` and `dist/QCBED-arm64.zip`.
ZIP size: **109,884,010 bytes**.
ZIP SHA-256: `233f97843a17cbb74a05c60c3af7efa3d5008974d1c2af304c830eb585d8133f`.
Build timestamp: **2026-10-05T22:07:00.655442+00:00**.
Source fingerprint (package Python plus packaging Python/spec):
`9a687a938ddd1ec94905ace540dc7bf2771fcd236d25b8baf54a0d9ffed4953d`.

Build versions: qcbed 0.1.0; pyinstaller 6.22.3; pyinstaller-hooks-contrib 2026.8; PySide6 6.11.2; numpy 2.4.6; scipy 1.17.1; h5py 3.16.0; matplotlib 3.11.2; psutil 7.2.2; PyYAML 6.0.3; ase 3.29.0.
The [machine-readable provenance](../packaging/macos/RELEASE_PROVENANCE.json)
records the exact accepted artifacts, dirty-tree base and acceptance results.
This is local acceptance, not a release publication or new commit/tag.

## Known limitations

- The packaged application is **Apple Silicon macOS only**, **ad-hoc signed** and
  **not notarized**. No Developer ID or Finder-quarantine distribution acceptance
  is claimed; no hosted release was published by this task.
- The FCC Al example is instructional. Experimental validity, mesh/slice convergence,
  parameter identifiability and physical convergence remain researcher judgments.
- Modern occupancy, atom-coordinate, reflection-multiplier and lattice/strain
  refinement remain unsupported. Vacancy/DILATION are deferred; linked-C absorption
  and absorption-enabled voltage refinement are unavailable on the modern path.
  See the catalogue for precise domains and legacy distinctions.
- Strict transverse periodicity remains required. Canonical modern execution uses
  batched scheduling, has no worker-count setting and checkpoints each evaluation.
- Project YAML references external inputs/cache; Save As is not a portable data bundle.
  A project-wide redistribution license has not been supplied.


## Shared refinement parameter UI clarification

The coefficient editor now uses **Share One Refinement Parameter**, **Make
Independent**, and **Parameter sharing** (`Independent` / `Shared: Group 1`).
Live counts distinguish physical targets from independent/shared optimizer
coordinates; confirmation explains the mode and previews the count reduction.
Results displays **Shared refinement parameter group**. Mandatory Friedel
relations remain separate. Physics, optimizer coordinates, YAML keys and stored
checkpoint identities are unchanged.

Validation for this presentation patch: **52 passed**, zero failures, 48 Qt
deprecation warnings (19.83 seconds). The previous **3,007-pass full suite** is the
preceding milestone baseline; it was not rerun for this UI patch. Frozen GUI
smoke **PASS**, including the packaged controls, Complex dimensionality 6 → 4,
compact YAML reload, Results and Resume, with **zero optimizer evaluations**.

Rebuilt app and ZIP: **PASS**. All **329 Mach-O binaries are arm64**; ad-hoc
signing, deep/strict signature verification and ZIP integrity **PASS**.
`git diff --check` and documentation/link/source hygiene checks passed.
Local validation records are under `tmp/sharing_ui/`; artifact hashes are recorded
in `packaging/macos/RELEASE_PROVENANCE.json`.

Limitations remain: Apple Silicon macOS only, ad-hoc signed and **not notarized**;
FCC Al is instructional, not a converged refinement. No refinement campaign was
run. Conventional 3D structure-factor refinement and automatic space-group
constraints remain deferred.


## Resolved refinement layout and dimension audit

The audited shared-group backend already allocated one coordinate per shared
Real/Magnitude component and two per Complex group. Duplicate optimizer coordinates
were not reproduced. Ambiguous **Active parameters** displays were replaced by
**Active physical targets** and **Optimizer parameters**, derived from compilation.
Validation now exposes each group's members, coordinates, starts, intersected bounds,
independent targets and mandatory Friedel relationships in a resolved layout.

A separate compiler defect was fixed: category-filtered coefficient-only selections
were rejected before coefficient coordinates could be appended to the model preview.
See [the implementation audit](shared_parameter_layout_audit.md) for the trace and
scope. Stored group/YAML/checkpoint identities and propagation mathematics are unchanged.

Focused regression: **119 passed**, zero failures, 200 Qt deprecation warnings,
48.79 seconds. This includes multiple groups, a three-member group, GUI validation,
coefficient-only dimensions, starts/bounds, writeback, YAML, checkpoint compatibility
and Results identities. The prior **3,007-pass full suite** remains a baseline;
it was not rerun for this focused change.

Frozen GUI smoke **PASS**, zero optimizer evaluations: resolved layout, A/B/C/D,
compact YAML, Results and Resume. Rebuilt app and ZIP **PASS**; all **329 Mach-O
binaries are arm64**, ad-hoc signed, deep/strict signature verification and ZIP
integrity passed. `git diff --check` and release hygiene checks passed. Local logs
are in `tmp/resolved_layout/`; current package hashes are in release provenance.

Limitations: Apple Silicon macOS only; **not notarized**. FCC Al remains instructional,
not a converged scientific refinement. Conventional 3D structure factors and automatic
space-group constraints remain deferred. No refinement campaign was run.


## Actual sharing-button validation fix

Fixed the missing coefficient commit in **Validate Refinement Targets**. Share and
Make Independent now invalidate refinement/C readiness explicitly and show compiled
dimension transitions. All visible optimizer counts use compiled dimensions.
Actual frozen widget flow: **PASS**, Real/Magnitude with depth **3 → 2 → 3**,
Complex **5 → 3 → 5**, zero optimizer evaluations. Final focused tests: **83 passed**,
176 Qt deprecation warnings, 43.64 seconds. Full suite was not rerun.

Rebuilt app/ZIP: **PASS**, **331 arm64 binaries**, ad-hoc signed, deep/strict and ZIP
verification passed. Not notarized; Apple Silicon only. FCC Al remains instructional.
No physics changes or refinement campaign. The earlier 3,007-pass full suite remains
historical baseline evidence. Current hashes are in release provenance.
See [the state-layer audit](sharing_gui_state_audit.md) for the reproduced defect,
exact compiled coordinates and receipts. The pre-fix diagnostic was intentionally
not counted as a passing smoke; final results are under `tmp/sharing_actual/after/`.


## GUI workflow and navigation usability

CIF loading now stays in **Structure**. Users inspect lattice/composition and
orientation, then explicitly continue to Site Groups. Apply stays in Site Groups;
local validation and Back/Next actions preserve g1/g2 and compatible view state.
Member sites are visible and read-only unless Advanced editing is enabled.

C retains the active sub-tab during sharing, Make Independent, mode/bound edits
and validation. Compatible drafts survive widget regeneration and upstream
revalidation using physical target identities. Clear/Revert deliberately discard
drafts. Coefficient edits explicitly invalidate C regardless of the visible
parent tab; this also fixes a pending-draft flag surviving Save.

The persistent A/B/C/D strip combines text, icons and muted state colors with a
Next step message. Each C page validates in place and offers explicit Next
navigation. Shared groups, component semantics and draft counts update in place;
validated counts and the resolved layout come from compilation. D shows readiness
and return links while keeping output validity separate from Resume availability.

Final frozen navigation smoke: **PASS**, including two groups, **5 physical
targets / 3 optimizer coordinates**, and Run readiness. Separate frozen
YAML/Results/Resume smoke: **PASS**. Both used **zero optimizer evaluations**.
Screenshots and instructions were refreshed in Quick Start, User Manual and the
FCC Al tutorial. No numerical behavior change or refinement campaign was needed.

Acceptance: **3,025 passed, zero failures**, 273 Qt deprecation warnings, 443.46
seconds on the final source/test tree. Focused GUI/state regression: **51 passed**;
the corrected member-site visibility/read-only assertions: **3 passed**. An initial
full run found one stale hidden-column assertion; after correcting it, the entire
suite was rerun to the green result above. Five provenance checks used private
maintainer fixtures in a byte-identical temporary source/test copy; no private
sources were restored to the shared tree or packaged app.

Rebuilt **QCBED.app** and **QCBED-arm64.zip**: PASS. All **331 Mach-O binaries** are
native arm64. Ad-hoc signing, deep/strict signature verification and ZIP integrity:
PASS. Source fingerprints match the tested tree and bundle. `git diff --check`,
local documentation links, generated-file and absolute-path hygiene checks: PASS.
Acceptance logs are local under `tmp/workflow_usability/`; artifact hashes are in
`packaging/macos/RELEASE_PROVENANCE.json`.

Limitations remain: Apple Silicon macOS only; ad-hoc signed, **not notarized**;
FCC Al is instructional, not a converged scientific refinement. Conventional 3D
structure-factor refinement, automatic space-group constraints and the unsupported
physical targets listed in the refinement catalogue remain deferred.


## C-tab layout cleanup and explicit sharing guidance

Released the old section scroll-area ownership before embedding C form bodies,
fixing content that remained constrained to 602 px despite window resizing.
Reproduced clipping, not a literal button-on-button collision. C now separates
scrolling content from responsive bottom validation/navigation actions. Technical
help, IAM/components, group bounds and the resolved layout are collapsed by default.
Sharing steps remain visible and distinguish Refine checkboxes from row selection;
live group/mode summaries and dimension changes remain immediate. The expandable
example explains Real-only 2 → 1 coordinates, or 3 → 2 with depth, separately from
mandatory Friedel constraints. No validation/scientific semantics changed.

**46 focused tests passed**, zero failures, 96 Qt deprecation warnings, 43.09 seconds.
The stricter geometry follow-up passed **2 tests**, 16 warnings, 5.70 seconds.
Final frozen A/B/C/D, actual sharing and responsive C-page smoke: **PASS**, at
960×720, 1250×850 and 1700×950, zero optimizer evaluations. The full suite was not
rerun for this layout-only pass; 3,025 passed is the preceding milestone baseline.

App/ZIP rebuilt: **PASS**. All **331 native binaries are arm64**; ad-hoc signing,
deep/strict signature verification and ZIP integrity passed. Not notarized.
User Manual, FCC Al tutorial, Quick Start and screenshots were updated.
Local logs: `tmp/c_tab_cleanup/`. Current artifact hashes are in release provenance.
Existing scientific/platform limitations remain unchanged; FCC Al is instructional.


## Final input-only cleanup and publication acceptance — 2026-10-06

The shared input directory now contains exactly **pattern.dat** (raw FCC Al CBED)
and **Al_FCC.cif** (crystal structure), plus its README. The FCC Al project starts
from these two inputs and regenerates preparation products; no hidden prepared
cache or saved sigma map is required. Ten historical numerical/parser/provenance
files moved byte-for-byte to `tests/fixtures/legacy_reference/`, with purpose and
SHA-256 manifest. The synthetic demo intensity moved to `tests/fixtures/synthetic/`.
Existing test-only geometry/catalogue fixtures remain documented. No tests were removed.

Obsolete fit/result snapshots, weights/potential dumps, Linux executables,
scheduler/crash records and preparation caches were archived locally in ignored
storage. Build/dist were already untracked; app/ZIP outputs remain local. Ignore
rules protect recurring caches/builds/logs and private source. Machine-specific
text paths and two old screenshots were removed from publication. Current README
Highlights cover projected coefficients and optional sharing; user/developer docs
and example fixture paths were reconciled. No scientific equations changed.

**Full suite: 3,027 passed, zero failures, zero skips**, 289 Qt deprecation warnings,
441.75 seconds, run directly from the cleaned working tree. The five original
private evidence checks passed against an external maintainer directory without
copying private source into the repository. Shared users without that external
setting see five explicit optional skips; numerical tests remain available.
Focused data/provenance checks: **413 passed**. Minimal automated GUI/state and
single compact YAML round-trip: **PASS**, zero objective evaluations. Exact user
configuration, shared groups, physical identities, optimizer dimension, checkpoint
identity, output validation and Results loading survived reload.

Final packaged **startup-only smoke: PASS**, zero evaluations. The user's full
manual workflow was not repeated. Rebuilt app/ZIP, **331 native arm64 binaries**,
ad-hoc signing, deep/strict verification and ZIP integrity: **PASS**. No private
source files are included. `git diff --check`, Markdown links, privacy and staged
artifact audits passed. Source/artifact fingerprints are in release provenance;
local detailed logs are under `tmp/final_cleanup/`.

Limitations: Apple Silicon macOS app only; ad-hoc signed, **not notarized**; FCC Al
remains instructional. Conventional 3D F_hkl refinement/mapping and automatic
space-group equality remain deferred. Occupancy, atomic-coordinate, lattice-strain
and vacancy refinement, voltage with absorption, and modern linked-C absorption
remain unsupported as documented in the target catalogue. No refinement campaign
or history rewrite was performed.
