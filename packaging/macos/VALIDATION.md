# Live Patch A frontend packaging closure

Build UTC: 2026-10-05T01:36:15.678421+00:00; native arm64; source/spec fingerprint
`ddef5a08d807ed8771cf029098b2b670f990847f68d639662517c9b058901959`.
`dist/QCBED.app` is ad-hoc signed, with deep/strict verification. ZIP refreshed.

Source and frozen acceptance reports:
`tmp/live_preparation_acceptance/source/acceptance.json`,
`tmp/live_preparation_acceptance/frozen-final/acceptance.json`.
The frozen run completed in 18.95 seconds via native Launch Services. It uses a
new synthetic raw pattern, not an archived catalogue: optional aperture threshold
preview, manual diameter, explicit ROI snapping, a spawned Patch A worker, nine
live beam ROIs, refined-origin confirmation, Pattern/Sigma/Weight, masking and
Validate A before Atomic Model setup. It also verifies the replay-mode transition.
No optimizer or refinement run is started.

Separate synthetic API execution completed actual A2 calibration, A3 correction,
A4 fitting and A5–A10 preparation (180×180 blurred three-disk scene with an intensity
ramp and seeded noise; existing quick PSF test budget), returning (1600,9) arrays.
A constant-background scene correctly fails A4's effective-fit-point requirement;
there is no automatic Poisson fallback after a calibration/noise failure.
The real pattern still explicitly fails live lattice detection; archived geometry
is never substituted. See the algorithm/workflow docs for integer-crop limitations.

Focused tests: 29 passed. GUI plus packaging: 106 passed (102 GUI, 4 packaging).
Full suite: **2874 passed** (343.14 seconds). Final GUI-state checks: **4 passed**.
`git diff --check` clean; source fingerprint matches the final build. No science kernels or refinement campaigns changed.
Ad-hoc signing is local integrity only, not notarized distribution.

---

# Compact modern YAML packaging closure — 2026-10-04

Build UTC: 2026-10-04T22:26:23.353148+00:00; native arm64, PyInstaller 6.22.3.
Source/spec SHA256: `0627fba86fcc884b0953d9902eddb32f342fbf3d80f7a3737f6c4fae076efdd3`.
Outputs: `dist/QCBED.app`, `dist/QCBED-arm64.zip` (ad-hoc signed, not notarized).

Frozen acceptance passed via Launch Services with CWD `/`, and with a clean
shell environment. `qcbed.modern_yaml` was loaded from the app's Frameworks,
not source-tree code. Reports: `tmp/compact_packaging_closure/frozen/acceptance.json`
and `tmp/compact_packaging_closure/frozen-clean/acceptance.json`.

Checks cover `examples/fcc_al.yml`, external arrays and pattern/provenance assets,
manual exclusion and origin restoration, nine beam annotations, A/B/C tabs,
atomic views/site groups/refinement controls, compact Save/Save As/reload with
exact canonical and array equality, verbose modern and legacy YAML opening.
The interactive check uses an explicitly prepared persisted-selection fixture;
the canonical example itself does not invent missing interactive metadata.

Canonical run initialization created output/checkpoint artifacts with **zero
evaluations**. The packaging harness intentionally intercepted stage dispatch
before optimizer entry. Its exception is recorded as `failed` by the unchanged
runner; this smoke checkpoint is not evidence of a completed refinement.
No ANNEAL/AMOEBA campaign was rerun. Physics and normal execution paths were unchanged.

Focused tests: **16 passed**. GUI tests: **98 passed**. Full suite: **2861 passed**
(343.33 seconds). `git diff --check` clean.
Deep/strict signature verification passes; source fingerprint matches the build.
External scientific assets must remain accessible at their resolved paths.
Native file-picker gestures and notarized distribution remain separate checks.

---

# Ordered A/B/C GUI workflow correction — 2026-10-04

Build UTC: 2026-10-04T22:00:00.957847+00:00; native arm64, existing PyInstaller path.
Source/spec SHA256: `c9368ddd69fdcbd91fae90bdcd81060990673bf2a1323fcf778f9a73b334282c`.
Outputs: `dist/QCBED.app`, refreshed `dist/QCBED-arm64.zip`.
Ad-hoc signed; deep/strict signature verification passes. Not notarized.

This build contains actual A/B/C tabs and ordered validation receipts, manual
symmetry-center selection using existing symmetry-center/A10 functions, separate Zone axis and
central-origin markers, experimental/conventional beam annotations, and retained
weight-only mask editing. Existing output-only paths are accepted after explicit
backup-and-replace confirmation at Fresh Run; canonical runner and CLI overwrite
rules and checkpoint compatibility remain unchanged.

Validation: **20 focused passed**, **101 GUI/packaging passed**, **2848 full-suite passed**
(336.59 seconds); `git diff --check` clean.

Short native Launch Services acceptance only:
`tmp/gui_workflow_closure/frozen-final.json` and
`tmp/gui_workflow_closure/frozen-final/acceptance.json`.
The smoke checks tabs, marker selection, nine labelled boxes, unchanged
signal/variance under masking, A→B→C gating, existing-output readiness without
mutating previous outputs, and Save/reload. It explicitly records
`refinement_executed: false`. No packaged Run/Resume refinement acceptance or
real-case scientific campaign was repeated. Normal regression tests remain separate.

The interactive pick is manual, not an automatic centrosymmetry fit. New origin
selection requires explicit existing Patch A geometry; unprepared raw-image
processing is not added. Lower Data Prep controls/overview remain scrollable.
Older prepared inputs need matching overview metadata and explicit origin
confirmation before A validation. Changes to canonical inputs remain subject to
the unmodified resume-identity checks.

# A/B/C packaging closure — 2026-10-04

Current outputs: `dist/QCBED.app` and `dist/QCBED-arm64.zip`.
Native arm64; ad-hoc signed, deep/strict verification. No notarization.
Build UTC: 2026-10-04T21:30:58.118681+00:00.
Source/spec SHA256: `b87ed34a86d147d2da8db9293e7e9f19f8e82b71ba5ebf4f66a8d4dc0e8f2802`.
Full provenance: `build/macos/build-environment.json`.

This closure extends the existing frozen acceptance harness for A/B/C and fixes
Data Prep canvas compression in the existing scrollable page (minimum canvas
sizes and layout minimum-size propagation). No numerical/configuration semantics
changed. A cancellation test no longer assumes a cooperative stop must take
multiple GUI event-loop turns; cancellation state/log/control checks remain.

Acceptance manifests and reports are under `tmp/abc_packaging_closure/`.
The final build uses `*-closure.json` manifests and `*-closure/acceptance.json`
reports. These are local generated artifacts, not committed fixtures.

- Native Cocoa launch through macOS Launch Services and reopening saved YAML.
- Real A/B/C inspection config: linked pattern/sigma/weight values and beam
  selection; exclude/restore/reset; exact intensity/variance preservation;
  Save, Save As and reload of a manual exclusion (7280 remaining active points).
- Full 256×256 pattern, nine boxes/labels, origin (127,126).
- Front/Top/Left projections, Al_fcc group/B editing, ASE extxyz/CIF resources.
- Dynamic two-target controls, signal mode, fresh-path selection and Resume state.
- Separate controlled run fixture retains baseline gradient support 7281
  (unit weights), YNORM 13447206.106343817, objective 1.1896911365294378.
  Short ANNEAL→AMOEBA, cooperative Cancel after completed evaluation,
  compatible Resume and incompatible ownership rejection; saved masked Results
  available after Cancel. Conventional smoke support remains 7974.
- Legacy open/compile, Run, Results, Resume/Cancel, migration/export.
- Runtime audit: 329 arm64 binaries, 103 ASE readers; no non-system absolute
  library dependencies or bundled pytest/PyInstaller/IPython/test modules.

Final validation: **2844 passed** in 324.12 seconds (full suite);
**97 passed** in GUI + packaging checks (94 GUI, 3 packaging).
The eight focused Data Prep/A-B-C/packaging checks passed.
The earlier post-layout run exposed the cancellation timing assertion described
above; final GUI and full-suite runs pass after its correction.
`git diff --check` is clean. ZIP checksum: `build/macos/QCBED-arm64.zip.sha256`.

Native screenshots now show readable linked plots; the overview and controls
below them are reached using the existing editor scroll area. Acceptance is
automated through the production GUI APIs; it is not a manual Finder/file-picker
usability certification. External project data/structure files remain external.
This is a local ad-hoc build, not a Developer-ID-signed/notarized distribution.

# Modern packaging closure — 2026-10-04

Final local bundle: `dist/QCBED.app`; refreshed archive: `dist/QCBED-arm64.zip`.
Build time: **2026-10-04T20:22:16.335496+00:00**. Version 0.1.0, native **arm64**, Python
3.11.16, ASE 3.29.0, PyInstaller 6.22.3 / hooks-contrib 2026.8.
Git revision `e4288b3aeab273df3ff02b118ccc4d8e086602e4` with an intentionally dirty working tree.
Source/spec fingerprint: `0778c65050e08a6b0c566f334ea6772e78fc22bef961113006bde7d64e194ff5`. Full version provenance is
in `build/macos/build-environment.json`; the Git revision alone is insufficient
for identifying this uncommitted build.

Changes are limited to package dependency/reader/resource collection, build
provenance, acceptance tooling and a GUI compatibility-report correction. The
latter now reports modern crystallographic export as unsupported without first
attempting unresolved pure compilation. No physical equations, optimizer algorithms,
checkpoint identities, masks or new parameter families changed.

## Measured final-bundle acceptance

- Final app explicitly re-signed ad hoc using `codesign --force --deep --sign -`;
  `codesign --verify --deep --strict` passed. No Developer ID or notarization.
- **329 Mach-O binaries**, all thin arm64; no non-system absolute
  dynamic-library dependencies. **103 ASE reader modules**
  collected; archive audit found no pytest, PyInstaller, IPython or package test
  modules. ASE extxyz and CIF loading (including spacegroup resources) passed.
- Native Cocoa launch via macOS Launch Services, working directory `/`: modern
  configuration and structure/projection views, original hkl fields, physical
  group table, B editing, signal modes, Save/Save As/reload passed.
- Actual real prepared gradient input: **7281 active points**, active weights
  **[1.0]**. First evaluation exactly reproduced YNORM **13447206.106343817** and
  objective **1.1896911365294378**.
- Short frozen ANNEAL→AMOEBA: cancelled after a completed evaluation; latest/best
  state, checkpoint and Results preserved; incompatible identity rejected;
  compatible Resume completed. Retained final, objective/parameter/YNORM histories,
  masked experimental/calculated/residual images and nine-beam diagnostics passed.
- Conventional smoke completed with **7974 active points** and distinct mode metadata.
- All **11 gradient evaluations** matched source execution exactly: physical vector,
  objective, YNORM, intermediate complex amplitudes, calculated/scaled intensities
  and support mask. This was a controlled smoke budget, not a new fitting campaign.
- Final legacy frozen run passed from `/private/tmp` under an empty environment
  (only standard PATH/HOME and test-only offscreen Qt): Open/compile, Save/reload,
  spawned Run, Resume, Cancel, Results/history/patterns and legacy import/export.
  Modern export remains blocked; no crystallographic provenance is fabricated.
- Fresh native app relaunch reopened the preceding session’s saved modern YAML
  with exact canonical equivalence.
- Source CLI smoke (`python -m qcbed.run --config ...`) completed after packaging
  changes. Source and editable workflows remain usable.

Acceptance artifacts and screenshots: `tmp/packaging-modern-closure/native-modern`,
`native-legacy`, `native-reopen`; source/frozen comparison:
`source-frozen-equivalence.json`. Earlier rehearsal folders are retained for
troubleshooting and are not the final build acceptance record.

## Source regression

- Prebuild focused modern backend + GUI: **151 passed**.
- GUI + packaging/export regression after the compatibility-report fix: **92 passed**.
- Prebuild full suite: **2838 passed**; new export regression passed separately.
- Final postbuild full suite: **2839 passed in 315.92 seconds**.
- Final `git diff --check`: clean; signature reverified after acceptance.

## Remaining limits

Ad-hoc local distribution only. No notarization, Gatekeeper/quarantined-download
acceptance, Intel/universal2 or other macOS-version certification. Literal Finder
mouse double-click and native file-picker clicking are not claimed: Launch Services
startup and the actual GUI file APIs are automated instead. The first sandboxed
Launch Services attempt returned `kLSNoExecutableErr`; the authorized normal macOS
launch succeeded, so this was not an executable/bundle defect.

Structure and experiment files remain external and document-relative; the reference
modern example requires its prepared gradient body. Other ASE formats have not been
individually frozen-tested. READY remains requested depth plus shared group B;
occupancy/absorption are read-only and other modern parameter families are deferred.
Build warnings about conda virtual `__osx` and optional platform imports remain in
the build log; tested workflows and the binary audit passed despite those warnings.

---

# Local Apple Silicon validation — 2026-09-27

Local artifact: `dist/QCBED.app`; symlink-preserving archive:
`dist/QCBED-arm64.zip` (103,408,024 bytes). Generated artifacts are ignored by git.

Environment: macOS 26.3.1, arm64, Python 3.11.16 (`qcbed-python`). PyInstaller
6.22.3 / hooks-contrib 2026.7; PySide6 6.11.2; NumPy 2.4.6; SciPy 1.17.1;
h5py 3.16.0; matplotlib 3.11.2; psutil 7.2.2; PyYAML 6.0.3. PyInstaller was
installed into a temporary build-tool directory, without changing conda packages.
The documented package extra provides the normal repeatable build command.

Measured checks:

- Thin arm64 main executable. Binary audit: **313 Mach-O binaries**, all arm64;
  no non-system absolute dynamic-library references (including no conda/Homebrew
  links). Archive inspection found no pytest, pandas, IPython, PyQt, PyInstaller
  developer modules or package test modules.
- `codesign --verify --deep --strict` passed. Signature is **ad hoc**, without
  Developer ID/team identity or notarization.
- Frozen app passed under an empty environment with standard PATH/HOME and
  offscreen Qt, starting from `/private/tmp`, with no conda activation/PYTHONPATH.
- Native Cocoa app launched through macOS Launch Services (`open -n -a`), with
  observed working directory `/`, and passed the same workflow using a small
  problem with an active depth parameter and a real optimizer loop.
- Open/reference compilation; edited YAML Save As and exact reload; spawned Run;
  progress/log delivery; checkpoint/report/output creation; compatible Resume;
  cooperative Cancel and control restoration; legacy import, FOR007 and bundle
  export all passed inside the frozen executable.
- Separate native app relaunch reopened the first session's saved YAML. The first
  relaunch harness attempt used invalid relative manifest paths and was rejected;
  the corrected absolute-path manifest passed. No application path fallback or
  semantics were changed to accommodate the test.
- Frozen versus normal Python execution: **exact equality** of final objective,
  driver evaluation count, active vector, raw intensities and scaled intensities.
- Normal Python full suite: **1,968 passed in 54.36 seconds**, including Qt tests.

Reports are local under `tmp/macos-acceptance/{frozen-1,native-1,native-reopen}`;
`build/macos/build-environment.json` and `binary-audit.json` record build/audit
metadata. The user's `examples/qcbed_config.yml` title remains `QCBED_test`.

Remaining manual acceptance: literal Finder mouse double-click and native Open /
experimental-file Browse dialog interaction have been requested from the user but
are not yet confirmed. Launch Services startup and file-handling APIs were tested;
those should not be misreported as manual interaction. No other Mac, older macOS,
quarantined download, Intel architecture, Developer ID signing or notarization was
validated.

Build warnings retained in the build log: conda's virtual `__osx` package has no
physical distribution; the SciPy hook references an unavailable optional
`scipy.special._cdflib`; Windows `user32` was not found. These did not prevent the
runtime checks above. No conda or Homebrew library paths remain in bundled Mach-O
load commands. A future dependency upgrade should repeat the same audit and
frozen acceptance checks rather than assuming these warnings remain harmless.

## Dynamic modern refinement acceptance — 2026-10-05

Full automated suite: **2923 passed** (225 existing Qt/Matplotlib deprecation warnings).
Initial focused adapter/Qt tests: **49 passed**; additional version-1 journal and
zero-start regressions are included in the final full suite. No convergence campaign ran.
The frozen editor smoke selected Debye–Waller plus a scan coordinate with depth
disabled, committed the real widgets, validated, saved and reloaded compact YAML:
**PASS**, frozen runtime, **0 optimizer evaluations**.

Native arm64 build: `dist/QCBED.app`; refreshed archive: `dist/QCBED-arm64.zip`.
Ad-hoc signing and `codesign --verify --deep --strict` passed. Not notarized.
Build timestamp: `2026-10-05T07:35:03.840755+00:00`.
Packaged source SHA-256: `94cbad1bfdf6f009e77007f3310c3663bcd17445d6b0c41bf379dcc36d92372b`.
Evidence: `tmp/dynamic_refinement_acceptance/` (generated, not source).
See [the parameter audit](../../docs/modern_refinement_catalogue.md) for the
explicitly unsupported physical mappings and optimizer-specific restrictions.

After the full suite, the final-output adapter was aligned with retained scan metadata.
The headerless output values are unchanged. All **14 dynamic-refinement tests**
passed after that change, followed by another frozen zero-evaluation PASS and
rebuild/ad-hoc deep/strict signature verification.

## GUI parameter-classification recovery — 2026-10-05

Recovered the existing dirty working tree without resetting or discarding any
earlier work. The previous full run had completed with 2,929 passes and three
outdated workflow-test failures; its five-test follow-up passed after adding the
required Simulation Setup validation. The earlier build and frozen usability
smoke also completed. No test, build, signing or refinement process remained
active when recovery was checked.

Preserved Basic / Physical, Advanced and Expert / Legacy classification and all
four C tabs: Simulation Setup, Refinement Targets, Optimizer, Execution / Output.
Fresh `0, 1` slice boundaries remain editable and incomplete. The UI explicitly
explains the internal-boundary requirement; validation still rejects one slab.
No numerical kernels, backend slicing rules, legacy sources or reference data
were changed during recovery. No expensive refinement campaign was rerun.

Validation from the recovered tree:

- Focused recovery checks: **36 passed** (45.98 s).
- Final boundary/control checks, including incomplete-draft recovery:
  **13 passed** (17.61 s).
- One full regression run after those checks: **2,933 passed** (414.55 s),
  including the existing numerical/reference comparisons; 257 Qt/Matplotlib
  `AA_UseHighDpiPixmaps` deprecation warnings, no failures.

Logs and the recovery checklist are in `tmp/usability_recovery/` (ignored local
artifacts). These results establish the tested tolerances; no new numerical-error
measurement or performance benchmark is claimed.

Rebuilt `dist/QCBED.app` as native arm64 at
`2026-10-05T15:55:01.111810+00:00`. Packaged source SHA-256:
`8d13c1cc1b17357415fa3827745636eff708008b6d92082dc5f76d3129822312`;
verified against the current source after acceptance. Ad-hoc signing and
`codesign --verify --deep --strict` passed. Refreshed
`dist/QCBED-arm64.zip` with `ditto` to preserve bundle symlinks.

Frozen usability acceptance: **PASS**, **zero optimizer evaluations**. Launched
from `/private/tmp` with an empty environment except standard HOME/PATH and
offscreen Qt, without conda activation or PYTHONPATH. The check prepared the real
pattern, loaded the CIF with explicit B, rejected the incomplete `0, 1` draft,
accepted the explicit FCC Al partition, selected physical targets, validated C,
saved/reloaded the project, and cancelled canonical initialization before any
evaluation. The resulting Optimizer screenshot was inspected. Evidence:
`tmp/usability_recovery/frozen/result.json`, `report.json`, `build.log` and
`sign.log` in the recovery directory.

The build retained the previously documented optional-import/conda warnings
(`__osx`, `scipy.special._cdflib`, Windows `user32`). Signing is ad hoc, not
Developer ID/notarized. This automated offscreen check does not establish manual
Finder/file-picker acceptance or compatibility on other machines/macOS versions.

## Modern GUI state and Run dashboard milestone — 2026-10-05

Scope: fresh-project requested thickness 1000 Å; preservation of independent
basis/orientation/thickness drafts; immediate configured hard bounds; and
D. Run & Results. No numerical kernels, legacy sources or reference data changed.
No expensive refinement campaigns ran.

The hkl reset was a presentation-state loss: Apply Site-Group Parameters built a
candidate from committed state and refreshed the widgets, dropping pending B
basis edits. It now commits those B fields with the site-group changes and
restores independent field drafts. Partial validation also preserves other
sections' drafts. The regression covers [2,0,0]/[0,2,0], B edits, Apply, Validate B,
tab switching and Save/reload, plus thickness drafts during CIF loading.

The hard-bound audit found specimen-specific example/test ranges, not universal
finite default intervals. Existing specific/category/relative rules are displayed
immediately, using the backend resolver. Missing canonical bounds are explicitly
identified. Unchanged inherited rules are not copied into per-target YAML entries;
user overrides persist. D1/D2/PLAMDA remain separate Expert settings.

Focused validation: initial compatibility set **42 passed**; state/action and
legacy Qt set **26 passed**; final focused set **15 passed** (33.50 s). The final
set includes six new state tests plus real small Run/Cancel/Open Results/Resume
button transitions. D uses the same controller and Results viewer, has no editable
configuration copy, and remains accessible while A/B/C are locked during a run.

One full regression suite after focused stabilization: **2,939 passed**, 225
Qt/Matplotlib deprecation warnings, no failures, in **423.09 s**. This includes
existing numerical/reference tests. No new numerical-error measurement or
performance benchmark is claimed for these frontend changes. Evidence is retained
under `tmp/gui_state_milestone/` (ignored local artifacts).

The first frozen smoke passed; screenshot inspection then found long external
paths clipping dashboard rows. The final presentation-only correction uses
concise successful input/output summaries, full-path tooltips and a scrollable
status area. **9 targeted state/action tests passed** afterward (25.57 s). The
single full-suite result above predates only this layout correction.

Final native arm64 rebuild: `2026-10-05T16:40:46.992238+00:00`.
Packaged source SHA-256:
`8f9c364be17b03732f616a3b35fe93fcac3c1f5361fbf5b3cf29175e2ffeedbe`.
The fingerprint matches the final working source. `dist/QCBED.app` was ad-hoc
signed; deep/strict signature verification passed. `dist/QCBED-arm64.zip` was
refreshed with symlink-preserving `ditto`.

Final frozen smoke: **PASS**, **zero optimizer evaluations**, launched from
`/private/tmp` without conda activation/PYTHONPATH. Verified fresh 1000 Å,
FCC Al g1/g2 preservation through Apply/Validate B, incomplete `0, 1` rejection,
four C sub-tabs, visible D actions/readiness, Save/reload, and canonical
initialization stopped before evaluation. Inspected `frozen-final/dashboard.png`:
status rows are scrollable and paths no longer clip successful summaries.
Evidence: `tmp/gui_state_milestone/report.json`, `full.log`,
`dashboard-final.log`, `build-final.log`, `sign-final.log`, and
`frozen-final/result.json`. Signing remains ad hoc, not notarized; automated
checks do not claim human Finder/file-picker acceptance.

## One-slice Simulation Setup correction — 2026-10-05

The explicit user request supersedes the earlier endpoint-only incomplete-draft
policy. `(0,1)` now selects one stored slice per proven simulation-cell period.
The strict transverse-periodicity proof, numerical kernels, Fourier conventions,
units and whole-slice depth accumulation remain unchanged. Legacy files and
reference data were read only. No expensive refinement campaign ran.

See [the full path/mesh audit](../../docs/one_slice_simulation.md). Required Slice
boundaries are always visible in Simulation Setup. Transverse simulation mesh
X/Y remain Advanced; compiled P/Q lengths, angle and lattice-direction sampling
steps are displayed when available. Non-[001] coverage verifies [100] on a
3×4×5 Å synthetic cell: transverse dimensions 4×5 Å, propagation period 3 Å,
mesh 8×6, array shape (6,8).

Focused acceptance: **75 passed**, 49 Qt/Matplotlib deprecation warnings,
**56.81 s**. Coverage includes one/two/four slices, atom ownership, absorptive and
elastic projected potentials/gratings, complex waves after every occurrence,
depth refinement/rebuild, compact Save/reload, one-slice checkpoint Resume and
rejection after changing to two slices. Resume checks cancelled before optimizer
evaluation; a few controlled forward evaluations checked depth changes/rebuild.

Measured numerical evidence:

| Comparison | Maximum absolute error | Scope/tolerance |
|---|---:|---|
| One-slice bridge vs explicit pre-sliced model, projected potential | 0 | complex128; with/without absorption |
| Same comparison, intermediate/final reciprocal waves | 0 | Each of 1/2/3 occurrences, two tilts; exact array equality |
| Actual compiled ONEIT MSEQ loop vs Python | 0 | Exact sequence/effective depth for 12 cases, including nextafter threshold and decimal accumulation |
| Saved Fortran elastic listing | 3.2861939374e-6 | Existing combined rtol=2e-7, atol=5e-7 |
| Saved Fortran absorptive listing | 5.3940549127e-8 | Existing combined rtol=3e-7, atol=5e-8 |
| Saved interaction constant | 3.9253656388e-9 | Existing atol=5e-9 |

The temporary Fortran harness compiles the actual sequence block extracted from
`private predecessor implementation`; no legacy file was edited. No historical saved one-slice
complex-wave fixture exists, so full Fortran-wave equivalence is not claimed.
Evidence: `tmp/one_slice_milestone/focused-final.log` and
`numerical-errors.json` (ignored local artifacts).

One full regression run after stabilization: **2,953 passed**, 225 existing
Qt/Matplotlib deprecation warnings, **426.18 s**, no failures or skips. The local
gfortran source-loop comparison executed in both focused and full runs.

Native arm64 rebuild completed at `2026-10-05T18:34:50.441794+00:00`.
Source SHA-256:
`4784adc29c7775b35409416616aace7932b95c50c1b6e526b40d539c3e1c19a0`,
verified against the final source after acceptance. Ad-hoc signing and
`codesign --verify --deep --strict` passed; the symlink-preserving
`dist/QCBED-arm64.zip` is refreshed.

Frozen smoke: **PASS**, **zero optimizer evaluations**, launched from
`/private/tmp` with standard HOME/PATH and offscreen Qt, no conda activation or
PYTHONPATH. It validated visible `(0,1)` with Advanced collapsed, fresh thickness
1000 Å, preserved g1/g2, physical target setup, D readiness, compact Save/reload
and initialization cancelled before evaluation. `frozen/simulation.png` was
inspected and shows Slice boundaries directly plus the compiled stored count 1.
Logs and complete provenance are in `tmp/one_slice_milestone/report.json`,
`full.log`, `build.log`, `sign.log`, and `frozen/result.json`.
The usual optional-import/conda build warnings remain; signing is ad hoc, not
notarized. No manual Finder/file-picker or other-machine acceptance is claimed.

## Output-directory workflow — 2026-10-05

Focused GUI/output/controller/project/packaging checks: **36 passed**.
Final display recheck: **16 passed**; lightweight directory recheck: **10 passed**.
The frozen app passed writable-directory validation, unrelated-file preservation,
separate Fresh Run/Resume readiness, and canonical initialization stopped before
evaluation (**0 optimizer evaluations**). Output, Run and Resume status are visible
at the top of D. No physics changes or refinement campaign were made.

QCBED.app was rebuilt, ad-hoc signed, and passed deep strict signature verification;
QCBED-arm64.zip was refreshed. Source fingerprint: `2de2b791cd20e3dd7635972c88cd3dbee391b6e462d12923908bf38078d9b487`.
Logs, screenshot and exact build metadata: `tmp/output_directory_milestone/`.

## C/D separation and Advanced execution audit — 2026-10-05

C contains only Simulation Setup, Refinement Targets and Optimizer. C validates
with a blank output directory; D alone owns writable-directory validation and
execution readiness. D Advanced exposes canonical serial/batched mode, batch size
and run safety limits. Modern grouped cadence remains fixed at every evaluation;
canonical worker-count configuration is unsupported. Existing YAML settings are retained.

**105 focused tests passed**, including 36 exact batch-invariance cases. Final
layout checks: **18 passed**. Batched/64 persistence and visible-region recheck:
**1 passed**, including preservation of non-default legacy cadence metadata.
The final frozen app exposes all four execution controls in its scrollable D
panel and passes save/reload and zero-evaluation initialization checks. No physics
kernels changed and no refinement campaign ran.

Rebuilt QCBED.app, ad-hoc signature and deep strict verification passed; ZIP refreshed.
Source fingerprint: `f292be54a4ae06381ce6ee4213793652a705214803b39ca53c8d7d026debb510`. Logs, screenshots and report:
`tmp/c_d_execution/` (final frozen artifacts under `frozen-verified/`).

## 2026-10-05 release cleanup and researcher documentation

Current closure: [release report](../../docs/RELEASE_REPORT.md) and
[provenance](RELEASE_PROVENANCE.json). Build source SHA-256
`c43e1e5e88442810599ce813a9b31d7970e6aec2d6839b6f3cadaa553f8cd106`. Full suite: 2964 passed / one obsolete C-only
readiness assertion failed in 410.38 s; corrected explicit-D initialization and
all 49 related follow-up tests passed. Full suite was run once, not relabeled as
an all-green run. No production source changed afterward.

Clean wheel GUI installation/dependency/import checks passed. Signed frozen smoke
opened the raw FCC Al template, prepared nine 38×38 beams, verified three C tabs
and D-only Advanced Batched/64, saved/reopened YAML, revalidated output and created
a compatible cancelled checkpoint with zero optimizer evaluations. Screenshots
were inspected. Deep/strict signature and ZIP integrity checks passed. No
refinement campaign or Finder/notarization acceptance was performed. See the
report for preserved references, cleanup and exact archive checksum.

## Final release acceptance closure — 2026-10-05T22:07:00.655442+00:00

**RELEASE ACCEPTANCE: PASS.** Complete suite: **2965 passed, 225 warnings,
422.22 s**, with no failures or skips. This all-green run supersedes the
previous stale readiness-assertion failure and focused-only follow-up record.

Current packaged smoke and rebuilt/signed frozen smoke passed with **zero optimizer
evaluations**. The extended harness verifies no-YAML new preparation, A/B/C/D,
editable one-stored-slice boundaries, site-group B/basis preservation, D validation
and Batched/64, compact single-YAML Save/reload, visible target bounds and explicit
optimizer overrides, fresh Resume messaging and initialization-only Open Results.
No production code, test source or physics changes were required.

The final build followed green regression/pre-build checks. Native arm64, ad-hoc
signing, deep/strict signature verification and ZIP integrity passed. Source hash:
`9a687a938ddd1ec94905ace540dc7bf2771fcd236d25b8baf54a0d9ffed4953d`. ZIP SHA-256:
`233f97843a17cbb74a05c60c3af7efa3d5008974d1c2af304c830eb585d8133f`.

See [the current release report](../../docs/RELEASE_REPORT.md) and
[accepted provenance](RELEASE_PROVENANCE.json) for versions and limitations.
Apple Silicon-only packaging, no notarization, instructional FCC Al results and
unsupported physical target families remain explicit. No expensive campaign or
Finder quarantine-distribution acceptance was performed. Evidence is retained in
ignored `tmp/release_acceptance/`.

## Projected elastic coefficient E_g milestone

Final complete suite: **2,984 passed, zero failures/skips, 225 existing warnings,
578.49 seconds**. The 19 focused E_g tests passed. Full validation used a temporary
byte-identical source/test copy with the private fixtures required by five source
provenance checks; private sources remain outside the shared tree and bundle.

Source GUI smoke and signed frozen GUI smoke: **PASS**, **zero optimizer evaluations**.
The frozen test selected original (200), slice [0,1], verified the Real-only IAM
default, selected Complex mode with separate bounds, saved/reopened the compact
YAML, and reconstructed E_g Results metadata from an initialization checkpoint.
Existing A/B/C/D, [0,1], B, g1/g2, D-only writable output, Batched/64, Resume and
Open Results checks also passed. Enable `projected_elastic: true` in the existing
usability smoke manifest to exercise these assertions.

The final app was rebuilt with the existing pinned packager dependency directory,
ad-hoc signed, and passed `codesign --verify --deep --strict` after smoke. All
**329 Mach-O binaries are arm64**. The refreshed symlink-preserving
`dist/QCBED-arm64.zip` passed `unzip -t`. Exact source and artifact hashes are in
[RELEASE_PROVENANCE.json](RELEASE_PROVENANCE.json). Local detailed logs are under
`tmp/projected_elastic_milestone/`; they are not release inputs.

Limitations: Apple Silicon macOS-only bundle, ad-hoc signed, **not notarized**.
FCC Al is instructional, not scientifically converged. E_g uses the source-scaled
projected convention (V·Å interpretation), not conventional 3D F_hkl. Space-group
constraint inference and conventional 3D mapping are deferred. (000) is read-only;
zero IAM has no magnitude-mode phase. Existing unsupported physical targets remain
listed in the [catalogue](../../docs/modern_refinement_catalogue.md).

Final E_g^(s) checklist closure: **PROJECTED ELASTIC COEFFICIENT REFINEMENT: PASS**.
The app was rebuilt after the exact-definition wording and component-details /
separate Complex-bound controls were added. The full suite remained at 2,984
passes, and the new frozen smoke passed at zero optimizer evaluations. Local
closure smoke records are in `tmp/projected_elastic_milestone/closure-frozen-smoke/`.
Ad-hoc signing, deep/strict verification and ZIP integrity were repeated for this
bundle; it is not notarized.

## Optional equality groups and termination closure

Final full suite: **3,007 passed**, zero failures/skips, 225 warnings,
420.77 seconds. Focused tests: **169 passed**. The isolated validation copy
used private fixtures for five provenance checks; the shared tree and bundle do
not include those fixtures. Source and frozen GUI smoke both **PASS** with zero
optimizer evaluations, including an optional linked Complex pair, compact YAML,
shared vector dimension, member Results rows, iteration/evaluation counts,
cancellation reason and Resume. Use `coefficient_groups: true` together with
`projected_elastic: true` in the usability manifest.

The app and ZIP were rebuilt. All **329 Mach-O binaries are arm64**. Ad-hoc signing,
`codesign --verify --deep --strict`, and `unzip -t` passed. This is **not notarized**.
See [release records](../../docs/RELEASE_REPORT.md) for limitations and the
[real-run audit](../../docs/25_iteration_audit.md) for the verified 12-iteration /
25-evaluation result. Exact current artifact hashes are in RELEASE_PROVENANCE.json.


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
in `RELEASE_PROVENANCE.json`.

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
See [the implementation audit](../../docs/shared_parameter_layout_audit.md) for the trace and
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
See [the state-layer audit](../../docs/sharing_gui_state_audit.md) for the reproduced defect,
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
