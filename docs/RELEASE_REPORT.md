# Final release acceptance — current QCBED milestone

**RELEASE ACCEPTANCE: PASS**

This closure verifies the final working tree after the earlier readiness-test
correction. It adds no features or physics changes and runs no refinement campaign.
The only executable-file change in this closure extends the opt-in packaging smoke
assertions; production Python and test sources remain unchanged.

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
