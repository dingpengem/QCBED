# Repository cleanup audit — 2026-10-06

This audit continued the existing working tree; it did not reset or discard earlier
scientific/GUI work. The inventory and link/syntax checks can be reproduced with
`python tools/release_audit.py`. The local detailed evidence is retained under
`tmp/release_cleanup/`; the distributable conclusions are in [the release report](RELEASE_REPORT.md).

| Area | Classification and disposition |
|---|---|
| `src/qcbed/` | Production: numerical kernels, data preparation, crystallography, canonical execution and Qt GUI. Keep compatibility entry points still referenced by tests/public callers. |
| `tests/` | Permanent scientific, preparation, crystallographic, dynamic-target, checkpoint and GUI regressions. No tests removed or consolidated. |
| `reference_data/` | Input-only: raw pattern and CIF. Derived snapshots moved byte-for-byte to tests/fixtures; caches and obsolete artifacts archived locally. |
| Private predecessor material | Preserved outside the shared repository; source listings and detailed inventories are not distributed. |
| `examples/` | Runnable demonstrations and historical diagnostic/replay cases. Keep the diagnostic scripts because reports/tests cite them. Replace only the main tutorial YAML with portable raw inputs. |
| `benchmarks/` | Permanent measured-performance tools. Keep; no new benchmark claims. |
| `tools/` | Maintained helper scripts, including the read-only release audit. |
| `docs/` | Current user/developer references plus historical milestone evidence. Archive the old long entry-point/specification texts; mark historical GUI reports. |
| `packaging/macos/` | Production build/signing instructions and frozen acceptance harnesses. Keep source and validation provenance. |
| `build/` | Generated build intermediates and isolated packager dependency installation. Work files may be regenerated; keep the packager environment and build fingerprint. |
| `dist/` | Generated app/ZIP. Keep the final signed app and archive; discard the redundant PyInstaller COLLECT copy after verification. |
| `tmp/` | Mixed generated validation, debugging, scientific campaign/replay results and project assets. Ignore in Git; retain scientific/provenance-referenced material. No blanket deletion. |

## Removed or relocated after reference inspection

- Removed reproducible root `pattern_dat_inspection.png`; its plotting script remains.
- Removed `examples/data/fcc_al.prepared.json`, which depended on a machine-local
  `tmp` campaign body. A backup remains in ignored local cleanup evidence. The raw
  tutorial now generates its own descriptor from real reference inputs.
- Relocated generated `examples/al_out.dat` and `examples/out.dat.backup-*` into
  ignored `tmp/release_cleanup/preserved_outputs/`; they are not tutorial inputs.
- `examples/out.dat` was already deleted in the incoming working tree; this cleanup
  did not recreate or independently delete it.
- Added ignores for run containers, backup files, generated example outputs and
  the reproducible inspection plot. Existing scientific reference files remain tracked.

## Small code/dependency cleanup

Removed unused imports from GUI modules and an unused document-validation local;
removed unreachable C Execution/All-fields branches and obsolete unused fresh-output
replacement text. Kept public reexports and compatibility wrappers with callers.
No numerical kernels, scientific defaults or optimizer equations were changed.
No large-module split was justified by this release task.

The GUI install extra now declares SciPy/ASE/image support used by real preparation
and structure loading. The environment YAML installs that extra. The frozen smoke
now opens the shipped raw template and reloads the saved prepared YAML, checking
physical values, slicing, Batched/64, D output validation and separate Resume.

The legacy user-editable `examples/qcbed_config.yml` also used absolute paths to
the developer checkout. They now resolve to the same inputs/outputs through
relative paths; the generated-body prerequisite remains explicit in its CLI guide.
The original local configuration is preserved in ignored cleanup evidence.

Permanent fixture purposes and SHA-256 hashes are recorded in
`tests/fixtures/legacy_reference/README.md` and its manifest. Generated binaries
and old scheduler/crash outputs were archived locally, not retained in source history.
Historical commits retain their original non-source artifacts; no history rewrite
is performed in this cleanup.


## Final input-only policy

The retained user-input files are exactly `reference_data/pattern.dat` and
`reference_data/Al_FCC.cif`. The accompanying README describes regeneration.
Ten historical files moved unchanged into `tests/fixtures/legacy_reference/`;
the small synthetic demonstration intensity moved into `tests/fixtures/synthetic/`.
Each retained fixture has a documented parser, numerical or provenance role.
Existing `tests/data/` geometry/catalogue fixtures remain test-only evidence.

Unused fit parameter/result snapshots, weights, potential dumps, empty outputs,
Linux executables, scheduler/crash records and ignored preparation caches were
archived under local ignored `tmp/final_cleanup/local_archive/`. None is needed by
the raw-input tutorial. Build/dist were already untracked and remain local.
Two old documentation screenshots containing machine-specific paths were archived;
current portable C-page screenshots remain intentional documentation assets.

No tests were removed. Five optional private evidence checks now take an external
maintainer root and retain their original scientific assertions. The final full
suite uses that external evidence directly, without restoring source in this tree.
