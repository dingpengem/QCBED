# QCBED.app — native Apple Silicon build

This is an additive PyInstaller **onedir, windowed arm64** bundle. Normal editable
installation, `python -m qcbed.gui`, and `python -m qcbed.run --config … [--resume]`
remain unchanged. No runtime interpreter or conda activation is needed to use the
built app. PyInstaller is only in the optional `package` extra.

## Build

Use native arm64 Python on the target macOS generation, not Rosetta. From a
checkout, with the qcbed-python environment selected:

```sh
python -m pip install -e '.[package]'
python packaging/macos/build.py
```

The script anchors every build path to this checkout, so it can also be invoked
by absolute path from another working directory. It runs `--clean --noconfirm`,
replacing generated QCBED build artifacts only. Outputs are `dist/QCBED.app` and
`build/macos/build-environment.json`; both locations are gitignored. No test suites,
reference datasets, conda activation scripts or developer tools belong in the app.
The spec excludes unrelated optional scientific packages, alternative Qt bindings,
and testing tools. PyInstaller hooks collect Qt Cocoa/platform plugins, NumPy/BLAS,
SciPy, HDF5/h5py, matplotlib Agg/QtAgg/fonts, psutil, PyYAML and transitive runtime libraries.
Distribution metadata and available third-party license notices are included.

The package extra pins the tested PyInstaller release. The environment report
records the actual Python/platform and main package versions for each build.
For repeatable release builds, retain that report and an explicit conda environment
export or pip lock alongside the release; the project's permissive runtime version
ranges are not a byte-for-byte lock. Build on the oldest macOS you intend to
support and validate there; this build does not establish compatibility with older
macOS versions or Intel Macs.

Archive with Apple's `ditto` to preserve bundle symlinks and permissions:

```sh
ditto -c -k --sequesterRsrc --keepParent dist/QCBED.app dist/QCBED-arm64.zip
```

Do not commit app/zip/build outputs. Distribute the entire `.app`, never just its
`Contents/MacOS/QCBED` executable. Users can unzip it, copy it to Applications (or a
writable user folder), and double-click QCBED. Use Open to select an existing YAML.
Experimental files remain external; moving the app does not change document-relative
path semantics. Save As preserves the canonical serializer's resolved paths.

## Frozen startup and workers

`launcher.py` calls `multiprocessing.freeze_support()` before Qt, qcbed, argparse
or smoke-test dispatch. PyInstaller's override dispatches spawned worker and
resource-tracker processes to the correct code instead of relaunching the GUI.
No process implementation, checkpoint schema, cancellation logic, numerical kernel
or configuration semantics are replaced. In-app Run/Resume still uses canonical
immutable snapshots; Cancel still uses the existing safe interruption boundary.

`argv_emulation=False` avoids a second event loop conflicting with Qt. This stage
uses the GUI's Open dialog; Finder YAML file associations/Open Document events
are not registered. Double-clicking the app itself needs no arguments or shell.

## Local, signed, and notarized distribution

The initial build has an **ad-hoc signature** applied by PyInstaller, required for
arm64 code. It has no Developer ID certificate/team identity and is not notarized;
this is the local/unsigned-distribution tier, not trusted public distribution.
Gatekeeper may block a downloaded/quarantined copy. Do not advertise this artifact
as signed/notarized or advise globally disabling Gatekeeper.

A later Developer ID build supplies a real signing identity to PyInstaller and
appropriate reviewed entitlements/hardened-runtime settings. Signing establishes
publisher identity; notarization is a separate Apple service submission, followed
by stapling and Gatekeeper validation of the distributed archive/app. This stage
does not request credentials, access signing keys, submit to Apple, or publish.
Use a separate release procedure after the local application is proven.

## Acceptance testing

The optional `--smoke-test /absolute/manifest.json` switch runs packaged Qt acceptance
checks using explicit external inputs. It is not used by normal startup and does
not import pytest or execute external scripts. The manifest has absolute paths:

```json
{
  "reference_yaml": "/project/examples/qcbed_config.yml",
  "run_yaml": "/scratch/small-fast-canonical-test.yml",
  "for007": "/project/reference_data/FOR007.DAT",
  "experimental": "/project/reference_data/in.dat",
  "output": "/scratch/NEW-acceptance-directory"
}
```

Output must not exist. Use a tiny valid numerical problem for `run_yaml`; the
harness replaces only run output destinations. It opens/compiles the reference,
edits/saves/reloads YAML, runs and resumes through the GUI's real spawned worker,
checks progress/artifacts, cancels another run, and imports/exports legacy files.
It writes `acceptance.json`, GUI logs, screenshot, runtime probe, and run artifacts.
A second launch with `reopen_only: true` and `reference_yaml` pointing to the first
saved YAML validates reopening in a new app process. Native file chooser interaction
and literal mouse double-click still require manual UI acceptance when automation
is unavailable.

```sh
open -n -a /absolute/dist/QCBED.app --args --smoke-test /absolute/manifest.json
codesign --verify --deep --strict --verbose=2 dist/QCBED.app
file dist/QCBED.app/Contents/MacOS/QCBED
python -m pytest
```

For launch independence, also run the executable from an unrelated CWD with an
empty environment except standard macOS PATH/HOME (and a test-only Qt platform
selection if needed). Check all acceptance results, app relaunch, new output paths,
checkpoint/report generation, and a real small numerical result, not merely that
the dock icon appears.

References: [PyInstaller multiprocessing pitfalls](https://pyinstaller.org/en/stable/common-issues-and-pitfalls.html#multi-processing)
and [macOS build/distribution notes](https://pyinstaller.org/en/stable/feature-notes.html#macos-multi-arch-support).

See [the measured local validation record](VALIDATION.md) for the actual build
versions, frozen-run checks, numerical comparison, signature audit and remaining
manual acceptance items.
