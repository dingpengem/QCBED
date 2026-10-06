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
  "for007": "/project/tests/fixtures/legacy_reference/FOR007.DAT",
  "experimental": "/project/tests/fixtures/legacy_reference/in.dat",
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

## Modern crystallography closure

The package extra now includes ASE. The same spec explicitly collects ASE's
runtime-selected readers, its spacegroup data, crystallography modules and modern
worker/Results adapters. extxyz and CIF are covered by frozen acceptance; support
for other ASE formats is not established by those tests. Structure files and
experimental datasets remain external document inputs, never bundle-relative.
Open `examples/fcc_al_grouped_gui.yml` with its matching prepared gradient body;
see [modern workflow and capability limits](../../docs/modern_grouped_gui.md).

For an isolated build-tool install (without changing conda packages):

```sh
python -m pip install --target build/macos/packager 'PyInstaller==6.22.3'
PYTHONPATH="$PWD/build/macos/packager" python packaging/macos/build.py
codesign --force --deep --sign - dist/QCBED.app
codesign --verify --deep --strict --verbose=2 dist/QCBED.app
```

The environment report includes UTC build time, Git revision/dirty status, Python,
architecture, dependency versions and a source/spec SHA-256 fingerprint. Retain it
with the test report; the revision alone does not identify an uncommitted build.

The existing smoke entry also accepts `modern: true`, `reference_yaml`, `run_yaml`,
`conventional_yaml` and a new absolute `output` directory. Modern acceptance uses
real prepared gradient/conventional inputs with a small explicit optimizer budget.
It exercises structure/projection displays, dynamic extxyz/CIF readers, B editing,
Save/Save As/reload, incompatible resume rejection, spawned ANNEAL/AMOEBA, cancel
after a completed evaluation, retained Results, history and masked patterns. The
legacy manifest/harness remains supported. Neither harness is normal startup code.
Modern export must report unsupported provenance without trying to fabricate it.

Ad-hoc signing provides local bundle integrity, not Developer ID identification or
notarization. Literal Finder clicks/native file-picker interaction and deployment
on other macOS versions remain separate manual/distribution checks.

## Compact-schema initialization-only acceptance

For a lightweight check, use a JSON smoke manifest with `modern: true`,
`compact_only: true`, absolute `reference_yaml`, `interactive_yaml`, `verbose_yaml`,
`legacy_yaml`, and a new absolute `output` directory. The interactive input must
have recorded manual exclusions and origin/overview metadata. Launch with:

```sh
open -n -a "$PWD/dist/QCBED.app" --args --smoke-test /absolute/manifest.json
```

This checks compact load, external assets, A/B/C displays, Save/Save As/reload,
and older input compatibility. It calls the canonical runner but intercepts stage
dispatch before optimizer entry. The runner writes a genuine initialization
checkpoint, records the deliberate smoke exception as `failed`, and releases its
lock. This is a zero-evaluation setup test, not a completed or resumable refinement
claim. The interception exists only in the opt-in packaging harness.

## Live Data Prep smoke

Use `{"modern": true, "live_only": true, "output": "/absolute/new/directory"}`
with the existing `--smoke-test` launcher. This constructs a fresh synthetic raw
pattern and drives the real spawned preparation worker, not a prepared fixture.
It checks aperture geometry preview, manual geometry, ROI confirmation, live
beam output, masking and Validate A. No refinement executes. Detector/noise
scientific tests and real-pattern failure checks run in the source test suite.


For real live A/B/C initialization acceptance, the opt-in smoke manifest supports
`modern: true`, `live_abc_only: true`, absolute `pattern`, `structure`, and new
`output` paths, plus explicit `radius_px` and `symmetry_roi`. It builds a fresh
preparation, confirms/masks A, applies FCC Al B/C controls, validates all stages,
checks replay, and deliberately stops before optimizer entry (zero evaluations).
No historical prepared manifest is accepted as input to this check.

For final new-user interaction acceptance add `interactive: true` and explicit
`diameter_endpoints` (two detector x/y points) to the real A/B/C smoke manifest.
The harness drives diameter/ROI Accept controls, the spawned preparation worker,
linked beam displays and weight-only masking, then a separate zero-evaluation
initialization worker. This is automated Qt event acceptance, not a refinement run.

## Accepted milestone artifacts

The parameter-classification usability smoke uses `usability: true`,
`require_frozen: true`, absolute `pattern`, `structure` and new `output` paths,
plus `diameter_endpoints` and `symmetry_roi` in detector x/y pixels. Run it with
the existing `--smoke-test /absolute/manifest.json` entry point. It checks live
A preparation, explicit B input, the three C tabs plus D-only Advanced execution settings, visible and valid
`0, 1` one-slice boundaries, physical refinement targets, Save/reload, and canonical
initialization cancelled before any optimizer evaluation. This is automated
widget acceptance, not a refinement campaign or manual Finder check.

See [release closure](../../docs/reconstruction_milestone.md) and
[the accepted build provenance](RELEASE_PROVENANCE.json). The provenance records
the original dirty build base and exact source/artifact fingerprints; documentation
closure does not require a rebuild. Distribution remains ad-hoc signed, not notarized.

The usability smoke also checks the fresh 1000 Å thickness, the exact FCC Al
basis-edit → site-group Apply → Validate B sequence, and D's visible Run/Resume/
Cancel states. It saves `dashboard.png` as well as the optimizer view. Actual
small Run/Cancel/Resume button transitions and retained Results are exercised by
`tests/test_gui_modern_qt.py`; the packaged smoke still stops before evaluation.
The one-slice check also saves `simulation.png` showing required Slice boundaries
outside the collapsed Advanced mesh/cutoff panel.
## Current researcher release

Use the [installation guide](../../docs/INSTALLATION.md) and
[current release report](../../docs/RELEASE_REPORT.md) as the distribution entry
points. Earlier validation sections are historical measurements, not instructions
to use old machine-local campaign YAMLs.

The usability manifest optionally accepts `example_project` pointing to the raw
`examples/fcc_al.yml` and `study_depth: 1000`. The harness copies the template to
its new output workspace, prepares the real pattern, validates A/B/C, verifies
D-only Advanced Batched/64 visibility, saves and reopens the project, revalidates D,
and initializes a cancellable checkpoint with zero optimizer evaluations. It must
use a new output workspace per invocation. The final app/ZIP is signed and verified
after the build; no Finder quarantine/notarization acceptance is implied.
