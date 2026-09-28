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
