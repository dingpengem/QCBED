# Installation

[Documentation index](README.md) · [Quick Start](QUICK_START.md)

## Recommended: the macOS application

The maintained package is native **Apple Silicon / arm64 macOS**. The recorded
release was built and checked on macOS 26.3.1. That is the tested platform, not a
claim of compatibility with every earlier macOS version. There is no Intel app
or Windows/Linux binary in this release.

1. Obtain `QCBED-arm64.zip` and its SHA-256 checksum from the maintainer. In a source
   checkout, the build creates it in `dist/`; the task producing this release does
   not publish to a download server. Obtain the matching source/reference-data
   archive as well to follow the FCC Al study: those inputs are external to the app.
2. Compare the ZIP checksum, for example with `shasum -a 256 QCBED-arm64.zip`.
3. Double-click the ZIP to extract `QCBED.app`. Move it to Applications or another
   writable application folder and open it.
4. This app is **ad-hoc signed, not signed with an identified Developer ID and not
   notarized by Apple**. macOS may report that the developer cannot be verified or
   that Apple cannot check the app. For a copy you have verified and trust, first
   attempt to open it, then use **System Settings → Privacy & Security → Open Anyway**
   if offered, and confirm the named app. Follow [Apple's current guidance](https://support.apple.com/en-us/102445).
   Do not disable Gatekeeper globally. If macOS reports damage/malware rather than
   an unidentified developer, stop and obtain a verified copy instead of bypassing it.
5. The blank window opens at **A. Data Prep**. Follow the [Quick Start](QUICK_START.md).

The bundle supplies its Python/Qt/numerical dependencies. **No Python installation
or Fortran compiler is needed for normal app use.** Keep project files and results
outside the app bundle. Updating the app does not bundle or migrate your datasets.

## Source/developer installation

Install Git and a native Python distribution (Miniforge/Conda is suitable on
Apple Silicon). Use an authorized clone URL if the repository is access-controlled.
From a fresh terminal:

```sh
git clone https://github.com/dingpengem/QCBED.git
cd QCBED
conda create -n qcbed-python python=3.11 pip
conda activate qcbed-python
python -m pip install -e '.[gui,test]'
python -m qcbed.gui
```

The `gui` extra includes Qt, YAML, plotting, SciPy preparation, ASE structure I/O
and image reading; it is the complete modern researcher workflow. `test` adds
pytest. `pip install -e .` alone installs the NumPy core, **not** the full GUI.
Smaller optional profiles are `yaml`, `crystallography` and `visualization`.
Dependencies are declared in [pyproject.toml](../pyproject.toml). Lower bounds are
not a lockfile: exact tested build versions are recorded in the release report.

Alternatively, the repository's full development environment is:

```sh
conda env create -f environment.yml
conda activate qcbed-python
python -m qcbed.gui
```

Do not run both environment-creation commands with the same environment name.
The environment file installs the project editable with GUI support. There are
no CUDA, CuPy, JAX or other GPU dependencies.

### GUI, command-line interface and tests

```sh
python -m qcbed.gui examples/fcc_al.yml
python -m qcbed.run --help
python -m pytest -q
```

`examples/fcc_al.yml` starts interactive preparation from a raw CBED image; it is
not yet a prepared CLI run. After saving a fully prepared project, CLI syntax is:

```sh
python -m qcbed.run --config saved_prepared_project.yml
python -m qcbed.run --config saved_prepared_project.yml --resume
```

The CLI uses the exact configured run directory and expects an unused path for a
fresh run. The GUI instead creates a unique managed subfolder inside its reusable
output directory. Do not treat a GUI output container as a CLI checkpoint directory.
`--resume` requires the matching checkpoint and configuration/data identity.
A safety-budget interruption can return a nonzero exit status with usable partial
results. These commands run the configured budget; they are not validation-only.

For automated GUI tests without a display:

```sh
QT_QPA_PLATFORM=offscreen python -m pytest -q
```

Reference files are included in the checkout and must not be regenerated to make
a failing scientific regression pass. A Fortran compiler is only needed for
optional developer work that explicitly rebuilds legacy executables; the normal
Python regression suite uses recorded references.

### Build the app (native macOS arm64 only)

```sh
python -m pip install -e '.[package,test]'
python packaging/macos/build.py
codesign --force --deep --sign - dist/QCBED.app
codesign --verify --deep --strict --verbose=2 dist/QCBED.app
ditto -c -k --sequesterRsrc --keepParent dist/QCBED.app dist/QCBED-arm64.zip
shasum -a 256 dist/QCBED-arm64.zip
```

The build checks architecture, freezes multiprocessing before importing Qt, and
writes exact versions/source fingerprint to `build/macos/build-environment.json`.
The ZIP is made **after** signing. These commands do not notarize or publish it.
See [packaging instructions](../packaging/macos/README.md) for release smoke tests.

## Installation troubleshooting

- **Missing PySide6 / ASE / SciPy:** activate the intended environment, check
  `python -m pip --version`, then install `'.[gui,test]'` from the repository root.
- **Wrong architecture:** use native arm64 Python to build; do not mix an Intel
  environment with the arm64 app specification.
- **Qt display error over SSH/CI:** use the offscreen test command; for interactive
  use, launch in a desktop session.
- **Missing input/cache after moving a project:** a YAML is an entry point, not a
  self-contained data archive. Restore its referenced assets or prepare again.
- **App starts but an example input is missing:** obtain the matching repository
  reference data; scientific inputs are deliberately not embedded in the app.
