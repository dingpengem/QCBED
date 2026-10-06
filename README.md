# QCBED

QCBED is a fully Python, GUI-driven platform for quantitative CBED data preparation,
plane-wave multislice simulation, and flexible parameter refinement.

## Highlights

- **Automated experimental data preparation** — Integrated CBED preprocessing, disk geometry, reciprocal-lattice/beam identification, noise and variance handling, masking, and validation, with user review of geometry and physical inputs.
- **Plane-wave multislice QCBED simulation** — Dynamical diffraction with crystallographic orientation, projected simulation cells, configurable slicing, and physically indexed beams.
- **CPU-accelerated computation** — Batched CPU execution for modern refinement, plus process-based parallelism across independent tilts through the simulation APIs. No GPU backend is required or implemented.
- **Modular optimization framework** — ANNEAL, AMOEBA, PARAMIN, and SIMPARA, supported multi-stage workflows, dynamic refinement parameter vectors, hard bounds, checkpointing, and compatible-run Resume.
- **Visual refinement-parameter selection** — Select slice-resolved projected elastic coefficients $E_g^{(s)}$ (displayed by original crystallographic hkl), requested depth, site-group Debye–Waller B, enabled absorption coefficients, scan coordinates, or voltage without absorption; inspect starting values, hard bounds, and activation status. See the [supported target catalogue](docs/modern_refinement_catalogue.md).
- **Optional shared refinement parameters** — Tie compatible coefficient targets to common optimizer coordinates: one for Real-only or Magnitude-only, two for Complex. Independent refinement remains available; automatic Friedel constraints stay separate.
- **Interactive refinement monitoring and results** — Live updates, parameter histories, objective/convergence plots, per-beam diagnostics, and experimental-versus-simulated pattern matching; available latest, best, and retained states remain inspectable after cancellation.
- **Fully Python implementation** — The modern simulation, refinement, GUI, configuration, and numerical workflow runs in Python with NumPy/SciPy and Qt. Normal users need no Fortran compiler.
- **User-friendly with expert control** — Basic / Advanced / Expert organization offers a streamlined workflow and recommended defaults while retaining numerical, optimizer, and execution controls for experienced users.

## Start here

- [Installation](docs/INSTALLATION.md): Apple Silicon macOS app or Python source.
- [Quick Start](docs/QUICK_START.md): the normal A → B → C → D workflow.
- [User Manual](docs/user_manual.md): preparation, model, refinement, execution and results.
- [FCC Al example study](docs/EXAMPLE_FCC_AL.md): real repository inputs, reproducible selections, instructional settings.
- [Documentation index](docs/README.md): configuration, numerical conventions and developer references.

The GUI separates **A. Data Prep**, **B. Atomic Model**, **C. Main QCBED Parameters**
(Simulation Setup / Refinement Targets / Optimizer), and **D. Run & Results**.
Output selection and Advanced execution settings belong only to D. Fresh modern
projects retain Batched execution with batch size 64; opening Advanced is optional.

## Source installation

```sh
git clone https://github.com/dingpengem/QCBED.git
cd QCBED
conda create -n qcbed-python python=3.11 pip
conda activate qcbed-python
python -m pip install -e '.[gui,test]'
python -m qcbed.gui
```

See [installation](docs/INSTALLATION.md) for dependency profiles, app installation,
CLI usage, tests and packaging. Normal use does not require a Fortran compiler.

Final cleanup acceptance passed **3,027 tests**, including five optional external
maintainer-provenance checks. Minimal GUI/state, compact YAML round-trip and
packaged startup checks passed with zero optimizer evaluations. See the
[release record](docs/RELEASE_REPORT.md) for scope and packaging limitations.

## Scientific scope and limitations

There are no GPU dependencies. Modern grouped refinement requires the batched
backend and strict transverse periodicity. Occupancy is a fixed physical input;
modern occupancy, atom-coordinate, legacy reflection-multiplier and lattice-strain refinement
are not implemented. QCBED currently refines the slice-resolved projected elastic
Fourier coefficients used directly by the multislice calculation, not conventional
3D unit-cell structure factors. See the [target catalogue](docs/modern_refinement_catalogue.md).
Numerical mesh/slice convergence and a scientifically adequate refinement budget
remain the researcher's responsibility. The teaching example is not a converged fit.

## Validation and distribution

The automated suite includes complex-wave, Fortran-reference, data-preparation,
crystallography, optimizer, checkpoint and GUI regressions. Run it with
`python -m pytest -q`. The [release report](docs/RELEASE_REPORT.md) records the exact
current test/build results and limitations. Before private-source separation, the acceptance suite passed **2,965 tests**; the rebuilt signed-app
GUI and compact-YAML round-trip smoke passed with zero optimizer evaluations. Five optional source-integrity checks require external maintainer evidence; without
`QCBED_PRIVATE_REFERENCE_ROOT`, they are explicitly skipped in a shared checkout. See
[source provenance](docs/source_provenance.md).
[Historical measurements](docs/history/DEVELOPMENT_LOG.md)
retain their original configurations and dates.

The maintained app build is **macOS arm64 (Apple Silicon)**. It is **ad-hoc signed,
not notarized**. Current build products are `dist/QCBED.app` and
`dist/QCBED-arm64.zip`; obtain the app ZIP and checksum from the maintainer.
This repository does not assert that a hosted binary release has been published.
Intel-macOS, Windows and Linux frozen packages are not provided.

Scientific reference data remain available for numerical comparison. The predecessor
implementation and detailed source inventories are maintained privately and are
not distributed with this repository.

Developer entry points: [architecture](docs/ARCHITECTURE.md),
[project specification](PROJECT_SPEC.md), [repository audit](docs/REPOSITORY_AUDIT.md),
and [source provenance](docs/source_provenance.md). A project-wide redistribution
license has not yet been supplied; third-party bundle notices are retained.
