# QCBED

QCBED is a CPU implementation of multislice quantitative convergent-beam electron
diffraction simulation and refinement, with a modern Qt workflow for TEM researchers.
It combines experimental-data preparation, an explicit crystallographic model,
physical refinement targets, resumable execution, and recorded-results inspection.
The original Fortran sources and scientific reference data remain available for
numerical comparison and provenance.

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

## Capabilities and boundaries

- Raw CBED preparation; optional detector correction; explicit uncertainty and mask layers.
- CIF/ASE structure input, zone-axis recasting, physical beam labels and site groups.
- Thickness, Debye–Waller B, supported absorption, voltage and scan-coordinate targets;
  explicit hard bounds; existing ANNEAL, AMOEBA, PARAMIN and SIMPARA engines.
- Batched CPU tilt evaluation, cooperative cancellation, compatible-checkpoint Resume,
  and recorded parameter/intensity comparisons.

There are no GPU dependencies. Modern grouped refinement requires the batched
backend and strict transverse periodicity. Occupancy is a fixed physical input;
modern occupancy, atom-coordinate, reflection-factor and lattice-strain refinement
are not implemented. See the [target catalogue](docs/modern_refinement_catalogue.md).
Numerical mesh/slice convergence and a scientifically adequate refinement budget
remain the researcher's responsibility. The teaching example is not a converged fit.

## Validation and distribution

The automated suite includes complex-wave, Fortran-reference, data-preparation,
crystallography, optimizer, checkpoint and GUI regressions. Run it with
`python -m pytest -q`. The [release report](docs/RELEASE_REPORT.md) records the exact
current test/build results and limitations. The final acceptance suite passed **2,965 tests**; the rebuilt signed-app
GUI and compact-YAML round-trip smoke passed with zero optimizer evaluations.
[Historical measurements](docs/history/DEVELOPMENT_LOG.md)
retain their original configurations and dates.

The maintained app build is **macOS arm64 (Apple Silicon)**. It is **ad-hoc signed,
not notarized**. Current build products are `dist/QCBED.app` and
`dist/QCBED-arm64.zip`; obtain the app ZIP and checksum from the maintainer.
This repository does not assert that a hosted binary release has been published.
Intel-macOS, Windows and Linux frozen packages are not provided.

Developer entry points: [architecture](docs/ARCHITECTURE.md),
[project specification](PROJECT_SPEC.md), [repository audit](docs/REPOSITORY_AUDIT.md),
and [source provenance](docs/source_provenance.md). A project-wide redistribution
license has not yet been supplied; third-party bundle notices are retained.
