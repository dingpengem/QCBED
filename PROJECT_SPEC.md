# QCBED project specification

QCBED is a CPU Python implementation of the legacy multislice calculation with
experimental preparation, crystallographic model construction, physical refinement
and a Qt researcher workflow. Scientific conventions remain governed by executable
reference tests and [source provenance](docs/source_provenance.md).

## Current contract

- A. Data Prep owns raw measurements, calibration, extraction, uncertainty, mask
  layers and the confirmed origin. Prepared inputs retain immutable provenance.
- B. Atomic Model owns the structure, original zone/basis indices, physical site
  groups, Debye–Waller B, fixed occupancy and explicit absorption.
- C. Main QCBED Parameters has exactly Simulation Setup, Refinement Targets and
  Optimizer. Validation of C is independent of output selection.
- D. Run & Results owns output validation, execution settings and Run/Resume/Cancel.
  A writable existing or creatable directory is sufficient. Fresh modern runs use
  unique managed children; existing unrelated files remain untouched. Resume checks
  a compatible checkpoint separately from directory validity.
- New GUI projects use Batched execution with batch size 64. Generic configuration
  defaults remain backward compatible. Canonical modern execution has no worker
  count control; the standalone process API is a separate scheduling interface.
- Slice boundaries [0, 1] describe one stored slice; [0, 0.5, 1] describe two.
  Strict transverse periodicity and whole-slice depth semantics remain enforced.
- Mesh X/Y samples the transverse simulation cell, not detector pixels.
- Schema-version-2 YAML stores editable inputs and references, not validation
  receipts, generated matrices or optimizer history. Raw projects must complete A.
- Supported targets and optimizer capabilities are defined by the
  [refinement catalogue](docs/modern_refinement_catalogue.md).

## Engineering constraints

Use NumPy/SciPy CPU kernels and preserve coordinates, units, FFT signs, array
orientation and normalization. No GPU dependencies or speculative numerical
changes belong in GUI/release work. Keep kernels separate from scheduling and I/O.
Legacy FOR007/pre-sliced compatibility remains supported within its documented
limits. Do not infer missing physical inputs or reinterpret legacy vacancy as
chemical occupancy. Preserve original scientific inputs and meaningful regressions.

The repository layout and module boundaries are in [Architecture](docs/ARCHITECTURE.md).
Follow [AGENTS.md](AGENTS.md) for numerical changes. Release changes require focused
checks, one complete regression suite, installed/frozen workflow validation and
signed artifact verification; report limitations rather than claiming convergence.
Normal release validation must not launch an expensive refinement campaign.

## Documentation and status

[Documentation index](docs/README.md), [release report](docs/RELEASE_REPORT.md),
and [packaging history](packaging/macos/VALIDATION.md) distinguish current instructions
from measured historical evidence. [The former specification](docs/history/PROJECT_SPEC_HISTORY.md)
is retained for provenance; its dated deferred-feature statements are not current
capability claims. The present macOS distribution is arm64 and ad-hoc signed,
not Developer ID signed or notarized. There is no claim of a public hosted release.
