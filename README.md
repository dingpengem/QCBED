# QCBED Python

Python rewrite of a legacy Fortran multislice QCBED simulation program.

## Current status

The project is currently in the source-analysis and numerical-validation
stage. The Python implementation is not yet validated for scientific use.

## Directories

- `private predecessor sources`: original Fortran implementation
- `reference_data/`: validation inputs and outputs
- `src/qcbed/`: Python package
- `tests/`: automated tests
- `benchmarks/`: performance benchmarks
- `docs/`: technical documentation

## Current target

- NumPy/SciPy serial reference implementation
- CPU process-based parallelism over independent incident tilts
- future optional GPU backend