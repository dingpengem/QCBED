# QCBED Python

Python rewrite of a legacy Fortran multislice QCBED simulation program.

## Current status

The first legacy input-parser stage is implemented. It supports the formatted
`FOR007.DAT` title and run header, beam voltage and zone axis, mesh and slice
geometry, species and atom records, and the output and experimental-intensity
filenames. Later refinement, scan, structure-factor, and numerical propagation
sections remain deferred; the parser does not calculate wavelength or physics.

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