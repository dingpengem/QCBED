# Multislice QCBED Python Rewrite

## Objective

Rewrite the existing Fortran multislice QCBED program in Python.

Current priorities:

1. Reproduce the Fortran numerical results.
2. Create a readable and modular implementation.
3. Parallelize independent incident tilts on CPU.
4. Optimize FFT-based propagation.
5. Preserve a path for future GPU acceleration.

## Known Fortran entry point

- Main source file: `private predecessor implementation`
- Main program name: `REPLACE_WITH_PROGRAM_NAME`
- Primary input file: `REPLACE_WITH_INPUT_FILENAME`
- Primary output file: `REPLACE_WITH_OUTPUT_FILENAME`

## Current execution procedure

Describe the existing procedure here:

1. ...
2. ...
3. ...

## Known physical workflow

The expected high-level workflow is:

1. Read microscope and specimen parameters.
2. Read atomic structure or potential parameters.
3. Construct the projected potential.
4. Construct transmission functions.
5. Construct reciprocal-space propagators.
6. Generate the incident wave for each tilt.
7. Propagate through successive slices.
8. Calculate the exit-wave diffraction intensity.
9. Integrate or save QCBED intensities.

Codex must verify this workflow against the source code.

## Parallelization target

Each incident tilt is assumed to be independent.

The intended CPU design is:

- one process handles one or several tilts;
- each process runs sequential multislice propagation;
- large tilt-independent arrays should be read-only and shared or loaded
  once per worker;
- internal FFT/BLAS thread count should be limited to one.

Codex must verify that no mutable state couples different tilts before
parallelizing them.

## Numerical validation

Reference cases are stored under `reference_data/`.

Important comparisons include:

- projected potential;
- transmission functions;
- wavefunction after selected slices;
- exit wave;
- diffraction intensity;
- integrated disk intensities.

## Known conventions

Fill in everything currently known:

- Electron energy:
- Wavelength convention:
- Coordinate convention:
- Tilt convention:
- FFT sign:
- FFT normalization:
- Array order:
- Boundary conditions:
- Absorptive potential:
- Units: