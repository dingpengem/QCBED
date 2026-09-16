# QCBED Python

Python rewrite of a legacy Fortran multislice QCBED simulation program.

## Current status

The first legacy input-parser stage is implemented. It supports the formatted
`FOR007.DAT` title and run header, beam voltage and zone axis, mesh and slice
geometry, species and atom records, and the output and experimental-intensity
filenames. Later refinement, scan, structure-factor, and numerical propagation
sections remain deferred; the parser does not calculate wavelength or physics.

Standalone electron-wavelength and grid utilities are also implemented.
`electron_wavelength_angstrom(voltage_kv)` accepts positive finite voltage in
kV and returns angstroms using float64-compatible arithmetic; convert volts
explicitly. Exact executable equivalence remains unestablished because
default-kind Fortran literal rounding may depend on the original compiler flags.

## Coordinate and atom API

- `SliceGeometry.x_shift_fractional` and `y_shift_fractional` are fractional
  direct-lattice translations, used by legacy PROPG as `H*X0 + K*Y0`.
  Cell lengths and thickness remain in angstroms; lattice gamma is in degrees.
- Exported `AtomSite` records expose `species_index_one_based`, `x_fractional`,
  `y_fractional`, `starting_occupancy`, and `record_type` (`S` or `G`). Positions
  and occupancy are input values, before gas randomization or vacancy effects.
  Use `atom.species_index_zero_based(len(parsed.species))` before Python species
  lookup; the helper and parser reject out-of-range species identifiers.
- Grid arrays use `(mesh_y, mesh_x)` with x along columns and y along rows.
  Generic PMAP conversions support odd and even dimensions and preserve dtype.
  Generic reciprocal indices for size 5 are `[-2, -1, 0, 1, 2]` centered and
  `[0, 1, 2, -2, -1]` unshifted. Legacy `oshift` requires both dimensions even.

The former angstrom shift names and ambiguous atom field names have been
removed; no compatibility aliases are provided. `SliceGeometry` and `AtomSite`
remain available from the package root. See [source provenance](docs/source_provenance.md).

## Legacy parser syntax

`parse_for007` follows each READ through the experimental-intensity filename.
It supports record continuation, comma/whitespace separators, E/D exponents,
explicit repetitions, and quoted list-directed S/G values. Numeric inputs are
finite; integer fields are parsed strictly, with field-specific bounds and
mask/flag domains. Errors include filename, section, record, and offending input.

Raw title, absorption and filename records are retained in `*_raw` fields;
`Species.label_raw` retains the full record while `label` follows FORMAT(A2).
Quotes in these formatted fields are literal. `AtomSite.record_type_raw`
preserves the list-directed token; `record_type` is its case-sensitive S/G value.
See [parser contract](docs/legacy_parser.md) for unsupported forms and the
limits of compatibility with the Fortran runtime.

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
