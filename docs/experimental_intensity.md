# Experimental intensity file parsing

## Source evidence

Readable executable statements are authoritative; no private predecessor implementation-only
declaration or fixed capacity is used by this parser.

| Source | Executable evidence |
| --- | --- |
| `private predecessor implementation` | READ MBOUT and ordered HIN/KIN pairs from FOR007; convert signed reciprocal indices to one-based PMAP sampling positions. |
| `private predecessor implementation` | READ EXPINT, OPEN unit 12, READ NUMP/NUML/GRADIENTOPT; compute counts; READ all XPINTS, all VARIANCE, then all WEIGHTS; sum weights and CLOSE. |
| `private predecessor implementation`, `private predecessor implementation` | Tilt row loop outside column loop; M increments once per scan position. |
| `private predecessor implementation`, `private predecessor implementation` | Each beam samples the same tilt list; optional GRADIENT call; beam-outer/tilt-inner flattening into THEORINTS. |
| `private predecessor implementation` | Later objective accesses XPINTS, VARIANCE, WEIGHTS and THEORINTS at the same flat index; YNORM scales theory. |

Searching the readable `.f` files found one experimental READ sequence and
no tilt-major input alternative. Both ONEIT and ONEITPARA use the same stored
experimental arrays. The supplied `private predecessor implementation` and `private predecessor implementation` are ELF binary
artifacts, not readable Fortran. The optional GRADIENT implementation cannot
be established from their filenames or related declarations.

## File structure and counts

EXPINT is a formatted filename record in FOR007 (currently `in.dat`). The
experimental file contains the scan header and all three data blocks:

```
NUMP NUML GRADIENTOPT
XPINTS:   MBOUT*NUML logical rows, NUMP real values each
VARIANCE: MBOUT*NUML logical rows, NUMP real values each
WEIGHTS:  MBOUT*NUML logical rows, NUMP real values each
```

There are no beam labels or beam counts in this header. MBOUT and beam order
come from the FOR007 HIN/KIN list. They are not independently specified by
the experimental file. The API therefore requires `beam_indices_hk` in that
external order, rather than guessing count from file length. A caller can
still supply incorrect labels of the correct count: this cannot be detected
from the experimental file alone. Duplicate H/K requests retain separate
beam columns and may have different supplied observations.

NUMP means scan columns, NUML means scan rows. RES2=NUMP*NUML is the point
count per beam, TOTALLINES=NUML*MBOUT is the logical row count per block, and
TOTALPOINTS=RES2*MBOUT is the value count per block. The source allocates the
same rectangular scan for every beam; no ragged or variable-length beam
dataset is read. GRADIENTOPT is preserved as an integer. The later source
branches only when it equals 1; this parser does not run that operation or
invent a restrictive enum for other integer values.

Each READ is list-directed and advancing. A logical row can continue across
physical records; leftover values on its final record are discarded. There
is no global flattening of all file tokens. The existing numeric lexer is
reused: decimal forms, leading/exponent signs, E/D exponents, comma/whitespace
separators, explicit repetitions and blank-record continuation are supported.
Finite values are required. Existing explicit restrictions on nulls, slash
termination, quoted numerics and unsupported syntax remain in force.

Incomplete blocks report filename, section, beam/row, physical record and
expected value. Nonblank records after all three blocks are rejected as an
explicit Python completeness check requested for this API; the legacy CLOSE
would ignore them. Trailing blank records are accepted. Ordinary surplus
tokens on a consumed row are discarded as the source permits. Missing a
physical line can cause continuation into subsequent lines, just as in
Fortran; with enough replacement data the parser cannot infer the author's
intended row boundaries. Counts describe READ values, not mandatory physical
line lengths.

## Ordering and public representation

Within each block, beam varies slowest, scan row next, scan column fastest.
For zero-based beam b, row r and column c:

```
t = r*NUMP + c
legacy_flat_index = b*(NUMP*NUML) + t + 1  # one-based Fortran position
python_values[t, b] = legacy_block[legacy_flat_index]
```

The original beam ordinal b+1 is preserved by column position; H/K themselves
are signed reciprocal indices, not one-based beam numbers. The driver adds
one to wrapped H/K only when constructing its PMAP sampling address. No such
offset is applied to the metadata in this data model.

This is the same t ordering produced by `generate_scan_tilts(scan)`. Duplicate
scan positions retain distinct entries. If `scan` is provided, the parser
checks both NUMP and NUML, not merely their product. Singleton experimental
dimensions can be represented because their READ loops are well-defined;
the existing scan generator separately rejects singleton interpolation.
Matching dimensions proves positional compatibility, not acquisition
calibration or equality to later refined scan coordinates.

```python
from qcbed import parse_for007_scan, generate_scan_tilts, parse_experimental_intensities

scan = parse_for007_scan("reference_data/FOR007.DAT")
tilts = generate_scan_tilts(scan)
# Exact ordered HIN/KIN list from the supplied FOR007 file:
beams = ((0, 0), (1, 1), (1, -1), (-1, -1), (-1, 1),
         (2, 0), (0, -2), (-2, 0), (0, 2))
data = parse_experimental_intensities("reference_data/in.dat",
                                    beam_indices_hk=beams, scan=scan)
assert data.intensities.shape == (len(tilts), len(beams))
```

`ExperimentalIntensityData` is frozen and owns bytes-backed immutable copies
of `intensities`, `variances`, and `weights`, each float64 `(tilt,beam)`.
It also stores `scan_columns`, `scan_rows`, `gradient_option`, and the ordered
`beam_indices_hk`. Constructor inputs are validated and snapshotted. Array
names describe source roles, not new weighting or uncertainty operations.
Only this model and `parse_experimental_intensities` are newly exported.

## Values, masks and unresolved preprocessing

**Source-proven:** XPINTS values are assigned directly by READ. There is no
clipping, normalization, rescaling, squaring or absolute-value conversion at
input. The only immediate computation in the reader is summing WEIGHTS for
reporting; arrays are not modified. No sentinel test, skipped point or separate
validity mask is present. Zero and negative intensities are accepted and
retained. Weights are read as real values, not Boolean flags. The objective
later multiplies terms by WEIGHTS and divides by VARIANCE; these operations
are documented as evidence only, not implemented here. In particular the
parser does not take square roots or reinterpret VARIANCE as sigma.

**Observed data:** all three experimental files contain signed intensities,
zero variances, and 0/1 weights. Zero weights can suppress contributions in
the later expressions, but points remain in the arrays and the source does
not short-circuit divisions at those positions. No values are removed or
treated as missing. Finite negative/zero variances or nonbinary/signed weights
are not rejected by this ingestion API because the source READ imposes no
such constraints; future objective validation may need stronger requirements.

**Unresolved:** physical detector units, external background correction,
normalization and other preprocessing are not proved by these READs. "Raw"
here means exactly the file values, not guaranteed raw detector counts.
The reference header has GRADIENTOPT=1. Both objective paths call GRADIENT
after calculating BEAM and before flattening theory/objective evaluation.
Without its readable implementation, changes it might make to shared arrays
and their eventual comparison semantics remain unknown. This does not change
the proven one-to-one mapping of the file positions to the initial scan.
No transformed-data compatibility or historical executable equivalence is
claimed. `pattern.dat`, `sigma.dat` and headerless `weights.dat` are not
alternative input paths established by the readable driver; their names do
not establish preprocessing rules.

## Reference validation

`in.dat`, `startingin.dat` and `bestin.dat` are tested independently. The latter
two have the matching experimental structure, but the current FOR007 files
explicitly name in.dat; their snapshot origin is not inferred from the name.
Each file has header `39 39 1`, 351 nonblank rows per data block and 39 values
per row. With the nine recorded FOR007 beams this is:

- 1521 scan positions, exactly matching `generate_scan_tilts`;
- 13689 intensity values, 13689 variance values and 13689 weight values;
- 41067 numeric data values in total, excluding the three header integers;
- 6100 negative intensities, 117 zero variances and 6140 zero weights;
- weight sum 7549, also printed by the driver's reference output.

Regression tests compare every value against independent physical-row indexing,
as well as hard-coded edge samples (intensities -490.594 and -18066;
variances 26233100 and 3608530). Synthetic tests separate all three blocks,
beam and tilt order, duplicates, signed/zero values, continuation, malformed
tokens, count mismatches, trailing records, scan compatibility and immutability.
There is no tilt-major input parser because no readable source path supports it.

Validation used exactly
`C:\Users\peng863\AppData\Local\miniconda3\envs\qcbed-python\python.exe`:
focused tests **39 passed in 2.65 s**; complete suite **1273 passed in 24.68 s**.
Every reference array matched the independent raw-row float64 values exactly.
`git diff --check` passed; protected source/reference files were unchanged.

No normalization, residuals, weighting, uncertainty propagation, objective,
refinement, optimization or forward-physics changes are part of this module.
