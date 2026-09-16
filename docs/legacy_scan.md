# Initial legacy scan definition

## Evidence and scope

The implementation uses readable executable statements, not private predecessor implementation
declarations. The latter remains a related-version declaration reference.

| Source | Evidence |
| --- | --- |
| `private predecessor implementation` | EXPINT filename; unit 12 READ of NUMP,NUML,GRADIENTOPT; RES2=NUMP*NUML; experimental records of NUMP points. |
| `private predecessor implementation` | One list-directed READ of AY,BEE,NOTHING,CEE,DEE,NAUGHT; scan-basis normalization and its inverse. |
| `private predecessor implementation` | Separate READs of A(1:2), B(1:2), C(1:2), followed by refinement masks. |
| `private predecessor implementation` | GAMMA for scan normalization is the mean reciprocal-cell angle, derived from the slice cells. |
| `private predecessor implementation`, `private predecessor implementation` | Scale second corner components by NORMFACTOR; divide corners by per-type dilation; interpolate with NUMP-1/NUML-1 denominators, columns innermost. |
| `private predecessor implementation`, `private predecessor implementation` | Undo normalization on component 2, transform into cell H/K, pass each type's TILT to PROPG. |
| `private predecessor implementation`, `private predecessor implementation` | All beams use the same tilt index; store BEAM(beam,tilt), then flatten with beam outermost. |
| `private predecessor implementation` | UPDPARS can replace A/B/C from refinement parameters; parsed initial corners are not necessarily later fitted corners. |
| `private predecessor implementation` and later objective calls | ONEITPARA takes a parameter vector; scan state is shared, not a tilt-array argument. Both objective variants construct their own TILTARRAY before PROPG calls. |

This task reads the initial geometry only. No new masks, objective options,
intensities, refinement parameters or optimizer behavior are implemented.
The existing prefix reader still validates its previously supported records.

## Records and units

After EXPINT, FOR007 contains one READ of six numeric values:

```
AY BEE NOTHING CEE DEE NAUGHT
A1 A2
B1 B2
C1 C2
```

The six-value READ may span records (as in both supplied reference inputs).
Each corner has its own READ. Extra ordinary values on the last record of
one READ do not become values for the next READ. Existing list-directed
syntax and contextual errors are reused, including commas, D exponents and
explicit repetitions. The parser stops after C, before the three MASKT
records and everything that follows.

The two in-plane scan basis vectors have cell coefficients `(AY,BEE)` and
`(CEE,DEE)`. NOTHING and NAUGHT are consumed and retained as the respective
third input components, but have no role in any found executable scan
transformation. No meaning such as out-of-plane tilt is assigned to them.
All six input components must be finite; neither in-plane basis vector may
be zero, since that makes source normalization singular. Collinear nonzero
basis vectors are allowed: the source does not invert this transformation.

Each A/B/C pair is dimensionless coefficients `(u,v)` of the supplied scan
basis. They are not angles or inverse-angstrom Cartesian coordinates.
There is no implicit H/K transpose: the cell transform is

```
h = AY*u + CEE*v
k = BEE*u + DEE*v
```

The first list-directed READ of the experimental file supplies NUMP (scan
columns), NUML (scan rows), and GRADIENTOPT. Only that header is read; the
third integer is consumed without interpreting its objective behavior.
Counts are strict positive integers and independent of the real-space wave
mesh MESHX/MESHY. No count is guessed from file length, corners, output images
or mesh size. Header continuation is supported without consuming subsequent
intensity records. Missing headers produce contextual errors.

## Generation and ordering

For zero-based row r and column c, the scan-basis coordinates are

```
(u,v) = A + r*(C-A)/(NUML-1) + c*(B-A)/(NUMP-1)
tilt_index = r*NUMP + c
```

A is the first point; B is the last point of the first row; C is the first
point of the last row; B+C-A is the final point. A sheared parallelogram is
allowed. Columns vary fastest. Neither scan dimension need align with H or K;
both components can change along either scan direction. Endpoints are
included algebraically for nonsingleton dimensions, subject to float64
rounding. Zero is not guaranteed even if coordinate ranges cross it.
Repeated corners or collinear steps can create duplicate tilts; no
deduplication, clipping, sorting or recentering is performed.

**Singleton dimensions:** the readable source divides by NUMP-1 and NUML-1
unconditionally, including when the corresponding corner difference is zero.
No finite singleton rule can be recovered. Parsing/representation preserves
positive counts of one, but generation raises a clear ValueError if either
dimension is one. Choosing A, a midpoint, or an endpoint would be an extension,
not source-proven behavior. H-only/K-only scans can still be represented with
both counts greater than one and a collapsed direction; duplicates remain.

The source first scales the second coordinate by NORMFACTOR and divides
both coordinates by APPLIEDDILATION(type). After interpolation it multiplies
the second coordinate by INVNORMFACTOR, then performs the linear H/K basis
transform. For finite nonzero source normalization these operations give
`(h0/d, k0/d)`, where h0,k0 are obtained directly from the input corners and
basis above. This module algebraically cancels normalization and returns
**undilated** h0,k0, leaving the existing simulation API to divide by d once.
No gamma, lattice lengths, wavelength, 2*pi or dilation is applied here.
It preserves the source's addition order separately for the two components
but does not claim the same floating-point rounding as the uncancelled path.

Every beam uses the same generated list. Legacy BEAM has `(beam,tilt)` order;
THEORINTS and experimental storage group NUMP-column, NUML-row images by beam.
The forward Python APIs return `(tilt,beam)`. No implicit image reshape is
added. If the caller reshapes, the established tilt axes are `(NUML,NUMP)`
with C-order flattening; no axis swap is required.

## Public interface

```python
from qcbed import parse_for007_scan, generate_scan_tilts

scan = parse_for007_scan("reference_data/FOR007.DAT")
tilts = generate_scan_tilts(scan)
# Pass tilts unchanged to serial, batched or parallel tilt-series execution.
```

`LegacyScanDefinition` is frozen and snapshots basis/corner sequences as
tuples. Fields are `basis_1`, `basis_2` (triples), `corner_a`, `corner_b`,
`corner_c` (pairs), `scan_columns`, and `scan_rows`. It can also be constructed
explicitly without files. Generated arrays are new float64 arrays of shape
`(scan_rows*scan_columns,2)`, with H in column 0 and K in column 1.

`parse_for007_scan(path, *, experimental_header_path=None)` reuses the
existing FOR007 prefix cursor rather than finding records by line number or
filename text. The default EXPINT path is relative to the FOR007 directory.
That is an explicit Python convenience policy: legacy OPEN instead resolves
relative paths from its process working directory. An explicit override is
accepted and resolved relative to the caller's working directory. Filename
characters otherwise follow the existing raw/formatted parser contract.
`parse_for007` itself retains its original prefix-only behavior and return
model; it does not require an experimental file or scan records.

## Regression evidence and limits

FOR007.DAT and STARTINGFOR007.DAT both name in.dat, whose header is `39 39 1`.
Both give basis vectors `(1,1,0)` and `(1,-1,0)`. The initial corners are
also printed in FOR0013.DAT: A=(-.508156,.407641), B=(.407641,.508156),
C=(-.407641,-.508156). The resulting 1521-pair list has cell-basis corners:

| Location | H | K |
| --- | ---: | ---: |
| A, index 0 | -.100515 | -.915797 |
| B, index 38 | .915797 | -.100515 |
| C, index 1482 | -.915797 | .100515 |
| B+C-A, index 1520 | .100515 | .915797 |
| Center, index 760 | 0 | 0 |

Tests use these independently calculated values (absolute tolerance 3e-16),
explicit integer-valued normal/sheared/one-direction/duplicate scans, and a
scalar reproduction of the source normalization/dilation flow with a
nontrivial basis and reciprocal angle (absolute tolerance 9e-16). They also
check record continuation, ignored trailing records, malformed data, finite
values, singleton rejection and immutable representation.

Validation with
`C:\Users\peng863\AppData\Local\miniconda3\envs\qcbed-python\python.exe`:
the focused scan and existing parser suites passed **309 tests in 4.25 s**;
the complete suite passed **1234 tests in 25.71 s**. `git diff --check` passed.

Remaining limits: singleton semantics are undefined; the two ignored third
basis components have no recovered scan role; normalization overflow or
underflow in extreme legacy inputs is not emulated by algebraic cancellation.
This initial scan is not a reconstruction of fitted/refined scan coordinates
or gradient-transformed experimental data. Exact executable/literal-rounding
equivalence has not been established. No protected source or reference file
is modified.
