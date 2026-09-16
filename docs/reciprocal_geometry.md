# Two-dimensional reciprocal geometry

## Authoritative source equations

The current readable equations agree across the driver, SFMAP, SFMAP2 and
PROPG. No declaration from the related-version private predecessor implementation is needed to
derive this geometry. In `private predecessor implementation`:

```text
input CELL2D(:,1:3) = (a, b, gamma_degrees)
CELL2D(:,3) = cos(3.14159265358 * gamma_degrees / 180.0)
VOL2D = a*b*sqrt(1 - CELL2D(:,3)*CELL2D(:,3))
CELL2R(:,1) = b/VOL2D
CELL2R(:,2) = a/VOL2D
CELL2R(:,3) = -CELL2D(:,3)
```

Thus CELL2D's third entry changes meaning after input: it stores the direct
angle cosine, not degrees. GAMMAREAL retains the input angle. The Python
representation keeps these meanings in separate named fields. VOL2D is area
in angstrom squared; a and b are lengths in angstroms. CELL2R's first two
entries are reciprocal basis lengths in inverse angstroms; its third entry
is the dimensionless reciprocal angle cosine. For a valid direct cell,
gamma-star is the supplementary angle, not the direct angle.

There is **no 2*pi factor** in the reciprocal lattice. Direct/reciprocal dot
products are a.a-star=b.b-star=1 and a.b-star=b.a-star=0. Phase expressions
introduce their own explicit 2*pi elsewhere; this unit adds none.

SFMAP (`private predecessor implementation`) and SFMAP2 (`private predecessor implementation`) use:

```text
d = APPLIEDDILATION
g_squared = (H*H*a_star*a_star + K*K*b_star*b_star
             + 2*H*K*a_star*b_star*cos_gamma_star) / (d*d)
S2 = 0.25*g_squared
STL = SATM = sqrt(S2)
```

PROPG (`private predecessor implementation`) uses the same metric through A2, B2 and ADOTB,
each divided by d squared. This module computes only that shared geometry,
not the excitation error or propagator. Dilation is dimensionless and scales
direct lengths by d, reciprocal vectors by 1/d.

STL is `sin(theta)/lambda = |g|/2`, in inverse angstroms. DSAP is the maximum
included STL; SFMAP/SFMAP2 exclude values when `sqrt(S2)>DSAP`. DSMAX is the
STL limit used by SFLIST: it computes `DS=0.5*sqrt(g_squared)` and compares
that with DSMAX (`private predecessor implementation`). Neither cutoff is applied here.

## Cartesian frame and API

The Fortran metric does not fix an absolute laboratory orientation. Python
uses the explicit right-handed frame with direct a on +x and b in the upper
half-plane. Write c=cos(gamma), t=sqrt(1-c*c):

```text
a_vector = (a, 0)
b_vector = (b*c, b*t)
a_star_vector = (a_star*t, -a_star*c)
b_star_vector = (0, b_star)
g = (H*a_star_vector + K*b_star_vector) / dilation
```

This is a Cartesian embedding of the source metric, not an additional
orientation inferred from the input. In an oblique lattice gx/gy are not
simply H*a-star and K*b-star. Array row/column axes still label K/H; they do
not become Cartesian vector-component axes.

The package exports four frozen dataclasses:

- `DirectLattice2D(a_angstrom, b_angstrom, gamma_degrees)` exposes direct
  vectors, `cos_gamma`, `sin_gamma`, and `area_angstrom_squared`.
- `ReciprocalLattice2D(direct)` derives `a_star_inverse_angstrom`,
  `b_star_inverse_angstrom`, `cos_gamma_star`, and reciprocal Cartesian vectors.
- `ReciprocalIndexGrid(mesh_x, mesh_y, centered=True)` contains separate
  dimensionless integer-valued arrays `h` and `k`, stored as float64 to agree
  with the existing helper. They are never overwritten by physical values.
- `ReciprocalVectorGrid(lattice, indices, dilation=1.0)` contains distinct
  float64 arrays `gx_inverse_angstrom`, `gy_inverse_angstrom`,
  `magnitude_inverse_angstrom` and `stl_inverse_angstrom`. The latter two are
  typed scalar-field representations belonging to the immutable grid.

```python
from qcbed import DirectLattice2D, ReciprocalLattice2D
from qcbed import ReciprocalIndexGrid, ReciprocalVectorGrid

cell = DirectLattice2D(2.0, 3.0, 60.0)
grid = ReciprocalVectorGrid(ReciprocalLattice2D(cell), ReciprocalIndexGrid(6, 4))
stl = grid.stl_inverse_angstrom  # shape (4, 6), inverse angstroms
```

Arrays use immutable byte-backed storage; element writes and re-enabling the
write flag fail. No caller-supplied writable arrays are retained. Derived
fields cannot be supplied independently and made inconsistent with the cell.
Arrays are allocated from the validated mesh, not legacy fixed capacities.

## Ordering and validation

Every array has shape `(mesh_y, mesh_x)`: axis 0=row/K and axis 1=column/H.
The source PMAP mapping is `I+(J-1)*MESHX` for one-based x index I and y
index J. Python uses `array[J-1,I-1]`, with zero-based flattened offset
`(J-1)*mesh_x+(I-1)` in C order. No transpose occurs.

The existing `reciprocal_index_grid` helper supplies ordering: centered H is
column-mesh_x//2, K is row-mesh_y//2. Odd size 5 is [-2,-1,0,1,2], or
[0,1,2,-2,-1] unshifted. Origin is at the central indices or [0,0],
respectively. This generic full grid does not apply the source's omitted
first row/column in SFMAP/PROPG. Existing legacy OSHIFT remains even-only.

Lengths must be positive and finite; gamma must be finite and strictly between
0 and 180 degrees. Geometry that degenerates or overflows in float64 is
rejected. Dilation and its square must be positive and finite. Mesh validation
is delegated to the existing helper. Bool/non-real lattice inputs and
non-boolean ordering flags are rejected. A nonrepresentable reciprocal metric
raises an error rather than clipping a negative squared magnitude.

## Validation and remaining limitations

Independent analytic tests cover square/rectangular orthogonal cells and
60/120-degree oblique cells, dual-basis identities, signed directions, origin,
odd/even rectangular grids, PMAP flattening, immutability, and invalid inputs.
For a=2, b=3, gamma=60 degrees, the ideal g(1,1) is
`(1/2, sqrt(3)/18)` with magnitude `sqrt(7/27)`.
For a=20 angstroms, H=1/K=0 gives |g|=0.05 and STL=0.025 inverse angstroms;
passing it to the Al elastic factor is checked against the independently
evaluated value 275.340445007875, with relative tolerance 2e-12.

Analytic geometry comparisons allow absolute error 2e-11 because the readable
source uses 3.14159265358 for pi. That literal is deliberately retained: an
input of exactly 90 degrees therefore has a small nonzero computed cosine,
which is not snapped to zero. Internal dual-basis and metric/vector identities
are tested at approximately float64 roundoff. The reference input lattice
a=b=2.8630754 gives reciprocal length 0.34927476936164517 inverse angstroms,
independently checked with Decimal arithmetic; this also agrees with the
rounded 0.34927 in FOR0013.DAT within its printed precision.

Exact executable-Fortran equivalence remains unestablished. Original compiler
literal kinds/flags and the reference build are unverified; the rounded
reference log cannot validate full precision. The Cartesian frame is an
explicit Python choice consistent with the metric. No laboratory orientation,
cutoff policy, scattering-domain filtering, structure-factor summation,
Debye-Waller factors, absorption, phase grating or propagation is implemented.
