# Projected-potential and phase-grating reference

## Decision: category B

The mathematical direct Fourier sum and PGMAP exponential are specified by
the available source. The exact legacy FFT operation is not fully established
for the historical executable; readable CDFT now confirms the stated signs.
The serial direct-sum reference remains the default, with separate
source-backed masks and PGMAP exponentiation. An explicit NumPy FFT path is
now validated against that Python mathematical reference; see fft_backend.md.
It does not emulate CDFT or claim complete legacy executable equivalence.

Primary evidence is the current readable `.f` code, its mathematical transform
specification, and saved input/output. private predecessor implementation remains a related-version
declaration reference and supplies no new numerical assumption or capacity.

## Source definitions, calls and intermediate arrays

There is one readable PGMAP definition, `private predecessor implementation`. Both callers use:

```text
ONEIT     (private predecessor implementation):     SFMAP2 -> OSHIFT -> FFT(-1) -> PGMAP
ONEITPARA (private predecessor implementation):   SFMAP2 -> OSHIFT -> FFT(-1) -> PGMAP
```

The only readable SFMAP and SFMAP2 definitions are `private predecessor implementation` and
`private predecessor implementation`. SFMAP is called at `private predecessor implementation`; SFMAP2 is additionally
called by the driver at lines 579,632,737,842,926,989,1129,1271,1406,1542.
Those diagnostic/refinement calls do not themselves invoke PGMAP. SFMAP
initializes atoms and applied state and stores separate SFS components;
SFMAP2 builds the combined PMAP, with later reflection refinement multipliers.
The raw elastic/absorptive equations are documented in structure_factors.md
and absorption.md and reused without changes.

SFMAP2 sets PMAP to zero (`private predecessor implementation`) before filling its active
reciprocal coefficients. At line 561 it writes:

```text
SFS[...,1:4] = SFRE, SFIM, AFRE, AFIM
C = PMAP = (SFRE-AFIM) + i*(SFIM+AFRE) = elastic + i*absorptive
```

It copies this reciprocal PMAP to PGARRAY at lines 565-567. Then the callers
transform PMAP in place and PGMAP overwrites both PMAP and PGARRAY with the
real-space transmission grating. PGARRAY's meaning therefore changes during
construction; the pre-transform copy is not already a phase grating.
PGMAP receives the already combined real-space quantity, not separate
elastic and absorptive arrays. No second PGMAP formula/call convention was found.

## Fourier derivation and its limit

`private predecessor implementation`, FFT, actually calls CDFT2D rather than the FOURN named
in its stale introductory comment. Its executable packing is:

Private implementation excerpt omitted from shared documentation.

Unpacking reverses those exact expressions, with no transpose, conjugation,
scale multiplication, or pixel permutation. The phase-grating callers pass
ISIGN=-1. The readable `private predecessor implementation` CDFT2D loops over its second
dimension and calls `cdft(n1,isgn,a(0,j),...)`, then CDFT2D_SUB.
CDFT2D_SUB (`private predecessor implementation`) gathers each first-dimension complex pair
into t(2*j),t(2*j+1), calls `cdft(2*n2,isgn,t,...)`, and copies back to the
same positions. Thus both axes receive the same sign and the wrapper applies
no normalization. The mathematical specification at `private predecessor implementation`
explicitly associates -1 with a negative-exponent sum and says the inverse
pair requires an extra 1/(n1*n2) multiplication.

The restored private predecessor implementation now supplies CDFT/MAKEWT and executable butterflies.
Its negative-exponent, unnormalized behavior agrees with this potential
construction. Original linkage, compiler options and historical executable
agreement remain unverified. See [FFT correction](fft_sign_correction.md).

For the direct reference, after interpreting centered H/K:

```text
x = column/mesh_x, column = 0,...,mesh_x-1
y = row/mesh_y,    row    = 0,...,mesh_y-1
V[row,column] = sum_HK C[H,K]*exp(-2*pi*i*(H*x+K*y))
```

There is no 1/(mesh_x*mesh_y) normalization. The mathematical expression
would be `numpy.fft.fft2(ifftshift(C))` for centered even grids, or
`numpy.fft.fft2(C)` for already-unshifted coefficients, under NumPy's default
normalization. This mapping is now implemented and tested against the Python
direct reference, but is not certified as a replacement for the legacy FFT.
Using the default IFFT would reverse the sign and add an unwanted scale.

OSHIFT (`private predecessor implementation`) is called **before** FFT(-1). Its half-grid swaps
move centered origin coefficients to index zero, separately in x/y, with no
transpose. The direct reference instead reads the grid's explicit signed H/K
and needs no OSHIFT operation. Centered and unshifted storage produce the same
spatial pixels. No final fftshift is applied. The DFT kernel uses mathematical
2*pi (`math.tau`), not the rounded 6.283185307 used in atomic position phases.
No readable transform twiddle constant is available to emulate.

## Coordinates, masks and normalization

All complex arrays are complex128 and shaped `(mesh_y,mesh_x)`, with y/K on
axis 0 and x/H on axis 1. PMAP flattening is row*mesh_x+column (zero-based).
Samples include fractional zero and exclude fractional one: the last x is
(mesh_x-1)/mesh_x. In an oblique lattice the physical point is x*a_vector
+ y*b_vector; the grid is a parallelogram, not a Cartesian rectangle.
The dual-lattice identity reduces g.r to H*x+K*y, so no gamma or Cartesian
component is inserted again in the Fourier phase. Geometry and dilation
already affect the coefficients and physical cutoff through STL.

The source's masks before PGMAP are distinct:

- **Physical cutoff:** SFMAP/SFMAP2 compare `sqrt(S2)>DSAP` and set all four
  components to zero when true (`private predecessor implementation`, `private predecessor implementation`). Thus
  STL<=DSAP is retained, with equality included. DSAP is in inverse angstroms,
  not |g| units. g=0 remains for nonnegative DSAP; the API rejects negative
  or nonfinite DSAP as explicit validation policy.
- **Nyquist boundary:** SFMAP2 starts with zero PMAP, uses
  MX=MESHX/2-1, MY=MESHY/2-1, IX=MESHX+1 and loops K=-MY..MY,
  H=-MX..MX with H-fastest insertion. On even meshes this fills centered
  `[1:,1:]`, leaving H=-MESHX/2 and K=-MESHY/2 zero. The mask identifies
  those frequencies in either ordering. It does not zero real-space pixels.
- DSMAX limits the SFLIST report, not phase-grating reconstruction. No RMAX,
  detector aperture, or extra circular cutoff was found in this construction
  path. PROPG has separate propagation boundaries that are not applied here.

Both masks apply before the direct sum in the composing API. The underlying
direct sum applies neither mask implicitly. Generic odd/non-power-of-two
grids are valid for the mathematical direct sum; the legacy Nyquist mask and
composing API require even dimensions. Odd legacy traversal is not inferred.
CDFT2D's documented power-of-two restriction is not a direct-sum restriction.

No Hermitian completion is performed. All supplied complex coefficients are
used. Elastic Friedel pairs give a real potential up to roundoff; combined
elastic+i*absorptive coefficients generally do not form a Hermitian array.
g=0 contributes its unmodified complex value uniformly to every pixel.
RMT=1/MT is applied later in the propagation loop, not in PGMAP construction.

## PGMAP, interaction constant and thickness

Executable `private predecessor implementation` gives:

```text
PR = real(V); PI = imag(V)
TEMP = -SIGMA*PR
amplitude = exp(SIGMA*PI)
T = amplitude*(cos(TEMP) + i*sin(TEMP)) = exp(-i*SIGMA*V)
```

The elastic phase sign is negative. For PI<0 and positive SIGMA, amplitude
decreases; for purely real V, |T|=1. Positive PI produces amplification by
the written equation. No sign clipping or guarantee that a truncated Fourier
potential is everywhere absorptive is imposed. Zero V gives exactly 1+0i.
The output is dimensionless.

The driver (`private predecessor implementation`) and UPDPARS (`private predecessor implementation`) calculate:

```text
WAVEL = 0.3878314 / sqrt(KV*(1+0.97846707e-3*KV))
SIGMA = 0.020886416 * WAVEL * (1+1.9569341e-3*KV)
```

The existing wavelength helper supplies WAVEL in angstroms from voltage in
kV. No separate wavelength override is accepted; the derived value is checked
for finite positive representability. The SIGMA unit 1/(V angstrom) follows
the projected-potential interpretation. As in previous modules, interpreting
the source-scaled V as projected potential in V angstrom is a dimensional
inference; no physical unit conversion beyond the visible constants is added.

Cell area/dilation and atom occupancies are already in the existing elastic
and absorptive coefficients (with their differing area scaling). Slice
thickness **has not been multiplied into those factors**, nor does PGMAP
multiply by thickness: the fixed atom content represents a projected slice.
Adding a thickness factor would change the visible source equation. Slice
thickness and fractional translations are used by PROPG instead. The new
APIs therefore have no thickness parameter; existing SliceGeometry continues
to validate thickness. Tests change finite thickness/shift values while
holding atom content and lattice fixed and confirm unchanged gratings.

## Public functions

```python
reciprocal_cutoff_mask(grid, dsap_inverse_angstrom) -> bool array
legacy_nyquist_mask(grid) -> bool array
projected_potential_direct(coefficients, grid) -> complex128 array
projected_potential_fft(coefficients, grid) -> complex128 array
electron_interaction_constant(voltage_kv) -> float
phase_grating_from_potential(projected_potential, *, voltage_kv) -> complex128 array
phase_grating_from_structure_factors(
    elastic, absorptive, grid, *, voltage_kv, dsap_inverse_angstrom, method="direct"
) -> complex128 array
```

The default method remains "direct"; "fft" explicitly selects the NumPy
backend and unknown method names are rejected. Both paths share combination,
masking, interaction constant and exponentiation. All are exported from qcbed.
`grid` is the existing ReciprocalVectorGrid;
arrays are required to be finite, nonempty 2D NumPy complex128 arrays, with
shapes matching the grid where supplied. There is no implicit transpose,
broadcasting or narrowing. Caller arrays and immutable grid data are unchanged.
Voltage is finite positive kV, bool/non-real scalar inputs are rejected,
and nonrepresentable sums/exponentiation raise ValueError. Exponential
underflow to zero attenuation is permitted. The separate mask outputs are
Boolean arrays; numerical inputs and outputs use float64/complex128.

For atom input, first call existing elastic_structure_factors and
absorptive_structure_factors with identical species, sites, grid, vacancy and
applied B; pass their results to phase_grating_from_structure_factors. That
function reuses combine_elastic_absorptive_structure_factors, then applies
both masks, reconstructs V, and exponentiates. No atomic or reciprocal metric
calculation is duplicated. The raw atomic helpers still validate their whole
input grid's scattering domain before this later masking stage; this wrapper
does not extend their domains or implement refinement.

## Numerical validation and historical comparison

Tests include explicit roots of unity on a 4x2 mesh, one signed H or K
coefficient, Friedel pairs, rectangular and odd direct grids, both reciprocal
orderings, fractional sample endpoints, periodic shifts, both masks, oblique
metric dependence, existing atomic API composition, dtype and immutability.
A separate scalar nested-loop reference evaluates every pixel of 4x4, 6x4,
and 5x3 cases. Additional tests now compare NumPy FFT with this reference.
The direct implementation is a serial loop
over coefficients with broadcast real-space trigonometric arrays, O(M^2)
work for M pixels; no performance or scalability claim is made.

Independent 65-digit Decimal constants and exponentials (Taylor sine/cosine)
give these hard-coded tests for uniform V=2-0.5i:

| kV | SIGMA | T |
| --- | --- | --- |
| 100 | .0009243913556182865 | .9995362029095982 - .0018479273565984993i |
| 200 | .0007288360746343612 | .9996345863381000 - .001457140527998951i |
| 300 | .0006526122090303823 | .9996728956001657 - .0013047982143670319i |

With Al at the origin, area=200, B=.86, INELA/B/C=(0,1,1), vacancy=0,
and DSAP=0, the retained DC coefficient is
1.409767975005-9.292395257301284i. At 200 kV the independently evaluated
uniform grating is .9932497254026199-.0010205542781868584i.

Measured maximum direct-versus-nested absolute residual is 4.45e-16
(4x4 unshifted, row=1,column=2; tied cases exist); maximum relative residual
is 3.73e-16 (6x4 unshifted case). Tests allow 8e-15 absolute for these sums.
The three independent uniform gratings differ by at most 1.12e-16 absolute
and relative; their absolute tolerance is 3e-16.

No saved projected-potential or PGARRAY pixel fixture was identified. No
WRITE of PMAP/PGARRAY occurs in the readable construction path. VGVSH.DAT
records depth and SFFAC refinement values (`private predecessor implementation` and driver WRITEs),
not phase-grating pixels; sigma.dat is a large data array, not the scalar
interaction constant. The saved FOR0013.DAT logs
SIGMA=0.72884E-03 at 200 kV (E12.5). Python gives
0.0007288360746343613: absolute difference 3.92537e-9, relative 5.38578e-6,
within the printed last-place half-step 5e-9. This scalar has no pixel
location and does not validate phase-grating output. Prior structure-factor
saved-output tests remain in the full suite.

Original 1D FFT code/linkage, literal kinds/compiler flags, and full complex
phase-grating reference output remain missing. The mathematical NumPy FFT
candidate is not yet a validated legacy replacement. Exact executable
equivalence is not claimed. PROPG, free-space propagation, multislice,
incident probes, detector extraction, tilt scans, parallelism and refinement
remain outside this stage.
