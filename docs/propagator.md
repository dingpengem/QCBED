# Reciprocal-space PROPG kernel

## Decision and source

**Category B:** the physical kernel is specified by readable executable code;
historical FFT normalization and linkage remain unverified. This module
evaluates that kernel on the documented Python reciprocal grid, with explicit
ordering and optional source boundary suppression. It performs no wave
propagation or Fourier transform. private predecessor implementation remains only a related-version
declaration reference; no include capacities or constants were copied.

There is one readable definition: `private predecessor implementation`, PROPG. Its two calls
are `private predecessor implementation` (ONEIT) and `private predecessor implementation` (ONEITPARA), both:

Private implementation excerpt omitted from shared documentation.

The call arguments and physical equations do not differ between these paths.
The single readable OSHIFT is `private predecessor implementation`; PROPG calls it at line 718.
All readable `.f` files were searched for PROPG, PRARRAY, PMAP, RMT, shifts
and tilt transfer. Binary artifacts are not alternative readable definitions.

## Complete source equation

With a*,b*,c*=CELL2R(:,1:3), d=APPLIEDDILATION, lambda=WAVL, z=DZ,
th=TILT(1), tk=TILT(2), the executable statements are:

```text
TWOPI = 6.2831853
TWOKR = 0.5*lambda
A2 = a*^2/d^2
B2 = b*^2/d^2
ADOTB = a* b* c*/d^2
A1 = 2*th + H
B1 = 2*tk + K
SG = -0.5*lambda * (A1*H*A2 + B1*K*B2 + (A1*K+B1*H)*ADOTB)
TEMP = -TWOPI*(SG*z*d + H*X0 + K*Y0)
P = cos(TEMP) + i*sin(TEMP)
```

The local printed constant is deliberately retained; it differs from the
atomic phase constant 6.283185307 and mathematical 2*pi in the Fourier sum.
The introductory PROPG comment's positive `exp(2*pi*i*Sg*Z)` does not override
the executable negative sign multiplying SG. Because SG itself is negative
at zero tilt, the resulting zero-tilt propagation phase is **positive**.
The separate shift phase is **negative**.

An equivalent geometric interpretation is:

```text
g = (H*a_star_vector + K*b_star_vector)/d
t = (th*a_star_vector + tk*b_star_vector)/d
SG = -lambda/2 * (|g|^2 + 2*t.g)
phi_propagation = -6.2831853*SG*z*d
phi_shift = -6.2831853*(H*X0+K*Y0)
```

The implementation evaluates the source A1/B1 metric grouping using the
existing reciprocal-lattice values. It does not reconstruct a second lattice
or replace the expression with a different Fresnel sign convention. Separating
the two phases before adding them can cause last-bit differences from the
source's single grouped TEMP expression; tests bound this effect.

Units: lambda and undilated z are angstroms; d is dimensionless and the
effective propagation distance is z*d. Direct lattice lengths are angstroms;
reciprocal lengths and g,t are inverse angstroms without a 2*pi factor.
SG is inverse angstroms; both phases are radians; P is dimensionless complex128.
The oblique cross term uses cos(gamma*)=-cos(gamma), not the direct angle
cosine. Both components of internal tilt can affect either axis in an oblique
cell. Holding internal th/tk fixed, the metric contributes 1/d^2 while
distance contributes d, giving a net propagation phase proportional to 1/d.
The caller also changes tilt with dilation; that prior operation is separate.

## Tilt audit and explicit API boundary

TILT(1/2) are dimensionless coefficients along the reciprocal cell a*/b*
basis, aligned with integer H/K respectively. This follows from executable
`2*TILT(1)+H` and `2*TILT(2)+K` and the common reciprocal metric, corroborated
by driver comments at `private predecessor implementation`. They are not radians, milliradians,
tangents or direction cosines. The physical transverse reciprocal displacement
is t above; no external angle-to-displacement conversion is introduced.

The input/caller trace is:

1. The driver reads the zone indices with voltage. It later reads
   AY,BEE,NOTHING,CEE,DEE,NAUGHT (`private predecessor implementation`), establishing two scan
   basis vectors in the cell basis, then A(1:2),B(1:2),C(1:2) scan corners
   at lines 435-437 and their refinement masks. These corners are not direct
   inputs to PROPG. The source comments mention ISCTL, but that alone does
   not establish an implemented angular-input branch.
2. Lines 403-430 compute basis lengths from AY/BEE/CEE/DEE and lattice gamma,
   converting the **lattice angle**, not tilt, with 3.1415926535/180.
   NORMFACTOR=LENGTHY/LENGTHX and INVNORMFACTOR=LENGTHX/LENGTHY.
3. ONEIT (`private predecessor implementation`) scales second scan coordinates by NORMFACTOR,
   divides both corner coordinates by each slice's APPLIEDDILATION, and
   interpolates TILTARRAY using NUML/NUMP and endpoints divided by N-1.
   ONEITPARA (`private predecessor implementation`) uses the same equations.
4. For each tilt and slice, the caller takes u=TILTARRAY(...,1),
   v=INVNORMFACTOR*TILTARRAY(...,2), then computes
   th=AY*u+CEE*v and tk=BEE*u+DEE*v (`private predecessor implementation`,
   `private predecessor implementation`) immediately before PROPG.

No degree/mrad/tangent conversion of the scan corners was found in this path.
The kernel API accepts **already transformed internal th/tk**. It neither
repeats the caller's division by dilation nor exposes an unverified user-facing
angle conversion. The existing generic TiltCoordinate model is not silently
reinterpreted or changed; the new parameter names state the precise basis.

## Slice shifts and zero thickness

The driver geometry READ at `private predecessor implementation` reads DELTAZ,XSHIFT,YSHIFT.
The calls pass these unchanged as DZ,X0,Y0. They are the existing
SliceGeometry thickness_angstrom, x_shift_fractional and y_shift_fractional.
Shifts are fractional **direct-cell translations**, entering H*X0+K*Y0;
they are not multiplied by thickness or dilation. No shift is added to
structure factors or phase gratings by this task.

Zero z gives zero physical propagation phase, but nonzero shifts still give
a translation phase. With zero shifts the unmasked zero-thickness kernel is
exactly 1+0i. The source-required masked kernel still has boundary zeros.
At H=K=0 both phases vanish for every tilt and shift, so P(0,0)=1.

A full fractional-cell shift is only approximately periodic with the source
constant: 6.2831853 differs from 2*pi by about -7.179586e-9. Thus a one-cell
x translation changes a coefficient by up to about |H|*7.18e-9, beyond ordinary
machine epsilon. Tests retain this known literal effect; they do not wrap
shifts modulo one, substitute mathematical pi, or hide the discrepancy.

## Boundary, ordering and PRARRAY

PROPG first zeros all PMAP (`private predecessor implementation`), then fills J=2..MESHY and
I=2..MESHX, with K=J-MESHY/2-1, H=I-MESHX/2-1. On even meshes the first
centered row/column are therefore the suppressed K=-Ny/2 and H=-Nx/2 terms.
No DSAP, DSMAX, RMAX, circular bandwidth mask or detector aperture appears
inside PROPG. The boundary matches the existing legacy_nyquist_mask and is
reused through propagation_boundary_mask without duplicating its logic.

The source then calls OSHIFT and copies PMAP to PRARRAY(IS,:). The even-axis
half swaps are already validated against the Python ordering helpers. In
unshifted order the zero row moves to row=Ny/2 and the zero column to
column=Nx/2 (zero-based), with each rectangular dimension handled separately.
PRARRAY is the unshifted reciprocal propagation multiplier for the current
tilt and slice; PMAP is scratch space during construction and later the
wavefunction working array. No transform occurs inside PROPG itself.

The new phase functions return arrays in the input grid's explicit ordering.
The final kernel requires an explicit output ordering. When a conversion is
needed it reuses centered_to_unshifted/unshifted_to_centered; ordering never
zeros values. `boundary="legacy"` is the default and requires even dimensions.
`boundary="none"` evaluates the scalar physical formula on all grid points,
including odd and mixed-parity grids. Generic odd FFT support is not restricted
by legacy OSHIFT. No exact odd-mesh legacy traversal is claimed.

All numerical outputs are `(mesh_y,mesh_x)`, x/H on columns, y/K on rows.
Phase arrays are float64; kernels are new C-contiguous complex128 arrays;
boundary masks are Boolean. Every unmasked finite-phase coefficient has unit
modulus up to rounding; source-suppressed coefficients are exactly 0+0i.

## RMT, reuse and changing parameters

`private predecessor implementation` sets RMT=1/DBLE(MT), MT=Nx*Ny. PROPG never reads RMT.
The caller later performs `PMAP=RMT*PMAP*PRARRAY` after FFT(-1), and follows
with FFT(+1) (`private predecessor implementation`, `private predecessor implementation`). This is a separate
transform-pair normalization, consistent with the readable CDFT2D mathematical
specification; it is not physical attenuation in PRARRAY. The missing readable
1D CDFT still prevents validating the exact historical transform implementation.

Python's reciprocal_to_real_sum is unnormalized and its mathematical inverse
already divides by Nx*Ny. RMT must not be embedded in this kernel or added
again to that inverse pair. The new API evaluates the literal source physical
phase at the project H/K labels; it does not certify how a future full
wave-propagation sequence maps every historical intermediate onto that pair.
No such sequence or normalization wrapper is implemented here.

Variation/reuse established by the callers:

| Scope | Source state |
| --- | --- |
| Simulation | Mesh, initial voltage/wavelength, zone and scan basis; RMT depends only on mesh. |
| Region | Depth and repeated slice types; vacancy/surface models affect applied dilation and scan geometry indirectly. |
| Slice | CELL2R, DELTAZ, XSHIFT/YSHIFT, APPLIEDDILATION and slice-specific TILTARRAY. |
| Tilt | Internal th/tk change; PROPG is called for every slice type inside every tilt iteration. |
| Refinement/objective evaluation | UPDPARS can change voltage/WAVEL, scan corners, depth and applied dilation; callers rebuild scan grids, gratings and propagators. |

MSEQ is populated with slice IDs while accumulating DELTAZ*d until region
depth is reached; three parameter/sequence updates precede each objective
calculation. PRARRAY for a slice type is reused at repeated appearances of
that ID during the current tilt. It is overwritten for the next tilt; it is
not a permanently tilt-independent cache. The previously documented unused
post-final MSEQ(MS+1) read remains outside this task.

## Public API and validation

```python
PropagationParameters(
    wavelength_angstrom, slice_thickness_angstrom,
    tilt_h_fractional=0.0, tilt_k_fractional=0.0,
    x_shift_fractional=0.0, y_shift_fractional=0.0,
)
propagation_phase_radians(reciprocal_grid, parameters)
fractional_slice_shift_phase_radians(
    reciprocal_grid, *, x_shift_fractional=0.0, y_shift_fractional=0.0
)
propagation_boundary_mask(reciprocal_grid)
propagation_kernel(
    reciprocal_grid, parameters, *, ordering, boundary="legacy"
)
```

All are exported from qcbed. Parameters are frozen; inputs must be finite
real scalars, excluding bool and strings. Wavelength is positive, thickness
nonnegative, and signed tilts/shifts are allowed. Existing SliceGeometry
requires strictly positive thickness; zero is deliberately allowed here to
test the kernel identity without changing that geometry model. Wavelength
can be obtained from the existing electron_wavelength_angstrom helper.
Invalid ordering/boundary names and nonrepresentable phases raise errors.
Input arrays and models are not modified.

Tests cover zero thickness, signed H/K, single/both internal tilts, independent
x/y shifts, full-cell translations, square/rectangular/60/120-degree cells,
dilation, both orderings, Nyquist zeros, odd generic grids, unit modulus,
invalid data, dtype and immutability. Hand cases use a=2,b=4,lambda=.025,DZ=2:
for H=1,K=0 at zero tilt, phase=.039269908125 and
P=.9992290362424846+.03925981571423079i. A quarter x shift at zero thickness
gives P(1,0)=1.7948965149208059e-9-i.

Independent scalar Cartesian dot-product sums check every pixel on a 6x4
oblique grid, including both tilts and shifts. Using ideal 60/120-degree
geometry allows 1e-11 absolute tolerance for the known rounded geometry-pi
constant; using that printed angle in the independent dual vectors allows
3e-15. Maximum residuals (absolute and relative, since |P|=1):

- Ideal oblique geometry: 8.159e-12 at gamma=120,H=-3,K=-2.
- Printed-angle Cartesian calculation: 8.887e-16 at gamma=120,H=2,K=-2.
- Orthogonal hard-coded hand cases: 1.923e-13 at H=K=1 (tiny nonzero source
  cosine at input gamma=90); tolerance 3e-12.

One-cell shift tolerance is 2.2e-8 for this mesh, derived from the printed
TWOPI residual and maximum |H|=3. No standalone saved PROPG/PRARRAY grid was
identified in reference records or readable WRITE statements, so there is
no historical pixel comparison. Literal kinds, original include/linkage and
CDFT implementation remain unresolved. Exact executable equivalence is not
claimed. Full propagation, multislice, probes, detectors, tilt orchestration,
parallelism and refinement remain deferred.
