# Elastic atomic structure factors

## Readable source and call paths

The only readable definitions are `SFMAP` in `private predecessor implementation`
and `SFMAP2` in `private predecessor implementation`. The complete readable call set is:

- SFMAP: `private predecessor implementation` (PROGRAM MULTIS).
- SFMAP2: `private predecessor implementation`;
  `private predecessor implementation` (ONEIT); `private predecessor implementation` (ONEITPARA).

Some `.f` filenames contain binary artifacts. Their linkage is unverified;
there is no conflicting readable definition. private predecessor implementation remains only the
related-version declaration reference described in source_provenance.md.
No capacity, equation, or constant is copied from it here.

The shared raw elastic kernel is `private predecessor implementation` and
`private predecessor implementation`. Both call `ATOMICSCATFACT(ZT(M),SATM,2)/VOL` for each
species, then accumulate atoms in their input order. SFMAP first reads atoms
and initializes vacancy/dilation/thermal state; SFMAP2 consumes applied state
and subsequently applies reflection refinement multipliers SFFAC with SYMOPT
selection. These later multipliers are not part of this module. SFMAP stores
elastic real/imaginary components separately in SFS slots 1/2; absorption is
in slots 3/4. SFMAP2 combines elastic and absorptive components as
`DCMPLX(SFRE-AFIM,SFIM+AFRE)` into PMAP and copies to PGARRAY.

## Equation and units

With a,b in angstroms, c=cos(gamma), dilation d=APPLIEDDILATION:

```text
VOL = a*b*sqrt(1-c*c)*d*d                    [angstrom^2]
S2 = (H^2*a_star^2 + K^2*b_star^2
      + 2*H*K*a_star*b_star*cos_gamma_star)/(4*d*d)
STL = sqrt(S2) = |g|/2                      [inverse angstrom]
AF[t] = ATOMICSCATFACT(Z[t], STL, 2)/VOL
phi[m] = 6.283185307*(H*x[m] + K*y[m])       [radians]
OCC[m] = STARTOCC[m]*(1-VAC)                 [dimensionless]
thermal[m] = exp(-APPLIEDDW[t]*S2)
F[H,K] = sum_m AF[t] * (cos(phi[m])+i*sin(phi[m])) * thermal[m] * OCC[m]
```

`t=ITYPE(IS,m)` is one-based in Fortran; Python calls the existing checked
`species_index_zero_based` helper before lookup. Atom x/y are fractional
direct-lattice coordinates; they are not Cartesian angstrom positions.
Negative/out-of-cell finite coordinates are accepted without wrapping: wrapping
would change results slightly with the printed, rounded phase constant.

The phase sign is **positive**, based on actual SFRE/SFIM additions, regardless
of generic negative-sign Fourier comments elsewhere. SFMAP assigns
`TWOPI=DBLE(6.283185307)`; SFMAP2 assigns `TWOPI=6.283185307` to REAL*8.
Python uses that printed decimal in binary64, not `2*math.pi`.

The return is the **source-scaled elastic coefficient**, not an unscaled sum
of electron scattering lengths. The atomic helper already includes 47.878009
or 1.145896 and its source interpolation/Mott conversion; see scattering.md.
There is exactly one area division. There is no thickness, wavelength,
interaction constant, relativistic multiplier, extra 2*pi, or beam-count
factor in this sum. The name VOL does not mean 3D volume here.
Interpreting the atomic factor as V angstrom^3 makes the result a projected
potential coefficient in V angstrom. That physical unit label remains an
inference from scaling and caller use; the source-scaled number is the contract.

APPLIEDDW has angstrom^2 units inferred from its dimensionless product with
S2. This is the applied B parameter, not an RMS displacement or a value to
multiply by 8*pi^2. The exponent is negative and already uses |g|^2/4.
SFMAP initializes it from regional DW and can alter it with vacancy/surface
expressions (`private predecessor implementation`); UPDPARS also updates it
(`private predecessor implementation`). The caller must supply already-applied B values.
No region, surface, or refinement state is guessed from the parser subset.

For fixed S atoms the source occupancy equation is explicit at
`private predecessor implementation` and `private predecessor implementation`. The API requires the slice vacancy
fraction explicitly, after any upstream vacancy model/clamping. It does not
calculate a vacancy profile or associated dilation. Dilation is supplied
through ReciprocalVectorGrid and consistently scales both STL and area.
This API describes fixed input S sites with that occupancy equation, not
arbitrarily refined effective atom occupancies.

## Slice shifts, gas atoms and symmetry

There is **no slice translation phase in SFMAP/SFMAP2**. The atomic-sum
translation factor is identity. XSHIFT/YSHIFT are read by the driver and
passed to PROPG (`private predecessor implementation`, `private predecessor implementation`). PROPG alone uses
`TEMP=-TWOPI*(SG*DZ*APPLIEDDILATION(IS)+H*X0+K*Y0)` at `private predecessor implementation`, with
its own printed `TWOPI=6.2831853`. Neither this translation nor the propagation
term is implemented here. SliceGeometry is deliberately not an argument:
only its lattice fields are needed to construct the reciprocal grid.
Tests use nonzero slice shifts and changed thickness and confirm that the
atomic sum is unchanged. No shift is folded into atom coordinates.

G records cause RANDO calls for coordinates and STARTOCC in
`private predecessor implementation`. Input G sites are unresolved, so this API rejects
them. It provides no randomization or materialized-gas bypass.

Both signs of H/K are calculated directly. No conjugate reflection is filled
from another reflection. Real scattering/occupancy/thermal factors imply
F(-g)=conj(F(g)) for this raw sum. SFMAP2's later symmetry-selected refinement
multipliers are a separate operation. g=0 uses fe(Z,0) with unit atomic and
thermal phase/factors; no special subtraction or zeroing occurs.

## Public API and grid boundary

```python
elastic_structure_factors(
    species: Sequence[Species],
    atoms: Sequence[AtomSite],
    reciprocal_grid: ReciprocalVectorGrid,
    *,
    applied_debye_waller_angstrom_squared: Sequence[float],
    vacancy_fraction: float,
) -> NDArray[np.complex128]
```

Exported from qcbed. One B per species is required, even for unused species.
Species must be nonempty; an empty atom sequence returns a zero grid. B and
starting occupancy must be finite and nonnegative; vacancy must be finite in
[0,1]. Nonnegative B/occupancy are explicit physical-domain validation policy,
not recovered Fortran runtime checks. Occupancies above one are not clipped.
Bool/non-real data, invalid species indices, unsupported atomic numbers and
scattering coordinates are rejected using the existing atomic helper's domain.
All species are evaluated at each grid point as in the source, even when
their effective occupancy is zero. Nonrepresentable intermediate phases,
thermal exponents, area, or final sums raise ValueError.

Output is a new C-contiguous complex128 array `(mesh_y,mesh_x)`. Row/axis 0
is K/y, column/axis 1 is H/x. Dimensionless H/K and physical STL remain
distinct. Real and imaginary terms are added separately in atom input order,
with the visible multiply ordering AF*trig*thermal*OCC. Float64 STL is squared
to recover S2 from the existing grid API; this can differ by a last-bit
rounding from retaining the source's original S2 intermediate.
Input models, lists and immutable reciprocal arrays are not mutated.

The function evaluates the **raw formula at every supplied grid point**, in
centered or unshifted order, including the existing generic odd-grid convention.
This explicitly extends evaluation of the scalar expression, not the legacy
mesh traversal. For even meshes SFMAP/SFMAP2 use MX=MESHX/2-1, MY=MESHY/2-1,
K outer/H inner, IX=MESHX+1 then preincrement: centered array `[1:,1:]`.
They omit the negative Nyquist row/column and zero points with STL>DSAP.
Neither mask is applied here. The full grid must not be passed off as the
complete SFMAP result. Odd legacy traversal equivalence is not claimed.
PMAP index is `1+row*mesh_x+column`; there is no hidden transpose or FFT shift.

## Validation and remaining limits

Hard-coded hand cases use independent 60-digit Decimal arithmetic: direct
Lagrange evaluation of the source Al/H table entries, Taylor sine/cosine of
the printed phase constant, and Decimal.exp. No expectation calls the module.
For a=20,b=10,gamma=90,area=200 and (H,K)=(1,0):

| Case | Expected complex coefficient |
| --- | --- |
| Al at origin | 1.376702225039375 |
| Al at x=1/4 | 6.180927559250617e-11 + 1.376702225039375i |
| Two origin Al sites | 2.75340445007875 |
| Al at x=0 and x=1/2 | 5.55005502248115e-21 + 1.2361855118501233e-10i |
| Starting occupancy 1/4 | 0.34417555625984375 |
| Vacancy 1/4 | 1.03252666877953125 |
| Applied B=1600 angstrom^2 (exp(-1)) | 0.5064604452069665 |
| Al + H at origin | 1.500825339960299475 |

At g=0, Al alone gives 1.409767975005; Al+H gives 1.53640530881.
For a=b=20,gamma=60 (using the printed geometry pi), the independent area
is 346.4101615131225761 and STL(1,0)=0.0288675134595356951; Al coefficients
are 0.8139299198655844 at zero and 0.7886596996082075 at (1,0).
The Al stencil here is (.01,.02,.03,.04) -> (5.867,5.800,5.692,5.547).
Hand-case tolerance is 3e-12 absolute, zero relative. Tests also check species
lookup, separate B values, vacancy, dilation, Friedel pairs, flattening,
orientation, shifts, invalid data, complex128 and input immutability.

The supplied FOR0013.DAT has a slice elastic amplitude/phase listing at lines
900-926. A separate test compares all 25 complex coefficients reconstructed
from that listing with the parsed FOR007 atoms and cell, explicit B=.86,
vacancy=0, dilation=1. These settings are recorded in the reference case.
SFLIST in private predecessor implementation thresholds small amplitudes to zero and formats amplitudes
with E14.7 and phases with F9.4. The comparison uses relative 2e-7 and absolute
5e-7 tolerance to accommodate recorded precision and default-literal differences.
It checks the elastic columns only, not absorption. It is a limited saved-output
comparison, not a rerun or proof of executable equivalence.

Validation measured a maximum absolute hand-case residual of 4.45e-16 and
a maximum saved-output complex residual of 3.2862e-6 (the origin coefficient).
The focused suite passes 49 tests; the complete package suite passes 460.
Both were run with the qcbed-python interpreter specified for this project.

Exact executable Fortran equivalence remains unestablished. Compiler flags,
default-kind literal rounding (including DBLE of an already-rounded literal),
and original include/linkage remain unknown. Absorption, gas materialization,
vacancy/surface model construction, refinement, masks, phase gratings, FFT,
PROPG and multislice remain deferred. No intermediate propagated wavefunction
is produced by this unit, and no propagation equivalence is claimed.
