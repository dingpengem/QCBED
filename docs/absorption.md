# Absorptive atomic and structure factors

## Decision and source audit

**Category A: source-backed and implementable on the explicitly supported
domain below.** The readable executable equations are complete and agree.
This does not establish an empirical calibration domain or executable
equivalence. The zero-Debye-Waller case is singular in the written expression
and is rejected, not assigned a speculative limiting value.

The relevant readable source is:

- `private predecessor implementation`: formatted absorption flag READ, per-species INELA,
  INELB, INELC and masks, and descriptions of the phenomenological model.
- `private predecessor implementation`, SFMAP: atom initialization, vacancy/surface corrections,
  elastic terms, and FAATM/AFRE/AFIM at lines 1901-1911.
- `private predecessor implementation`, SFMAP2: the same FAATM/AFRE/AFIM equations at 339-349;
  reflection ABSFAC multipliers later; PMAP assembly at 561.
- `private predecessor implementation`, SFLIST: prints separate absorptive amplitude/phase,
  thresholds small amplitudes, and formats E14.7/F9.4 output.
- `private predecessor implementation`, ARRANGEPARS, and `private predecessor implementation`, UPDPARS: transfer
  refinement parameters, including INELA/B/C and reflection ABSFAC. UPDPARS
  at 923-944 uses absolute values for refined coefficients and ties C=B
  when MASKINELC=-1. These operations are outside the new module.
- `private predecessor implementation`, PGMAP: consumes the transformed complex potential;
  inspected to establish signs, not translated.

ABSFAC is a **reflection/region refinement array**, not an atomic-scattering
routine. The driver reads it with HK and SFFAC at `private predecessor implementation`; UPDPARS
updates it at `private predecessor implementation`. No readable ABSFAC function or alternative
absorptive atomic routine exists. FAATM is an inline scalar expression.
The search covered readable `.f` files by names and variables (ABSORB,
INELA/B/C, FAATM, AFRE/AFIM, ABSFAC, imaginary and diffuse-scattering terms).
The source's Bird and King description identifies motivation for this fitted
form, not a separate numerical thermal-diffuse-scattering integration. Its
Nakashima citation has placeholder journal/page details; no external formula
has been substituted. Binary artifacts with `.f` names do not resolve linkage.
private predecessor implementation supplies no authoritative numerical dependency for this module.

SFMAP is called at `private predecessor implementation`. SFMAP2 is called at
`private predecessor implementation`,
`private predecessor implementation` (ONEIT), and `private predecessor implementation` (ONEITPARA).
The latter two paths immediately call OSHIFT, FFT(-1), then PGMAP.
SFMAP reads/materializes atoms and sets applied parameters; SFMAP2 consumes
that state and applies later symmetry-selected refinement multipliers.
The raw absorptive equations are identical in both.

## Exact mathematical expression

Let Z=ZT(ITYPE), s=SATM=STL=|g|/2, B=APPLIEDDW, and a,b,c denote the
applied INELA, INELB, INELC weights (not cell lengths or the thermal B).
With `V=1000*voltage_kv`:

```text
q = exp(-0.0041225*sqrt(V))
p = 0.64579/(1 + 0.079647*Z^0.40196)
E = 0.0012822*(0.5 + q/(1-q))
              *(1-exp(-0.86245*Z^0.377))*Z^1.6191
D = B^p
T0 = a*exp(-1.0e12*s^2)
T1 = b*exp(-2.7964*s^1.83/(1+0.18088*Z^0.48919))
T2 = c*(0.27183*Z^0.42387 + 0.45264*Z^(-0.2099))*s^2
        / ((1+0.60054*Z^0.63377+3.5729*s^4)*D)
FAATM = E*D*(T0+T1+T2)
```

These are all constants in the atomic equation. There are no coefficient
tables, interpolation boundaries, integration variables/bounds, quadrature,
or convergence controls. INELA/B/C come from input or applied refinement
state, not atomic-number lookups. The function contains no wavelength.
The driver explicitly sets VOLTS=KV*1000 at `private predecessor implementation`.
Python uses kV at the API boundary and performs that conversion once.

At s=0 the ordinary finite expression gives FAATM=E*D*(a+b); T2=0.
The source contains no special zero-g branch. The first term is described in
source output as a Dirac delta, but is actually a finite-width Gaussian:
at s=1e-6 its attenuation is exp(-1), not zero. Python preserves this term.
At B=0, D=0 and the written T2 divides by zero, even for c=0. Python requires
B>0 and does not cancel D to manufacture a zero-temperature extension.

For atom m, the source sum is:

```text
area = cell_a*cell_b*sqrt(1-cos(gamma)^2)*dilation*dilation
phi_m = 6.283185307*(H*x_m + K*y_m)
OCC_m = starting_occupancy_m*(1-vacancy_fraction)
A(H,K) = AFRE+i*AFIM
       = -sum_m (cos(phi_m)+i*sin(phi_m))
                 *area*FAATM_m*exp(-B_m*s^2)*OCC_m
```

The area **multiplies** absorption, whereas it divides the elastic factor.
No inverse-area correction, additional elastic conversion constant, thickness,
or wavelength is introduced. Dilation affects area and s through the existing
reciprocal grid. B affects both FAATM and the external thermal damping;
omitting either is incorrect. Occupancy/vacancy enter outside FAATM, once.
Thermal and coefficient parameters are species-specific; occupancy and phase
are atom-specific. Atoms are summed in input order with separate real/imaginary
subtractions, retaining CN*area*FAATM*damping*OCC ordering.

For nonnegative weights, FAATM is nonnegative real. Signed unrefined weights
are accepted by the actual READ and may produce negative FAATM. The API keeps
them signed rather than applying the refinement-only DABS operations.
The stored sum A is generally complex; at the origin with nonnegative
weights/occupancies it is nonpositive real. Both reciprocal signs are
evaluated directly, yielding A(-g)=conj(A(g)) for real weights.

SFMAP stores AFRE/AFIM in SFS slots 3/4. SFMAP2 combines them as:

```text
PMAP = DCMPLX(SFRE-AFIM, SFIM+AFRE) = elastic + i*A
```

After OSHIFT and FFT(-1), PGMAP uses PR=real(PMAP), PI=imag(PMAP),
`exp(SIGMA*PI)*(cos(-SIGMA*PR)+i*sin(-SIGMA*PR))`. This explains the loss
sign for negative real-space PI. Neither the transform nor exponentiation
is implemented. The combined reciprocal array need not itself satisfy
Friedel conjugacy. Slice translations belong to PROPG, not either atomic sum.

## Units and numerical domain

Confirmed input conventions: voltage in kV (converted to volts), s in inverse
angstroms, fractional x/y, area in angstrom^2, dimensionless vacancy/occupancy.
B has angstrom^2 units inferred from exp(-B*s^2). The fitted constants assume
these numeric units; their individual dimensional units are not specified.
INELA/B/C are applied numeric weights in the source parameterization.
The exact return contract is **source FAATM scaling**, not a conventional
electron scattering length. If the elastic coefficient's projected-potential
interpretation is V angstrom, A has those same units and FAATM must have
V/angstrom before multiplying area. This is a dimensional inference from
assembly, not an explicit unit label in the source. No unit repair is made.

The API deliberately limits itself to the existing joint SFMAP domain:
integer Z=1..98, s=0..6 inverse angstroms (hydrogen:0..1.5). These bounds
come from the elastic function called alongside absorption, **not** from an
absorption table or a documented calibration range of this fitted expression.
The absorption formula itself states no finite upper coordinate, atomic-number,
energy or B calibration bounds. No empirical accuracy claim outside the saved
case is implied by arithmetic tests over the chosen domain.

Voltage and applied B must be finite and positive; vacancy must be in [0,1];
starting occupancy must be finite and nonnegative (values above one are not
clipped). Coordinates and weights may be signed but must be finite. Bool and
wrong scalar types are rejected. Nonrepresentable arithmetic raises ValueError.
Ordinary exponential underflow in vanishing tails is permitted. The energy
denominator uses `-expm1(t)` instead of `1-exp(t)` to avoid cancellation at
small positive voltage; the numerator and mathematical expression are unchanged.
Printed constants are interpreted as binary64, without emulating unknown
default-kind Fortran rounding. Recovered S2 uses the grid's float64 STL squared,
with the same last-bit limitation documented for the elastic unit.

## Public API and integration boundary

```python
AbsorptionCoefficients(inela: float, inelb: float, inelc: float)

absorptive_atomic_scattering_factor(
    atomic_number: int, stl_inverse_angstrom: float, *, voltage_kv: float,
    applied_debye_waller_angstrom_squared: float,
    coefficients: AbsorptionCoefficients,
) -> float

absorptive_structure_factors(
    species: Sequence[Species], atoms: Sequence[AtomSite],
    reciprocal_grid: ReciprocalVectorGrid, *, voltage_kv: float,
    applied_debye_waller_angstrom_squared: Sequence[float],
    coefficients: Sequence[AbsorptionCoefficients],
    vacancy_fraction: float, absorption_raw: str,
) -> NDArray[np.complex128]

combine_elastic_absorptive_structure_factors(
    elastic: NDArray[np.complex128], absorptive: NDArray[np.complex128],
) -> NDArray[np.complex128]
```

All are exported from qcbed. The frozen coefficient model preserves applied
values; no masks or coefficient tying are interpreted. The existing parser
consumes but does not retain INELA/B/C and B, so callers must supply them
explicitly. No parser changes or guessed defaults were introduced.
Use the same grid, atom order, vacancy and applied B with the elastic API.
Both sums reject unresolved G atoms because SFMAP randomizes their input
positions and occupancies. Species indices use the existing checked conversion.

The absorption flag follows `ABSORB(1:2).EQ.'YE'` literally, also used by
LegacyInput.absorption_enabled. `YES_ABSORPTION` and `YEanything` enable it;
`yes`, ` YE`, empty text and `NO` do not. No strip/uppercase or include-only
character width is imposed. The summation validates required input parameters
even when disabled, then returns a zero complex128 grid.

Arrays have shape `(mesh_y,mesh_x)` and C order, K/y on axis 0, H/x on axis 1.
Both full centered and unshifted grids (including generic odd grids) are
supported without transpose. H/K are dimensionless, separate from physical s.
PMAP flattening is row*mesh_x+column. This raw formula does not apply DSAP,
negative-Nyquist exclusions or ABSFAC/SYMOPT refinement, and is not complete
SFMAP/SFMAP2. The combination API requires identical nonempty 2D complex128
NumPy arrays; finite values and shapes are checked, with no implicit dtype
conversion or broadcasting. Callers must ensure both grids have the same
physical coordinates and ordering. Inputs are never mutated.

## Validation and unresolved issues

Independent references use 65-digit Decimal arithmetic. The evaluator uses
energy=(1+q)/(2*(1-q)) and distributes D through the bracket so it cancels
from the third term; production retains the written source grouping.
Complex phases use independent Decimal Taylor sine/cosine. Tests contain
hard-coded expected values, not calls through another implementation path.

With 200 kV, B=.86 and weights (2,.5,1.5):

| Z | s=0 | s=1e-6 | s=.25 | upper joint-domain s |
| --- | --- | --- | --- | --- |
| 1 | .0011644736797559489 | .0005756034771997381 | .00021455287057936757 | .00006491914736578576 (1.5) |
| 13 | .11615494071626605 | .057415799886858502 | .021531909702254684 | .0006275195646752409 (6) |
| 98 | 3.4351027768917547 | 1.6979835073148182 | .6567876701842111 | .03531580129877369 (6) |

Atomic references use relative tolerance 3e-14. Structure cases use absolute
5e-13: for Al, area=200, s=.025, B=.86, weights (0,1,1), occupancy=1,
A=-9.270486844482137. At x=1/4 this becomes
-4.1621351794569144e-10 - 9.270486844482137i. Fractional occupancy, vacancy,
species-specific weights, zero g, positive/negative H/K, rectangular layouts,
dilation, oblique geometry, flags, input validation and immutability are tested.
There are no integration-convergence or interpolation-boundary tests because
neither operation exists in this source expression.

The saved FOR0013.DAT listing at lines 900-926 supplies 25 absorptive complex
coefficients through amplitude and phase. FOR007/FOR0013 confirm Al,
200 kV, weights (0,1,1), applied B=.86, vacancy=0 and dilation=1. Tests compare
those absorptive columns (not only intensity) using relative 3e-7 and absolute
5e-8 tolerance, allowing printed precision, literal rounding and SFLIST's
thresholded zeros. This is a saved-output comparison, not executable rerunning.

Measured errors for the table above: maximum absolute 2.22e-16 at Z=98,s=1e-6;
maximum relative 3.46e-16 at Z=13,s=6. For the four primary complex sum
cases the maxima are 1.80e-15 absolute and 1.94e-16 relative, at x=1/4.
For saved output, maxima are 5.40e-8 absolute and 1.94e-7 relative at H=0,K=-2
(a tied symmetry case). Relative error excludes zero reference coefficients;
their maximum absolute residual is 3.51e-11.

One surrounding source issue remains: UPDPARS updates KV/WAVEL/SIGMA at
`private predecessor implementation` but does not update VOLTS there. Only the initial driver
assignment to VOLTS was found. Consequently historical energy-refinement
behavior may use stale voltage in absorption. This stateless API uses its
explicit voltage; it does not implement or silently repair that workflow.
Empirical fit range, exact dimensional calibration, original linkage/include,
and compiler literal rounding remain unverified. Exact executable equivalence
is not claimed. Masks, refinement, gas randomization, phase gratings, FFT,
propagation, detectors and parallelism remain outside this implementation.
