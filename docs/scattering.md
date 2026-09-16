# Elastic atomic scattering factor

## Definition and call path

The only readable definition of `ATOMICSCATFACT(ATNUM,STL,OPTION)` is in
`private predecessor implementation`, lines 1–1638. All readable `.f` files were searched.
Despite its filename, `private predecessor implementation` defines SIMPARA, not ATOMICSCATFACT.
Some other `.f` files are binary artifacts; they do not supply readable
competing definitions. Object/executable linkage has not been reconstructed.
There is no conflicting readable active definition.

The source call path is:

- PROGRAM MULTIS calls SFMAP in `private predecessor implementation`; its line 1887 calls
  `ATOMICSCATFACT(ZT(M),SATM,2)/VOL`.
- The driver and ONEIT/ONEITPARA call SFMAP2 in `private predecessor implementation`; its line 323
  uses the same function arguments and division.

Both paths calculate `S2=0.25*|g|^2` and `SATM=sqrt(S2)` before the call.
ATNUM is the neutral atom's atomic number. STL is
`sin(theta)/lambda = |g|/2`, in inverse angstroms, not the reciprocal-vector
magnitude itself or its square. OPTION=1 is the source X-ray factor branch;
OPTION=2 is the electron branch used by both current callers.

Absorption is calculated separately in the callers' `ABSORB`/`FAATM` blocks.
Gas positions and occupancies are randomized separately in SFMAP. Neither
absorption nor random-number behavior belongs to ATOMICSCATFACT.

## Public API and domain

```python
elastic_atomic_scattering_factor(
    atomic_number: int, stl_inverse_angstrom: float, *, option: int = 2
) -> float
```

The function is exported from `qcbed`. It returns a Python float using
binary64-compatible arithmetic. The source-supported atomic range is 1–98
(H through Cf). The documented STL domain is 0–6 inverse angstroms, except
hydrogen is restricted to 0–1.5 inverse angstroms. Negative STL is rejected;
the reflected interpolation point near zero does not extend the public input
domain to negative coordinates. Bool, non-real/non-integer inputs, non-finite
coordinates, and out-of-range inputs are rejected with named-argument errors.

Only OPTION=2 is supported by this elastic-electron API. OPTION=1 is confirmed
as X-ray behavior but intentionally not exposed here; all other options are
also rejected. The original function's fallback for arbitrary OPTION values
is not a supported API contract.

The output retains the source's potential-factor scaling, including 47.878009
and 1.145896 below. It is **not** an unscaled electron scattering length.
Interpreting 47.878009 as the scattering-length-to-potential conversion gives
V angstrom^3; the caller divides by the two-dimensional cell area, giving a
projected-potential contribution. This dimensional interpretation is an
inference from the scaling and caller use, not an explicit output-unit label
in the function. The source-scaled numeric return is the API contract; no
volume/area division or additional unit conversion occurs in this unit.

## Tables and mathematical calculation

The function contains tabulated data, not a Gaussian coefficient fit:

- STLARRAY: 62 ordered sampling coordinates.
- EFARRAY and XFARRAY: 98 rows of 62 electron and X-ray values.
- EMARRAY: 98 method selectors (1: direct RHF, 2: Mott conversion, 3: hydrogen HF).

All values were copied literally from current readable `private predecessor implementation`, including
negative sentinels: -1 denotes missing entries and -2 denotes entries rejected
by the source authors. Both are excluded from interpolation. Source comments
attribute the data to International Tables for Crystallography, chapters 4.3
and 6.1, and Doyle & Turner (1968), Acta Cryst. A24, 390–397. They also describe
manual corrections to suspected table/OCR errors. The Python copy preserves
those existing choices; it does not substitute external tables or repair data.
No table or coefficient comes from the related-version private predecessor implementation.

Let P(s) be the unique cubic through the four source-selected usable points:

```text
P(s) = A*s^3 + B*s^2 + C*s + D
```

Point selection follows lines 1503–1585: the second point is the last usable
sample at or below s; neighboring usable points bracket it. Near zero, the
first point is reflected using the source's one-based `PP=ABS(PT(1))+2`
rule. At the upper table boundary, the final four usable samples are used.
This is a local four-point interpolating polynomial, despite the source's
"cubic spline" terminology; it is not a global spline fit.

For method 1 or 3, interpolate EFARRAY. For method 2, interpolate XFARRAY:

```text
direct:        result = 47.878009 * P(s)
Mott, s > 0:   result = 1.145896 * (Z - P(s)) / s^2
Mott, s = 0:   result = 47.878009 * EFARRAY[Z, first sample]
```

The production implementation preserves the source's F1–F4, D13–D33,
DT1–DT3, D21/D22 and A/B/C/D coefficient evaluation ordering and grouping
(lines 1591–1634). Tables are immutable tuples indexed by Z-1 in Python.
They are finite scientific lookup data, not legacy simulation array capacities.
No SciPy spline or physics outside this function is introduced.

## Independent numerical references

Expected values were evaluated independently using exact rational arithmetic
on the printed decimal table entries and the Lagrange interpolation formula.
The reference evaluator did not call the implementation or reuse its cubic
coefficient algebra. Tests store the resulting values as literals.

| Z | STL (inverse angstroms) | Expected source-scaled return |
| --- | --- | --- |
| 1 | 0 | 25.327466761 |
| 1 | 0.025 | 24.824622984184895 |
| 1 | 0.375 | 6.886204263203125 |
| 1 | 1.5 | 0.526658099 |
| 13 | 0 | 281.953595001 |
| 13 | 0.025 | 275.340445007875 |
| 13 | 0.375 | 48.5632630038125 |
| 13 | 6 | 0.430902081 |
| 98 | 0 | 806.313549569 |
| 98 | 0.025 | 802.5855584 |
| 98 | 0.375 | 246.79875153955555 |
| 98 | 6 | 2.9426290975555554 |

Reference stencils, with table values in coordinate order:

| Z | Test STL | Coordinates | Values |
| --- | --- | --- | --- |
| 1 | 0, 0.025 | -0.04, 0, 0.04, 0.05 | 0.510, 0.529, 0.510, 0.505 |
| 1 | 0.375 | 0.35, 0.36, 0.38, 0.40 | 0.160, 0.153, 0.141, 0.130 |
| 1 | 1.5 | 1.2, 1.3, 1.4, 1.5 | 0.017, 0.014, 0.012, 0.011 |
| 13 | 0 | -0.01, 0, 0.01, 0.02 | 5.867, 5.889, 5.867, 5.800 |
| 13 | 0.025 | 0.01, 0.02, 0.03, 0.04 | 5.867, 5.800, 5.692, 5.547 |
| 13 | 0.375 | 0.35, 0.36, 0.38, 0.40 | 1.111, 1.070, 0.997, 0.932 |
| 13 | 6 | 3.5, 4, 5, 6 | 0.024, 0.019, 0.012, 0.009 |
| 98 | 0.025 | 0.01, 0.02, 0.03, 0.04 | 97.929, 97.718, 97.375, 96.912 |
| 98 | 0.375 | 0.35, 0.36, 0.38, 0.40 | 69.707, 68.898, 67.325, 65.810 |
| 98 | 6 | 3.5, 4, 5, 6 | 11.695, 9.933, 7.049, 5.553 |

For Z=98 at zero, the independent expectation uses EFARRAY(98,1)=16.841
and the explicit zero-coordinate branch. Aluminium (Z=13) is present in
reference_data/FOR007.DAT. Hydrogen exercises missing entries and reflected
origin interpolation; Cf exercises the Mott branch and upper atomic bound.

Tests use relative tolerance 2e-12 and absolute tolerance 1e-12, allowing
roundoff in the source coefficient algebra and Mott subtraction. Measured
maximum errors across these 12 cases were 1.10e-11 absolute and 1.36e-14
relative. An additional provenance test checks every copied table entry,
including ordering and sentinels, against readable source.

## Remaining validation limits

Exact executable-Fortran equivalence is not established. The original compiler
flags may change default-kind literal rounding before assignment to REAL*8;
these tests instead interpret the printed decimals directly as binary64.
No standalone ATOMICSCATFACT executable-reference values have been identified
in the supplied reference outputs. The tests do not certify full SFMAP output.
Near-zero Mott subtraction/division retains the source's cancellation and
underflow risks; no analytic limit or numerical stabilization was substituted.
Absorption, gas randomization, structure-factor summation, phase gratings,
FFT and propagation remain outside this implementation.
