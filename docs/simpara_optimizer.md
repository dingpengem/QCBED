# Source-derived SIMPARA

> Clean-import provenance update: the [canonical source inventory](legacy_source_inventory.md)
> supersedes source-location, binary-only, duplicate and missing-provider claims below.
> Driver: private predecessor implementation; ONEIT: private predecessor implementation; ONEITPARA: private predecessor implementation; optimizers now
> occupy their matching .f files; RNG: private predecessor implementation; ATOMICSCATFACT:
> private predecessor implementation (all tables confirmed). The four recently corrected core
> file bodies are unchanged. Historical test results below remain historical.


`run_legacy_simpara(state, initial_vector, experimental, *, uncertainties,
max_evaluations, high_temperature, low_temperature, cooling, seed,
random_order=False, phase_grating_method="fft")` returns immutable
`LegacySimparaResult`. Coordinates and uncertainties are physical active-vector
values; uncertainties are absolute widths in each parameter's own units.
All objective evaluations use `evaluate_refinement_trial` in physical mode.

## Source and scope

The complete readable routine is `private predecessor implementation`.
No second readable SIMPARA definition was found. It is a simultaneous stochastic
search with parabolic proposals, sign-direction restrictions and adaptive
temperature. All coordinates participate, although zero widths or zero
direction components can freeze individual coordinates. It is not coordinate
descent, a Euclidean orthogonal-direction optimizer, or Metropolis annealing.

| Executable lines | Evidence |
| --- | --- |
| 48–59 | argument types, 2000-entry local arrays, RNG declarations |
| 73–104 | initialization, first evaluation, cycle condition |
| 125–145 | listed or duplicate-rejection random coordinate ordering |
| 169–225 | sign-direction rules and simultaneous plus/minus proposals |
| 231–281 | two evaluations, coordinate parabola, conditional third evaluation |
| 286–319 | strict local acceptance and no-improvement counter |
| 387–415 | periodic reevaluation and direction reset |
| 417–430 | temperature adaptation, widths, C0 and NFUNC |

The sole driver call is `private predecessor implementation`:
`SIMPARA(ONEITPARA,MAXFSIMP,MAXTEMP,MINTEMP,COOLRATE,UNCERTAINTIES,P,PARMIN,NPAR,PARAOPT,NFUNCEVAL,LOUD,13)`.
At 506–513 PLAMDA values become absolute UNCERTAINTIES for SPAR. The call
is inside SPAR's LOOPS loop (583–605); it is not part of the enumerated hybrid
methods. The driver first evaluates ONEITPARA(P), then restores PARMIN through
UPDPARS/ARRANGEPARS/UPDATE afterward. This API implements one invocation,
excluding that preliminary driver evaluation and all driver sequencing/output.

## Explicit startup policy

UNITVEC is read at line 171 and SINCECHANGE may be read at 189 before either
is initialized in the readable source. The user explicitly approved initializing
both to zero in Python. This is a deterministic startup policy, **not evidence
that historical Fortran locals were zero-initialized**. No other startup values
are inferred from uninitialized storage.

## Proposal and RNG order

OUTVAR initially copies the supplied vector; CENTRAL copies the evaluated
OUTVAR. PERTSIZE[j]=CURRENTTEMP*UNCERT[j], initially HITEMP*UNCERT[j].
POPT=0 visits coordinates in packed order. POPT=1 draws `int(uniform*N)` and
rejects duplicates until every coordinate is selected. This only changes the
assignment of random draws; objective samples still change coordinates together.

UNITVEC is approximately a sign vector, not normalized to Euclidean length.
After an improvement, in parameter order:

```
UNITVEC[j] = (OUTVAR[j]-CENTRAL[j]) / (abs(OUTVAR[j]-CENTRAL[j])+1e-30)
PERTVEC[j] = UNITVEC[j] * abs(normal(0, PERTSIZE[j]))
CENTRAL[j] = OUTVAR[j]
```

Without improvement, zero direction or an odd SINCECHANGE draws unrestricted
normal components, repeating the entire vector until
`abs(sum(UNITVEC[j]*PERTVEC[j]/(abs(PERTVEC[j])+1e-20))) <= 1.5`.
With zero direction the first draw always passes. Even SINCECHANGE with nonzero
direction uses the sign-aligned absolute normal rule. This is an approximate
sign-direction condition, not exact perpendicularity. There is no artificial
rejection-attempt cap or extra objective evaluation during direction sampling.

The two full vectors are CENTRAL-PERTVEC and CENTRAL+PERTVEC; evaluate them
in that order. Perturbations are additive absolute widths. Negative and zero
parameter values are supported without multiplicative scaling or a near-zero
floor. Zero uncertainties are permitted and can produce degenerate parabolas.
D1/D2 and SCAL play no role.

A private seeded PCG64 stream replaces clock-dependent RANDO. Normal sampling
reuses the existing source polar-normal transform from `anneal._normal`, based
on readable `private predecessor implementation`: two uniforms, rejection outside the unit disk or
at zero radius, keep the first normal and discard the second. No URANDO call
occurs. No existing optimizer implementation was changed.

The logical SEED counter starts at zero, increments for each ordering draw,
each normal request (including rejected direction vectors), then each parabola
coordinate. Its parity selects fallback signs. Internal polar rejection draws
do not advance this outer counter. Historical RANDO also mutates its seed using
clock values (`private predecessor implementation`); Python intentionally omits those clock mutations.
Consequently historical fallback parity as well as random samples may differ.

## Parabolic proposal

For coordinate i, write x=CENTRAL[i], l=INVTEMPVAR[i], r=TEMPVAR[i]. Every
coordinate uses the same total objective samples C0, C1, C2:

```
AA = ((C1-C0)*(r-x) - (C2-C0)*(l-x)) /
     ((l*l-x*x)*(r-x) - (r*r-x*x)*(l-x))
BB = ((C2-C0) - AA*(r*r-x*x)) / (r-x)
NEWVAR[i] = -BB/(2*AA)
```

BB and the minimum are used only when AA>0. The source's `NEWVAR/NEWVAR==1`
test rejects zero, NaN and infinity. Rejected minima, or nonpositive/NaN
curvature, use `x + (-1)**SEED*(r-x)/2`. Each valid nonzero minimum increments
PARAFLAG. Inside the positive-curvature branch only, values above 10000 are
capped at 10000, including fallback values. There is no lower cap or general
clipping. A fallback from the nonpositive-curvature branch is not capped.
No near-zero curvature tolerance is introduced.

Evaluate the combined NEWVAR only if PARAFLAG>=1; otherwise set C3=1e20 without
evaluating. Thus a mixed vector may contain both parabolic and fallback
components. Python permits IEEE invalid/division results locally to reach the
source fallback; it does not introduce an alternative interpolation formula.

## Acceptance and callback aliasing

SIMPARA checks C1, C2, C3 in order against the updated CBEST using strict `<`.
It copies the winning vector to OUTVAR and increments its improvement counter
once for an improving cycle. There is no tolerance helper, acceptance
probability or deliberate uphill acceptance.

The actual driver passes PARMIN as OUTVAR. ONEITPARA (`private predecessor implementation`)
also writes PARMIN on `score <= CHISMIN`, including ties, **before** SIMPARA's
strict comparisons. Python explicitly reproduces that write to its internal
OUTVAR copy. A tie can therefore move the returned vector without setting
CHANGEFLAG or updating CENTRAL. CBEST and callback CHISMIN are separate;
periodic CBEST reevaluation does not reset callback CHISMIN. This differs from
a pure standalone callback with no alias, and from ANNEAL/PARAMIN resets.
The fresh Python invocation starts callback CHISMIN at 1e20; it does not import
an incoming driver's prior global minimum.

ONEITPARA's UPDPARS can mutate PT for vacancy, absorption and Debye-Waller
categories (`private predecessor implementation`). Python copies back the values from
the existing mapped trial state for those categories before interpolation;
it does not repack the entire vector. Unsupported vacancy/free-surface/gas
paths still fail through the existing evaluator, without penalties. Other
physical transformations do not overwrite the raw optimizer coordinates.
No physics or parameter transformation equations are duplicated.

## Temperature, resets and stopping

Every NPARAM cycles, reevaluate OUTVAR, set CBEST to that score, copy
CENTRAL=OUTVAR and zero UNITVEC through their difference. CHANGEFLAG is then
zeroed. Consequently that cycle **always cools**, even if it improved.
Otherwise improvement heats by dividing CURRENTTEMP by COOLING; failure
multiplies it by COOLING. Update all widths to UNCERT*CURRENTTEMP afterward,
then set C0=CBEST. HITEMP is not an upper bound. There is no independent
per-coordinate temperature adaptation; widths differ through UNCERT only.
Resetting direction does not randomize the retained starting point.

The loop condition is `TOTSTEPS <= MAXITS` and
`CURRENTTEMP >= LOTEMP*COOLING`. Initialization counts as one evaluation.
Each cycle adds two, plus one if a parabola qualifies, plus one every NPARAM
cycles. Checks occur only at cycle boundaries: equality permits another full
cycle, so up to four evaluations can be added at the nominal limit. MAXITS=0
still evaluates initialization, but runs no cycle. For NPARAM=1 every cycle
resets and cools. There are no FTOL, FMIN, objective-spread, parameter-spread,
or no-improvement stopping tests. NFUNC is the complete callback count.

No unconditional final reevaluation occurs. last_vector/objective can describe
a rejected trial, while best_vector/objective represents OUTVAR/CBEST. The
result also exposes CENTRAL, UNITVEC, next perturbation widths, temperature,
cycles, evaluations, strict improvement count and stop_reason. If both limits
fail, Python labels the reason max_evaluations. A temperature stop is not proof
of convergence. Periodic reevaluation can replace CBEST with a worse score;
there is no independent all-time archive. Results and inputs are immutable
snapshots, with no mutation of the caller's refinement state.

## Capacity and limitations

SIMPARA's own arrays have capacity 2000 with no source guard. The actual
ONEITPARA callback has PT(1000), so this integrated API rejects counts outside
1..1000 and allocates dynamically. Zero parameters would also make the source
MOD(CYCLE,NPARAM) invalid. Finite positive temperatures, 0<COOLING<1, finite
nonnegative uncertainties, nonnegative integer limit/seed excluding bool,
and finite float64 vectors are explicit Python validation requirements.
Nonfinite objectives or objectives above the source 1e20 sentinel fail
explicitly. Other unsafe arithmetic propagates an error, not partial success.

SIMPARA declares REAL*8 RNG results, but its default-real 0.0 actual argument
to NRANDO conflicts with the readable REAL*8 dummy unless compiler settings
reconcile them. Historical build/link provenance, literal rounding, seed/time
behavior, uninitialized startup storage and floating-point trap settings remain
unresolved. The approved startup and explicit float64/RNG choices do not
establish executable-Fortran equivalence. No compiled optimizer trace is used.

| Optimizer | Proposal | Retention and adaptation |
| --- | --- | --- |
| ANNEAL | simultaneous multiplicative normals around retained best | ties accepted; per-coordinate fractional amplitudes and scalar cooling/reset |
| PARAMIN | three sequential samples for one coordinate | callback ties; strict local score updates; relative-range adaptation; random ordering optional |
| SIMPARA | simultaneous additive plus/minus normals and optional combined parabola | strict local updates plus driver callback ties; sign-direction rules; global heating/cooling |

The refinement driver, serialization, CLI and other numerical kernels remain
out of scope. Tests use controlled draws, exact branch examples, a separate
scalar symmetric-interpolation oracle and a real trial-evaluator integration;
they validate source-derived Python behavior rather than historical execution.

Validation interpreter:
`C:\Users\peng863\AppData\Local\miniconda3\envs\qcbed-python\python.exe`.
Focused suite: **24 passed in 1.05 s**. Complete suite:
**1484 passed in 38.72 s**. `git diff --cached --check` passed with nothing
staged. Protected source/reference files and existing optimizer implementations
were not modified.
