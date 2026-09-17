# Source-derived PARAMIN

> Clean-import provenance update: the [canonical source inventory](legacy_source_inventory.md)
> supersedes source-location, binary-only, duplicate and missing-provider claims below.
> Driver: private predecessor implementation; ONEIT: private predecessor implementation; ONEITPARA: private predecessor implementation; optimizers now
> occupy their matching .f files; RNG: private predecessor implementation; ATOMICSCATFACT:
> private predecessor implementation (all tables confirmed). The four recently corrected core
> file bodies are unchanged. Historical test results below remain historical.


`run_legacy_paramin(state, initial_vector, experimental, *, max_cycles,
tolerance, random_order=False, seed=None, phase_grating_method="fft")`
returns immutable `LegacyParaminResult`. Inputs and outputs use float64
physical active-parameter vectors. Every objective uses the existing
`evaluate_refinement_trial(..., coordinate_mode="physical")` path.
No physics, SCAL, YNORM or objective formula is duplicated.

## Source evidence

The complete readable definition is `private predecessor implementation`,
despite its artifact extension. `private predecessor implementation` is binary, not a competing readable
algorithm. Executable statements take precedence over the routine's comments
and the related-version private predecessor implementation declarations.

| Lines in private predecessor implementation | Behavior |
| --- | --- |
| 1–14 | PARAMIN(BOT,TOP,GOF,TOLER), capacities and dependencies |
| 19–33 | UNCERT, initial projection/evaluation, counters |
| 35–65 | cycle loop, fixed or rejection-sampled random order |
| 72–91 | one coordinate, first and second projected samples |
| 93–120 | literal parabola, fallback and fraction limits |
| 125–151 | third sample, callback-retained point, strict local score updates |
| 153–167 | average old range, three independent adaptation conditions |
| 169–190 | printed count, mean displacement and mandatory cycle-end reevaluation |

All driver calls are `PARAMIN(LOW,HIGH,ONEITPARA,FTOL)` at
`private predecessor implementation`. Immediately preceding each call,
HIGH=D1 and LOW=D2, so **BOT=D2 and TOP=D1**. The reference pairs 1.05,0.95
therefore sample at approximately 95% and 105% of the current value.
Driver FTOL and MAXFPARA are read at 484 and 483. ARRANGEPARS
(`private predecessor implementation`) supplies only active parameters, in the existing
mapper's order. UPDPARS is in `private predecessor implementation`.
`private predecessor implementation` supplies ONEITPARA's objective and shared best update.

## Algorithm, coordinates and traversal

This is coordinate-wise parabolic-section search with adaptive relative
ranges. It is neither simultaneous multidimensional interpolation nor a
bracketed line minimizer. PARAOPT=0 visits active parameters in packed order.
PARAOPT=1 constructs a new order each cycle: draw `int(uniform*N)`, discard
already-used indices, and continue until all N are selected. There is no
random proposal displacement or random uphill acceptance.

The physical initial vector is passed through existing UPDPARS/ARRANGEPARS
mapping before evaluation. Each cycle copies the retained vector into TEMPVAR.
For each parameter m, record STARTVAL, then perform these assignments in order,
projecting the entire vector and evaluating after each assignment:

```
TEMPVAR[m] = BOT[m] * TEMPVAR[m]          # C1
TEMPVAR[m] = TOP[m] * TEMPVAR[m] / BOT[m] # C2
TEMPVAR[m] = NEWFRAC * TEMPVAR[m] / TOP[m]# C0TEMP
```

These are sequential assignments using the **projected previous sample**,
not three independent multiples of STARTVAL. Projection can alter the effective
sample positions. With identity projection, the interpolation abscissae are
BOT, 1, TOP; they need not be evenly spaced. No SCAL/simplex conversion occurs.

## Literal interpolation and safeguards

With b=BOT[m], t=TOP[m], C0 the local retained score:

```
AA = ((C1-C0)/b**2 + (C2-C0)/(b**2*(t-1))
      - (C2-C0)/(b*(t-1))) / (1-(t+1)/b+t/b**2)
BB = (C2-C0)/(t-1) - AA*(t+1)
CC = C0-AA-BB
```

CC is unused and need not be stored. If AA>0, NEWFRAC=-BB/(2*AA).
Otherwise NEWFRAC=2+BB/(2*AA), including negative curvature. No epsilon or
near-zero-curvature tolerance is present.

The source tests `NEWFRAC/NEWFRAC == 1`. Zero, NaN and infinity fail that
IEEE test and use `1 + (-1)**(cycle+m)*(1-b)/2`, with one-based m and cycle.
Then values above 2 are capped at 2; negative values use the same fallback.
The fallback is not subjected to an additional final clipping pass. Python
locally permits IEEE division/invalid results during interpolation to reach
this explicit fallback, rather than substituting a modern line search. Equal
range factors or TOP=1 may consequently use the degenerate fallback.
Historical behavior with floating-point traps enabled is unknown.

BOT and TOP are relative factors, not absolute bounds. No clipping to the
sampled interval occurs: a valid prediction outside it is evaluated. Physical
sign folding and near-zero replacement remain delegated to the existing
mapper. Python explicitly rejects zero BOT/TOP divisors and undefined zero
STARTVAL displacement; the supported public path projects near-zero starting
values to +1e-10. Other unsafe arithmetic raises instead of returning success.

## Acceptance, adaptation and state

ONEITPARA retains each evaluated point on `score <= CHISMIN`, including ties.
Thus either side sample can be retained even if the parabolic proposal fails.
After the third sample PARAMIN sets only TEMPVAR[m]=PARMIN[m]. The separate
local BEST/C0 update checks C0TEMP, C1, C2 in that order using strict `<`.
UNCERT[m] is set to abs(1-BOT[m]) on such an improvement. UNCERT is not a
statistical confidence interval and is distinct from the final range.

Let change=abs(STARTVAL-PARMIN[m])/abs(STARTVAL) and d=abs(1-BOT[m]).
Add d to AVEDISP **before** range adaptation, then execute independent IFs:

1. If change>=d/2 and change<=0.5: delta=change.
2. If change>0.5: delta=0.5.
3. If change<d/2: delta=d/2.

Finally BOT=1-delta, TOP=1+delta. The last condition can override the apparent
0.5 cap for unusual starting ranges; the implementation preserves that order.
With no improvement and ordinary ranges, the relative half-width halves.

Every cycle ends by dividing AVEDISP by N, resetting CHISMIN to 1e20,
projecting PARMIN, and reevaluating it to reset BEST/C0. This is a reevaluation,
not a random restart. For a nondeterministic objective a worse reset score
can replace the previous retained score; there is no all-time archive.
The result records best_vector/objective, last_vector/objective, final BOT/TOP,
UNCERT, old-range average, cycles, coordinate steps, evaluations and converged.
The last evaluation is always initialization or that cycle-end reevaluation.
All result arrays are immutable; baseline state and input vector are untouched.

## Stopping and counts

AVEDISP starts at 1. The loop runs while `AVEDISP > TOLER` and
`CYCLE <= MAXFPARA`, incrementing CYCLE inside it. Hence max_cycles=0 permits
one cycle, and up to max_cycles+1 cycles can run. tolerance>=1 permits only
initial evaluation. Convergence means average **pre-adaptation relative range**
is <= tolerance, not objective convergence or a proven minimum.
There is no FMIN/MINCH test, separate evaluation budget, per-coordinate early
exit, or objective-change test in PARAMIN.

For N active parameters and C completed cycles, actual evaluations are
`1 + C*(3*N+1)`. The source print `TOTSTEPS*3+1` omits cycle-end evaluations.
Fallback uses the same third evaluation, with no extra callback. All cycles
finish completely before another stopping check.

## Python guards, RNG and capacity

Python allocates from N and explicitly supports 1..1000 active parameters.
The source's TEMPVAR, BOT, TOP, PARAMORDER and PARAMSWITCH are 1000 entries,
as are the callback parameter arrays; related shared 2000-entry declarations
do not extend that safe domain. There is no corresponding source guard.
N=0 is rejected because the source divides AVEDISP by N.

API validation requires finite nonnegative tolerance, nonnegative integer
max_cycles excluding bool, bool random_order, complete finite range factors,
float64 finite vector and an explicit nonnegative integer seed for random
traversal. Objective NaN/Inf or values above the 1e20 CHISMIN sentinel fail
explicitly, with no artificial penalties or partial successful result.
Unsupported refinement mappings propagate the existing evaluator's exceptions.
These guards are documented restrictions on unsafe source inputs.

Random traversal uses a private PCG64 stream and the source's duplicate-rejection
ordering procedure, not NumPy permutation or hidden global RNG state. Fixed
order uses no RNG. Source RANDO is clock-dependent (`private predecessor implementation`), and its
REAL*8 definition conflicts with PARAMIN's default-REAL caller unless original
compiler settings reconcile them. Python does not reseed from cycle*N or
reproduce historical timing/ABI behavior. Numeric calculations use float64.

## Driver and equivalence limits

PARAMIN appears standalone (PARA), before/after ANNEAL (PAAN/ANPA), and
before/after AMOEBA (PASI/SIPA), inside the driver's LOOPS repetitions.
The driver restores PARMIN through UPDPARS, ARRANGEPARS and UPDATE after calls.
This API implements one PARAMIN invocation, not those combinations or shared
global side effects. No other optimizer is called inside PARAMIN.

Source-derived tests cover the readable algorithm and a real Python trial
integration. They are not historical executable equivalence tests. Linked
artifact provenance, original literal kinds, RNG ABI/sequences and compiler
floating-point trap settings remain unresolved. No optimizer reference trace
has been used as a compiled-Fortran oracle.

Validation using exactly
`C:\Users\peng863\AppData\Local\miniconda3\envs\qcbed-python\python.exe`:
focused tests **26 passed in 0.83 s**, full suite **1460 passed in 27.55 s**.
Tests include an independent divided-difference scalar oracle, exact sample
ordering, degeneracy/fallbacks, non-improvement, range adaptation, resets,
seeded order generation, constraints and real trial integration.
`git diff --cached --check` passed with no staged changes. Protected source
and reference files were not edited.
