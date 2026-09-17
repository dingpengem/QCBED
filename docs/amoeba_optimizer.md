# Source-derived AMOEBA

The production implementation follows readable `private predecessor implementation`,
not SciPy. This artifact contains complete Fortran text despite its extension.
The routine identifies itself as Numerical Recipes (1986), modified 1990.
It is modified Nelder–Mead with canonical coefficients but specific legacy
acceptance and termination rules. Historical executable equivalence has not
been established.

## Evidence

| Source | Executable evidence |
| --- | --- |
| private predecessor implementation| signature, dimensions, coefficients, NUMF initialization |
| private predecessor implementation| exact vertex ranking |
| private predecessor implementation| ordered stopping tests |
| private predecessor implementation| centroid, reflection and counting |
| private predecessor implementation| expansion and replacement |
| private predecessor implementation| contraction and shrink |
| private predecessor implementation| reflection-only replacement and loop |
| private predecessor implementation| LOOPS, MAXFSIMP, FTOL, CHIMIN and PLAMDA reads |
| private predecessor implementation| SCAL, simplex construction, initial evaluations, call, postprocessing |
| private predecessor implementation| same construction/callback in hybrid stages |
| private predecessor implementation| ONEIT saves simplex coordinates, maps to physical, updates and restores input |
| private predecessor implementation| objective, shared best bookkeeping and return |

All five driver calls (894,1052,1239,1329,1510) pass
`AMOEBA(SIMP,YP,1001,1000,NPAR,ONEIT,FTOL,MAXFSIMP,CHIMIN,NUMF,SCAL)`.
There is no RNG or tolerance comparison helper in this body.

## Coordinates and construction

`build_legacy_simplex(state, initial_vector, plamda)` separates driver setup
from `run_legacy_amoeba(state, initial_simplex, experimental, *, ftol,
max_evaluations, fmin, phase_grating_method="fft")`.

The builder's initial vector is physical. SCAL is anchored to the baseline
state's packed starting values p0, using existing mapping code:
`a=1/(D2-D1)/p0`, `b=D1*p0`, `q=(p-b)*a`. The inverse is `p=q/a+b`.
D1/D2 are relative range factors, not enforced bounds. SCAL is outside the
driver's repetition loop; the caller must preserve that baseline between stages.
Packing's near-zero replacement is +1e-10; singular range factors are rejected.
There is no invented step for zero parameters or clipping to [0,1].

For N parameters, the immutable float64 simplex has shape `(N+1,N)`:
all rows initially equal q, and row j+1 adds PLAMDA[j] only to column j.
PLAMDA is an independent displacement in simplex coordinates. Signed or zero
displacements, including duplicate vertices, are preserved. The driver's
`P(NPAR+1)=0.5` is outside the coordinates copied into simplex vertices.
Python rows correspond to Fortran's first index, without a transpose.

Every initial vertex and subsequent trial is evaluated through
`evaluate_refinement_trial(..., coordinate_mode="legacy_simplex")`.
This reuses SCAL/UPDPARS, forward simulation, YNORM and the objective. It does
not call ONEITPARA or duplicate its physical-coordinate path. Unsupported
refinement categories propagate exceptions without penalties or partial results.

## Ranking and operations

Best starts at row 0 and updates on strict `<`: the first minimum wins.
Worst/second-worst start (0,1) if y[0]>y[1], otherwise (1,0).
Scan vertices in order: strict `>` replaces worst, moving its old index into
second-worst; otherwise strict `>` replaces second-worst if not the worst.
All ties preserve these source choices; there is no sorting or tolerance helper.

The centroid is freshly accumulated in vertex order excluding worst, divided
by N. There is no PSUM incremental update in this implementation's source.
Reflection alpha=1, expansion gamma=2, contraction beta=0.5 and shrink=0.5:

1. Reflect `r=2*c-w` and evaluate.
2. If fr<=best, expand `e=2*r-c`. Accept expansion if fe<**old best**,
   otherwise reflection. Expansion can be worse than reflection. This differs
   from the usual modern comparison against the reflected score.
3. Otherwise, if fr>=second-worst, first replace worst with reflection when
   fr<worst. Contract `t=0.5*w+0.5*c` using this possibly replaced worst.
   This gives outside or inside contraction without separate named branches.
   Accept only if ft<the current worst score; equality triggers shrink.
4. Shrink each nonbest vertex halfway toward the original best, in row order,
   evaluating exactly N vertices. Otherwise accept reflection alone.

## Stopping and counters

At each loop head, after ranking, test in this order:

1. `2*abs(worst-best)/(abs(worst)+abs(best)) < FTOL`.
2. `NUMF > MXF`.
3. `0.5*(worst+best) <= FMIN` (driver CHIMIN).

Only objective spread is tested; there is no parameter-spread test. The source
has no special zero-denominator branch. Python float64 0/0 yields NaN for this
test, so the comparison is false and subsequent tests still run. Finite input
objectives are required explicitly in Python; source NaN/Inf behavior is not
defined as an intentional policy. Arithmetic overflow fails explicitly.
Python also validates finite nonnegative FTOL, finite FMIN and nonnegative
integer MXF excluding bool. These are API guards, not extra stopping tests.

NUMF starts at zero after the driver's N+1 initial evaluations. Each reflection,
expansion, contraction or individual shrink evaluation increments it. There
is no mid-cycle budget check: the limit can be exceeded by a complete cycle
(up to N+2 evaluations), and equality permits another cycle. `search_evaluations`
is NUMF; `evaluations` includes the initial N+1; `iterations` counts completed
cycles. `converged` means either spread or target termination, not proof of a
minimum. `stop_reason` distinguishes these from the evaluation limit.

## Final state and capacity

AMOEBA returns without swapping the best vertex to row 0 or reevaluating it.
The immutable result exposes final simplex/objectives, best_index and
best_vector/best_objective in **simplex coordinates**. It separately exposes
best_evaluated_vector/objective using the callback's `<=` tie rule: an evaluated
reflection can remain the callback's best even if discarded by the simplex.
This diagnostic covers this call's evaluations, not an incoming legacy shared
CHISMIN or the driver's preliminary ONEITPARA evaluation. Python initializes
this diagnostic from the first finite score rather than the historical 1e20
sentinel. It does not emulate unrelated shared-state side effects.

ONEIT restores PT before returning but mutates global physical state and best
bookkeeping. The driver subsequently calls UPDPARS(PARMIN), ARRANGEPARS and
UPDATE. Python's evaluator is immutable: no baseline state changes, and no
implicit final reevaluation occurs. To obtain the retained physical result,
evaluate best_evaluated_vector explicitly with legacy_simplex mode. Repeated
driver stages, shared-state restoration and LOOPS are out of scope.

Source work vectors have NMAX=2000; driver SIMP is (2001,2000), YP(2000),
while calls pass MP=1001,NP=1000 to explicit-shape P(MP,NP),Y(MP).
This is a real leading-dimension mismatch, not merely an unresolved declaration.
ONEIT has PT/STOREPT(1000). Python dynamically allocates arrays and restricts
N to 1..1000, avoiding both empty-simplex Y(2) access and callback capacity
overflow. This does not certify any safe historical multidimensional execution
with the mismatched stride. It implements the visible indexed algorithm rather
than reproducing memory corruption. Historical link/build provenance, default
literal rounding and executable comparisons remain unresolved.

## Validation

Focused tests check exact branch trial coordinates (including expansion worse
than reflection), all contractions and shrink, ranking ties, zero denominator,
source soft budget, a separate scalar control-flow oracle, construction and
SCAL anchoring, capacity and validation, immutable inputs/results, deterministic
runs and real trial-evaluator integration. No SciPy optimizer is an oracle.

Validation interpreter:
`C:\Users\peng863\AppData\Local\miniconda3\envs\qcbed-python\python.exe`.
Focused suite: **21 passed in 6.00 s**. Complete suite:
**1434 passed in 38.10 s**. `git diff --cached --check` passed with no staged
changes. Protected source and reference files were not edited.
