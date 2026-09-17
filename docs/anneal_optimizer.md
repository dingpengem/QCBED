# Source-derived ANNEAL

> Clean-import provenance update: the [canonical source inventory](legacy_source_inventory.md)
> supersedes source-location, binary-only, duplicate and missing-provider claims below.
> Driver: private predecessor implementation; ONEIT: private predecessor implementation; ONEITPARA: private predecessor implementation; optimizers now
> occupy their matching .f files; RNG: private predecessor implementation; ATOMICSCATFACT:
> private predecessor implementation (all tables confirmed). The four recently corrected core
> file bodies are unchanged. Historical test results below remain historical.


`run_legacy_anneal(state, initial_vector, experimental, *, max_iterations,
seed, maximum_temperature, minimum_temperature, cooling_rate,
phase_grating_method="fft")` returns immutable `LegacyAnnealResult`.
Vectors are float64 active **physical** parameters in the existing refinement
layout order, not simplex-scaled coordinates. Every objective evaluation goes
through `evaluate_refinement_trial`; no forward physics or objective is copied.

## Executable evidence

The current readable `private predecessor implementation` is authoritative:

| Lines | Behavior |
| --- | --- |
| 1–13 | `ANNEAL(GOF,MAXITS)`, local arrays, declarations and include |
| 15–29 | CHISMIN sentinel, UPDPARS/ARRANGEPARS projection, initial GOF and temperatures |
| 31–40 | loop condition, increment, normal proposals and fractional amplitudes |
| 42–54 | projection, evaluation and literal `FIT <= BESTFIT` adaptation |
| 71–81 | count small amplitudes before cooling failed proposals |
| 82–94 | scalar cooling, amplitude reset and retained-point reevaluation |
| 95–103 | independent reevaluation every tenth iteration |

`private predecessor implementation` calls
`ANNEAL(ONEITPARA,MAXFSIMP)`. The callback implementation in
`private predecessor implementation` updates CHISMIN/PARMIN on `CSQ <= CHISMIN`, then
unpacks/repacks the retained state; line 250 returns CSQ. The existing mapper
provides the UPDPARS/ARRANGEPARS transformations. Related-version private predecessor implementation
does not override these executable statements.

The readable artifact `private predecessor implementation` defines NRANDO's polar-normal transform;
lines 52–77 define its clock-dependent uniform generator. Its file extension
does not describe its actual readable contents.

## Algorithm and states

This is adaptive Gaussian random search around a retained best point, with
greedy acceptance including ties. It is not standard simulated annealing:
there is no Metropolis probability, tolerance helper, or uphill acceptance.
Gaussian proposals support heuristic global exploration, but offer no global
minimum guarantee and do not form a walk through worse states.

For each coordinate, before parameter projection:

```
proposal[m] = best[m] + (best[m] * amplitude[m]) * normal[m]
trial_amplitude[m] = 1.48258 * abs((proposal[m] - best[m]) / best[m])
```

All coordinates are proposed each iteration. The proposal is passed through
UPDPARS/ARRANGEPARS-equivalent mapping before evaluation. On `FIT <= BESTFIT`,
BESTFIT and all amplitudes are replaced; otherwise amplitudes are multiplied
by COOLRATE. Retaining a vector is the callback's separate CHISMIN comparison,
emulated explicitly by the optimizer because the Python evaluator is stateless.
The original acceptance branch does not itself copy the best vector.

`best_vector`/`best_objective` represent PARMIN/CHISMIN after the most recent
reset. `adaptation_objective` represents local BESTFIT. `last_vector` and
`last_objective` describe the last evaluation, which may be a rejected proposal
or a reset reevaluation. There is no separate all-time archive: resetting
CHISMIN to 1e20 permits a changed reevaluation of the retained point to replace
its previous score. No random restart point is introduced. Result arrays are
immutable snapshots; the input vector and baseline state are not mutated.

## Temperatures and stopping

Scalar temperature and all amplitudes start at MAXTEMP. The number of
amplitudes strictly below MINTEMP is counted **before** failure cooling. When
all are below it, scalar temperature is multiplied by COOLRATE, all amplitudes
are reset to that scalar, CHISMIN is reset, and the retained point is projected
and reevaluated. Every tenth iteration independently resets CHISMIN and
reevaluates that point; both branches can evaluate on the same iteration.

The loop starts while `steps <= MAXITS` and scalar temperature `> MINTEMP`,
then increments steps. Thus `max_iterations=0` permits one proposal and a
nonnegative MAXITS permits up to MAXITS+1 proposals. Initialization always
evaluates, even if MAXTEMP <= MINTEMP. Evaluations also count every reset;
they are not the iteration count. Empty active layouts follow the source's
vacuous all-amplitudes-small condition. There is no objective convergence test.

## RNG, domain and limitations

An explicit required integer seed creates a private NumPy PCG64 generator.
Pairs of uniform draws are mapped to (-1,1); reject squared radius zero or
at least one, then use `sqrt(-2*log(radius)/radius)` times the first coordinate.
The second normal is discarded as in readable NRANDO. Negative parameter
centers retain the source's signed width. Existing parameter packing replaces
near-zero parameters with positive 1e-10, avoiding the relative-amplitude
division by zero on the supported path. Arithmetic errors are not penalties.

Python uses float64 throughout, without the source's default-REAL casts.
ANNEAL increments KER and supplies it as SEED for each coordinate; historical
RANDO also depends on clocks. Python intentionally uses one reproducible
stream, not clock reseeding. Default-REAL caller arguments/results conflict
with readable REAL*8 RNG definitions unless historical build settings reconcile
them. Compiler flags, linked RNG and runtime behavior remain unresolved.
Neither historical RNG sequence nor executable equivalence is claimed.

Arrays are dynamically allocated. This API explicitly rejects more than 1000
active parameters, matching the local source capacity without emulating its
unguarded overflow. As explicit Python validation policies, controls require
nonnegative integer iteration limit/seed (excluding bool), finite positive
temperatures and `0 < cooling_rate < 1`. Objective values must be finite and
at most the source CHISMIN sentinel 1e20; otherwise no successful result is
returned. These guards restrict unsafe source cases rather than add acceptance
rules. Mapping and trial-evaluator constraints remain authoritative; unresolved
vacancy/free-surface and other unsupported states propagate as exceptions,
without artificial penalties or partial successful results.

No parameter clipping or additional optimizer bounds are introduced. Source
sign folding, dilation handling and near-zero projection are delegated to the
existing mapper. ANNEAL calls no AMOEBA, PARAMIN or SIMPARA; the legacy driver
combines stages externally. Those stages and the driver remain unimplemented.

## Validation

Tests include controlled normal draws, signed widths, literal-source control
flow on independent scalar objectives, ties and failures, pre-cooling counts,
both reset branches, evaluation counts, threshold and MAXITS boundaries,
seed reproducibility and isolation, unsupported-state propagation, immutable
inputs/results and a real small refinement-trial integration. These establish
Python source-derived behavior, not a compiled-Fortran comparison.

Validation with
`C:\Users\peng863\AppData\Local\miniconda3\envs\qcbed-python\python.exe`:
focused tests **22 passed in 6.55 s**; complete suite **1413 passed in 35.97 s**.
`git diff --cached --check` passed with nothing staged. The protected ANNEAL
source remained unchanged during implementation.
