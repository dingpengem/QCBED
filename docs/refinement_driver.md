# Legacy refinement driver

`parse_legacy_refinement(path, *, seed, experimental_header_path=None)` returns
the existing `RefinementState` and immutable `LegacyRefinementConfig`.
`run_legacy_refinement(state, experimental, config, *, phase_grating_method="fft")`
returns immutable `LegacyRefinementResult`, including stage history and the
final physical state. It reuses all four existing optimizers and the existing
trial evaluator. No optimizer algorithm or forward physics was changed.

## Executable evidence and selections

All line numbers below refer to readable `private predecessor implementation`.
Lines 475–528 read controls and vectors; ARRANGEPARS is at 2094–2296.
The guard at 531 selects optimization only for LOOPS>0.

| METHOD first four characters | Stages per loop | Branch / optimizer calls |
| --- | --- | --- |
| ANNE | ANNEAL | 533–581 / 549 |
| SPAR | SIMPARA | 583–634 / 602 |
| PAAN | PARAMIN, ANNEAL | 636–739 / 666,707 |
| ANPA | ANNEAL, PARAMIN | 741–844 / 766,812 |
| SIMP | AMOEBA | 846–929 / 894–895 |
| PARA | PARAMIN | 931–993 / 959 |
| SIPA | AMOEBA, PARAMIN | 996–1133 / 1052–1053,1099 |
| PASI | PARAMIN, AMOEBA | 1136–1275 / 1170,1239–1240 |
| SIAN | AMOEBA, ANNEAL | 1278–1410 / 1329–1330,1376 |
| ANSI | ANNEAL, AMOEBA | 1413–1547 / 1441,1510–1511 |

The source accepts these exact uppercase or lowercase prefixes, including
longer strings with a matching prefix. Arbitrary mixed case is not accepted.
Unknown positive-LOOPS modes are rejected explicitly in Python; the source has
no defined fallback for them. LOOPS<=0 takes the single ONEITPARA evaluation
branch (1548–1560), even for an otherwise unknown METHOD. Python preserves
this behavior for zero and negative integer loop counts.

DO R=1,LOOPS repeats the **whole** selected sequence inclusively. Each stage
runs once per loop; hybrids are not grouped into separate blocks of repeated
first/second stages. No SPAR hybrid exists in the source.

## State transfer and objective bookkeeping

Before every stage, the driver resets CHISMIN=1e20 and calls ONEITPARA(P).
This includes intermediate hybrid boundaries, not just the initial loop.
Examples are 647–648 / 693–694, 752–753 / 793–794, and
1152–1153 / 1197–1198. Python performs exactly one physical trial evaluation
at each boundary. Each optimizer then performs its own documented internal
evaluations. The first boundary evaluation supplies initial_objective; no
additional initial diagnostic call is inserted.

After each stage, UPDPARS(PARMIN), ARRANGEPARS and UPDATE restore the retained
physical state and prepare P for the next stage (e.g. 550–552, 667–669,
897–899, 1055–1057, 1242–1244, 1513–1515).

| Optimizer | Coordinates | Retained handoff |
| --- | --- | --- |
| ANNEAL | physical / ONEITPARA | result.best_vector, best_objective |
| PARAMIN | physical / ONEITPARA | result.best_vector, best_objective |
| SIMPARA | physical / ONEITPARA | result.best_vector, best_objective, including implemented callback alias policy |
| AMOEBA | legacy_simplex / ONEIT | best_evaluated_vector/objective, converted back to raw physical coordinates |

For AMOEBA, the lowest final simplex vertex is insufficient: a better evaluated
reflection may have been discarded by the simplex algorithm. Furthermore its
callback CHISMIN includes the pre-stage evaluation. Python compares the
optimizer's best evaluated score to that boundary score using `<=`; if all
simplex evaluations were worse, retain the pre-stage vector. This is source
callback retention within a stage, not an added global-best archive.

Retained PARMIN and next-stage P remain distinct. Existing scalar UPDPARS
mapping and packing reconstruct P from the retained vector; packing can apply
the source near-zero floor. History and final_vector preserve the raw retained
physical coordinates. Thus a retained zero can appear as +1e-10 in the next
stage's input without rewriting the recorded PARMIN value.

The baseline is never replaced: inactive input values remain fixed and SCAL
always uses the original packed starting values. Existing coordinate-conversion
code was extracted into a shared private helper, with the same arithmetic,
so the driver does not duplicate inverse SCAL or sign-folding equations.
The existing parser similarly exposes already-consumed control values through
a private shared reader; `parse_refinement_state` retains its original behavior.

The source does not maintain an independent all-time best across stages.
Reset scores and retained states can supersede earlier stage results. Python
does not choose the lowest score in history as the final state.

## AMOEBA reconstruction

SCAL is initialized before LOOPS for SIMP (847–850), SIPA (1005–1008),
PASI (1146–1149), SIAN (1282–1285), ANSI (1417–1420). It is not reanchored
after ANNEAL, PARAMIN or a completed loop. Each AMOEBA occurrence reconstructs
a fresh simplex from current PARMIN using the existing `build_legacy_simplex`:
row zero is the affine-transformed retained point; row j+1 adds PLAMDA[j] to
coordinate j. The source performs N+1 initial ONEIT evaluations, which the
existing AMOEBA API includes in its evaluation count. The original Fortran
leading-dimension mismatch remains an executable-equivalence limitation;
Python uses the already validated logical `(N+1,N)` convention.

## Controls

| Config field / source | Use |
| --- | --- |
| method / METHOD | exact dispatch above |
| loops / LOOPS | whole-sequence repetition, nonpositive means one evaluation |
| maxfsimp / MAXFSIMP | ANNEAL max_iterations; AMOEBA and SIMPARA max_evaluations |
| maxfpara / MAXFPARA | PARAMIN max_cycles |
| ftol / FTOL | AMOEBA objective spread; PARAMIN relative-range tolerance |
| chimin / CHIMIN | AMOEBA fmin only |
| maxtemp, mintemp, coolrate | ANNEAL/SIMPARA temperature controls |
| plamda / PLAMDA | AMOEBA simplex displacement; SPAR absolute uncertainties |
| paraopt / PARAOPT | random parameter order for PARAMIN/SIMPARA |
| loud / LOUD | parsed and validated; no legacy unit-file logging implemented |
| seed | explicit Python policy, not a legacy input record |

Each optimizer retains its own inclusive/soft budget rules; the driver neither
normalizes their limits nor truncates a stage. Total result.evaluations counts
boundary evaluations plus all optimizer evaluations. It is not a claim that
the legacy printed NFUNCEVAL is consistently maintained by every optimizer.
Original parameter ranges are reused when PARAMIN starts, as LOW/HIGH are
reloaded from D2/D1 before each source call; adapted ranges are not carried
into another invocation. Temperature and amplitude state likewise restart
according to each optimizer's existing initialization.

For stochastic stages, a deterministic stage seed is derived using
`SeedSequence([config.seed, one_based_loop_index, zero_based_stage_position])`.
History records that seed. It reproduces Python runs without hidden global
RNG state; it does not emulate legacy clock/seed continuity.

## Final state, supported scope and limits

Lines 1564–1569 restore PARMIN with UPDPARS and write the retained parameters,
YNORMBEST and CHISMIN. They do **not** call an objective again. Python returns
the retained objective and mapped final_state without adding a final callback.
For nonpositive LOOPS, the sole evaluation supplies both initial and final
scores. No serialization, legacy output files, CLI or extra optimizer is added.

History entries provide optimizer name, one-based loop index, input physical
vector/score, retained physical vector/score, original immutable optimizer
result and applicable seed. Termination details stay in optimizer_result,
without fabricated uniform convergence labels. Arrays are immutable snapshots.

The driver supports the existing deterministic trial subset. Unsupported
vacancy/free-surface/gas or other mappings fail through the evaluator; no
penalties or partial successful result are returned. Optimization requires
1..1000 active parameters; zero active parameters are permitted only in the
evaluation-only branch. PLAMDA must match active count. Unknown/mixed-case
positive-loop methods, invalid controls, nonfinite/sentinel-exceeding boundary
objectives and optimizer-specific unsupported configurations are rejected.
Finite positive temperatures, 0<COOLRATE<1 and nonnegative optimizer limits
are Python validation restrictions even if a selected stage does not use them.

The adapter relies on deterministic objective reevaluation in the supported
trial subset when reinitializing standalone optimizer callback state. It does
not reproduce arbitrary stateful callbacks, gas randomness, unrelated shared
COMMON mutations or historical output side effects. SIMPARA retains its
approved zero-startup policy; historical RNG ABI/sequences and build provenance
remain unresolved. Source-derived sequencing does not prove numerical
equivalence to the historical executable.

## Validation scope

Tests cover all ten modes and two-loop ordering, exact handoffs, fixed SCAL
and simplex reconstruction, pre-stage callback retention, raw PARMIN versus
packed P, evaluation-only runs, rejection policies, no invented global archive,
immutable baseline/results and a real deterministic hybrid integration.
Both FOR007.DAT and STARTINGFOR007.DAT retain 18 active parameters, select SPAR,
and supply expected controls to a stubbed first invocation. No long full
reference optimization or executable-equivalence comparison is claimed.

The inspected reference controls are SPAR, PARAOPT=1, LOUD=1, LOOPS=1,
MAXFSIMP=1000, MAXFPARA=1, FTOL=CHIMIN=1e-6, MAXTEMP=100, MINTEMP=1,
COOLRATE=0.9. Their PLAMDA begins 0.5, 0.005, 0.005, 0.0002, followed by
fourteen 0.0005 values; all 18 D1/D2 pairs are 1.05/0.95.

Validation interpreter:
`C:\Users\peng863\AppData\Local\miniconda3\envs\qcbed-python\python.exe`.
Focused driver/mapping tests: **61 passed in 1.69 s**. Full suite:
**1511 passed in 45.49 s**. An initial sandboxed focused run encountered a
Windows temporary-directory permission error in an existing parser fixture;
the rerun with required filesystem access passed. `git diff --cached --check`
passed with nothing staged. Protected source/reference files were not edited.
