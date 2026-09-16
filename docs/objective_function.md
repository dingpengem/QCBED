# Legacy objective evaluation

## Executable evidence

The sole readable CHISQCAL definition is in `private predecessor implementation`.
Lines 736-770 compute YNORM; line 741 initializes NORMTOTAL, line 772
initializes CHI2. Lines 773-781 visit every beam (outer loop) and tilt (inner
loop), and line 782 replaces CHI2 with DSQRT(CHI2/NORMTOTAL).
private predecessor implementation is only a related-version declaration reference; no include-only
declaration determines this translation.

With c=THEORINTS, x=XPINTS, v=VARIANCE, w=WEIGHTS and s=YNORM:

```text
term = s*c - x                         private predecessor implementation
NORMTOTAL += ((w*x)*x)/v               private predecessor implementation
CHI2      += ((w*term)*term)/v         private predecessor implementation
objective = sqrt(CHI2/NORMTOTAL)       private predecessor implementation
```

Thus the returned quantity is a normalized weighted root squared error,
not conventional chi-square, reduced chi-square, or unnormalized weighted
SSE, despite the source's CHI2 name and printed “Chi squared” label.
The normalizer is the weighted experimental squared norm, not observation
count, degrees of freedom, beam/tilt count or total weight. There is no
additional normalization, penalty term or parameter-constraint contribution
in this routine or visible callback return path.

The residual is scaled calculated minus experimental. Negative observations
retain their sign in that subtraction; no absolute value or clipping is
applied here. Squaring removes the residual sign in ordinary finite arithmetic,
but the literal operation is `(w*term)*term`, preserving its floating-point
order. Weight multiplies once, before the second residual multiplication;
variance divides once, without a square root. YNORM multiplies calculated
intensity before subtraction and is recomputed using the existing scaling API.
IREF changes only that scale: positive=fixed, zero=weighted least-squares
scale, negative=weighted absolute-count scale. There is no residual branch.

## Callbacks and optimizer paths

ONEIT (`private predecessor implementation`) optionally calls GRADIENT, copies BEAM to THEORINTS
with beam outermost and tilt innermost, then calls CHISQCAL(CSQ). Its best-fit
bookkeeping stores CSQ/YNORM and prints scaled intensities; line 261 returns
ONEIT=CSQ without a visible subsequent adjustment to CSQ.
ONEITPARA (`private predecessor implementation`) performs the same conversion and calls the
same CHISQCAL, returning ONEITPARA=CSQ at line 250. Despite its filename,
private predecessor implementation is this callback, not a second random residual definition.

The readable driver is `private predecessor implementation`, despite its filename:

| Path | Executable callback evidence |
| --- | --- |
| Initial evaluation | 536: STARTINGCHISQ=ONEITPARA(P) |
| AMOEBA | 886,891: YP from ONEIT; 894: CALL AMOEBA(...,ONEIT,...) |
| PARAMIN | 666: CALL PARAMIN(LOW,HIGH,ONEITPARA,FTOL) |
| ANNEAL | 549: CALL ANNEAL(ONEITPARA,MAXFSIMP) |
| SIMPARA | 602: CALL SIMPARA(ONEITPARA,...) |

Later driver branches repeat these callback choices. Readable SIMPARA is in
`private predecessor implementation`: GOF evaluations at 93,231,234,277,395 receive the callback
value. These establish the common objective at the callback interface, not
identical optimizer internals. No readable AMOEBA optimizer definition is
available; private predecessor implementation instead defines ONEIT and CHISQCAL. The supplied private predecessor implementation
and private predecessor implementation are ELF artifacts. Their internal acceptance rules, constraints
and transformations cannot be reconstructed from those files.

GRADIENTOPT=1 invokes GRADIENT before CHISQCAL in both callbacks; it does not
branch inside CHISQCAL. private predecessor implementation is also an ELF artifact, so preprocessing
of shared arrays remains unresolved. This API evaluates the supplied arrays
without claiming that raw simulated/parsed intensities are post-GRADIENT data.

## Array and arithmetic contract

Public `compute_legacy_objective(calculated, experimental, *, iref)` returns
the frozen `LegacyObjectiveResult` with np.float64 fields `objective`,
`scaling`, `weighted_squared_error` (CHI2 before the square root) and
`weighted_experimental_norm` (NORMTOTAL). It reuses both existing YNORM APIs;
their finite float64 input, exact matching shape and real non-bool IREF
validation applies. Input arrays are never modified.

Python `[tilt, beam]` maps to legacy M=beam*n_tilts+tilt+1. An internal helper
produces float64 numerator/normalizer arrays in this shape. Scalar sums use
beam outer, tilt inner loops, preserving sequential addition rather than a
NumPy pairwise reduction. Beam labels, scan positions and GRADIENTOPT do not
exclude any points in this routine.

No readable guard rejects negative variance/weight or masks zero weights.
Zero weights contribute zero for ordinary finite terms and nonzero variance,
but do not bypass division: zero weight with zero variance still evaluates
0/0. A zero normalization denominator is likewise unguarded. Negative
variance/weight can yield negative sums and a negative ratio whose square root
is NaN. Finite input validation is a Python interface policy; nonfinite values
generated by arithmetic are retained, including scaling overflow and 0*Inf.

NumPy IEEE float64 arithmetic is explicitly nontrapping, with warnings locally
suppressed even if the caller selected `np.seterr(all="raise")`. NaN/Inf are
not repaired. Historical compiler trap settings, reassociation and runtime
behavior are unknown. No historical executable equivalence is claimed.

## Validation scope

Tests use hand-calculable signed, nonuniform data, all three IREF branches,
one/multiple beams and tilts, a separate flat-M literal-loop oracle, explicit
pre-scaling, immutable results and unchanged inputs. Cancellation distinguishes
beam-major from row-major summation; overflow distinguishes w*t*t from w*(t*t).
Singular cases cover zero variance, zero weights, zero normalization, negative
ratios and arithmetic-generated NaN/Inf.

All three reference experimental files (in.dat, startingin.dat, bestin.dat)
are evaluated separately with deterministic all-one calculated arrays, for
each IREF branch. They contain 117 zero variances and produce NaN objectives
under this unguarded policy. These are reference-input regressions, not fits
or evidence of the historical executable's output. Their negative-branch
scaling is checked against independent decimal evaluation, 12975.673278983971
with absolute tolerance 1e-9. Literal float64 objective-oracle comparisons are
exact, with NaN compared explicitly. No optimizer or parameter mapping is added.

Validation used exactly
`C:\Users\peng863\AppData\Local\miniconda3\envs\qcbed-python\python.exe`:
48 focused tests passed in 9.70 s; the complete suite passed 1357 tests in
33.39 s. `git diff --check` passed. Protected source/reference files were
unchanged; no files were staged or committed.
