# Legacy global intensity scaling (YNORM)

> Clean-import provenance update: the [canonical source inventory](legacy_source_inventory.md)
> supersedes source-location, binary-only, duplicate and missing-provider claims below.
> Driver: private predecessor implementation; ONEIT: private predecessor implementation; ONEITPARA: private predecessor implementation; optimizers now
> occupy their matching .f files; RNG: private predecessor implementation; ATOMICSCATFACT:
> private predecessor implementation (all tables confirmed). The four recently corrected core
> file bodies are unchanged. Historical test results below remain historical.


## Source and data flow

`private predecessor implementation`, inside CHISQCAL, contains three YNORM branches.
`private predecessor implementation` reads and prints IREF; its formatted output and the
assignment of positive IREF directly to YNORM establish a real-valued scale
selector, not a beam index. The reference FOR007/STARTINGFOR007 use -1.
No include-only declaration is needed; private predecessor implementation is a related-version
reference, not authority over these executable statements.

ONEIT (`private predecessor implementation`) and ONEITPARA (`private predecessor implementation`) optionally
call GRADIENT, flatten BEAM into THEORINTS with beam outermost and tilt
innermost, then call CHISQCAL. Consequently data-dependent YNORM is recomputed
for each objective evaluation, including initial evaluations. It is not a
one-time calibration. The positive fixed branch assigns the same input value
each time. Saved YNORMBEST (`private predecessor implementation`, `private predecessor implementation`) records a best
evaluation; it is not a separate normalization formula.

## Exact formulas

Let c=THEORINTS (calculated), x=XPINTS (experimental), v=VARIANCE and w=WEIGHTS.
There is **one global scalar**, common to all beams and all scan points:

| Selector | YNORM | Source |
| --- | --- | --- |
| IREF > 0 | IREF | `private predecessor implementation` |
| IREF = 0 | sum(w*c*x/v) / sum(w*c*c/v) | `private predecessor implementation` |
| IREF < 0 | sum(w*abs(x)) / sum(w*abs(c)) | `private predecessor implementation` |

Both sums cover every tilt of every beam. VARIANCE participates only in the
zero branch; WEIGHTS participates directly in both data-dependent branches.
The zero branch is a scalar least-squares projection when w/v defines valid
least-squares weights, as also indicated by its source comment. The code
does not validate that condition. The negative branch is called equal-count
normalization but sums weighted absolute values, not signed counts. The
magnitude of negative IREF does not affect the answer.

The zero branch multiplies `(w*c)*x` and `(w*c)*c` before division by v.
Python retains these operations and sequential accumulation, rather than
precomputing w/v or using a pairwise NumPy reduction. Legacy flat index
`M=b*n_tilts+t+1` maps to Python `[t,b]`; the implementation loops over b
outside t. A cancellation test distinguishes this from row-major accumulation.
There is no per-beam fit. Beam labels, wavelength and scan geometry are not
used directly; geometry affects only the supplied calculated values.

Negative observations retain their sign in the zero branch and can produce
a negative YNORM; there is no clipping. Both calculated and experimental
values enter absolute values in the negative branch. No square root, offset,
background subtraction, threshold or normalization to unit intensity occurs
in YNORM itself. Finite signed weights and variances are not prohibited by
this translation.

## Singular arithmetic: evidence versus execution policy

The readable code performs divisions unconditionally. It does not reject
zero variance, mask zero weights before division, add epsilon, substitute a
factor, or guard a zero denominator. Zero variance is unused in fixed and
equal-count scaling. In the zero branch, even w=0,v=0 evaluates 0/0 and can
contaminate both sums. The reference data contain 117 zero variances.

The Python implementation explicitly uses NumPy IEEE float64 arithmetic:
nonzero/zero produces signed infinity; 0/0 and infinity/infinity produce NaN.
Overflow and subsequent invalid operations propagate. Floating-point warnings
are suppressed locally, independently of caller `np.seterr`; no value is
repaired. Finite source inputs can therefore yield nonfinite scaling/output.
Tests cover these cases and multiplication-before-division overflow.

**Unresolved:** readable Fortran proves the absence of guards, but not the
historical compiler's trap settings, reassociation or exact evaluation order.
An executable could trap instead of returning IEEE nonfinite values. Python
reproduces the unguarded arithmetic under an explicit nontrapping float64
policy; it does not claim the historical runtime's singular-case behavior has
been observed or establish executable equivalence.

## Public API and mutation

```python
factor = compute_legacy_intensity_scaling(calculated, experimental, iref=-1.0)
scaled = apply_legacy_intensity_scaling(calculated, factor)
```

`experimental` is the existing ExperimentalIntensityData. `calculated` must
be finite float64 with exactly the same `(n_tilts,n_beams)` shape. No beam or
tilt alignment is guessed from labels. IREF must be a finite real scalar,
excluding bool. All three experimental fields retain their existing finite
data contract, although the selected formula may not use variance/weights.

Compute returns an immutable `np.float64` scalar. Apply returns a new,
bytes-backed immutable float64 `(tilt,beam)` array with every value multiplied
by that scalar. Apply accepts NaN/Inf factors from the compute step and
propagates them; finite calculated inputs and a real non-bool scalar are
otherwise required. No source data, experimental arrays or raw simulated
array is modified. The functions introduce no global state or cached factor.

This matches the source's use of YNORM at output
(`private predecessor implementation`, `private predecessor implementation`, `private predecessor implementation` and other output paths).
The source does not overwrite THEORINTS or XPINTS in the scaling branches.
Its later residual term is `YNORM*THEORINTS(M)-XPINTS(M)` at `private predecessor implementation`;
lines 778-783 then aggregate with WEIGHTS/VARIANCE and take a square root.
Those residual/objective calculations are recorded for the next stage and
are **not implemented** here. The square root is not part of YNORM.

GRADIENTOPT does not select a YNORM branch; IREF does. However, GRADIENT is
called before CHISQCAL when GRADIENTOPT=1. The supplied private predecessor implementation is an ELF
artifact, not readable source, so its transformation and effects on shared
arrays remain unresolved. These APIs scale the arrays supplied by the caller;
they do not pretend that raw simulation arrays are equivalent to the
historical post-GRADIENT arrays.

## Validation

Hand-calculable cases exercise all branches, multiple tilts/beams, negative
and zero values, nonuniform/zero weights, zero variance and denominators,
shape/dtype checks, input/output immutability and scaling recomputation.
A separate flat-M loop oracle agrees exactly with the beam/tilt implementation.
No vectorized reduction is introduced.

Each of in.dat, startingin.dat and bestin.dat is tested separately using
synthetic all-one calculated data, isolating scaling from unresolved forward
and GRADIENT equivalence. Independent Decimal evaluation of the stored tokens
gives weighted absolute sum 97953357.58305 and weight sum 7549, hence factor
12975.67327898397138693866737. Sequential float64 accumulation yields
12975.673278983528 (about 4.44e-10 absolute difference); the hard-coded decimal
comparison permits 1e-9, while the direct float64 loop oracle matches exactly.
Doubling calculated values halves the two fitted factors; fixed IREF is
unchanged. On these files, zero-branch scaling is NaN because of unguarded
zero variances. This is a domain observation, not an objective result or an
assertion that the historical executable produced NaN.

Validation used exactly
`C:\Users\peng863\AppData\Local\miniconda3\envs\qcbed-python\python.exe`:
**36 focused tests passed in 1.31 s**; **1309 complete-suite tests passed in
24.57 s**. `git diff --check` passed. No protected source/reference files were
modified, and no files were staged or committed.

No fitting pipeline, chi-square, residual aggregation, optimization, new
parser/refinement fields, or forward-physics changes are included.

## Corrected-source reconciliation

The restored `private predecessor implementation` and `private predecessor implementation` confirm that only calculated BEAM
values undergo optional GRADIENT, before YNORM/objective evaluation. Input
experimental intensity, variance and weight are unchanged. Legacy divisions
retain their IEEE zero-variance behavior; masked mode remains a separate
non-legacy policy. ONEIT takes simplex coordinates and applies inverse SCAL;
ONEITPARA takes physical coordinates. See the
[four-source reconciliation](four_source_reconciliation.md) for selective
UPDPARS best-state writeback, source lines and validation results.
