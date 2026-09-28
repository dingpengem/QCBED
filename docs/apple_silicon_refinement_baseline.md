# Apple Silicon refinement reuse baseline

Stage 1 only: instrumentation and exact identity counts, without caching or
changes to production kernels, parameter mapping, objectives, or optimizers.

## Reproduction and scope

From the repository root, activate `qcbed-python` and run:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
python benchmarks/profile_refinement_reuse.py
```

The script reads the current editable full-refinement example's configuration:
depth plus six scan coordinates, existing bounds and rejection policy, ANSI,
seed 2026, one loop, MAXFSIMP=20, masked objective, batched NumPy with batch 64.
It keeps the complete reference mesh, scan, beams and specimen. The configuration,
parameter layout, environment, and SHA256 source/input fingerprints are saved.
Output goes to a unique `tmp/refinement_reuse/<timestamp>/results.json`;
`--output` selects a new directory. No normal production run is written.

The driver runs its complete configured stages with original stopping rules;
the durable runner's outer budgets, previews, checkpointing and independent
verification calls are not part of this optimizer census. Two uninstrumented
initial trials and one instrumented initial trial are recorded separately.
Neither these controls nor reporting create extra optimizer evaluations.

## Instrumentation and interpretation

Only the benchmark imports `TrialProfiler`. With it disabled, no instrumentation
is installed and production has no new overhead. Scoped wrappers invoke the
original functions; temporary line tracing times the inline sequence loop and
the existing batched chunk's kernel construction, propagation and beam extraction.
No mathematical operation is replaced or reordered. Existing trace/debugger
sessions are rejected explicitly. Hooks are restored on exit and exceptions.

Disjoint timing fields are mapping, applied-state setup, exact sequence, base
sample preparation, tilt generation, propagator construction, propagation,
complex beam extraction, intensity calculation, GRADIENT, and objective/YNORM.
`other_seconds` is the remaining trial time, including validation, assembly,
snapshots and wrapper overhead. `forward_inclusive` and
`sequence_setup_inclusive` are nested cross-check timers: do not add them again.
Propagator construction includes chunk setup, allocation, ordering and masking.
Sample preparation times the original preparation function; additional nonunit
reflection assembly, if present in another workload, belongs to the remainder.

All detailed times are instrumented wall times. Tracing overhead is distributed
among measured regions, not confined to the remainder. Uninstrumented controls
show its approximate impact but are not a calibrated correction. Signature
construction and output writing occur outside individual trial timings; the
whole driver wall time includes them. RSS is sampled after each trial, including
the census storage; its maximum is an observed endpoint maximum, not peak RSS.
No multiprocessing or allocation tracer is introduced.

## Exact reuse opportunities

Identity uses full typed keys with IEEE float bytes and array shape/dtype/bytes;
Python sets check key equality, not just a digest. No rounding or tolerance is
used. Signed zero is conservatively distinct. Mapped states, rather than raw
optimizer coordinates, are used. Each count is the number of trials after the
first occurrence of an identity, not the number of repeated pairs.

- Mapped identity includes the complete mapped dataclass; unchanged configuration
  metadata is included conservatively.
- Scan identity includes the complete mapped scan definition (corners, bases,
  dimensions); propagator identity includes actual generated ordered H/K tilts.
- Prepared identity includes mapped sample, applied slice parameters, absorption,
  DSAP, reflection records/symmetry/region mapping, grating method and conventions.
  Region depth and scan are excluded; their physical effects enter other keys.
- Propagator identity includes mesh, stored geometry, applied voltage/dilation,
  tilts, boundary and reciprocal ordering.
- Forward identity combines preparation, propagators, the exact ordered sequence,
  ordered beams, and execution configuration. It excludes experimental values,
  IREF and objective policy because only raw forward results could be reused.

These are unlimited-history opportunity estimates, not implemented cache hit
rates. Rejected trials are counted in total/mapped/scan statistics where a mapped
state exists, but do not enter forward-stage denominators. The report supplies
both repeats per eligible stage invocation and repeats per all trials. It does
not construct sequences for bound-rejected proposals just to increase counts.
No invalid result is treated as a reusable forward calculation.

Prepared and propagator reuse can occur together, while a forward-result hit
already avoids both. Their savings must not be added independently. A finite
cache, validation, lookup/copy costs and verification bypass would reduce actual
benefits. This audit does not authorize Stage 2.

## Measurements: 2026-09-27

Native arm64, macOS 26.3.1; `/opt/miniconda3/envs/qcbed-python/bin/python`,
Python 3.11.16, NumPy 2.4.6. Thread environment variables in the command above
were all 1. This records requested limits, not independent verification of every
library's thread count. No test suite ran concurrently with these measurements.

Raw report: `tmp/refinement_reuse/20260927_160403_606891/results.json`.
Reference dimensions: 16-by-16 mesh, 1521 tilts, nine beams, one stored slice type.
The initial sequence has 415 occurrences; successful refinement trials covered
395 through 507 occurrences. The exact ordered sequence, not this count alone,
determines identity.

Two uninstrumented initial trials took **4.388520 s** and **4.175147 s**.
The instrumented identical trial took **4.711812 s**, about 10.0% above their
mean (7.4–12.9% above the individual controls). This is a small sequential sample,
not a statistical estimate of tracing overhead. Complex amplitudes, raw and
post-GRADIENT intensities were bitwise equal; maximum amplitude and raw-intensity
errors were **0**. Objective and YNORM were exactly equal within this Mac run.
No cross-platform bitwise or compiled-Fortran equivalence is claimed.

| Disjoint component | Initial traced trial (s) | Whole traced refinement (s) |
| --- | ---: | ---: |
| Parameter mapping | 0.000265 | 0.012355 |
| Applied-state setup | 0.000036 | 0.001013 |
| Exact sequence construction | 0.000733 | 0.028596 |
| Sample preparation | 0.011480 | 0.411039 |
| Tilt generation | 0.000059 | 0.001716 |
| Propagator construction and chunk setup | 0.058063 | 2.060746 |
| Wave propagation | 4.633627 | 173.345629 |
| Complex beam extraction | 0.000140 | 0.005550 |
| Intensity calculation | 0.000038 | 0.001436 |
| GRADIENT | 0.000166 | 0.004496 |
| Objective and normalization | 0.006128 | 0.217408 |
| Remaining trial work | 0.001078 | 0.043124 |
| **Total trial time** | **4.711812** | **176.133107** |

The complete driver took **176.523507 s**, including benchmark bookkeeping and
incremental report writes. Signature work summed to 0.084826 s. Successful
trials averaged 5.032182 s (median 4.974270 s). Propagation consumed **98.4%** of
summed trial time; preparation consumed 0.23%, propagator/chunk construction 1.17%.
These fractions apply to this traced workload, not an uninstrumented prediction.

The configured ANSI sequence made **56 evaluations**, matching the profiler's
trial count: **35 successful forward evaluations and 21 bound rejections**.
Initial/final retained objectives were 1.4154022826663109 and 0.7635082330287782.
This describes the existing driver's output, not scientific convergence.

| Exact identity | Eligible trials | Unique | Repeats after first | Potential hit rate among eligible | Repeats / all 56 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Mapped physical state | 56 | 51 | 5 | 8.93% | 8.93% |
| Exact slice sequence | 35 | 21 | 14 | 40.00% | 25.00% |
| Mapped scan definition | 56 | 50 | 6 | 10.71% | 10.71% |
| Prepared state | 35 | 1 | 34 | 97.14% | 60.71% |
| Propagators | 35 | 29 | 6 | 17.14% | 10.71% |
| Complete forward problem | 35 | 30 | 5 | 14.29% | 8.93% |

The preparation hit opportunity is frequent but cheap: repeated preparation
accounts for only **0.398976 s**. Repeated propagator construction accounts for
**0.343900 s**. The five repeated complete forward problems account for
**23.062086 s** of inclusive forward time. Eliminating that last cost for free
would imply an approximately **1.15x** whole-driver ceiling in this instrumented
run, not a measured caching speedup. Those savings overlap and must not be added.
Identical sequences by themselves substantially overstate reusable forward work.

Post-trial RSS observations ranged from **61.94 to 75.50 MiB** during refinement.
They include the process, retained optimizer state and exact-key census, and do
not measure transient peaks or a future result cache. For scale, one retained
amplitude/intensity pair would occupy 328,536 bytes (about 321 KiB) of raw arrays;
one all-tilt complex grid or one-type propagator table would occupy 6,230,016 bytes
(5.94 MiB), excluding metadata/copies. These latter figures are storage estimates.

## Validation and next decision

Focused tests compare instrumented/uninstrumented numerical results, optimizer
results for ANSI/SPAR/PARA, ordered callbacks and evaluation counts; they exercise
mapped sign folding, exact sequence boundaries, scan/thermal invalidation,
bound rejection, full-key identity and hook restoration after failure.

On this Mac, the focused suite passed **7 tests in 0.38 s** and the complete suite
passed **1752 tests in 20.13 s**. `git diff --check` passed. No production Python,
legacy Fortran, or reference-data files were changed for this task.

The measured priority remains exact forward-result reuse, supported by a correct
prepared-state boundary. Prepared-only reuse cannot materially accelerate this
workload. A later bounded cache must retain independent uncached verification,
and its actual speedup must be measured. **Stage 2 has not begun.**
