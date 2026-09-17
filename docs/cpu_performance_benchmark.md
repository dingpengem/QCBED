# CPU performance benchmark: corrected reference case

Benchmark-only work: no production backend, loop, physics, objective or optimizer changes.
No package installations, persistent pool, caching, shared memory or acceleration was added.

## Environment and method

- Interpreter: `C:\Users\peng863\AppData\Local\miniconda3\envs\qcbed-python\python.exe`; Python `3.11.16 | packaged by Anaconda, Inc. | (main, Aug 27 2026, 14:36:16) [MSC v.1942 64 bit (AMD64)]`.
- CPU: 13th Gen Intel(R) Core(TM) i7-1365U; 10 physical / 12 logical CPUs.
- NumPy 2.4.6; SciPy 1.17.1; Windows-10-10.0.26200-SP0.
- Input: FOR007.DAT + in.dat; depth 840.0 angstrom; mesh [16, 16]; 1521 tilts; 9 beams; 415 occurrences; 1 stored slice type.
- Fixed physical values; only refinement flags are reduced to one depth parameter.
- The scoped-timer baseline calls evaluate_refinement_trial with physical coordinates,
  phase_grating_method="fft" and objective_mode="masked". Preparation and actual final
  gratings are captured and checked equal; every sweep uses that same prepared sample.
- Forward tables exclude objective/GRADIENT/report generation and preparation. Sampling
  order, all 1521 tilts and nine measured beams are preserved. No optimizer is run.
- Sequential wall-clock measurements, no deliberately concurrent benchmarks. Single
  sweep samples are followed by two winner repeats; thermals/OS scheduling are not controlled.
- 100-ms RSS sampling is approximate. Process RSS sums double-count shared pages and
  may miss short allocation peaks. CPU affinity is not pinned.
- Raw evidence: tmp/cpu_performance/results.json, baseline.npz, batched.prof and
  batched_profile.txt. Reproduce with the exact interpreter and benchmarks/cpu_performance.py.
  --resume reuses completed stages in that directory; use only for the same code/input/environment.

## Threads and backend evidence

Inherited environment: `{'OMP_NUM_THREADS': None, 'MKL_NUM_THREADS': None, 'OPENBLAS_NUM_THREADS': None, 'NUMEXPR_NUM_THREADS': None}`.
The benchmark sets OMP_NUM_THREADS, MKL_NUM_THREADS, OPENBLAS_NUM_THREADS and
NUMEXPR_NUM_THREADS to 1 before NumPy/SciPy import, only in its own process/children.
SciPy FFT workers=2/4 are explicit single-process comparisons. Process-worker tests
use NumPy single-thread PocketFFT; FFT and process parallelism are not combined.
NumPy loads _pocketfft_umath; SciPy loads _pocketfft/pypocketfft. BLAS/LAPACK
are MKL-backed: conda libblas/libcblas/liblapack builds end in mkl; loaded DLLs
include mkl_rt.3, mkl_core.3, mkl_intel_thread.3 and libiomp5md. This does NOT
make np.fft/scipy.fft MKL FFTs. mkl_fft and pyFFTW are absent; neither was installed.
Full config text, package builds and loaded library paths are in the JSON.

## Corrected numerical baseline and objective cost

The numerical baseline is the current production scalar trial result, not the historical
Fortran runtime claim. Previous direct-DFT and firstout correctness results are not
newly measured executable equivalence here. All listed residuals use complex128
amplitudes and raw calculated intensities against this corrected production baseline.

| Scoped production trial component | Seconds |
|---|---:|
| mapping | 0.000838 |
| sequence_and_applied_state | 0.000432 |
| sample_preparation | 0.019481 |
| scan_generation | 0.000205 |
| propagator_construction | 0.212813 |
| tilt_multislice | 192.396803 |
| gradient | 0.000445 |
| ynorm_objective | 0.011039 |
| Full scoped trial wall | 192.763361 |
| NumPy FFT calls (nested inside above) | 56.326412 |

Masked YNORM=4723684.37705874; objective=1.41540228266631. FFT time is inclusive
inside multislice/preparation and must not be added again to the component totals.
Minimally instrumented full serial forward: 189.755449 s.
It excludes the very small preparation/objective stages and scoped per-FFT timer
overhead. Different timing conditions are reported, not silently combined.

## Full-reference batch sweep

| Batch | Wall s | ms/tilt | Approx peak RSS MiB | Increment MiB | Max amplitude error | Max intensity error | Bitwise equal |
|---:|---:|---:|---:|---:|---:|---:|---|
| 1 | 76.672710 | 50.4094 | 62.74 | 0.04 | 0 | 0 | True |
| 4 | 35.483939 | 23.3293 | 63.29 | 0.22 | 0 | 0 | True |
| 8 | 23.397189 | 15.3828 | 63.29 | 0.02 | 0 | 0 | True |
| 16 | 18.834325 | 12.3829 | 63.39 | 0.02 | 0 | 0 | True |
| 32 | 14.805142 | 9.7338 | 63.41 | 0.05 | 0 | 0 | True |
| 64 | 14.758088 | 9.7029 | 63.93 | 1.84 | 0 | 0 | True |
| 128 | 15.848954 | 10.4201 | 65.38 | 3.12 | 0 | 0 | True |
| 256 | 20.245718 | 13.3108 | 67.77 | 5.59 | 0 | 0 | True |
| 512 | 25.466294 | 16.7431 | 73.32 | 10.08 | 0 | 0 | True |

Best observed sweep batch: **64**. Winner repeat wall times: 14.439718 s, 15.362765 s.
This is the best measured configuration in this sweep, not a universal hardware optimum.
RSS increments are affected by allocator reuse across sequential stages; absolute RSS
includes Python/SciPy import state. No tested batch exceeded available memory.

## Full propagation with scoped FFT backend substitutions

All configurations use batch 64, identical inputs and existing production _chunk.
Only np.fft.fft2/ifft2 callables are temporarily substituted; axes and norm="forward"
are retained. Patches are restored after each run. All four paths include FFT timers.
NumPy FFT totals exclude validation/multiplication and include both transform axes.

| Backend | FFT workers | Full wall s | Forward FFT s | Inverse FFT s | Max amplitude error | Max intensity error | Relative intensity L2 |
|---|---:|---:|---:|---:|---:|---:|---:|
| numpy | 1 | 15.614322 | 5.091969 | 4.775731 | 0 | 0 | 0 |
| scipy | 1 | 10.109318 | 3.203224 | 3.207651 | 4.11e-15 | 5.163e-15 | 3.564e-15 |
| scipy | 2 | 17.041834 | 5.493787 | 5.459891 | 4.11e-15 | 5.163e-15 | 3.564e-15 |
| scipy | 4 | 10.749409 | 3.573232 | 3.252497 | 4.11e-15 | 5.163e-15 | 3.564e-15 |

Microbenchmarks: random complex128 C-order production-sized batches; three timed
repeats, median per-call time, after warmup. Both directions use norm="forward".
All measured calls return separate output storage; no in-place optimization is used.

| Backend/workers | Direction | Shape | us/call | calls/s | FFT-only max error |
|---|---|---|---:|---:|---:|
| numpy/1 | forward | [64, 16, 16] | 294.747 | 3392.7 | 0 |
| numpy/1 | inverse | [64, 16, 16] | 266.727 | 3749.2 | 0 |
| scipy/1 | forward | [64, 16, 16] | 197.602 | 5060.7 | 7.38e-17 |
| scipy/1 | inverse | [64, 16, 16] | 220.294 | 4539.4 | 2.14e-14 |
| scipy/2 | forward | [64, 16, 16] | 388.586 | 2573.4 | 7.38e-17 |
| scipy/2 | inverse | [64, 16, 16] | 376.830 | 2653.7 | 2.14e-14 |
| scipy/4 | forward | [64, 16, 16] | 309.246 | 3233.7 | 7.38e-17 |
| scipy/4 | inverse | [64, 16, 16] | 394.048 | 2537.8 | 2.14e-14 |

## Slice-loop timing

The following timers trace actual _chunk source lines; the loop itself is unchanged.
Dispatch includes validation/context management plus the separately timed FFT. These
are nested values, not disjoint buckets. Loop-control timing includes tracing overhead.

| Existing production operation | Traced seconds |
|---|---:|
| kernel_setup_and_initialization | 0.207226 |
| loop_control | 0.032235 |
| forward_dispatch_inclusive | 5.494929 |
| propagator_multiplication | 1.936537 |
| post_product_validation | 0.442113 |
| inverse_dispatch_inclusive | 5.230075 |
| phase_grating_multiplication | 0.639765 |
| beam_extraction | 0.000674 |
| Unaccounted outside traced lines | 0.173522 |
| Full traced wall | 14.157076 |
| Raw FFT (nested) | 8.704090 |

Tracing and cProfile change runtime. Use the untraced sweep/repeat totals for speed
comparisons. cProfile full-batch wall=12.870084 s; allocation-traced
wall=14.104726 s. Their overhead is not production cost.

## Allocation/copy audit

| Operation | Classification / evidence |
|---|---|
| initial np.broadcast_to(grating, shape) | zero-stride view; no wave data allocation |
| later wave * grating | out-of-place complex128 wave allocation each occurrence |
| np.fft.fft2/ifft2 | out-of-place; two 1-D axis stages in NumPy; backend temporaries implementation-dependent |
| reciprocal *= kernels[index] | in-place; kernel is not modified |
| _batch_array: np.isfinite + np.all | finite-check boolean temporary and reduction, input array returned unchanged |
| kernel batch np.empty + per-tilt kernel calculation | allocated once per batch/type, reused through all occurrences |
| centered kernel[permutation] | advanced-index copy once per tilt at setup |
| reciprocal[:, beam_rows, beam_cols].copy() | advanced-index beam-sized allocation followed by explicit copy, once per tilt batch |
| TiltSeriesResult _snapshot | immutable bytes copies for tilts, amplitudes, intensities, outside slice loop |
| scalar _product(...).copy(order="C") | out-of-place product then an extra copy; three finite checks |
| scalar centered_to_unshifted/unshifted_to_centered | fresh frequency grids/index arrays and copied/reordered wave every transform |
| scalar _centered_bins astype(np.intp) | new index arrays every call; no reuse in scalar path |
| batched FFT result layout | no explicit astype/ascontiguousarray/copy in FFT wrappers; backend output allocation remains |

tracemalloc full-run peak: 1.468 MiB; retained: 0.338 MiB.
This is live traced allocation, not cumulative traffic or a bandwidth measurement.
Per 16x16 complex128 plane is 4096 bytes. Across 1521 tilts x 415 occurrences:
phase multiplication allocates 629,694 planes (2.402 GiB cumulative payload).
NumPy two-axis forward/return FFT outputs account for approximately 9.620 GiB
additional cumulative payload. Five finite masks per interior occurrence add about
0.752 GiB. These operation-count estimates exclude backend workspace, allocator
metadata, kernel construction, copies, and cache reuse. They are NOT measured DRAM
traffic. Live batched storage is small and repeatedly reused by the allocator/cache;
large cumulative allocation does not by itself establish a memory-bandwidth bottleneck.

## Current process parallelism

The production API constructs and shuts down a new spawn ProcessPoolExecutor per
simulation. A benchmark-only executor subclass records submission/result/shutdown
times; top-level wrappers time the unchanged initializer and worker chunk. Worker
file-write markers add small diagnostic overhead. No persistent pool is introduced.
Process chunk size equals selected batch size for these comparisons.

| Mode | Full wall s | FFT batch | Workers | Approx summed peak RSS MiB | Max amplitude / intensity error |
|---|---:|---:|---:|---:|---|
| batched | 14.758088 | 64 | 1 | 63.93 | 0 / 0 |
| process | 13.978934 | 64 | 2 | 151.60 | 0 / 0 |
| process | 7.833471 | 64 | 4 | 239.06 | 0 / 0 |

| Process timing component | 2 workers s | 4 workers s |
|---|---:|---:|
| executor_construction_seconds | 0.003695 | 0.001694 |
| submission_seconds | 8.273623 | 4.945560 |
| future_result_seconds | 0.000178 | 0.000176 |
| shutdown_seconds | 0.111748 | 0.091590 |
| start_to_first_initializer_seconds | 7.473275 | 1.714714 |
| initializer_body_sum_seconds | 0.001308 | 0.003438 |
| worker_computation_sum_seconds | 11.577494 | 15.685334 |
| worker_compute_span_seconds | 6.398873 | 6.040341 |
| worker_cpu_sum_seconds | 11.375000 | 15.281250 |

Prepared initializer pickle: 19378 bytes; isolated mean dump 0.00010228 s, load 8.498e-05 s (100 repetitions).
Start-to-first-initializer includes process launch, imports, argument transfer and
scheduling; initializer-body time excludes deserialization. These were not separately
identifiable without changing the production executor. Submission includes pool
startup work and queues asynchronous serialization, so it is not pure task transfer.
Worker computation sums overlap across processes; do not add them to wall time.
Future.result timing is post-wait retrieval; deserialization/collection also occurs in
executor management threads. Compute span and total wall bound the noncompute tail/
startup, not a clean disjoint IPC decomposition.

## Supplementary scalar and process checks

The 32-real-tilt scalar prefix took 4.0947 s unprofiled and 6.7519 s with cProfile
and scoped product/FFT timers. This is **not** a full-run timing extrapolation.
The profiled prefix made 4,819,011 calls; the complete batched profile made
1,425,464 calls. Scalar cumulative timings include _array 1.6126 s (159,328 calls),
_centered_bins 1.3505 s (26,528), centered-to-unshifted 1.1244 s and the reverse
1.2242 s. These overlap and must not be added. Product wrappers (validation and
copies included) took 0.5913 s for phase gratings and 0.6698 s for propagators;
raw FFT calls took 1.8803 s. This locates repeated validation/index construction
and array API overhead around the scalar FFTs, rather than an expensive objective.
See tmp/cpu_performance/scalar_prefix_profile.txt for the complete profile.

Both NumPy and SciPy forward outputs at (64,16,16) are C-contiguous with byte
strides (4096,256,16), matching the test input. Thus a backend-layout difference
was not demonstrated in this probe. No layout conversion should be proposed on
that assumption. MKL_Get_Max_Threads returned **1** in the benchmark process.
Microbenchmark FFT times are warm, standalone calls and are lower than some full
interleaved-run per-call times; use measured full totals instead of extrapolating
microbenchmark rates to claim an end-to-end runtime.

| Minimally instrumented repeat | Wall s | Approx summed peak RSS MiB | Max amplitude / intensity error |
|---|---:|---:|---|
| 2 process workers | 7.194758 | 149.67 | 0 / 0 |
| 4 process workers | 7.081543 | 237.39 | 0 / 0 |

The first two-worker launch was substantially slower than its repeat. The cause
of launch variability was not isolated; avoid attributing the whole difference to
profiling instrumentation or to worker scaling. Each run creates a fresh pool.
All NumPy sweep, plain serial, worker, line-traced, allocation-traced, cProfile,
repeat and scalar-prefix results were bitwise equal to their matching baseline.
SciPy full-run differences were identical across worker settings and at ordinary
floating-point FFT roundoff scale; no unexplained numerical discrepancy occurred.

## Conclusions and ranked next steps

The current trial API still runs the scalar engine. The isolated existing batched
engine is about **12.9x faster** (189.76 / 14.76) on this case with identical
amplitudes. Batch 32 and 64 differ by only 0.3% in the sweep; the batch-64 repeats
span 14.44-15.36 s, so treat 32-64 as a useful plateau, not a precise optimum.
Batch 1 already removes much scalar API/reordering overhead (76.67 s); combining
independent tilts reduces remaining call overhead further without skipping slices.

The newly measured scalar times (189.76 s forward, 192.76 s scoped objective)
are slower than the earlier 166-174 s report. FFT time is still similar: 56.33 s
versus about 55 s. Hardware power/thermal state and prior-run instrumentation are
not controlled enough to assign the change to one cause. Preparation remains
0.0195 s, propagators 0.2128 s, GRADIENT 0.00045 s and masked objective 0.0110 s.
Those stages cannot explain the gap. Even eliminating every scalar FFT for free
would leave approximately **136.44 s** in the scoped trial, above the nominal
60-second target. FFT replacement alone is therefore insufficient.

The batch-64 traced run partitions into raw FFT 8.7041 s; phase multiplication
0.6398 s; propagator multiplication 1.9365 s; FFT-wrapper work excluding raw FFT
2.0209 s; post-product validation 0.4421 s; loop/control 0.0322 s; setup 0.2072 s;
and unaccounted time 0.1735 s. These sum to 14.1571 s without double-counting.
The tiny loop-control bucket does not justify starting with a compiled loop.
Native operations and their validation/dispatch matter much more than the bare for loop.

The nominal 60-second target is met by the isolated batched engine and process
runs. The historical Fortran number is user-reported, on unspecified build/hardware
and iteration semantics; it is not a measured comparator and does not establish
that Python is faster than the historical executable.

| Rank / next step | Measured evidence and speedup estimate | Complexity / numerical risk | Scope |
|---|---|---|---|
| 1. Let trial evaluation select the existing batched engine, initially batch 32-64 | Scalar 189.76 s vs batch 14.44-15.36 s: about 12-13x measured forward benefit. Preparation/objective costs are tiny. | Moderate integration; low risk if trial-specific reflection-scaled gratings, ordering, unsupported-case rejection and objective flow are preserved. No optimizer/physics redesign. | General supported refinement parameters |
| 2. Offer SciPy FFT workers=1 through the backend interface | Full timed NumPy 15.61 s vs SciPy 10.11 s: 1.54x in this comparison. Microbenchmarks also favor SciPy/1. Workers 2/4 did not improve it. | Low/moderate; roundoff changes of 4.11e-15 require numerical and callback/threshold regression checks. Keep backend choice explicit. | General |
| 3. Expose the existing process engine to trial execution, prefer testing 2 workers first | Plain repeated NumPy process totals 7.19 s (2), 7.08 s (4), about 2x vs batched single process; four uses much more RSS with little repeated-run gain. | Moderate integration; preserve immutable per-trial sample and deterministic ordering; startup variance must be visible. | General; most useful when enough tilt work exists |
| 4. Reduce internal temporary/copy/validation overhead with tested buffers | Batched wrappers/validation and multiplications consume several seconds; peak live traced allocations only 1.47 MiB. No isolated buffer optimization was measured. Eliminating all non-FFT work gives an unrealistic broad ceiling of 1.63x for the traced run; actual buffer-only benefit is smaller and unknown. | Moderate/high; aliasing, read-only samples, overflow checks and complex128 must remain correct. Do not simply disable scientific validation. | General |
| 5. Consider persistent pools for repeated evaluations after integration | Fresh-pool first-initializer latency 1.71-7.47 s; pickle only 19,378 bytes and ~0.1 ms. Wall/compute-span ratios give hypothetical overhead-free ceilings 1.30x (four-worker run), 2.18x (cold two-worker run); not measured persistent-pool gains, and warm gains may be much smaller. | Moderate/high lifecycle and trial-state update complexity; stale prepared samples are a numerical risk. | General repeated trials, provided sample updates are explicit |

These speedups are separate comparisons; they must not be multiplied into a claim
about an unbenchmarked combined backend/process implementation. **Implement step 1 first.**
No production optimization was implemented here.

Thickness checkpoints remain a separate, potentially valuable narrow optimization:
415 stored occurrences versus one extra occurrence suggests roughly 99.8% fewer
propagation steps for a valid one-slice extension, but no checkpoint runtime speedup
was measured. Validity depends on unchanged gratings, propagators and scan, not just
the name of the refined parameter. They do not address general parameter changes.
A pure compiled loop is lower priority: the measured bare loop/control fraction is
about 0.23%; reducing the broader function/array-operation cost requires more than
translating that for loop. Memory bandwidth itself was not measured or proven to be
the bottleneck; high cumulative small-array allocation and repeated API work are
confirmed, while cache effects/allocator/native computation are not fully separable.


## Validation and working tree

Every benchmark records amplitude/intensity residuals. Backend runs assert
allclose(rtol=2e-12, atol=2e-12) on amplitudes; hard equality is reported separately.
No unexplained numerical differences are accepted as a speed improvement.
Final pytest, diff-check and status are recorded below.

Complete suite with the exact interpreter: **1585 passed in 81.23 s**.
Benchmark residual/shape/dtype sanity checks passed for all 24 recorded full/prefix
comparison configurations (plus independent FFT microbenchmark checks).
`git diff --check`: passed. New untracked files were separately checked for
trailing whitespace. `git diff --exit-code -- src/qcbed private predecessor sources reference_data`:
passed with no changes. Nothing staged or committed; bestin.dm4 was ignored.

Reproduction:

```powershell
& "C:\Users\peng863\AppData\Local\miniconda3\envs\qcbed-python\python.exe" benchmarks/cpu_performance.py
```

Final `git status --short`:

```text
?? benchmarks/cpu_performance.py
?? docs/cpu_performance_benchmark.md
?? examples/validate_fft_correction.py
?? reference_data/bestin.dm4
```

The latter two paths were already untracked before this task and were not modified.

## Integrated refinement trial backend

The follow-up integration adds `ForwardExecutionConfig` and
`LegacyRefinementConfig.execution`. Serial remains the default reference;
explicit batched execution calls the existing `simulate_tilt_series_batched`
with a default maximum batch size of 64. This replaces only the tilt execution
inside a trial, not any optimizer or objective algorithm. Parallel refinement
integration and SciPy FFT remain deferred.

Source/control-flow checks used `private predecessor implementationâ€“101` (SFMAP2 and phase-grating
preparation), `111â€“127` (per-tilt propagators), `129â€“163` (ordered slice propagation
and beam extraction), and `164â€“176` (GRADIENT, beam-major assembly, CHISQCAL).
The existing Python coordinate, FFT, sequence and objective conventions remain
unchanged. The trial's final reflection-adjusted gratings are placed in immutable
prepared-sample snapshots for batching; base gratings are never substituted for
nonunit elastic/absorptive multipliers. There is one preparation implementation.

All pre-stage, optimizer callback, stage handoff, local-runner final diagnostic,
and validation reevaluation calls receive the same execution config. The four
optimizer modules changed only to forward that config. Tests compare exact
evaluation histories and all returned optimizer state for ANNEAL, AMOEBA,
PARAMIN, SIMPARA, ANSI and PAAN with fixed seeds/order and both objective modes.
Additional trial coverage includes optional GRADIENT, three IREF branches,
nonbinary weights and zero-variance masked points, nonzero scan coordinates,
partial batches, multiple stored slices, absorption, nonunit reflection
multipliers and direct/FFT phase-grating construction.

The reproducible follow-up benchmark is `benchmarks/refinement_execution.py`:

```powershell
& "C:\Users\peng863\AppData\Local\miniconda3\envs\qcbed-python\python.exe" benchmarks/refinement_execution.py
```

It times complete production `evaluate_refinement_trial` calls for serial and
batch sizes 32/64 on the physical reference case (840 angstrom, 1521 tilts,
nine beams, 415 occurrences), with only the active flags reduced to depth as in
the previous benchmark. No optimizer or shortened scan/sequence is used.
Each call independently prepares its sample. Scoped forward timers include
propagator construction and propagation; preparation/assembly are in the
remaining total. GRADIENT and combined YNORM/objective have separate timers.
Native thread environment limits are one, and timings are single observations,
not a statistical performance guarantee.

After all timed calls, the benchmark uses `CbedComparison`, `_simulated` and
the existing `plot_cbed_comparison` to produce separate serial/batched galleries
for all nine beams. It verifies the actual Matplotlib image data, weight masks
and color limits agree exactly. Experimental and YNORM-scaled simulated maps
share color limits; residuals use symmetric limits about zero. No independent
display normalization is applied. The existing four-panel layout repeats the
evaluated trial in its two simulation slots, explicitly relabeled; this is a
diagnostic trial preview, not evidence of refinement improvement. Plotting,
compression and report generation are outside every measured trial scope.

### Measured integrated trial results (2026-09-16)

Exact interpreter: `C:\Users\peng863\AppData\Local\miniconda3\envs\qcbed-python\python.exe`. NumPy 2.4.6.
The reference uses GRADIENT=1 and the existing mask-before-YNORM policy.

| Backend | Complete trial (s) | Forward (s) | GRADIENT (s) | YNORM/objective (s) | Other/preparation (s) | Speedup |
|---|---:|---:|---:|---:|---:|---:|
| serial | 333.206213 | 333.139277 | 0.000415 | 0.006682 | 0.059838 | 1.00x |
| batched32 | 7.990145 | 7.969735 | 0.000231 | 0.007684 | 0.012494 | 41.70x |
| batched64 | 7.309462 | 7.294363 | 0.000235 | 0.005692 | 0.009171 | 45.59x |

All three calls returned objective **1.4154022826663093** and YNORM
**4723684.377058742**. Max complex-amplitude, raw-intensity, post-GRADIENT,
YNORM and objective residuals were **exactly zero**. Mapped physical states
also matched exactly. No objective or optimizer policy was modified.

This invocation differs substantially from the earlier 192.76/14.76-second
measurements: serial was slower and batching faster. The 41.70x/45.59x ratios
are observations from this invocation, not a universal or statistically established
speedup. CPU frequency, thermal/power state and background load were not controlled;
no single cause for the difference is established. Batch 64 was about 8.5% faster
in wall time than batch 32 here. The configurable default remains 64, not an
automatic tuner or a promise that 64 wins on every mesh/machine.

Artifacts (ignored benchmark output; not committed):

- Gallery index: `tmp/refinement_execution/20260916_233916_909590/index.html`.
- Serial preview: `tmp/refinement_execution/20260916_233916_909590/serial.html`.
- Batched preview: `tmp/refinement_execution/20260916_233916_909590/batched64.html`.
- Each preview contains all nine beams as separate labeled PNGs; the HTML and PNGs
  should be kept together. Results JSON and three NPZ numerical snapshots are in
  the same timestamped directory.
- Actual plotted data, masks and color limits were checked with exact equality
  for every panel and beam. The generated origin-beam PNG was visually inspected.

Production files touched: new `execution.py`; `refinement_trial.py`,
`refinement_driver.py`, the four optimizer callback wrappers (`anneal.py`,
`amoeba.py`, `paramin.py`, `simpara.py`), `local_runner.py`, `validation.py`,
`validation_report.py`, and package exports. New focused tests:
`tests/test_refinement_execution.py`. Physics kernels, objective implementations,
GRADIENT, scan generation, parameter mapping and optimizer algorithms were not edited.

### Integration validation

- Focused refinement/backend/driver/local-runner suite: **137 passed in 22.80 s**.
- Complete suite: **1648 passed in 62.80 s**, using the exact qcbed-python interpreter.
- Serial/batched optimizer decisions: all returned optimizer fields, ordered
  objective-evaluation histories, retained vectors and stage handoffs matched
  exactly in the deterministic/seeded short runs. No tolerance-based branching
  differences were accepted.
- `git diff --check`: passed. New untracked implementation/test/benchmark files
  separately passed trailing-whitespace checks.
- `git diff --exit-code -- private predecessor sources reference_data`: passed; protected
  tracked files unchanged. The two user-owned untracked files were not modified.
- No stage or commit performed. The measured reference run used masked mode;
  finite legacy-objective equivalence is tested with synthetic data. Historical
  executable equivalence and a full real-data optimization remain outside this
  performance integration validation.

Final working tree for this task:

```text
 M README.md
 M docs/cpu_performance_benchmark.md
 M docs/local_optimization.md
 M src/qcbed/__init__.py
 M src/qcbed/amoeba.py
 M src/qcbed/anneal.py
 M src/qcbed/local_runner.py
 M src/qcbed/paramin.py
 M src/qcbed/refinement_driver.py
 M src/qcbed/refinement_trial.py
 M src/qcbed/simpara.py
 M src/qcbed/validation.py
 M src/qcbed/validation_report.py
?? benchmarks/refinement_execution.py
?? examples/validate_fft_correction.py
?? reference_data/bestin.dm4
?? src/qcbed/execution.py
?? tests/test_refinement_execution.py
```
