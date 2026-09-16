# CPU process execution

## API and unchanged physics

```python
simulate_tilt_series_parallel(
    prepared_sample, tilts,
    max_workers=2, process_chunk_size=32, batch_size=16,
    start_method="spawn",
)
```

The hierarchy is explicit tilts -> contiguous process tasks -> bounded FFT
batches inside each task. Every task calls the existing validated
`simulate_tilt_series_batched`; no propagation equations are duplicated here.
Slices remain sequential. Only independent plane-wave experiments run in
different processes. There is no coherent sum, cross-tilt FFT, convergent
probe, extra normalization or change to the source-backed tilt convention.
The independence audit is in [tilt series](tilt_series.md), based on readable
ONEIT/ONEITPARA loops, not the related-version private predecessor implementation declarations.
Historical CDFT/executable equivalence remains unresolved (category B).

`tilts` retains the finite float64 `(n_tilts,2)` contract: undilated fractional
reciprocal-cell H then K, not angles. `max_workers`, `process_chunk_size` and
`batch_size` are positive integers excluding bool. The process chunk size is
the number of tilts submitted in one task; the internal batch size limits
simultaneous spatial FFTs, independently of the process chunk size. Final
incomplete chunks, duplicates, noncontiguous input arrays and chunks/batches
larger than the input are supported. Empty tilts are rejected.

Only `start_method="spawn"` is supported. Even `max_workers=1` creates a
spawned worker; there is no automatic batched/serial fallback. For direct
one-process execution call `simulate_tilt_series_batched` explicitly.

The returned `TiltSeriesResult` retains immutable float64 tilts and
intensities, complex128 amplitudes, requested beam H/K and duplicate/order
preservation. Shapes remain `(n_tilts,n_beams)` for beam data. Its
`execution_mode` is `parallel`; `batch_size` is the requested internal batch.
Optional frozen `parallel: ParallelExecutionMetadata` contains:

- `requested_workers`: requested max_workers;
- `worker_limit`: configured executor limit, capped by the number of tasks
  and, on Windows, ProcessPoolExecutor's 61-worker wait-handle limit;
- `actual_workers`: distinct processes that completed at least one task;
- `process_chunk_size` and `start_method`.

Participation can be less than the configured limit when early workers finish
tiny tasks before others start. Worker PIDs are used only to count distinct
participants; neither IDs nor completion order are exposed in the result.
Serial/batched results retain their existing behavior with `parallel=None`.

## Windows spawn and data ownership

Workers and their initializer are top-level importable functions. Each call
owns a new `ProcessPoolExecutor` using `multiprocessing.get_context("spawn")`.
Windows user scripts must be importable and guard pool-starting code:

```python
import numpy as np
from qcbed import simulate_tilt_series_parallel


def main():
    # Build prepared_sample using prepare_multislice_sample and explicit state.
    prepared_sample = build_my_prepared_sample()
    tilts = np.array([[0., 0.], [.1, -.2]], dtype=np.float64)
    return simulate_tilt_series_parallel(
        prepared_sample, tilts, max_workers=2,
        process_chunk_size=1, batch_size=1,
    )


if __name__ == "__main__":
    result = main()
```

`build_my_prepared_sample` is the caller's preparation function; see
[single-tilt preparation](single_tilt_simulation.md) for its required state.
Interactive stdin/notebook execution is not promised by this spawn API.
Importing qcbed or the benchmark does not start workers.

The prepared sample is serialized in executor initializer arguments once per
worker, not once per task. Tasks contain only start/stop indices and a tilt
chunk. The internal batch size is also installed by the initializer. Worker
state is private to that new process and is never installed in the parent.
Later executors start fresh workers; no public global cache or persistent
pool lifecycle is introduced.

Prepared dataclasses and arrays are picklable under spawn. NumPy pickle may
restore writable arrays even when the source was bytes-backed/read-only.
The initializer recursively restores read-only flags on prepared arrays.
These flags are defensive, not the original immutable bytes backing: an
owning array can deliberately have its flag re-enabled. The validated batched
engine never writes to prepared arrays. Tests verify read-only worker flags,
unchanged values after tasks, reuse of the same object within a worker, and
unchanged parent state. The parent is not frozen or mutated by the initializer.

## Scheduling, assembly and failure

All tilt/type/shape/settings validation happens in the parent before creating
workers. Existing prepared-sample compatibility checks are rerun without
reconstructing reciprocal grids or gratings. Finite inputs that overflow a
later physical calculation can still fail in a worker; no extreme-value
clipping or guessed numerical repair is introduced.

Tasks cover zero-based half-open `[start:stop]` ranges, contiguous and
nonoverlapping. At most `2*worker_limit` futures are pending. Completion can
be out of order, but each result is placed into its original index interval.
Workers return only the interval, beam amplitudes/intensities and private
participation metadata; no waves or propagators are transferred back.
Final intensities are derived from assembled complex amplitudes using the
existing result invariant.

Task/submission errors raise RuntimeError identifying the failed range, with
the original exception as `__cause__`. Pending futures are cancelled where
possible, and executor shutdown waits for already-running tasks. Cancellation
does not forcibly terminate a currently executing kernel. No partial result
or silent serial fallback is returned. Initializer/process failures surface
through their affected future; pool construction failures propagate directly
because no tilt task has yet been assigned.

## Threads and memory

The package never changes process-wide environment variables and adds no
thread-control dependency. Set these before launching Python so spawned
workers inherit them:

```powershell
$env:OMP_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
$env:OPENBLAS_NUM_THREADS='1'
$env:NUMEXPR_NUM_THREADS='1'
& 'C:\Users\peng863\AppData\Local\miniconda3\envs\qcbed-python\python.exe' benchmarks/benchmark_parallel.py
```

The installed NumPy/FFT/BLAS stack determines whether and how these settings
are honored. More processes combined with internal library threads can
oversubscribe CPU resources. Tilt computation is process-based; the standard
executor's internal IPC-management machinery is not a thread-based physics
backend. No assumption about universal NumPy FFT threading is made.

Each worker owns a full prepared-sample copy plus local batched waves and
propagators. For S stored types, M mesh points and B internal batch entries,
working storage grows roughly as `16*M*B*(S+constant)` bytes, in addition to
the prepared sample and Python/NumPy runtime. Process-result buffering is
bounded by the task window. The parent necessarily stores all requested
beam results and tilt coordinates; there is no all-tilt wave allocation.
No shared memory or automatic memory-based worker selection is provided.

## Validation and performance

`tests/test_parallel.py` includes real Windows spawn tests and lightweight
scheduler tests. Coverage includes one/two workers, fewer tilts than workers,
both reciprocal orders, absorption on/off, two stored types with A-B-A-B,
shifts/dilation, duplicate tilts/beams, internal batch sizes above/below task
sizes, partial chunks, pickling, invalid inputs before spawn, out-of-order
completion, queue bounds, worker failures, cancellation and immutable results.
Changing one tilt leaves other complex amplitudes and intensities exactly
unchanged. Amplitude tolerance is 3e-15 absolute and intensity tolerance 4e-15
for the small comparisons against serial and batched references.

Final validation with the required interpreter: **149 passed in 5.94 s** for
`tests/test_parallel.py tests/test_tilt_series.py` (62 new parallel tests),
and **1197 passed in 26.02 s** for the complete suite, including existing
reference-input checks. The sandbox denied Windows named-pipe creation; the
spawn tests and benchmark were rerun successfully with execution outside that
sandbox. `git diff --check` passed. No executable-Fortran equivalence is claimed.

`benchmarks/benchmark_parallel.py` uses a guarded entry point and reports
single end-to-end calls including pool startup, prepared transfer and
shutdown. It compares serial, one-process batch sizes 8/16/32, and requested
2/4 workers with process/batch pairs 16/8,32/16,64/32. Both workloads use two
elastic-only types, five beams, oblique cells, shifts, unequal dilation and
FFT-prepared gratings. The small workload is 16x16,32 tilts,20 occurrences;
the larger is 64x64,256 tilts,256 occurrences. Preparation is timed separately.

The approximate per-worker array budget uses serialized prepared size plus
`16*M*B*(S+6)` and task input/output arrays. It is not measured RSS and excludes
interpreter/library overhead, opaque FFT scratch and IPC copies. Multiply by
the worker limit for a rough array-only pool budget; actual process memory is
higher. Benchmark speedups include startup, and small tasks may use fewer
workers than requested. No timing assertions are placed in pytest.

An initial synthetic long absorptive benchmark overflowed its float64
intensities, so the timing workload uses elastic-only gratings. That numerical
domain issue was not changed by this scheduling task; short absorptive spawn
tests pass. Neither timing data nor Python agreement establishes historical
Fortran executable equivalence.

### Measured timings

Measured on 2026-09-16 with the exact qcbed-python interpreter above:
Python 3.11.16, NumPy 2.4.6, 12 reported logical CPUs,
`platform.platform() = Windows-10-10.0.26200-SP0`. All four documented thread
environment variables were explicitly set to 1 in the launching PowerShell
process. These are individual local end-to-end measurements, not portable
performance guarantees. Speedups use the serial and fastest measured
single-process batched row of the same workload.

Small workload preparation: **0.012630 s**; prepared pickle: **0.030 MiB**.

| Execution (requested/participating workers) | Process chunk / batch | Seconds | ms/tilt | vs serial | vs best batch | Array budget MiB/worker |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Serial | -- | 0.082462 | 2.5769 | 1.000 | 0.127 | 0.06 |
| Batched | -- / 8 | 0.013943 | 0.4357 | 5.914 | 0.751 | 0.28 |
| Batched | -- / 16 | 0.013012 | 0.4066 | 6.338 | 0.804 | 0.53 |
| Batched | -- / 32 | 0.010464 | 0.3270 | 7.881 | 1.000 | 1.03 |
| Parallel 2/1 | 16 / 8 | 1.078888 | 33.7153 | 0.076 | 0.010 | 0.28 |
| Parallel 2/1 | 32 / 16 | 0.476486 | 14.8902 | 0.173 | 0.022 | 0.53 |
| Parallel 2/1 | 64 / 32 | 0.356860 | 11.1519 | 0.231 | 0.029 | 1.03 |
| Parallel 4/1 | 16 / 8 | 0.632999 | 19.7812 | 0.130 | 0.017 | 0.28 |
| Parallel 4/1 | 32 / 16 | 0.360608 | 11.2690 | 0.229 | 0.029 | 0.53 |
| Parallel 4/1 | 64 / 32 | 0.398996 | 12.4686 | 0.207 | 0.026 | 1.03 |

There are only one or two process tasks in this small case. The pool limit
is capped accordingly; one worker finished all tasks in each measured run.

Larger workload preparation: **0.093683 s**; prepared pickle: **0.440 MiB**.

| Execution (requested/participating workers) | Process chunk / batch | Seconds | ms/tilt | vs serial | vs best batch | Array budget MiB/worker |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Serial | -- | 20.836934 | 81.3943 | 1.000 | 0.648 | 0.97 |
| Batched | -- / 8 | 13.499655 | 52.7330 | 1.544 | 1.000 | 4.47 |
| Batched | -- / 16 | 16.781307 | 65.5520 | 1.242 | 0.804 | 8.47 |
| Batched | -- / 32 | 20.302601 | 79.3070 | 1.026 | 0.665 | 16.47 |
| Parallel 2/2 | 16 / 8 | 9.232947 | 36.0662 | 2.257 | 1.462 | 4.44 |
| Parallel 2/2 | 32 / 16 | 11.151700 | 43.5613 | 1.868 | 1.211 | 8.44 |
| Parallel 2/2 | 64 / 32 | 12.015540 | 46.9357 | 1.734 | 1.124 | 16.45 |
| Parallel 4/4 | 16 / 8 | 7.052545 | 27.5490 | 2.955 | 1.914 | 4.44 |
| Parallel 4/4 | 32 / 16 | 9.775254 | 38.1846 | 2.132 | 1.381 | 8.44 |
| Parallel 4/4 | 64 / 32 | 9.678145 | 37.8053 | 2.153 | 1.395 | 16.45 |

The maximum measured serial/batched/parallel amplitude and intensity residuals
were **0** in both benchmark workloads. Separate spawn correctness tests use
the tolerances above, including absorption and both orderings; no-cross-talk
checks require unchanged neighboring results exactly.

**Recommendation is workload-dependent:** use one-process batching for the
small case; for the larger measured case use four workers, process chunks of
16 and internal batches of eight. That setting was 1.914 times faster than
the best measured one-process batched setting and 2.955 times faster than
serial, including startup. Larger internal batches were slower here, so do
not select a batch size merely because more memory is available.

Deferred: shared-memory preparation, persistent pools, distributed computing,
thread-based tilt execution, GPU, coherent probes, scan/experimental parsing,
YNORM, chi-square/refinement and CLI. Shared read-only storage or a future GPU
backend could reduce duplication; the current physics/backend boundary is
unchanged.
