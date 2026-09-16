"""End-to-end CPU process benchmark; launch with numerical thread limits set.

Every parallel timing includes a new executor, spawn, sample transfer and
shutdown. Single measured calls, not timing assertions or warm-pool timings.
"""

import os
import pickle
import platform
import sys
from time import perf_counter

import numpy as np

from qcbed import (
    AtomSite, Species, MeshConfig, SliceGeometry,
    SliceAppliedParameters, SamplePreparationParameters, prepare_multislice_sample,
    simulate_tilt_series_serial, simulate_tilt_series_batched, simulate_tilt_series_parallel,
)
from qcbed.legacy_io import LegacyInput


def prepare(mesh_size, occurrences):
    model = LegacyInput("process benchmark", 200, False, (0, 0, 1), MeshConfig(mesh_size, mesh_size),
        (SliceGeometry(2, .125, -.2, 80, 60, 60), SliceGeometry(3, -.1, .3, 80, 60, 60)),
        (Species("Al", 13),), ((AtomSite(1, .25, .125, .5, "S"),),
                             (AtomSite(1, .375, .25, .75, "S"),)),
        "out.dat", "in.dat", title_raw="process benchmark", absorption_raw="NO")
    parameters = SamplePreparationParameters(200,
        (SliceAppliedParameters((.86,), .1, 1.1), SliceAppliedParameters((1.2,), .2, .95)),
        .5, None, (0, 1)*(occurrences//2),
        ((0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)), "centered", "fft", "legacy", "unity")
    return prepare_multislice_sample(model, parameters)


def working_mib(prepared_bytes, mesh, batch, chunk, beams=5):
    # Two kernel stacks plus an allowance of six complex grids per batch entry.
    # Excludes Python/NumPy runtime, opaque FFT scratch and IPC buffers.
    return (prepared_bytes+16*mesh*mesh*batch*(2+6)+chunk*(24*beams+16))/2**20


def benchmark(mesh, count, occurrences):
    start = perf_counter()
    prepared = prepare(mesh, occurrences)
    preparation_time = perf_counter()-start
    prepared_bytes = len(pickle.dumps(prepared, protocol=pickle.HIGHEST_PROTOCOL))
    t = np.arange(count, dtype=np.float64)
    tilts = np.column_stack((.007*(t-count/2), .013*((t*7) % 19-9)))
    # Warm imports and parent numerical dispatch only, outside timed runs.
    simulate_tilt_series_batched(prepared, tilts[:1], batch_size=1)
    start = perf_counter()
    reference = simulate_tilt_series_serial(prepared, tilts)
    serial_time = perf_counter()-start
    rows = [("serial", serial_time, 1, count, 1, 0., 0.)]
    batch_times = {}
    for batch in (8, 16, 32):
        start = perf_counter()
        actual = simulate_tilt_series_batched(prepared, tilts, batch_size=batch)
        elapsed = perf_counter()-start
        batch_times[batch] = elapsed
        rows.append((f"batched {batch}", elapsed, batch, count, 1,
                     float(np.max(abs(actual.amplitudes-reference.amplitudes))),
                     float(np.max(abs(actual.intensities-reference.intensities)))))
    print(f"\nmesh={mesh}x{mesh}; tilts={count}; stored types=2; occurrences={occurrences}; beams=5; elastic only", flush=True)
    print(f"preparation={preparation_time:.6f}s; prepared pickle={prepared_bytes/2**20:.3f} MiB", flush=True)
    print("execution | seconds | ms/tilt | serial speedup | best batched speedup | approximate MiB/worker | amplitude residual | intensity residual", flush=True)
    baseline = min(batch_times.values())
    for name, elapsed, batch, chunk, workers, ae, ie in rows:
        print(f"{name} | {elapsed:.6f} | {elapsed/count*1000:.4f} | {serial_time/elapsed:.3f} | {baseline/elapsed:.3f} | {working_mib(prepared_bytes,mesh,min(batch,count),chunk):.2f} | {ae:.3e} | {ie:.3e}", flush=True)
    for workers in (2, 4):
        if (os.cpu_count() or 1) < workers:
            print(f"Skipping {workers} workers: fewer reported logical CPUs", flush=True)
            continue
        for chunk, batch in ((16, 8), (32, 16), (64, 32)):
            start = perf_counter()
            actual = simulate_tilt_series_parallel(prepared, tilts, max_workers=workers,
                process_chunk_size=chunk, batch_size=batch)
            elapsed = perf_counter()-start
            ae = float(np.max(abs(actual.amplitudes-reference.amplitudes)))
            ie = float(np.max(abs(actual.intensities-reference.intensities)))
            np.testing.assert_allclose(actual.amplitudes, reference.amplitudes, atol=3e-12, rtol=0)
            np.testing.assert_allclose(actual.intensities, reference.intensities, atol=3e-12, rtol=0)
            memory = working_mib(prepared_bytes, mesh, min(batch,chunk,count), min(chunk,count))
            print(f"parallel {workers}, chunk {chunk}, batch {batch} (used {actual.parallel.actual_workers}) | {elapsed:.6f} | {elapsed/count*1000:.4f} | {serial_time/elapsed:.3f} | {baseline/elapsed:.3f} | {memory:.2f} | {ae:.3e} | {ie:.3e}", flush=True)


def main():
    print(f"Platform: {platform.platform()}; Python: {sys.version.split()[0]}; NumPy: {np.__version__}")
    print(f"Interpreter: {sys.executable}; logical CPUs: {os.cpu_count()}")
    print("Thread environment:", {name: os.environ.get(name, "unset") for name in
        ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS")})
    print("All process timings include spawn, transfer and shutdown. Memory is an array-budget estimate, not RSS.")
    benchmark(16, 32, 20)
    benchmark(64, 256, 256)


if __name__ == "__main__":
    main()
