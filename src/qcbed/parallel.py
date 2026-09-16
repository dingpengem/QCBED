"""Spawned CPU processes running independent chunks through the batched engine."""

from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from dataclasses import fields, is_dataclass, replace
import multiprocessing
import os
import sys

import numpy as np
from numpy.typing import NDArray

from .fft_backend import _array
from .multislice import _integer, _snapshot
from .simulation import PreparedMultisliceSample
from .tilt_series import (
    ParallelExecutionMetadata, TiltSeriesResult, _inputs, _batch_size,
    simulate_tilt_series_batched,
)


# Set only by a spawned executor's initializer, never by the public parent API.
_worker_prepared: PreparedMultisliceSample | None = None
_worker_batch_size: int | None = None


def _readonly_arrays(value) -> None:
    """Restore defensive NumPy flags lost in pickle, only on worker-owned data."""
    if isinstance(value, np.ndarray):
        value.setflags(write=False)
    elif is_dataclass(value):
        for entry in fields(value):
            _readonly_arrays(getattr(value, entry.name))
    elif isinstance(value, tuple):
        for entry in value:
            _readonly_arrays(entry)


def _initialize_worker(prepared: PreparedMultisliceSample, batch_size: int) -> None:
    global _worker_prepared, _worker_batch_size
    _readonly_arrays(prepared)
    _worker_prepared = prepared
    _worker_batch_size = batch_size


def _worker_chunk(start: int, stop: int, tilts: NDArray[np.float64]):
    if _worker_prepared is None or _worker_batch_size is None:
        raise RuntimeError("tilt worker was not initialized")
    result = simulate_tilt_series_batched(_worker_prepared, tilts, batch_size=_worker_batch_size)
    # PID is used only to count participating workers; it is not exposed.
    return start, stop, result.amplitudes, result.intensities, os.getpid()


def _positive_integer(value: int, name: str) -> int:
    value = _integer(value, name)
    if value <= 0:
        raise ValueError(f"{name} must be positive")
    return value


def _ranges(count: int, size: int):
    for start in range(0, count, size):
        yield start, min(start+size, count)


def simulate_tilt_series_parallel(
    prepared_sample: PreparedMultisliceSample, tilts: NDArray[np.float64], *,
    max_workers: int, process_chunk_size: int, batch_size: int, start_method: str = "spawn",
) -> TiltSeriesResult:
    """Execute contiguous tilt tasks with bounded submission and indexed assembly.

    process_chunk_size controls task payloads; batch_size controls simultaneous
    spatial FFTs inside each task. Only spawn is supported, including for one
    worker. Call from an importable script behind an if __name__ == '__main__'
    guard. Numerical-library thread limits must be set by the caller before
    launch; this module never changes the environment. No serial fallback.
    """
    values = _inputs(prepared_sample, tilts)
    requested = _positive_integer(max_workers, "max_workers")
    chunk_size = _positive_integer(process_chunk_size, "process_chunk_size")
    batch_size = _batch_size(batch_size)
    if not isinstance(start_method, str) or start_method != "spawn":
        raise ValueError("start_method must be 'spawn'")
    # Reuse the existing compatibility checks without reconstructing grids or
    # gratings. This also validates a normally unpickled prepared object.
    prepared = replace(prepared_sample)
    for i, stored in enumerate(prepared.slice_types):
        grating = _array(stored.phase_grating, f"slice type {i} phase_grating")
        if grating.shape != (prepared.mesh.mesh_y, prepared.mesh.mesh_x):
            raise ValueError(f"slice type {i} phase_grating must match sample mesh")
    values = _snapshot(values)
    task_count = (len(values)+chunk_size-1)//chunk_size
    worker_limit = min(requested, task_count)
    # ProcessPoolExecutor's Windows wait-handle limit, not a physics capacity.
    if sys.platform == "win32":
        worker_limit = min(worker_limit, 61)
    beams = prepared.parameters.beam_indices_hk
    amplitudes = np.empty((len(values), len(beams)), dtype=np.complex128)
    ranges = iter(_ranges(len(values), chunk_size))
    pending = {}
    worker_pids = set()
    executor = ProcessPoolExecutor(max_workers=worker_limit,
        mp_context=multiprocessing.get_context(start_method),
        initializer=_initialize_worker, initargs=(prepared, batch_size))
    try:
        exhausted = False
        while pending or not exhausted:
            # Bound queued tilt copies and returned chunk buffers, independently
            # of the total number of requested experiments.
            while not exhausted and len(pending) < 2*worker_limit:
                indices = next(ranges, None)
                if indices is None:
                    exhausted = True
                    break
                start, stop = indices
                try:
                    future = executor.submit(_worker_chunk, start, stop, values[start:stop])
                except Exception as exc:
                    raise RuntimeError(f"parallel tilt chunk [{start}:{stop}] submission failed") from exc
                pending[future] = indices
            if not pending:
                break
            done, _ = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                start, stop = pending.pop(future)
                try:
                    result_start, result_stop, chunk_amplitudes, chunk_intensities, pid = future.result()
                    if (result_start, result_stop) != (start, stop):
                        raise ValueError("worker returned incorrect tilt indices")
                    shape = (stop-start, len(beams))
                    if (chunk_amplitudes.shape != shape or chunk_amplitudes.dtype != np.complex128
                            or chunk_intensities.shape != shape or chunk_intensities.dtype != np.float64
                            or not np.all(np.isfinite(chunk_amplitudes))
                            or not np.all(np.isfinite(chunk_intensities))):
                        raise ValueError("worker returned invalid beam data")
                    amplitudes[start:stop] = chunk_amplitudes
                    worker_pids.add(pid)
                except Exception as exc:
                    raise RuntimeError(f"parallel tilt chunk [{start}:{stop}] failed") from exc
            # Drop completed futures/results before refilling the bounded queue.
            done.clear()
            del future, chunk_amplitudes, chunk_intensities
    finally:
        for future in pending:
            future.cancel()
        executor.shutdown(wait=True, cancel_futures=True)
    metadata = ParallelExecutionMetadata(requested, worker_limit, len(worker_pids), chunk_size, start_method)
    # Keep the existing result invariant: intensities are derived once from the
    # assembled complex amplitudes, independent of process completion order.
    return TiltSeriesResult(values, amplitudes, beams, "parallel", batch_size, metadata)
