from concurrent.futures import Future
from dataclasses import FrozenInstanceError, fields, is_dataclass, replace
import os
import pickle

import numpy as np
import pytest

from qcbed import (
    AtomSite, Species, MeshConfig, SliceGeometry, AbsorptionCoefficients,
    SliceAppliedParameters, SamplePreparationParameters, prepare_multislice_sample,
    ParallelExecutionMetadata, simulate_tilt_series_serial, simulate_tilt_series_batched,
    simulate_tilt_series_parallel,
)
from qcbed.legacy_io import LegacyInput
import qcbed.parallel as parallel


TILTS = np.array([[.2, -.1], [0., 0.], [-.3, .4], [.2, -.1], [.1, .6], [-.5, -.2], [0., -.1]])
_observed_prepared_id = None


def prepared_sample(ordering="centered", absorption=True):
    model = LegacyInput("parallel", 200, False, (0, 0, 1), MeshConfig(6, 4),
        (SliceGeometry(2, .125, -.25, 20, 10, 60), SliceGeometry(3, -.3, .2, 20, 10, 60)),
        (Species("Al", 13),), ((AtomSite(1, .25, .125, .5, "S"),),
                             (AtomSite(1, .375, .25, .75, "S"),)),
        "out.dat", "in.dat", title_raw="parallel", absorption_raw="YES" if absorption else "NO")
    parameters = SamplePreparationParameters(200,
        (SliceAppliedParameters((.86,), .1, 1.1), SliceAppliedParameters((1.2,), .2, .95)),
        .5, (AbsorptionCoefficients(0, 1, 1),) if absorption else None, (0, 1, 0, 1),
        ((1, 0), (0, -1), (0, 0), (-1, 0), (1, 0)), ordering, "fft", "legacy", "unity")
    return prepare_multislice_sample(model, parameters)


def arrays(value):
    if isinstance(value, np.ndarray):
        yield value
    elif is_dataclass(value):
        for field in fields(value):
            yield from arrays(getattr(value, field.name))
    elif isinstance(value, tuple):
        for entry in value:
            yield from arrays(entry)


def probe_worker(start, stop, tilts):
    """Importable spawn probe: verify read-only data and reuse within a worker."""
    global _observed_prepared_id
    prepared = parallel._worker_prepared
    assert prepared is not None
    assert all(not array.flags.writeable for array in arrays(prepared))
    if _observed_prepared_id is not None:
        assert _observed_prepared_id == id(prepared)
    _observed_prepared_id = id(prepared)
    before = pickle.dumps(prepared)
    result = simulate_tilt_series_batched(prepared, tilts, batch_size=parallel._worker_batch_size)
    assert pickle.dumps(prepared) == before
    return start, stop, result.amplitudes, result.intensities, os.getpid()


@pytest.mark.parametrize("ordering,absorption,batch_size", [("centered", True, 2), ("unshifted", False, 1)])
def test_two_spawned_workers_agree_with_serial_and_batch(ordering, absorption, batch_size):
    prepared = prepared_sample(ordering, absorption)
    snapshot = pickle.dumps(prepared)
    tilts = TILTS.copy()
    serial = simulate_tilt_series_serial(prepared, tilts)
    batch = simulate_tilt_series_batched(prepared, tilts, batch_size=batch_size)
    result = simulate_tilt_series_parallel(prepared, tilts, max_workers=2,
        process_chunk_size=3, batch_size=batch_size)
    for expected in (serial, batch):
        np.testing.assert_allclose(result.amplitudes, expected.amplitudes, atol=3e-15, rtol=0)
        np.testing.assert_allclose(result.intensities, expected.intensities, atol=4e-15, rtol=0)
    np.testing.assert_array_equal(result.tilts, tilts)
    np.testing.assert_array_equal(tilts, TILTS)
    np.testing.assert_array_equal(result.amplitudes[0], result.amplitudes[3])
    np.testing.assert_array_equal(result.amplitudes[:, 0], result.amplitudes[:, 4])
    assert result.beam_indices_hk == prepared.parameters.beam_indices_hk
    assert result.execution_mode == "parallel" and result.batch_size == batch_size
    assert result.parallel.requested_workers == result.parallel.worker_limit == 2
    assert 1 <= result.parallel.actual_workers <= 2
    assert result.parallel.start_method == "spawn" and result.parallel.process_chunk_size == 3
    assert pickle.dumps(prepared) == snapshot
    assert parallel._worker_prepared is parallel._worker_batch_size is None
    assert result.amplitudes.dtype == np.complex128 and result.intensities.dtype == np.float64
    for array in (result.tilts, result.amplitudes, result.intensities):
        with pytest.raises(ValueError):
            array.setflags(write=True)
    with pytest.raises(FrozenInstanceError):
        result.parallel.actual_workers = 4


def test_spawn_no_cross_talk():
    prepared = prepared_sample()
    tilts = TILTS[:3].copy()
    before = simulate_tilt_series_parallel(prepared, tilts, max_workers=2,
        process_chunk_size=2, batch_size=2)
    tilts[1] = [7, -9]
    after = simulate_tilt_series_parallel(prepared, tilts, max_workers=2,
        process_chunk_size=2, batch_size=2)
    np.testing.assert_array_equal(before.amplitudes[[0, 2]], after.amplitudes[[0, 2]])
    np.testing.assert_array_equal(before.intensities[[0, 2]], after.intensities[[0, 2]])
    assert np.max(abs(before.amplitudes[1]-after.amplitudes[1])) > 1e-8


def test_one_worker_pickle_readonly_and_reuse(monkeypatch):
    monkeypatch.setattr(parallel, "_worker_chunk", probe_worker)
    prepared = prepared_sample()
    restored = pickle.loads(pickle.dumps(prepared))
    for original, copy in zip(arrays(prepared), arrays(restored)):
        np.testing.assert_array_equal(original, copy)
    expected = simulate_tilt_series_batched(restored, TILTS, batch_size=2)
    result = simulate_tilt_series_parallel(restored, TILTS, max_workers=1,
        process_chunk_size=2, batch_size=2)
    np.testing.assert_allclose(result.amplitudes, expected.amplitudes, atol=3e-15, rtol=0)
    assert result.parallel.actual_workers == result.parallel.worker_limit == 1
    assert result.execution_mode == "parallel"  # One worker still uses spawn.
    assert parallel._worker_prepared is None


@pytest.mark.parametrize("count,chunk", [(1, 10), (2, 1)])
def test_fewer_tilts_than_workers_and_large_internal_batch(count, chunk):
    prepared = prepared_sample()
    tilts = TILTS[:count]
    expected = simulate_tilt_series_batched(prepared, tilts, batch_size=8)
    result = simulate_tilt_series_parallel(prepared, tilts, max_workers=4,
        process_chunk_size=chunk, batch_size=8)
    np.testing.assert_allclose(result.amplitudes, expected.amplitudes, atol=3e-15, rtol=0)
    assert result.parallel.requested_workers == 4
    assert result.parallel.worker_limit == count


def test_real_worker_failure_has_range_cause_and_no_partial_result():
    tilts = np.array([[0., 0.], [np.finfo(np.float64).max, 0.], [.1, .2]])
    # Valid finite inputs overflow the downstream physical phase in the worker.
    with pytest.raises(RuntimeError, match=r"tilt chunk \[1:2\] failed") as error:
        simulate_tilt_series_parallel(prepared_sample(), tilts, max_workers=2,
            process_chunk_size=1, batch_size=1)
    assert isinstance(error.value.__cause__, ValueError)
    assert "finite" in str(error.value.__cause__)
    assert parallel._worker_prepared is None


class ForbiddenExecutor:
    def __init__(self, **kwargs):
        raise AssertionError("invalid inputs must fail before creating the executor")


@pytest.mark.parametrize("name", ["max_workers", "process_chunk_size", "batch_size"])
@pytest.mark.parametrize("bad", [True, np.bool_(False), 0, -1, 1., .5, np.nan, np.inf, "2", None])
def test_invalid_integer_settings_before_spawn(monkeypatch, name, bad):
    monkeypatch.setattr(parallel, "ProcessPoolExecutor", ForbiddenExecutor)
    kwargs = dict(max_workers=2, process_chunk_size=3, batch_size=2)
    kwargs[name] = bad
    with pytest.raises((TypeError, ValueError), match=name):
        simulate_tilt_series_parallel(prepared_sample(), TILTS, **kwargs)


@pytest.mark.parametrize("bad", [[], np.zeros((0, 2)), np.zeros(2), np.zeros((2, 3)),
    np.zeros((2, 1, 2)), np.array([[np.nan, 0.]]), np.array([[0., np.inf]]),
    np.array([[0., -np.inf]]), np.zeros((2, 2), dtype=int), np.zeros((2, 2), dtype=bool)])
def test_invalid_tilts_before_spawn(monkeypatch, bad):
    monkeypatch.setattr(parallel, "ProcessPoolExecutor", ForbiddenExecutor)
    with pytest.raises((TypeError, ValueError), match="tilts"):
        simulate_tilt_series_parallel(prepared_sample(), bad, max_workers=2, process_chunk_size=3, batch_size=2)


@pytest.mark.parametrize("method", ["fork", "forkserver", "auto", None, True])
def test_invalid_start_method_before_spawn(monkeypatch, method):
    monkeypatch.setattr(parallel, "ProcessPoolExecutor", ForbiddenExecutor)
    with pytest.raises(ValueError, match="start_method"):
        simulate_tilt_series_parallel(prepared_sample(), TILTS, max_workers=2,
            process_chunk_size=3, batch_size=2, start_method=method)


def test_prepared_validation_before_spawn(monkeypatch):
    monkeypatch.setattr(parallel, "ProcessPoolExecutor", ForbiddenExecutor)
    with pytest.raises(TypeError, match="prepared_sample"):
        simulate_tilt_series_parallel(None, TILTS, max_workers=2, process_chunk_size=3, batch_size=2)
    restored = pickle.loads(pickle.dumps(prepared_sample()))
    object.__setattr__(restored, "mesh", MeshConfig(4, 4))
    with pytest.raises(ValueError, match="mesh shape"):
        simulate_tilt_series_parallel(restored, TILTS, max_workers=2, process_chunk_size=3, batch_size=2)
    restored = pickle.loads(pickle.dumps(prepared_sample()))
    object.__setattr__(restored.slice_types[0], "phase_grating",
                       np.full((4, 6), np.nan, dtype=np.complex128))
    with pytest.raises(ValueError, match="phase_grating.*finite"):
        simulate_tilt_series_parallel(restored, TILTS, max_workers=2, process_chunk_size=3, batch_size=2)


class RecordingExecutor:
    """Parent-only scheduler double; true spawn is exercised separately above."""
    instances = []
    fail = False

    def __init__(self, *, max_workers, mp_context, initializer, initargs):
        assert mp_context.get_start_method() == "spawn"
        assert initializer is parallel._initialize_worker
        self.prepared, self.batch_size = initargs
        self.max_workers = max_workers
        self.submissions, self.futures = [], []
        self.shutdown_args = None
        self.instances.append(self)

    def submit(self, function, start, stop, tilts):
        assert function is parallel._worker_chunk
        assert tilts.shape == (stop-start, 2)
        self.submissions.append((start, stop))
        future = Future()
        future.start = start
        self.futures.append(future)
        if self.fail:
            if start == 0:
                future.set_exception(ValueError("synthetic worker failure"))
        else:
            result = simulate_tilt_series_batched(self.prepared, tilts, batch_size=self.batch_size)
            future.set_result((start, stop, result.amplitudes, result.intensities, 123))
        return future

    def shutdown(self, *, wait, cancel_futures):
        self.shutdown_args = (wait, cancel_futures)


def reverse_completion(pending, *, return_when):
    executor = RecordingExecutor.instances[-1]
    assert len(pending) <= 2*executor.max_workers
    chosen = max(pending, key=future_start)
    return {chosen}, set(pending)-{chosen}


def future_start(future):
    return future.start


@pytest.mark.parametrize("count,chunk,batch,workers", [(1, 1, 1, 1), (2, 1, 8, 8),
    (7, 3, 2, 2), (7, 20, 1, 4), (7, 1, 2, 2)])
def test_partition_indexed_assembly_and_bounded_submission(monkeypatch, count, chunk, batch, workers):
    RecordingExecutor.instances.clear()
    monkeypatch.setattr(parallel, "ProcessPoolExecutor", RecordingExecutor)
    monkeypatch.setattr(parallel, "wait", reverse_completion)
    environment = dict(os.environ)
    prepared = prepared_sample()
    tilts = TILTS[::-1][:count]
    expected = simulate_tilt_series_serial(prepared, tilts)
    result = simulate_tilt_series_parallel(prepared, tilts, max_workers=np.int64(workers),
        process_chunk_size=np.int64(chunk), batch_size=np.int64(batch))
    np.testing.assert_allclose(result.amplitudes, expected.amplitudes, atol=3e-15, rtol=0)
    np.testing.assert_allclose(result.intensities, expected.intensities, atol=4e-15, rtol=0)
    np.testing.assert_array_equal(result.tilts, tilts)
    executor = RecordingExecutor.instances[-1]
    assert executor.submissions == [(i, min(i+chunk, count)) for i in range(0, count, chunk)]
    assert executor.shutdown_args == (True, True)
    assert dict(os.environ) == environment
    assert parallel._worker_prepared is None


def test_pending_cancellation_cleanup_and_no_fallback(monkeypatch):
    RecordingExecutor.instances.clear()
    monkeypatch.setattr(RecordingExecutor, "fail", True)
    monkeypatch.setattr(parallel, "ProcessPoolExecutor", RecordingExecutor)
    with pytest.raises(RuntimeError, match=r"\[0:1\] failed") as error:
        simulate_tilt_series_parallel(prepared_sample(), TILTS, max_workers=2,
            process_chunk_size=1, batch_size=1)
    assert isinstance(error.value.__cause__, ValueError)
    executor = RecordingExecutor.instances[-1]
    assert executor.shutdown_args == (True, True)
    assert len(executor.submissions) == 4
    assert all(future.cancelled() for future in executor.futures[1:])


class SubmissionFailureExecutor(RecordingExecutor):
    def submit(self, function, start, stop, tilts):
        if start:
            raise OSError("synthetic submission failure")
        return super().submit(function, start, stop, tilts)


class WrongIntervalExecutor(RecordingExecutor):
    def submit(self, function, start, stop, tilts):
        future = super().submit(function, start, stop, tilts)
        first, last, amplitudes, intensities, pid = future.result()
        invalid = Future()
        invalid.set_result((first+1, last+1, amplitudes, intensities, pid))
        return invalid


@pytest.mark.parametrize("executor_type,expected_cause,match", [
    (SubmissionFailureExecutor, OSError, r"\[1:2\] submission failed"),
    (WrongIntervalExecutor, ValueError, r"tilt chunk .* failed"),
])
def test_submission_and_bad_result_cleanup(monkeypatch, executor_type, expected_cause, match):
    RecordingExecutor.instances.clear()
    monkeypatch.setattr(parallel, "ProcessPoolExecutor", executor_type)
    with pytest.raises(RuntimeError, match=match) as error:
        simulate_tilt_series_parallel(prepared_sample(), TILTS, max_workers=2,
            process_chunk_size=1, batch_size=1)
    assert isinstance(error.value.__cause__, expected_cause)
    assert RecordingExecutor.instances[-1].shutdown_args == (True, True)


def test_metadata_validation_and_existing_result_contract():
    for changes in (dict(actual_workers=0), dict(actual_workers=3), dict(worker_limit=True),
                    dict(process_chunk_size=0), dict(start_method="fork")):
        fields = dict(requested_workers=2, worker_limit=2, actual_workers=1, process_chunk_size=3)
        fields.update(changes)
        with pytest.raises((ValueError, TypeError)):
            ParallelExecutionMetadata(**fields)
    result = simulate_tilt_series_batched(prepared_sample(), TILTS, batch_size=2)
    assert result.parallel is None
    with pytest.raises(TypeError, match="ParallelExecutionMetadata"):
        replace(result, execution_mode="parallel")
    with pytest.raises(ValueError, match="parallel metadata"):
        replace(result, parallel=ParallelExecutionMetadata(2, 2, 1, 3))
