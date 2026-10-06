from dataclasses import FrozenInstanceError, fields, is_dataclass, replace
import json
from pathlib import Path

import numpy as np
import pytest

from qcbed import ForwardExecutionConfig
from qcbed.refinement import pack_refinement_parameters, RefinementReflection
from qcbed.refinement_trial import evaluate_refinement_trial
from qcbed.refinement_driver import run_legacy_refinement
from qcbed.local_runner import run_local_optimization
from test_refinement import data
from test_refinement_driver import problem, config


def assert_exact(a, b):
    if is_dataclass(a):
        for f in fields(a):
            assert_exact(getattr(a, f.name), getattr(b, f.name))
    elif isinstance(a, (tuple, list)):
        assert len(a) == len(b)
        for x, y in zip(a, b):
            assert_exact(x, y)
    elif isinstance(a, np.ndarray):
        np.testing.assert_array_equal(a, b)
        assert a.dtype == b.dtype
    else:
        assert a == b


def test_execution_defaults_and_immutability():
    assert ForwardExecutionConfig().execution_backend == "serial"
    assert ForwardExecutionConfig().batch_size is None
    c = ForwardExecutionConfig("batched")
    assert c.batch_size == 64
    with pytest.raises(FrozenInstanceError):
        c.batch_size = 32


@pytest.mark.parametrize("kwargs", [dict(execution_backend="unknown"),
    dict(execution_backend="parallel"), dict(batch_size=64),
    *[dict(execution_backend="batched", batch_size=v) for v in (0, -1, True, 1.5)],
    dict(max_workers=2), dict(process_chunk_size=4)])
def test_invalid_execution(kwargs):
    with pytest.raises((ValueError, TypeError)):
        ForwardExecutionConfig(**kwargs)


@pytest.mark.parametrize("mode", ["legacy", "masked"])
@pytest.mark.parametrize("gradient", [0, 1])
@pytest.mark.parametrize("iref", [-1., 0., 1.])
@pytest.mark.parametrize("batch_size", [1, 3, 64])
def test_trial_exact_and_batch_invariance(mode, gradient, iref, batch_size):
    s = problem()
    # Nonzero tilts, rectangular wave mesh, multiple slice occurrences and a
    # nonunit reflection multiplier catch accidental reuse of base gratings.
    s = replace(s, iref=iref, scan=replace(s.scan, corner_a=(-.2, .1),
        corner_b=(.3, .1), corner_c=(-.2, .4)),
        reflections=((RefinementReflection(0, 0, 1.7, .8, "NN"),),))
    e = data(gradient)
    if mode == "masked":
        e = replace(e, weights=np.array([[0.], [2.], [.5], [1.]]),
                    variances=np.array([[0.], [1.], [2.], [3.]]))
    v = np.array([5.])
    before = repr(s)
    serial = evaluate_refinement_trial(s, v, e, objective_mode=mode)
    seen = []
    batch = evaluate_refinement_trial(s, v, e, objective_mode=mode,
        execution=ForwardExecutionConfig("batched", batch_size), evaluation_callback=seen.append)
    assert_exact(serial, batch)
    assert seen == [batch]
    assert repr(s) == before
    np.testing.assert_array_equal(v, [5.])
    for array in (batch.amplitudes, batch.calculated_intensities, batch.objective_intensities):
        with pytest.raises(ValueError):
            array.setflags(write=True)


def test_backend_dispatch_and_invalid_type(monkeypatch):
    import qcbed.refinement_trial as module
    original = module.simulate_tilt_series_batched
    calls = []
    def tracked(*args, **kwargs):
        calls.append(kwargs)
        return original(*args, **kwargs)
    monkeypatch.setattr(module, "simulate_tilt_series_batched", tracked)
    s = problem(); v = pack_refinement_parameters(s)
    evaluate_refinement_trial(s, v, data())
    assert calls == []
    evaluate_refinement_trial(s, v, data(), execution=ForwardExecutionConfig("batched"))
    assert calls == [dict(batch_size=64)]
    with pytest.raises(TypeError, match="execution"):
        evaluate_refinement_trial(s, v, data(), execution="batched")


@pytest.mark.parametrize("method", ["direct", "fft"])
def test_multiple_stored_slices_and_reflection_corrections(method):
    s = problem()
    geometry = replace(s.model.slices[0], thickness_angstrom=1., x_shift_fractional=.2)
    model = replace(s.model, slices=(s.model.slices[0], geometry),
        atoms_by_slice=(s.model.atoms_by_slice[0],
                       (replace(s.model.atoms_by_slice[0][0], x_fractional=.3),)))
    s = replace(s, model=model, regions=(replace(s.regions[0], slice_count=2),),
        reflections=((RefinementReflection(0, 0, 1.3, 1., "NN"),),))
    v = np.array([7.])
    a = evaluate_refinement_trial(s, v, data(1), phase_grating_method=method)
    b = evaluate_refinement_trial(s, v, data(1), phase_grating_method=method,
        execution=ForwardExecutionConfig("batched", 3))
    assert_exact(a, b)


def test_reference_absorption_and_nonunit_elastic_absorptive_multipliers():
    from qcbed.refinement_io import parse_refinement_state
    from qcbed.refinement import build_refinement_parameter_layout
    original = parse_refinement_state(Path(__file__).parents[1] / "tests/fixtures/legacy_reference" / "FOR007.DAT")
    s = replace(original, dsap=.1,
        regions=(replace(original.regions[0], depth_angstrom=6.074),),
        scan=replace(original.scan, scan_rows=2, scan_columns=2))
    v = pack_refinement_parameters(s).copy()
    count = 0
    for p in build_refinement_parameter_layout(s).parameters:
        if p.category == "reflection_multiplier":
            v[p.vector_index] = 1.3 if p.target[-1] == 0 else .7
            count += 1
    assert count > 0 and s.model.absorption_enabled
    a = evaluate_refinement_trial(s, v, data(1))
    b = evaluate_refinement_trial(s, v, data(1), execution=ForwardExecutionConfig("batched", 3))
    assert_exact(a, b)


@pytest.mark.parametrize("method", ["ANNE", "SIMP", "PARA", "SPAR", "ANSI", "PAAN"])
@pytest.mark.parametrize("mode", ["legacy", "masked"])
def test_optimizer_decisions_and_stage_handoff_exact(method, mode):
    s = problem(); e = data()
    c = config(method, maxfsimp=2, maxfpara=0)
    events = []
    serial = run_legacy_refinement(s, e, c, objective_mode=mode,
        progress_callback=lambda *args: events.append(args))
    actual = []
    batch = run_legacy_refinement(s, e, replace(c, execution=ForwardExecutionConfig("batched", 3)),
        objective_mode=mode, progress_callback=lambda *args: actual.append(args))
    # Includes every evaluated objective/vector in order and all optimizer
    # bookkeeping, acceptance/adaptation counters, retained states and seeds.
    assert_exact(serial, batch)
    assert_exact(events, actual)


def test_local_runner_metadata_and_observation(tmp_path):
    pytest.importorskip("matplotlib")
    s = problem(); e = data(); c = config("ANSI", maxfsimp=0)
    serial = run_local_optimization(s, e, c, objective_mode="masked")
    c = replace(c, execution=ForwardExecutionConfig("batched", 3))
    path = tmp_path / "run"
    batch = run_local_optimization(s, e, c, objective_mode="masked",
        snapshot_every=2, output_directory=path)
    plain = run_local_optimization(s, e, c, objective_mode="masked", record_progress=False)
    assert_exact(serial.validation, batch.validation)
    assert_exact(batch.validation, plain.validation)
    assert_exact(serial.progress, batch.progress)
    assert batch.execution == c.execution
    saved = json.loads((path / "summary.json").read_text())
    assert saved["execution_backend"] == "batched"
    assert saved["batch_size"] == 3 and saved["process_workers"] is None
    assert saved["objective_mode"] == "masked"
    html = (path / "report.html").read_text(encoding="utf-8")
    assert "Execution backend" in html and "batched" in html and "Batch size" in html
