from dataclasses import replace
import json
from pathlib import Path

import numpy as np
import pytest

from qcbed import ForwardExecutionConfig
from qcbed.full_runner import FullRefinementConfig, run_full_refinement, load_refinement_checkpoint
from qcbed.refinement_driver import run_legacy_refinement, LegacyRefinementRestart
from qcbed.refinement import pack_refinement_parameters
from qcbed.output_io import parse_legacy_out_dat
from test_refinement_driver import problem, config
from test_refinement import data


def settings(**kwargs):
    return FullRefinementConfig(config("ANSI", maxfsimp=0,
        execution=ForwardExecutionConfig("batched", 3)), snapshot_every=None,
        checkpoint_every=1, max_total_evaluations=100, **kwargs)


@pytest.fixture
def cheap_render(monkeypatch):
    # Most lifecycle tests exercise real numerical calls/I/O without repeated PNGs.
    import qcbed.full_runner as runner
    monkeypatch.setattr(runner, "_preview", lambda *a, **k: None)
    import qcbed.validation_report as report
    monkeypatch.setattr(report, "generate_refinement_validation_report",
        lambda path, *a, **k: Path(path).write_text("<html></html>"))
    from matplotlib.figure import Figure
    monkeypatch.setattr(Figure, "savefig", lambda self, path, **kw: Path(path).write_bytes(b"test"))


@pytest.mark.parametrize("changes", [dict(checkpoint_every=0), dict(snapshot_every=True),
    dict(max_total_evaluations=1), dict(max_wall_time=float("nan")), dict(objective_mode="other")])
def test_invalid_config(changes):
    with pytest.raises((ValueError, TypeError)):
        replace(settings(), **changes)


def test_complete_outputs_and_exact_driver(cheap_render, tmp_path):
    s, e, c = problem(), data(), settings()
    before = repr(s)
    direct = run_legacy_refinement(s, e, c.refinement, objective_mode="masked")
    folder = tmp_path / "run"
    seen = []
    def observe(row):
        if row["classification"] == "trial":
            assert (folder / "initial_snapshot.npz").exists()
            assert load_refinement_checkpoint(folder)["retained_objective"] is not None
        seen.append(row)
    result = run_full_refinement(s, e, c, output_directory=folder, progress_hook=observe)
    assert result.status == "completed"
    assert result.summary["final_objective"] == direct.final_objective
    np.testing.assert_array_equal(np.load(folder / "final_active_vector.npy"), direct.final_vector)
    assert repr(s) == before
    for name in ("run_config.json", "summary.json", "progress.jsonl", "initial_parameters.json",
        "final_parameters.json", "initial_active_vector.npy", "final_active_vector.npy",
        "initial_simulated_intensity.npy", "final_simulated_intensity.npy", "final_scaled_intensity.npy",
        "initial_snapshot.npz", "final_snapshot.npz", "final_out.dat", "convergence.png",
        "parameter_trajectories.png", "report.html"):
        assert (folder / name).exists(), name
    rows = [json.loads(x) for x in (folder / "progress.jsonl").read_text().splitlines()]
    assert rows == seen
    assert all(x["accepted"] is None for x in rows)
    assert {x["classification"] for x in rows} >= {"trial", "stage_retained", "initial_retained"}
    out = parse_legacy_out_dat(folder / "final_out.dat", scan=s.scan, beam_indices_hk=e.beam_indices_hk)
    np.testing.assert_array_equal(out.scaled_intensity, np.load(folder / "final_scaled_intensity.npy"))
    cp = load_refinement_checkpoint(folder)
    assert cp["completed"] and cp["rng_state"] is None and cp["evaluation_count"] == result.summary["evaluations"]
    assert cp["execution_backend"] == "batched" and cp["objective_mode"] == "masked"
    assert cp["retained_ynorm"] == result.summary["final_ynorm"]


def test_interrupt_and_stage_resume(cheap_render, tmp_path):
    s, e, c = problem(), data(), settings()
    folder = tmp_path / "run"
    def stop(row):
        if row["classification"] == "stage_retained":
            raise KeyboardInterrupt
    first = run_full_refinement(s, e, c, output_directory=folder, progress_hook=stop)
    assert first.status == "interrupted"
    cp = load_refinement_checkpoint(folder)
    assert len(cp["stages"]) == 1 and cp["stage_position"] == 1
    old_log = (folder / "progress.jsonl").read_text()
    final = run_full_refinement(s, e, c, output_directory=folder, resume=True)
    direct = run_legacy_refinement(s, e, c.refinement, objective_mode="masked")
    assert final.status == "completed" and final.summary["restart_count"] == 1
    assert final.summary["final_objective"] == direct.final_objective
    np.testing.assert_array_equal(np.load(folder / "final_active_vector.npy"), direct.final_vector)
    assert (folder / "progress.jsonl").read_text().startswith(old_log)
    assert len(load_refinement_checkpoint(folder)["stages"]) == 2


def test_midstage_trial_is_not_retained(cheap_render, tmp_path):
    c = settings(); s = problem(); e = data(); folder = tmp_path / "run"
    def stop(row):
        if row["classification"] == "trial" and row["stage_local_index"] >= 2:
            raise KeyboardInterrupt
    result = run_full_refinement(s, e, c, output_directory=folder, progress_hook=stop)
    assert result.status == "interrupted"
    cp = load_refinement_checkpoint(folder)
    assert cp["stage_position"] == 0 and cp["stages"] == []
    np.testing.assert_array_equal(cp["retained_vector"], pack_refinement_parameters(s))
    final = run_full_refinement(s, e, c, output_directory=folder, resume=True)
    assert final.status == "completed"


def test_budgets_and_resume_identity(cheap_render, tmp_path):
    c = replace(settings(), max_total_evaluations=3)
    folder = tmp_path / "run"
    limited = run_full_refinement(problem(), data(), c, output_directory=folder)
    assert limited.status == "limited" and limited.summary["evaluations"] == 3
    assert limited.summary["final_verification"]
    with pytest.raises(ValueError, match="fingerprint"):
        run_full_refinement(problem(), replace(data(), intensities=data().intensities*2), c,
                            output_directory=folder, resume=True)
    result = run_full_refinement(problem(), data(), replace(c, max_total_evaluations=100),
        output_directory=folder, resume=True)
    assert result.status == "completed"


def test_snapshot_cadence(cheap_render, tmp_path):
    c = replace(settings(), snapshot_every=2)
    run_full_refinement(problem(), data(), c, output_directory=tmp_path / "run")
    assert list((tmp_path / "run").glob("retained_*_snapshot.npz"))


def test_real_preview_and_report(tmp_path):
    pytest.importorskip("matplotlib")
    c = replace(settings(), refinement=replace(settings().refinement, loops=0))
    folder = tmp_path / "run"
    result = run_full_refinement(problem(), data(), c, output_directory=folder)
    assert result.status == "completed"
    assert (folder / "initial_preview" / "beam_01.png").exists()
    assert (folder / "final_preview" / "index.html").exists()
    text = (folder / "report.html").read_text(encoding="utf-8")
    assert "before YNORM computation" in text and "Production run / restart summary" in text


def test_restart_validation():
    with pytest.raises(ValueError):
        LegacyRefinementRestart(0, 0, np.ones(1))
    with pytest.raises(ValueError, match="outside"):
        run_legacy_refinement(problem(), data(), settings().refinement,
            restart=LegacyRefinementRestart(4, 0, np.ones(1)))


@pytest.mark.parametrize("method", ["ANNE", "SIMP", "PARA", "SPAR"])
def test_instrumentation_preserves_standalone_decisions(cheap_render, tmp_path, method):
    c = replace(settings(), refinement=replace(settings().refinement, method=method, maxfpara=0))
    direct = run_legacy_refinement(problem(), data(), c.refinement, objective_mode="masked")
    a = run_full_refinement(problem(), data(), c, output_directory=tmp_path / "a")
    b = run_full_refinement(problem(), data(), c, output_directory=tmp_path / "b")
    assert a.summary["final_objective"] == b.summary["final_objective"] == direct.final_objective
    assert a.summary["evaluations"] == b.summary["evaluations"]
    for folder in (tmp_path / "a", tmp_path / "b"):
        np.testing.assert_array_equal(np.load(folder / "final_active_vector.npy"), direct.final_vector)


def test_final_verification_failure_is_explicit(cheap_render, tmp_path, monkeypatch):
    import qcbed.full_runner as module
    original = module.evaluate_refinement_trial
    count = 0
    def evaluate(*args, **kwargs):
        nonlocal count
        count += 1
        value = original(*args, **kwargs)
        return replace(value, scaling=value.scaling+1) if count == 2 else value
    monkeypatch.setattr(module, "evaluate_refinement_trial", evaluate)
    folder = tmp_path / "run"
    with pytest.raises(ValueError, match="fresh retained"):
        run_full_refinement(problem(), data(), settings(), output_directory=folder)
    assert load_refinement_checkpoint(folder)["status"] == "failed"
    assert not json.loads((folder / "summary.json").read_text())["final_verification"]


def test_checkpoint_version_and_snapshot_loop(cheap_render, tmp_path):
    folder = tmp_path / "run"
    run_full_refinement(problem(), data(), settings(), output_directory=folder)
    meta = json.loads((folder / "final_snapshot.json").read_text())
    assert meta["loop_index"] == 1 and meta["stage"] == "AMOEBA"
    cp = load_refinement_checkpoint(folder)
    cp["version"] = 99
    (folder / "checkpoint_latest.json").write_text(json.dumps(cp))
    with pytest.raises(ValueError, match="version"):
        load_refinement_checkpoint(folder)


def test_pre_stage_loop_seed_preserved_on_restart(cheap_render, tmp_path):
    c = replace(settings(), refinement=replace(settings().refinement, loops=2))
    folder = tmp_path / "run"
    def stop(row):
        if row["classification"] == "stage_retained" and row["stage"] == "AMOEBA":
            raise KeyboardInterrupt
    run_full_refinement(problem(), data(), c, output_directory=folder, progress_hook=stop)
    assert load_refinement_checkpoint(folder)["loop_index"] == 2
    final = run_full_refinement(problem(), data(), c, output_directory=folder, resume=True)
    direct = run_legacy_refinement(problem(), data(), c.refinement, objective_mode="masked")
    assert final.summary["final_objective"] == direct.final_objective
    assert [h["seed"] for h in load_refinement_checkpoint(folder)["stages"]] == [h.seed for h in direct.stage_history]


def test_interrupt_before_initial_evaluation_can_resume(cheap_render, tmp_path, monkeypatch):
    import qcbed.full_runner as module
    original = module.evaluate_refinement_trial
    def interrupted(*a, **k):
        raise KeyboardInterrupt
    monkeypatch.setattr(module, "evaluate_refinement_trial", interrupted)
    folder = tmp_path / "run"
    result = run_full_refinement(problem(), data(), settings(), output_directory=folder)
    assert result.status == "interrupted"
    assert load_refinement_checkpoint(folder)["retained_objective"] is None
    assert (folder / "progress.jsonl").exists()
    monkeypatch.setattr(module, "evaluate_refinement_trial", original)
    assert run_full_refinement(problem(), data(), settings(), output_directory=folder, resume=True).status == "completed"


def test_resume_initial_preview_uses_mapped_scan(cheap_render, tmp_path, monkeypatch):
    import qcbed.full_runner as module
    s = replace(problem(), scan_masks=("YN", "NN", "NN"), ranges=((1.5, .5), (1.5, .5)))
    c = replace(settings(), refinement=replace(settings().refinement, loops=0, plamda=(.2, .2)))
    corners = []
    monkeypatch.setattr(module, "_preview", lambda path, e, initial, *args: corners.append(initial.scan.corner_a[0]))
    folder = tmp_path / "run"
    run_full_refinement(s, data(), c, output_directory=folder)
    run_full_refinement(s, data(), c, output_directory=folder, resume=True)
    # ARRANGEPARS floors the selected zero corner to 1e-10; the immutable
    # baseline still contains zero and must not supply resumed preview geometry.
    assert corners == [1e-10, 1e-10, 1e-10]
    assert s.scan.corner_a[0] == 0
