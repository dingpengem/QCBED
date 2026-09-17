from dataclasses import replace
import json

import numpy as np
import pytest

from qcbed.local_runner import run_local_optimization
from qcbed.refinement_driver import run_legacy_refinement
from qcbed.refinement import pack_refinement_parameters, build_refinement_parameter_layout
from qcbed.validation_report import generate_refinement_validation_report
from test_refinement_driver import problem, config
from test_refinement import data


@pytest.mark.parametrize("method",["ANNE","SIMP","PARA","SPAR","ANSI"])
def test_observation_changes_no_results(method):
    s=problem(); e=data(); original=repr(s); p=pack_refinement_parameters(s); before=p.copy()
    cfg=config(method,maxfsimp=0,maxfpara=0)
    direct=run_legacy_refinement(s,e,cfg,objective_mode="masked")
    a=run_local_optimization(s,e,cfg,objective_mode="masked",snapshot_every=2)
    b=run_local_optimization(s,e,cfg,objective_mode="masked",record_progress=False)
    for r in (a.validation.refinement,b.validation.refinement):
        np.testing.assert_array_equal(r.final_vector,direct.final_vector)
        assert r.final_objective==direct.final_objective and r.evaluations==direct.evaluations
        for actual,expected in zip(r.stage_history,direct.stage_history):
            np.testing.assert_array_equal(actual.retained_vector,expected.retained_vector)
            # Includes stochastic proposal/adaptation state, not just final score.
            for field in actual.optimizer_result.__dataclass_fields__:
                np.testing.assert_equal(getattr(actual.optimizer_result,field),
                                        getattr(expected.optimizer_result,field))
    evaluations=[r for r in a.progress if r.kind=="evaluation"]
    assert len(evaluations)==direct.evaluations
    assert [r.evaluation_index for r in evaluations]==list(range(1,direct.evaluations+1))
    expected=list(range(2,direct.evaluations+1,2))
    assert [x.progress.evaluation_index for x in a.snapshots[1:-1]]==expected
    assert [x.label for x in b.snapshots]==["initial","final"] and b.progress==()
    assert a.snapshots[0].label=="initial" and a.snapshots[-1].label=="final"
    np.testing.assert_array_equal(p,before); assert repr(s)==original
    with pytest.raises(ValueError): evaluations[0].parameter_vector.setflags(write=True)
    with pytest.raises(ValueError): a.snapshots[0].simulated.scaled_intensity.setflags(write=True)


@pytest.mark.parametrize("cadence",[None,1,3])
def test_report_and_checkpoint_files(tmp_path,cadence):
    pytest.importorskip("matplotlib")
    s=problem(); cfg=config("ANNE",maxfsimp=0)
    folder=tmp_path/"run"
    run=run_local_optimization(s,data(),cfg,objective_mode="masked",snapshot_every=cadence,
                               output_directory=folder,reference_input_name="synthetic")
    direct=run_legacy_refinement(s,data(),cfg,objective_mode="masked")
    np.testing.assert_array_equal(run.validation.refinement.final_vector,direct.final_vector)
    assert run.validation.refinement.final_objective==direct.final_objective
    from qcbed.visualization import plot_refinement_convergence
    assert plot_refinement_convergence(direct,objective_mode="masked").axes[0].get_ylabel().startswith("masked")
    text=(folder/"report.html").read_text(encoding="utf-8")
    assert "Objective mode: masked" in text and "Non-legacy masked objective:" in text
    assert "Included zero-variance points" in text and "Snapshot gallery" in text
    assert "initial beam" in text
    for snap in run.snapshots:
        assert (folder/f"{snap.label}.npz").exists()
        assert (folder/f"{snap.label}_beam_01.png").exists()
        if snap.label not in ("initial","final"): assert snap.label in text
    saved=json.loads((folder/"summary.json").read_text())
    assert saved['final_objective']==run.validation.refinement.final_objective
    assert (folder/"final_physical_parameters.json").exists()
    with pytest.raises(FileExistsError):
        run_local_optimization(s,data(),cfg,output_directory=folder)


def test_default_legacy_and_mode_report(tmp_path):
    pytest.importorskip("matplotlib")
    s=problem(); e=data(); cfg=config("ANNE",maxfsimp=0)
    run=run_local_optimization(s,e,cfg)
    direct=run_legacy_refinement(s,e,cfg)
    assert run.objective_mode=="legacy"
    assert run.validation.refinement.final_objective==direct.final_objective
    path=generate_refinement_validation_report(tmp_path/"legacy.html",run.validation,
         config=cfg,layout=build_refinement_parameter_layout(s),reference_input_name="test",local_run=run)
    assert "Objective mode: legacy" in path.read_text()
    e=replace(e,weights=np.array([[0.],[1.],[1.],[1.]]),variances=np.array([[0.],[1.],[1.],[1.]]))
    with pytest.raises(ValueError,match="CHISMIN"):
        run_local_optimization(s,e,cfg)
    masked=run_local_optimization(s,e,cfg,objective_mode="masked")
    assert np.isfinite(masked.validation.refinement.final_objective)
    assert np.all(np.isfinite(masked.validation.comparison.numerator_by_beam))


@pytest.mark.parametrize("cadence",[0,-1,1.5,True])
def test_invalid_cadence(cadence):
    with pytest.raises(ValueError,match="snapshot_every"):
        run_local_optimization(problem(),data(),config(),snapshot_every=cadence)
