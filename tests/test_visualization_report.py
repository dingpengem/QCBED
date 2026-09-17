import numpy as np
import pytest

pytest.importorskip("matplotlib")

from qcbed.visualization import plot_cbed_comparison, plot_refinement_convergence, plot_parameter_trajectories
from qcbed.validation_report import generate_refinement_validation_report
from qcbed.validation import validate_refinement
from qcbed.refinement import build_refinement_parameter_layout
from test_output_validation import comparison
from test_refinement_driver import problem,config
from test_refinement import data


def test_plot_orientation_and_common_scale_bad_pattern():
    from dataclasses import replace
    c=comparison()
    bad=replace(c.final_simulated,raw_intensity=np.full((6,2),1000.),objective_intensity=np.full((6,2),1000.))
    c=replace(c,final_simulated=bad)
    fig=plot_cbed_comparison(c,beam_index=1)
    axes=fig.axes[:4]
    images=[a.images[0] for a in axes]
    assert images[0].get_clim()==images[1].get_clim()==images[2].get_clim()
    assert images[2].get_clim()[1]==2000
    assert np.max(images[0].get_array())==12 and np.min(images[2].get_array())==2000
    assert all(im.origin=="lower" for im in images)
    np.testing.assert_array_equal(images[0].get_array(),[[2,4,6],[8,10,12]])
    np.testing.assert_array_equal(images[3].get_array(),c.residual[:,1].reshape(2,3))
    lo,hi=images[3].get_clim()
    assert lo==-hi and hi>1900
    assert "(0, 2)" in fig._suptitle.get_text()
    assert "Δ(H,K)" in axes[0].get_xlabel()


def test_report_all_beams_and_stage_level_figures(tmp_path):
    s=problem()
    controls=config("ANSI",maxfsimp=0)
    validation=validate_refinement(s,data(),controls)
    layout=build_refinement_parameter_layout(s)
    r=validation.refinement
    fig=plot_refinement_convergence(r)
    assert len(fig.axes[0].lines[0].get_ydata())==len(r.stage_history)+1
    assert "not individual evaluations" in fig.axes[0].get_xlabel()
    p=plot_parameter_trajectories(r,layout)
    np.testing.assert_array_equal(p.axes[0].lines[0].get_ydata(),
                                 [r.initial_vector[0]]+[h.retained_vector[0] for h in r.stage_history])
    path=generate_refinement_validation_report(tmp_path/"report.html",validation,layout=layout,
             config=controls,reference_input_name='<unsafe & input>')
    html=path.read_text(encoding="utf-8")
    assert '&lt;unsafe &amp; input&gt;' in html and 'data:image/png;base64,' in html
    assert "Initial YNORM" in html and "Numerator" in html and "ANNEAL" in html and "AMOEBA" in html
    assert "Display masking is separate from legacy objective arithmetic" in html
    assert "n_zero_variance_and_zero_weight" in html and "RMS: valid display" in html
    assert "weight, zero variance and overlap maps" in html
    for beam in validation.comparison.experimental.beam_indices_hk:
        assert f"Beam {beam}" in html
    # Separate two-beam fixture verifies no beam is dropped from the gallery.
    from qcbed.validation import RefinementValidation
    multiple=RefinementValidation(r,comparison())
    html=generate_refinement_validation_report(tmp_path/"multi.html",multiple,layout=layout,config=controls,
                                              reference_input_name="synthetic").read_text(encoding="utf-8")
    assert 'Beam (1, -1)' in html and 'Beam (0, 2)' in html
