from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from qcbed.display_mask import valid_display_mask, experimental_mask_statistics
from qcbed.experimental import parse_experimental_intensities
from qcbed.refinement_driver import parse_legacy_refinement
from qcbed.refinement_io import _CaptureRecords
from qcbed.legacy_io import _parse_for007_records
from qcbed.visualization import plot_cbed_comparison, plot_weight_diagnostics
from test_output_validation import comparison


def masked_comparison():
    c=comparison()
    w=np.array([[0,1],[.25,0],[1,1],[0,.5],[-.5,1],[0,1]],dtype=np.float64)
    return replace(c,experimental=replace(c.experimental,weights=w))


def test_exact_mask_fractional_weights_and_immutable_data():
    c=masked_comparison()
    expected=np.array([[False,True],[True,False],[True,True],
                       [False,True],[True,True],[False,True]])
    np.testing.assert_array_equal(valid_display_mask(c.experimental),expected)
    np.testing.assert_array_equal(c.valid_display,expected)
    stats=experimental_mask_statistics(c.experimental)
    assert stats['n_zero_weight']==4 and not stats['binary_weights']
    assert stats['n_negative_weights']==1 and stats['n_nonfinite_weights']==0
    assert stats['unique_weights']==[-.5,0,.25,.5,1]
    with pytest.raises(ValueError): c.valid_display.setflags(write=True)
    for b in range(2):
        assert c.rms_by_beam[b]==pytest.approx(np.sqrt(np.mean(c.residual[expected[:,b],b]**2)))
        assert c.rms_masked_by_beam[b]==pytest.approx(np.sqrt(np.mean(c.residual[~expected[:,b],b]**2)))
        assert c.rms_all_by_beam[b]==pytest.approx(np.sqrt(np.mean(c.residual[:,b]**2)))


def test_shared_mask_reconstruction_and_raw_mode():
    pytest.importorskip('matplotlib')
    c=masked_comparison()
    before=c.experimental.weights.copy(),c.residual.copy(),c.final_simulated.scaled_intensity.copy()
    images=[ax.images[0] for ax in plot_cbed_comparison(c).axes[:4]]
    expected=[[True,False,False],[True,False,True]]
    for image in images:
        np.testing.assert_array_equal(np.ma.getmaskarray(image.get_array()),expected)
        assert np.all(np.isnan(image.get_array().filled(np.nan)[expected]))
        assert image.cmap.get_bad()[3]==1  # explicitly visible grey mask
    np.testing.assert_array_equal(images[3].get_array().compressed(),c.residual[[1,2,4],0])
    assert images[0].get_clim()==images[1].get_clim()==images[2].get_clim()
    raw=plot_cbed_comparison(c,mask_zero_weight=False)
    for ax in raw.axes[:4]: assert not np.any(np.ma.getmaskarray(ax.images[0].get_array()))
    np.testing.assert_array_equal(raw.axes[3].images[0].get_array(),c.residual[:,0].reshape(2,3))
    for original,current in zip(before,(c.experimental.weights,c.residual,c.final_simulated.scaled_intensity)):
        np.testing.assert_array_equal(original,current)


def test_all_masked_is_explicit_and_does_not_repair_contributions():
    pytest.importorskip('matplotlib')
    c=comparison()
    c=replace(c,experimental=replace(c.experimental,weights=np.zeros((6,2)),variances=np.zeros((6,2))))
    assert np.all(np.isnan(c.rms_by_beam)) and np.all(np.isfinite(c.rms_masked_by_beam))
    assert np.all(np.isnan(c.numerator_by_beam)) and np.all(np.isnan(c.denominator_by_beam))
    fig=plot_cbed_comparison(c)
    assert all(np.all(np.ma.getmaskarray(ax.images[0].get_array())) for ax in fig.axes[:4])


def test_masked_outlier_does_not_control_visible_scale():
    pytest.importorskip('matplotlib')
    c=masked_comparison()
    outlier=c.final_simulated.objective_intensity.copy(); outlier[0,0]=1e12
    c=replace(c,final_simulated=replace(c.final_simulated,objective_intensity=outlier))
    masked=plot_cbed_comparison(c)
    raw=plot_cbed_comparison(c,mask_zero_weight=False)
    assert masked.axes[2].images[0].get_clim()[1]<100
    assert raw.axes[2].images[0].get_clim()[1]==2e12
    assert c.final_simulated.scaled_intensity[0,0]==2e12


def test_visible_zero_variance_location_and_diagnostic_maps():
    pytest.importorskip('matplotlib')
    c=masked_comparison()
    v=c.experimental.variances.copy(); v[1,0]=0; v[0,0]=0
    e=replace(c.experimental,variances=v)
    stat=experimental_mask_statistics(e)
    assert stat['visible_zero_variance']==[dict(beam_index=0,beam_hk=(1,-1),tilt_index=1,
        row=0,column=1,intensity=3.,variance=0.,weight=.25)]
    fig=plot_weight_diagnostics(e,c.final_simulated.scan)
    np.testing.assert_array_equal(fig.axes[0].images[0].get_array(),[[0,.25,1],[0,-.5,0]])
    np.testing.assert_array_equal(fig.axes[1].images[0].get_array(),[[1,1,0],[0,0,0]])
    np.testing.assert_array_equal(fig.axes[2].images[0].get_array(),[[1,0,0],[0,0,0]])


def test_reference_counts_and_spatial_overlap():
    p=Path(__file__).parents[1]/'reference_data'/'FOR007.DAT'
    state,_=parse_legacy_refinement(p,seed=0)
    records=_CaptureRecords(p); _parse_for007_records(records)
    row=records.values['output beam count and indices']
    beams=tuple(zip(row[1::2],row[2::2]))
    e=parse_experimental_intensities(p.parent/'in.dat',beam_indices_hk=beams,scan=state.scan)
    s=experimental_mask_statistics(e)
    assert [s[k] for k in ('n_total','n_zero_variance','n_zero_weight',
        'n_zero_variance_and_zero_weight','n_zero_variance_and_nonzero_weight',
        'n_nonzero_variance_and_zero_weight')]==[13689,117,6140,117,0,6023]
    assert s['unique_weights']==[0.,1.] and s['binary_weights']
    assert s['n_negative_weights']==s['n_nonfinite_weights']==0
    assert not s['same_mask_all_beams'] and s['visible_zero_variance']==[]
    assert [b['n_zero_weight'] for b in s['per_beam']]==[682,685,681,683,684,682,681,681,681]
    assert [b['n_zero_variance'] for b in s['per_beam']]==[0,0,0,0,0,0,78,39,0]
    expected=np.zeros((39,39,9),dtype=bool)
    expected[:2,:,6]=True; expected[:,0,7]=True
    np.testing.assert_array_equal((e.variances==0).reshape(39,39,9),expected)
    assert np.all(e.weights[e.variances==0]==0)
