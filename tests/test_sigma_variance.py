"""Uncertainty-unit boundary and independent-pixel A9 marginal variance algebra."""
import numpy as np
import pytest
from qcbed.data_preparation import (variance_map_from_sigma_map,load_legacy_sigma_map,
    sigma_map_from_noise_model,variance_map_from_noise_model,NoiseModel,
    construct_signal_fields,SignalMode,ArcBeamCatalogue,BootstrapSite,DiskGeometry,
    ExtractionConfig,prepare_experimental_qcbed)


def test_sigma_square_shape_ownership_and_model_convergence(tmp_path):
    sigma=np.array([[0.,2.,3.],[4.,5.,6.]])
    result=variance_map_from_sigma_map(sigma)
    np.testing.assert_array_equal(result,[[0,4,9],[16,25,36]])
    assert result.shape==sigma.shape
    sigma[:]=9
    assert result[0,0]==0
    with pytest.raises(ValueError):result.setflags(write=True)
    p=tmp_path/'sigma.txt';np.savetxt(p,np.array([[2.,3.],[4.,5.]]))
    loaded=load_legacy_sigma_map(p,expected_shape=(2,2))
    np.testing.assert_array_equal(loaded,[[2,3],[4,5]])
    with pytest.raises(ValueError):loaded.setflags(write=True)
    with pytest.raises(ValueError,match='shape'):load_legacy_sigma_map(p,expected_shape=(2,3))
    image=np.array([[-9.,0,9],[16,25,36]])
    model=NoiseModel(.2,2.,3.,'counts')
    expected=.2*np.abs(image)+2*np.sqrt(abs(image))+3
    np.testing.assert_array_equal(sigma_map_from_noise_model(image,model),expected)
    np.testing.assert_array_equal(variance_map_from_noise_model(image,model),
        variance_map_from_sigma_map(expected))


@pytest.mark.parametrize('value',[np.ones(3),np.ones((2,2,2)),np.empty((0,2)),
    np.array([[np.nan]]),np.array([[np.inf]]),np.array([[-1e-30]]),
    np.array([[1e300]]),np.ones((2,2),complex),np.ones((2,2),bool)])
def test_invalid_sigma_not_hidden(value):
    with pytest.raises((ValueError,TypeError)):variance_map_from_sigma_map(value)


def extraction(mode,**extra):
    image=np.arange(256,dtype=float).reshape(16,16)+10
    site=BootstrapSite((0,0),(8.,8.),None,None,'reference','synthetic',beam_id='test')
    return prepare_experimental_qcbed(image,ArcBeamCatalogue((site,),'test',False),DiskGeometry(4,4),
        selected_beam_ids=['test'],config=ExtractionConfig.from_disk_geometry(DiskGeometry(4,4),final_binning=2),
        signal_mode=mode,**extra)


@pytest.mark.parametrize('mode',list(SignalMode))
def test_hand_calculated_crop_and_two_by_two_variance(mode):
    v=np.arange(256,dtype=float).reshape(16,16)+1
    r=extraction(mode,variance_map=v,uncertainty_provenance='declared independent detector variance')
    baseline=extraction(mode)
    if mode==SignalMode.GRADIENT:
        offsets=[(-2,0),(2,0),(0,-2),(0,2),(-2,-2),(2,2),(-2,2),(2,-2)]
        def pixel(y,x): return v[y,x]+sum(v[y+dy,x+dx] for dx,dy in offsets)/64
    else:
        def pixel(y,x): return v[y,x]
    expected=np.array([[sum(pixel(4+2*y+dy,4+2*x+dx) for dy in (0,1) for dx in (0,1))/16
                        for x in range(4)] for y in range(4)])
    np.testing.assert_array_equal(r.beams[0].variance,expected)
    for name in ('signal','support_fraction','support_weight','final_weight'):
        np.testing.assert_array_equal(getattr(r.beams[0],name),getattr(baseline.beams[0],name))
    if mode==SignalMode.GRADIENT:
        assert r.fields.variance[8,8]==pixel(8,8)
    via_sigma=extraction(mode,sigma_map=np.sqrt(v),uncertainty_provenance='synthetic sigma')
    via_squared=extraction(mode,variance_map=np.sqrt(v)**2,uncertainty_provenance='same sigma squared')
    np.testing.assert_array_equal(via_sigma.beams[0].variance,via_squared.beams[0].variance)


def test_conventional_minimum_never_subtracted_from_variance():
    r=extraction(SignalMode.CONVENTIONAL,variance_map=np.full((16,16),4.),uncertainty_provenance='synthetic')
    np.testing.assert_array_equal(r.beams[0].variance,1.) # (4+4+4+4)/16, not signal minimum
    assert not np.array_equal(r.beams[0].signal,np.full((4,4),1.))


def test_out_of_frame_placeholders_and_zero_sigma_no_floor():
    site=BootstrapSite((0,0),(0.,0.),None,None,'reference','synthetic',beam_id='edge')
    r=prepare_experimental_qcbed(np.ones((16,16)),ArcBeamCatalogue((site,),'test',False),DiskGeometry(4,4),
        selected_beam_ids=['edge'],config=ExtractionConfig(8,8,final_binning=2),
        variance_map=np.ones((16,16)),uncertainty_provenance='synthetic')
    b=r.beams[0]
    assert np.all(b.signal[:2]==0) and np.all(b.variance[:2]==1) and np.all(b.final_weight==0)
    zeros=variance_map_from_sigma_map(np.zeros((16,16)))
    assert np.count_nonzero(zeros)==0


def test_synthetic_adapter_accepts_variance_without_objective():
    from test_prepared_dataset_integration import assembled
    from qcbed.data_preparation import adapt_prepared_dataset
    _,_,data=assembled(SignalMode.CONVENTIONAL,True)
    assignments=tuple((b.id,n) for b,n in zip(data.beams,data.lattice_coordinates))
    result=adapt_prepared_dataset(data,beam_hk_assignments=assignments,
        basis_1=(1,0,0),basis_2=(0,1,0),objective_mode='masked')
    np.testing.assert_array_equal(result.experimental.variances,data.variances)
    with np.errstate(divide='raise',invalid='raise'):
        np.divide(result.experimental.weights,result.experimental.variances,
            out=np.zeros_like(data.weights),where=result.experimental.weights>0)
