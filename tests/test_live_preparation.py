"""Live orchestration uses existing kernels; no historical geometry is loaded."""
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
import json
import numpy as np
import pytest
from qcbed.data_preparation.live import LiveSelection,resolve_geometry,prepare_signal,prepare_from_raw_cbed
from qcbed.data_preparation.models import DiskGeometry
from qcbed.gui.preparation_model import publish_preparation,new_preparation_document,load_image,save_raw_project


def scene():
    y,x=np.indices((180,180));a=np.full((180,180),10.)
    for cx,cy in [(90,90),(140,90),(90,140)]:
        a+=100/(1+np.exp(np.clip((np.hypot(x-cx,y-cy)-20)*3,-700,700)))
    return a+30*np.exp(-((x-90)**2+(y-90)**2)/10)


def selection():return LiveSelection(DiskGeometry(20.,20.),(87.,87.,93.,93.),40,40,symmetry_crop_half_size_px=(8,8))


def test_geometry_precedence_and_measured_aperture_ellipse():
    y,x=np.indices((80,120));image=100.*(((x-60)/12)**2+((y-40)/8)**2<1)
    g,source=resolve_geometry(aperture=image)
    assert g==DiskGeometry(11.5,7.5) and source=='aperture threshold geometry'
    override=DiskGeometry(3.25,4.5)
    assert resolve_geometry(geometry=override,aperture=image)[0]==override
    assert resolve_geometry(disk_profile=image,aperture=np.ones((1,1)))[0]==g
    with pytest.raises(ValueError,match='Draw and confirm'):resolve_geometry()
    with pytest.raises(ValueError):resolve_geometry(aperture=np.ones((80,120)))


def test_no_aperture_exact_poisson_fallback_no_correction(monkeypatch):
    import qcbed.data_preparation.live as live
    monkeypatch.setattr(live,'calibrate_psf_from_aperture',lambda *a,**k:pytest.fail('No calibration permitted'))
    image=np.array([[-1.,0.,4.],[9.,16.,25.]])
    corrected,sigma,variance,mask,g,meta=prepare_signal(image,geometry=DiskGeometry(2.,2.))
    np.testing.assert_array_equal(corrected,image)
    np.testing.assert_array_equal(variance,np.maximum(image,0))
    np.testing.assert_array_equal(sigma,np.sqrt(np.maximum(image,0)))
    np.testing.assert_array_equal(mask,image>0)
    assert meta['noise_source']=='Poisson / intensity fallback' and meta['correction']=='none'


def test_calibrated_branch_uses_distinct_backend_roles(monkeypatch):
    import qcbed.data_preparation.live as live
    from qcbed.data_preparation.models import NoiseModel
    calls=[];image=scene();aperture=np.zeros_like(image);aperture[80:100,70:110]=100
    # Orchestration spy: underlying PSF/noise numerical equivalence has its own tests.
    response=SimpleNamespace(provenance=('calibration record',))
    monkeypatch.setattr(live,'calibrate_psf_from_aperture',lambda a,**k:(calls.append('calibrate') or response))
    monkeypatch.setattr(live,'correct_cbed_with_otf',lambda a,r:(calls.append('correct') or a.copy()))
    from qcbed.data_preparation.noise import NoiseFitConfig
    fit=SimpleNamespace(model=NoiseModel(1.,0.,0.,'input intensity units'),noise_mask=np.ones(image.shape,bool),config=NoiseFitConfig())
    monkeypatch.setattr(live,'characterize_noise',lambda a,**k:(calls.append('noise') or fit))
    r=prepare_from_raw_cbed(image,selection(),calibration_aperture=aperture)
    assert calls==['calibrate','correct','noise']
    assert r.geometry==selection().geometry  # explicit override wins; calibration geometry never overwrites it
    assert r.provenance['noise_source']=='calibrated' and r.dataset.variances is not None


@pytest.fixture(scope='module')
def live_result():return prepare_from_raw_cbed(scene(),selection())


def test_live_catalogue_origin_extraction_and_input_immutable(live_result):
    r=live_result
    assert len(r.dataset.beams)==9
    assert r.symmetry.reference_kind=='direct' and r.symmetry.trial_count==36
    assert r.dataset.reciprocal_frame.zone_axis_origin_px==r.symmetry.symmetry_center_px
    assert r.dataset.scan_geometry.scan_rows==40
    np.testing.assert_array_equal(r.corrected_pattern,scene())
    assert set(np.unique(r.dataset.weights))<={0.,1.}


def test_publish_validate_mask_and_replay_without_B(tmp_path,live_result):
    from qcbed.gui.data_model import read_data,edit_mask,select_zone_axis
    from qcbed.gui.workflow import validate_data
    path=publish_preparation(live_result,tmp_path/'prep',{})
    c=new_preparation_document(path)
    with pytest.raises(ValueError,match='explicitly confirm'):validate_data(c,tmp_path)
    c=select_zone_axis(c,tmp_path,beam=0)
    assert 'Poisson' in validate_data(c,tmp_path)
    before=read_data(c,tmp_path)
    masked=edit_mask(c,tmp_path,beam=0,box=(15,15,16,16))
    after=read_data(masked,tmp_path)
    np.testing.assert_array_equal(before.intensities,after.intensities)
    np.testing.assert_array_equal(before.variances,after.variances)
    assert after.weights[15*40+15,0]==0
    assert Path(path).is_file()


def test_raw_document_no_prepared_input(tmp_path):
    from qcbed.yaml_config import read_yaml_document
    path=tmp_path/'raw.yml'
    save_raw_project({'schema_version':2,'data_prep':{'pattern':'p.dat'}},path,base_directory=tmp_path)
    assert read_yaml_document(path)['data_prep']=={'pattern':str(tmp_path/'p.dat')}
    with pytest.raises(FileExistsError):save_raw_project({'data_prep':{'pattern':'p.dat'}},path,base_directory=tmp_path)


@pytest.fixture(scope='module')
def real_live_result():
    pattern=Path(__file__).resolve().parents[1]/'reference_data/pattern.dat'
    settings=replace(selection(),geometry=DiskGeometry(38.,38.),symmetry_roi=(111.,105.,155.,151.),selection_width=76,selection_height=76,final_binning=2,symmetry_crop_half_size_px=(19,19))
    return prepare_from_raw_cbed(load_image(pattern),settings)


def test_real_pattern_live_arc_recovery(real_live_result):
    r=real_live_result
    assert r.geometry==DiskGeometry(38.,38.)
    d=r.provenance['lattice_detection']
    assert d['candidate_count']==3
    assert d['seed']['reference_center_px']==(127.,126.)
    assert d['seed']['basis_vector_1_px']==(81.,9.)
    assert d['seed']['basis_vector_2_px']==(-9.,82.)
    assert 'insufficient independent observations' in r.provenance['lattice_status']
    extracted=[b for b in r.provenance['extraction'] if b['status']=='extracted']
    assert {b['indices'] for b in extracted}=={(x,y) for x in (-1,0,1) for y in (-1,0,1)}
    assert all(b['coverage']==1. for b in extracted)
    assert r.dataset.intensities.shape==(38*38,9)
    assert any(s['status']=='missing' for s in r.provenance['lattice_sites'])


def test_sparse_arc_selection_is_order_independent():
    from itertools import permutations
    from qcbed.data_preparation.live import _live_arc_seed
    from qcbed.data_preparation.lattice_detection import detect_lattice_candidates
    image=load_image(Path(__file__).resolve().parents[1]/'reference_data/pattern.dat')
    geometry=DiskGeometry(38.,38.)
    candidates=detect_lattice_candidates(image,geometry,(127.,126.))
    seeds=[_live_arc_seed(p,(127.,126.),geometry) for p in permutations(candidates)]
    assert all(s==seeds[0] for s in seeds)
    # An unrelated third center must not be silently discarded to force a seed.
    bad=replace(candidates[0],center_px=(160.,170.))
    with pytest.raises(ValueError,match='consistent non-collinear'):
        _live_arc_seed((bad,*candidates[1:]),(127.,126.),geometry)


def test_aperture_calibration_failure_is_not_poisson_fallback(monkeypatch):
    import qcbed.data_preparation.live as live
    def failed(*a,**k):raise ValueError('PSF calibration failed')
    monkeypatch.setattr(live,'calibrate_psf_from_aperture',failed)
    with pytest.raises(ValueError,match='PSF calibration failed'):
        prepare_signal(scene(),calibration_aperture=scene(),geometry=DiskGeometry(20.,20.))


def test_profile_geometry_without_calibration_remains_poisson():
    image=np.zeros((20,30));image[4:12,5:17]=50
    *_,g,meta=prepare_signal(np.ones((24,32)),disk_profile=image)
    assert g==DiskGeometry(6.,4.)
    assert meta['noise_source']=='Poisson / intensity fallback'


@pytest.mark.parametrize('mode',['conventional','gradient'])
def test_recorded_mask_layers_reconstruct_objective_weight(mode):
    from qcbed.data_preparation.models import SignalMode
    result=prepare_from_raw_cbed(scene(),replace(selection(),signal_mode=SignalMode(mode)))
    for i,beam in enumerate(result.dataset.beams):
        layers=result.provenance['mask_layers'][beam.id]
        geometry=np.asarray(layers['geometric'])
        processing=np.asarray(layers['processing'])
        assert set(np.unique(geometry)) <= {0,1}
        assert set(np.unique(processing)) <= {0,1}
        np.testing.assert_array_equal((geometry*processing).ravel(),result.dataset.weights[:,i])
        if mode=='conventional':
            # Strictly positive synthetic detector: no intensity-validity holes.
            np.testing.assert_array_equal(geometry.ravel(),result.dataset.weights[:,i])
