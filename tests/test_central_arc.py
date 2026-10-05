"""Independent central geometry validation; canonical A5 remains untouched."""
from dataclasses import FrozenInstanceError
import numpy as np
import pytest
from qcbed.data_preparation import (DiskGeometry, CentralArcConfig, LatticeDetectionConfig,
                                    compare_central_disk_centers)
from qcbed.data_preparation.legacy_disk_localization import locate_central_disk_intensity as locate_central_disk
from qcbed.data_preparation.lattice_detection import arc_edge_support
from qcbed.data_preparation.central_arc import write_central_center_comparison


def scene(kind='clean',radii=(24.,24.)):
    y,x=np.indices((256,256));dx,dy=x-128.,y-128.
    r=np.hypot(dx/radii[0],dy/radii[1])
    profile=1/(1+np.exp(np.clip((r-1)*40,-700,700)))
    image=50*profile
    if kind=='asymmetric': image*=1+15/(1+np.exp(np.clip(-dx/4,-700,700)))
    if kind=='fringe': image+=500*np.exp(-dy**2/2)*np.exp(-(dx/14)**8)
    if kind=='partial':
        theta=np.arctan2(dy/radii[1],dx/radii[0]);local=(theta+np.pi/8)%(np.pi/4)-np.pi/8
        image*=1/(1+np.exp(np.clip((abs(local)-.16)*200,-700,700)))
    if kind=='wrong': image+=800/(1+np.exp(np.clip((np.hypot(dx,dy)-8)*3,-700,700)))
    return image,DiskGeometry(*radii)


@pytest.mark.parametrize('kind',['clean','asymmetric','fringe','partial','wrong'])
def test_true_boundary_center_despite_interior_and_missing_arcs(kind):
    image,g=scene(kind);before=image.copy()
    a5=locate_central_disk(image,g)
    result=compare_central_disk_centers(image,g)
    assert result.arc_center_px is not None
    np.testing.assert_allclose(result.arc_center_px,(128,128),atol=1.5)
    assert result.intensity_center_px==a5.central_disk_center_px
    assert result.canonical_center_method=='known_radius_partial_arcs'
    assert result.difference_vector_px==tuple(np.array(result.arc_center_px)-result.intensity_center_px)
    assert result.difference_magnitude_px==pytest.approx(np.linalg.norm(result.difference_vector_px))
    assert result.arc_support.supporting_edge_pixels>=12
    assert result.arc_support.radius_residual_rms_px<=2
    np.testing.assert_array_equal(image,before)
    assert result==compare_central_disk_centers(image,g)
    with pytest.raises(FrozenInstanceError): result.arc_center_px=(0,0)


def test_elliptical_geometry_and_polarity():
    image,g=scene(radii=(30,20))
    a=compare_central_disk_centers(image,g,intensity_center_px=(126,129))
    b=compare_central_disk_centers(-image,g,intensity_center_px=(126,129))
    np.testing.assert_allclose(a.arc_center_px,(128,128),atol=1)
    assert a.arc_center_px==b.arc_center_px
    assert a.arc_support.supporting_pixels==b.arc_support.supporting_pixels


def test_wrong_radius_alone_not_returned_as_central_boundary():
    y,x=np.indices((256,256));image=800/(1+np.exp(np.clip((np.hypot(x-128,y-128)-8)*3,-700,700)))
    result=compare_central_disk_centers(image,DiskGeometry(24,24),intensity_center_px=(128,128))
    assert result.status=='no_validated_arc'
    assert result.arc_center_px is None and result.difference_magnitude_px is None


def test_default_edge_helper_retains_a6_central_exclusion():
    image,g=scene();seed=(128,128)
    a,ad=arc_edge_support(image,g,seed)
    b,bd=arc_edge_support(image,g,seed,exclude_central=True)
    for key in a: np.testing.assert_array_equal(a[key],b[key])
    assert ad==bd
    y,x=np.indices(image.shape)
    assert not a['selected_edges'][np.hypot(x-128,y-128)<=24*1.2].any()
    central,_=arc_edge_support(image,g,seed,exclude_central=False)
    assert central['selected_edges'][np.hypot(x-128,y-128)<=24*1.2].any()


def test_diagnostics_and_no_overwrite(tmp_path):
    import json
    image,g=scene();result=compare_central_disk_centers(image,g)
    output=tmp_path/'comparison'
    path=write_central_center_comparison(output,image,result)
    assert path.is_file()
    data=json.loads((output/'central_disk_center_comparison.json').read_text())
    assert data['canonical_center_method']=='known_radius_partial_arcs'
    assert data['arc_center_px']==list(result.arc_center_px)
    with pytest.raises(FileExistsError): write_central_center_comparison(output,image,result)


@pytest.mark.parametrize('kwargs',[dict(search_radius_px=0),dict(search_step_px=True),dict(search_step_px=9),dict(search_step_px=.001)])
def test_invalid_search_config(kwargs):
    with pytest.raises(ValueError): CentralArcConfig(**kwargs)


def test_search_remains_bounded_and_marks_boundary_optimum():
    image,g=scene()
    result=compare_central_disk_centers(image,g,intensity_center_px=(126,128),
        search_config=CentralArcConfig(search_radius_px=1,search_step_px=.5))
    assert result.intensity_center_px==(126,128)
    assert result.status=='search_boundary'
    assert result.difference_magnitude_px<=1
    assert result.arc_center_px!=(128,128)


def test_asymmetric_interior_geometric_error_smaller_than_intensity_error():
    image,g=scene('asymmetric')
    result=compare_central_disk_centers(image,g)
    intensity_error=np.linalg.norm(np.array(result.intensity_center_px)-(128,128))
    arc_error=np.linalg.norm(np.array(result.arc_center_px)-(128,128))
    assert intensity_error>=2
    assert arc_error<=1.5 and arc_error<intensity_error


@pytest.mark.parametrize('kind',['clean','asymmetric','fringe','partial','wrong','noise','missing'])
def test_canonical_geometry_without_intensity_initializer(kind,monkeypatch):
    from qcbed.data_preparation import locate_central_disk as geometric
    import qcbed.data_preparation.legacy_disk_localization as legacy
    def forbidden(*a,**kw): raise AssertionError('Legacy intensity localization called')
    monkeypatch.setattr(legacy,'locate_central_disk_intensity',forbidden)
    image,g=scene(kind)
    if kind=='noise': image+=np.random.default_rng(21).normal(0,.8,image.shape)+np.indices(image.shape)[1]*.01
    if kind=='missing': image[:128,:128]=0
    before=image.copy()
    result=geometric(image,g)
    assert result.method=='known_radius_partial_arcs'
    assert result.search_prior_px==(127.5,127.5)
    np.testing.assert_allclose(result.central_disk_center_px,(128,128),atol=1.5)
    assert result.arc_support.radius_offset_px==0
    assert result.arc_support.radius_residual_rms_px<=2
    assert result==geometric(image,g)
    np.testing.assert_array_equal(image,before)
    with pytest.raises(FrozenInstanceError): result.central_disk_center_px=(0,0)


def test_canonical_ellipse_and_explicit_prior_not_preset():
    from qcbed.data_preparation import locate_central_disk as geometric
    image,g=scene(radii=(30,20))
    a=geometric(image,g,initial_center_px=(126,129))
    b=geometric(-image,g,initial_center_px=(126,129))
    np.testing.assert_allclose(a.central_disk_center_px,(128,128),atol=1)
    assert a.central_disk_center_px!=a.search_prior_px
    assert a==b


def test_canonical_wrong_radius_failure_has_no_intensity_fallback():
    from qcbed.data_preparation import locate_central_disk as geometric
    y,x=np.indices((256,256))
    image=800/(1+np.exp(np.clip((np.hypot(x-128,y-128)-8)*3,-700,700)))
    with pytest.raises(ValueError,match='No validated central disk boundary'):
        geometric(image,DiskGeometry(24,24))


def test_canonical_dynamic_frame_and_localization_artifact(tmp_path):
    from qcbed.data_preparation import locate_central_disk as geometric
    from qcbed.data_preparation.central_arc import write_central_arc_localization
    import json
    y,x=np.indices((121,163));r=np.hypot(x-81,y-60)
    image=50/(1+np.exp(np.clip((r-20)*2,-700,700)))
    result=geometric(image,DiskGeometry(20,20))
    assert result.central_disk_center_px==(81,60)
    path=write_central_arc_localization(tmp_path,image,result)
    assert path.endswith('central_disk_arc_localization.png')
    data=json.loads((tmp_path/'central_disk_arc_localization.json').read_text())
    assert data['central_disk_center_px']==[81,60]
    assert 'intensity_center_px' not in data
    with pytest.raises(FileExistsError): write_central_arc_localization(tmp_path,image,result)
