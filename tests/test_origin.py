"""A8 source-equation regression; no DM fractional-subarea equivalence claim."""
from dataclasses import replace,FrozenInstanceError
from types import SimpleNamespace as NS
import numpy as np
import pytest
from qcbed.data_preparation.origin import *
from qcbed.data_preparation import DiskGeometry


def problem():
    y,x=np.indices((100,120));image=(x+3*y).astype(float)
    # The local symmetric feature is off the beam center; the rest is asymmetric.
    rng=np.random.default_rng(42);crop=rng.random((20,20));crop+=crop[::-1,::-1]
    image[40:60,55:75]=crop
    reference=SymmetryReference('chosen',(60,50),(5,0),(1,0))
    return image,reference


@pytest.mark.parametrize('shape',[(3,3),(4,4),(3,6),(6,3)])
def test_rotation_exact_even_odd(shape):
    crop=np.arange(np.prod(shape)).reshape(shape)
    expected=np.array([[crop[shape[0]-1-y,shape[1]-1-x] for x in range(shape[1])] for y in range(shape[0])])
    np.testing.assert_array_equal(rotate180(crop),expected)
    np.testing.assert_array_equal(rotate180(rotate180(crop)),crop)


def test_legacy_loop_equivalence_and_distinct_centers():
    image,reference=problem();before=image.copy()
    m=measure_symmetry_origin(image,DiskGeometry(20,20),(30,50),reference)
    assert m.feature_seed_px==(65,50) and m.reference.beam_center_px==(60,50)
    assert m.symmetry_center_px==(65,50) and m.crop_pixel_midpoint_px==(64.5,49.5)
    assert m.zone_axis_origin_px==(50,50) and m.mismatch==0
    assert m.reference_kind=="diffracted"
    assert m.trial_count==100 and m.search_half_size_px==(5,5) and m.crop_half_size_px==(10,10)
    trials=[]
    for j in range(45,55):
        for i in range(60,70):
            cut=image[j-10:j+10,i-10:i+10]
            rotated=np.array([[cut[19-y,19-x] for x in range(20)] for y in range(20)])
            trials.append((i,j,float(abs(cut-rotated).sum())))
    assert m.trials==tuple(trials)
    assert m==measure_symmetry_origin(image,DiskGeometry(20,20),(30,50),reference)
    np.testing.assert_array_equal(image,before)
    with pytest.raises(FrozenInstanceError): m.mismatch=4


def test_selection_not_order_and_no_auto_selection():
    image,ref=problem()
    records=[NS(beam_id='other',integer_indices=(0,1),predicted_center_px=(30,50),status='observed'),
             NS(beam_id='unused',integer_indices=(1,1),predicted_center_px=(80,80),status='missing'),
             NS(beam_id='chosen',integer_indices=(1,0),predicted_center_px=(60,50),status='observed')]
    for order in (records,records[::-1]):
        selected=select_symmetry_reference(NS(sites=order),beam_id='chosen',feature_offset_px=(5,0))
        assert selected==ref
        assert measure_symmetry_origin(image,DiskGeometry(20,20),(30,50),selected).mismatch==0
    other=select_symmetry_reference(NS(sites=records),lattice_indices=(0,1),feature_offset_px=(5,0))
    assert measure_symmetry_origin(image,DiskGeometry(20,20),(30,50),other).mismatch>0
    with pytest.raises(ValueError,match='specified'): select_symmetry_reference(NS(sites=records),feature_offset_px=(0,0))
    with pytest.raises(ValueError,match='specified'): measure_symmetry_origin(image,DiskGeometry(20,20),(30,50))
    with pytest.raises(ValueError,match='reliable'): select_symmetry_reference(NS(sites=records),beam_id='unused',feature_offset_px=(0,0))


def test_offset_follows_geometry_and_odd_crop():
    image,ref=problem()
    moved=replace(ref,beam_center_px=(61,50))
    m=measure_symmetry_origin(image,DiskGeometry(20,20),(30,50),moved)
    assert m.feature_seed_px==(66,50) and m.symmetry_center_px==(65,50)
    assert m.zone_axis_origin_px==(49.5,50)
    # Half-integer DM trial coordinates yield integral bounds and odd crops.
    m=measure_symmetry_origin(image,DiskGeometry(20,20),(30,50),ref,
        config=SymmetrySearchConfig((1.5,1.5),(2.5,2.5)))
    x0,y0,x1,y1=m.crop_bounds_xyxy
    assert x1-x0==y1-y0==5


def test_consistency_never_silently_averages():
    image,ref=problem();m=measure_symmetry_origin(image,DiskGeometry(20,20),(30,50),ref)
    shifted=replace(m,zone_axis_origin_px=(51,50))
    assert combine_origin_measurements([m,shifted],tolerance_px=1).combined_origin_px is None
    result=combine_origin_measurements([m,shifted],tolerance_px=1,combine=True)
    assert result.consistent and result.combined_origin_px==(50.5,50)
    result=combine_origin_measurements([m,shifted],tolerance_px=.5,combine=True)
    assert not result.consistent and result.combined_origin_px is None and len(result.measurements)==2


def test_invalid_crop_and_flat_tie_diagnostics(tmp_path):
    image,ref=problem()
    with pytest.raises(ValueError,match='Fractional'): measure_symmetry_origin(image,DiskGeometry(38,38),(30,50),ref)
    with pytest.raises(ValueError,match='fit the frame'): measure_symmetry_origin(image,DiskGeometry(20,20),(30,50),replace(ref,beam_center_px=(1,1),feature_offset_px=(0,0)))
    m=measure_symmetry_origin(np.ones_like(image),DiskGeometry(20,20),(30,50),ref)
    assert m.tied_minimum_count==100 and m.boundary_limited and not m.crop_has_contrast
    assert m.symmetry_center_px==(60,45)  # legacy first tie, not invented quality
    path=write_symmetry_diagnostics(tmp_path/'diagnostic',np.ones_like(image),m)
    assert path.exists()
    with pytest.raises(FileExistsError): write_symmetry_diagnostics(tmp_path/'diagnostic',np.ones_like(image),m)


def test_formula_and_inconsistent_selectors():
    assert legacy_zone_axis_origin((10,30),(40,50),(45,62))==(30,52)
    catalogue=NS(sites=(NS(beam_id='a',integer_indices=(1,0),predicted_center_px=(60,50),status='observed'),))
    with pytest.raises(ValueError,match='conflicting'):
        select_symmetry_reference(catalogue,beam_id='a',lattice_indices=(0,1),feature_offset_px=(0,0))
    with pytest.raises(ValueError): SymmetrySearchConfig((0,2),(3,3))
    with pytest.raises(ValueError): combine_origin_measurements([],tolerance_px=1)


def test_feature_offset_is_required_and_may_not_leave_disk():
    with pytest.raises(TypeError): SymmetryReference('chosen',(60,50))
    image,ref=problem()
    with pytest.raises(ValueError,match='inside'):
        measure_symmetry_origin(image,DiskGeometry(20,20),(30,50),replace(ref,feature_offset_px=(21,0)))
