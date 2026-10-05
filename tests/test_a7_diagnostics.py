"""Small diagnostic fixtures; no expensive real-data geometry assertions."""
import csv
import json
from types import SimpleNamespace as NS
import numpy as np
import pytest
from qcbed.data_preparation import DiskGeometry, BeamObservation, BeamStatus, PredictedBeam
from qcbed.data_preparation.diagnostics import (load_pattern, beam_filename, crop_bounds,
    ellipse_overlay, beam_rows, write_beam_table, write_figures, validate_pattern)


def catalogue():
    geometry=DiskGeometry(3,2)
    beams=(PredictedBeam(BeamObservation('a',(10.,9.)),(0,0),(10.,9.),(10.5,8.75),True,BeamStatus.OBSERVED,7.,1.,1.),
           PredictedBeam(BeamObservation('b',(3.,2.)),(-1,1),(3.,2.),None,False,BeamStatus.MISSING,.01,1.,1.))
    return NS(geometry=geometry,beams=beams,lattice=NS(assignments=(NS(integer_indices=(0,0),residual_px=.3,inlier=True),)))


def test_text_shape_not_reshaped(tmp_path):
    path=tmp_path/'pattern.dat'; path.write_text('1 2 3\n4 5 6\n')
    result=load_pattern(path)
    assert result.shape==(2,3) and result.dtype==np.float64
    assert result[1,0]==4


def test_crop_geometry_and_overlay():
    g=DiskGeometry(3,2)
    assert crop_bounds((20,30),(10.,9.),g)==((6,15,6,13),False)
    assert crop_bounds((20,30),(3.,2.),g)==((0,8,0,6),True)
    assert ellipse_overlay((10.,9.),g)==dict(xy=(10.,9.),width=6.,height=4.)
    assert beam_filename(2,(-1,3))=='beam_002_nm1_3.png'
    with pytest.raises(ValueError):
        crop_bounds((20,30),(3,2),g,0)


def test_table_keeps_absent_observations(tmp_path):
    rows=write_beam_table(tmp_path/'table.csv',catalogue())
    assert rows[0]['predicted_x_px']==10
    assert rows[0]['observed_x_px']==10.5
    assert rows[0]['predicted_to_observed_shift_px']==pytest.approx(np.hypot(.5,.25))
    assert rows[0]['lattice_residual_px']==.3
    assert rows[1]['observed_x_px'] is None and rows[1]['lattice_residual_px'] is None
    with (tmp_path/'table.csv').open() as stream:
        serialized=list(csv.DictReader(stream))
    assert serialized[1]['observed_x_px']==''
    with pytest.raises(FileExistsError):
        write_beam_table(tmp_path/'table.csv',catalogue())


def test_figures_small_pattern_and_clipping(tmp_path):
    image=np.arange(600.).reshape(20,30); before=image.copy()
    seed=NS(assignments=(NS(candidate=NS(center_px=(10.2,8.8))),))
    records=write_figures(tmp_path,image,NS(central_disk_center_px=(10,9)),seed,catalogue())
    assert (tmp_path/'beam_catalogue_overview.png').stat().st_size>0
    assert len(list((tmp_path/'beams').glob('*.png')))==2
    assert records[1]['clipped'] and not records[0]['clipped']
    np.testing.assert_array_equal(image,before)


def test_missing_inputs_reported_without_geometry_execution(tmp_path,monkeypatch):
    import qcbed.data_preparation.diagnostics as module
    def fail(*a,**k):
        raise AssertionError('Geometry ran despite missing inputs')
    monkeypatch.setattr(module,'locate_central_disk',fail)
    path=tmp_path/'pattern.dat'; np.savetxt(path,np.ones((16,24)))
    output=tmp_path/'result'
    result=validate_pattern(path,output,expected_shape=(2048,2048))
    assert result['status']=='blocked' and len(result['blockers'])==2
    assert result['a5'] is None and result['a7'] is None
    assert json.loads((output/'validation_summary.json').read_text())['input_shape']==[16,24]
    assert not (output/'beam_catalogue.csv').exists()
    with pytest.raises(FileExistsError):
        validate_pattern(path,output)


def test_external_geometry_dynamic_shape_and_provenance(tmp_path,monkeypatch):
    import qcbed.data_preparation.diagnostics as module
    g=DiskGeometry(3.,2.)
    path=tmp_path/'pattern.dat'; np.savetxt(path,np.ones((20,30)))
    seen=[]
    def locate(image,geometry):
        seen.append((image.shape,geometry))
        return NS(central_disk_center_px=(15.,10.))
    def forbidden(*a,**kw):
        raise AssertionError('Manual geometry must not run profile estimation')
    monkeypatch.setattr(module,'estimate_disk_geometry',forbidden)
    monkeypatch.setattr(module,'locate_central_disk',locate)
    monkeypatch.setattr(module,'detect_lattice_candidates',lambda *a,**kw: ())
    monkeypatch.setattr(module,'write_candidate_diagnostics',lambda *a,**kw:
                        dict(candidate_count=0,overview='overview',candidate_table='table',
                             arc_overview='arc',arc_table='arc_table',arc_support_pixels='pixels'))
    result=validate_pattern(path,tmp_path/'manual',geometry=g,geometry_provenance='Independent visual estimate')
    assert seen==[((20,30),g)]
    assert result['expected_shape'] is None and result['input_shape']==[20,30]
    assert result['geometry_input']['source']=='manually_supplied'
    assert result['geometry_input']['radii']==(3.,2.)
    assert not result['a5']['geometry_estimated_by_a5']
    assert result['a6'] is None and result['a7'] is None
    assert 'Too few' in result['blockers'][0]


def test_geometry_provenance_required_and_exclusive(tmp_path):
    path=tmp_path/'pattern.dat'; np.savetxt(path,np.ones((20,30)))
    result=validate_pattern(path,tmp_path/'missing_provenance',geometry=DiskGeometry(3,2))
    assert result['status']=='blocked' and result['blockers']==['Explicit geometry provenance is required']
    result=validate_pattern(path,tmp_path/'conflict',geometry=DiskGeometry(3,2),radii=(3,2),geometry_provenance='manual')
    assert result['status']=='failed' and 'only one' in result['blockers'][0]


def test_empty_candidate_overview_is_not_a_beam_catalogue(tmp_path):
    from qcbed.data_preparation.diagnostics import write_candidate_diagnostics
    from qcbed.data_preparation import LatticeDetectionConfig
    image=np.arange(65536.).reshape(256,256); before=image.copy()
    record=write_candidate_diagnostics(tmp_path,image,DiskGeometry(34.5,34.5),
                                       (123.5,129.5),(),LatticeDetectionConfig())
    assert record['eligible_center_pixels']>966
    assert record['physical_disk_radii']==(34.5,34.5)
    assert record['candidate_count']==0
    assert (tmp_path/'candidate_detection_overview.png').is_file()
    assert not (tmp_path/'beam_catalogue.csv').exists()
    np.testing.assert_array_equal(image,before)


def test_validated_arc_review_table_and_pixel_membership(tmp_path):
    from qcbed.data_preparation import LatticeDetectionConfig,detect_lattice_candidates,arc_center_votes
    from qcbed.data_preparation.diagnostics import write_arc_validation
    y,x=np.indices((160,220));g=DiskGeometry(20,20)
    image=50/(1+np.exp(np.clip((np.hypot(x-70,y-70)-20)*3,-700,700)))
    candidates=detect_lattice_candidates(image,g,(170,120))
    assert candidates
    _,arrays,_=arc_center_votes(image,g,(170,120))
    before=image.copy()
    result=write_arc_validation(tmp_path,image,g,candidates,arrays)
    with open(result['arc_table']) as stream:rows=list(csv.DictReader(stream))
    assert len(rows)==len(candidates)
    c=candidates[0];r=rows[0]
    assert float(r['center_vote_rms_px'])==c.center_vote_rms_px
    assert float(r['mean_normal_residual_deg'])==c.mean_normal_residual_deg
    assert float(r['angular_coverage_deg'])==360*c.angular_coverage
    assert int(r['support_pixel_count'])==len(c.supporting_pixels)
    assert int(r['fragment_count'])==c.fragment_count
    with open(result['arc_support_pixels']) as stream:points=list(csv.DictReader(stream))
    assert [(int(r['x_px']),int(r['y_px'])) for r in points if r['candidate_id']==c.id]==list(c.supporting_pixels)
    assert result['retained_support_pixels']>0
    assert result['not_retained_gradient_pixels']>0
    assert (tmp_path/'arc_validation_overview.png').is_file()
    np.testing.assert_array_equal(image,before)
