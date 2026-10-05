from dataclasses import replace
import json
from pathlib import Path
import numpy as np
import pytest
from qcbed.data_preparation import (convert_scan_geometry, scan_from_catalogue,
    assemble_prepared_dataset, ExtractionConfig, ArcBeamCatalogue, BootstrapSite,
    DiskGeometry, prepare_experimental_qcbed)


def convert(v1=(1.,0.),v2=(0.,1.),**kw):
    args=dict(reference_point_px=(10.,20.),zone_axis_origin_px=(8.,17.),
        reference_point_1_px=tuple(np.array((8.,17.))+v1),reference_point_2_px=tuple(np.array((8.,17.))+v2),
        selection_width=8,selection_height=6,final_binning=2)
    args.update(kw)
    return convert_scan_geometry(**args)


def test_identity_xy_and_source_half_width():
    r=convert()
    assert r.point4==(2.,3.)
    assert r.points==((-1.,1.),(5.,1.),(-1.,5.))
    assert r.scan_geometry.corner_a==(-1.,1.)
    assert r.scan_geometry.corner_b==(5.,1.)
    assert r.scan_geometry.corner_c==(-1.,5.)
    assert (r.scan_geometry.scan_rows,r.scan_geometry.scan_columns)==(3,4)


@pytest.mark.parametrize('v1,v2',[((0.,2.),(-2.,0.)),((4.,1.),(2.,3.)),((-4.,-1.),(2.,3.)),((2.,3.),(4.,1.))])
def test_general_solve_against_literal_dm(v1,v2):
    r=convert(v1,v2)
    corners=(r.scan_geometry.corner_a,r.scan_geometry.corner_b,r.scan_geometry.corner_c)
    for (x,y),actual in zip(r.points,corners):
        a=(x*v2[1]-y*v2[0])/(v1[0]*v2[1]-v1[1]*v2[0])
        b=y/v2[1]-(v1[1]/v2[1])*a if v2[1] else np.linalg.solve(np.column_stack((v1,v2)),[x,y])[1]
        np.testing.assert_allclose(actual,[a,b],rtol=0,atol=1e-14)
        np.testing.assert_allclose(np.column_stack((v1,v2))@actual,[x,y],atol=1e-14)


@pytest.mark.parametrize('v2',[(2.,0.),(2.,1e-14)])
def test_singular_and_pathological(v2):
    with pytest.raises(ValueError,match='Singular/pathological'): convert((1.,0.),v2)



def catalogue():
    def s(id,indices,point): return BootstrapSite(indices,point,None,None,'missing','test',beam_id=id)
    return ArcBeamCatalogue((s('upstream_second',(0,1),(9.,20.)),s('origin',(0,0),(10.,10.)),
                             s('upstream_first',(1,0),(20.,11.))),'test',False)


@pytest.mark.parametrize('transform',[((1,0),(0,1)),((1,0),(0,-1)),((0,1),(1,0))])
def test_explicit_basis_change_and_noncentral_a8_origin(transform):
    c=catalogue(); cfg=ExtractionConfig(8,8,physical_selection_shape_px=(8,8))
    r=scan_from_catalogue(c,central_beam_id='origin',canonical_basis_beam_ids=('upstream_first','upstream_second'),
        zone_axis_origin_px=(11.,12.),extraction_config=cfg,basis_transform=transform)
    canonical=np.array([[10.,-1.],[1.,10.]])
    np.testing.assert_array_equal(np.column_stack((r.vec1,r.vec2)),canonical@transform)
    assert r.point4==(-1.,-2.)
    assert r.canonical_to_scan_basis==transform
    assert r.reference_point_1_px==tuple(np.array((11.,12.))+np.array(r.vec1))


def test_dataset_integration_and_missing_variance_readiness():
    c=catalogue(); cfg=ExtractionConfig(8,8,physical_selection_shape_px=(8,8))
    r=scan_from_catalogue(c,central_beam_id='origin',canonical_basis_beam_ids=('upstream_first','upstream_second'),
        zone_axis_origin_px=(11.,12.),extraction_config=cfg,basis_transform=((1,0),(0,-1)))
    kwargs=dict(selected_beam_ids=['origin'],config=cfg)
    p=np.arange(900,dtype=float).reshape(30,30)
    absent=prepare_experimental_qcbed(p,c,DiskGeometry(4,4),**kwargs)
    assert assemble_prepared_dataset(absent,r,intensity_units='counts').variances is None
    extraction=prepare_experimental_qcbed(p,c,DiskGeometry(4,4),variance_map=np.ones_like(p),uncertainty_provenance='synthetic',**kwargs)
    data=assemble_prepared_dataset(extraction,r,intensity_units='counts')
    for actual,expected in zip((data.intensities,data.variances,data.weights),extraction.tilt_beam_arrays()):
        np.testing.assert_array_equal(actual,expected)
    assert data.scan_geometry==r.scan_geometry
    assert data.reciprocal_frame.zone_axis_origin_px==(11.,12.)
    assert data.beams[0].indices is None


def test_engine_corner_mapping():
    from qcbed.scan import LegacyScanDefinition, generate_scan_tilts
    # Checked against production scan implementation via its public API below.
    r=convert(); s=r.scan_geometry
    scan=LegacyScanDefinition(basis_1=(1.,0.,0.),basis_2=(0.,1.,0.),corner_a=s.corner_a,
        corner_b=s.corner_b,corner_c=s.corner_c,scan_rows=s.scan_rows,scan_columns=s.scan_columns)
    values=generate_scan_tilts(scan)
    np.testing.assert_array_equal(values[0],s.corner_a)
    np.testing.assert_array_equal(values[s.scan_columns-1],s.corner_b)
    np.testing.assert_array_equal(values[-s.scan_columns],s.corner_c)


def test_expansion_does_not_enter_corners_and_geometry_mismatch_rejected():
    assert convert(expansion_factor=2).scan_geometry==convert(expansion_factor=5).scan_geometry
    with pytest.raises(ValueError,match='divisible'): convert(selection_width=7)


def test_invalid_change_of_basis_and_missing_ids():
    kwargs=dict(central_beam_id='origin',canonical_basis_beam_ids=('upstream_first','upstream_second'),
        zone_axis_origin_px=(10.,10.),extraction_config=ExtractionConfig(8,8,physical_selection_shape_px=(8,8)))
    for transform in (((1,0),(0,0)),((1.,0.),(0.,1.))):
        with pytest.raises(ValueError,match='basis_transform'):
            scan_from_catalogue(catalogue(),basis_transform=transform,**kwargs)
    with pytest.raises(ValueError,match='missing beam ID'):
        scan_from_catalogue(catalogue(),basis_transform=((1,0),(0,1)),**(kwargs|{'central_beam_id':'absent'}))




def test_current_geometry_binning_two_uses_actual_selection_dimensions():
    r=convert_scan_geometry(reference_point_px=(127.,126.),zone_axis_origin_px=(127.,126.),
        reference_point_1_px=(208.,135.),reference_point_2_px=(137.,45.),
        selection_width=64,selection_height=64,final_binning=2,expansion_factor=1)
    assert r.points==((-31.,-31.),(31.,-31.),(-31.,31.))
    assert r.scan_geometry.scan_rows==r.scan_geometry.scan_columns==32
    # Literal DM algebra, not inversion of historical output targets.
    expected=[]
    for x,y in r.points:
        a=(x*(-81)-y*10)/(-6651)
        b=y/(-81)-(9/(-81))*a
        expected.append((a,b))
    np.testing.assert_allclose([r.scan_geometry.corner_a,r.scan_geometry.corner_b,r.scan_geometry.corner_c],expected,rtol=0,atol=1e-15)
    assert convert_scan_geometry(reference_point_px=(127.,126.),zone_axis_origin_px=(127.,126.),
        reference_point_1_px=(208.,135.),reference_point_2_px=(137.,45.),
        selection_width=64,selection_height=64,final_binning=2,expansion_factor=1)==r


def test_current_cli_explicit_binning_comparison_and_no_legacy_read(tmp_path,monkeypatch):
    import runpy
    import sys
    import hashlib
    from dataclasses import asdict
    root=Path(__file__).resolve().parents[1]
    cat=tmp_path/'catalogue.csv'
    cat.write_text('beam_id,n1,n2,predicted_x,predicted_y,observed_x,observed_y,status,source\n'
        'center,0,0,127,126,,,reference,test\nright,1,0,208,135,,,observed,test\ndown,0,1,117,207,,,observed,test\n')
    a8=tmp_path/'a8.json'; a8.write_text(json.dumps(dict(zone_axis_origin_px=[127,126],central_disk_center_px=[127,126])))
    a9=tmp_path/'a9.json'; a9.write_text(json.dumps(dict(config=asdict(ExtractionConfig(64,64,final_binning=2,physical_selection_shape_px=(64,64))),
        provenance=dict(catalogue_sha256=hashlib.sha256(cat.read_bytes()).hexdigest()),beams=[],uncertainty='unavailable')))
    original=Path.open
    def guarded(path,*args,**kw):
        if path.name in ('inputparams.dat','inputparas'): raise AssertionError('Legacy runtime dependency')
        return original(path,*args,**kw)
    monkeypatch.setattr(Path,'open',guarded)
    monkeypatch.setattr(sys,'argv',['a10','--a8-measurement',str(a8),'--a9-summary',str(a9),'--catalogue',str(cat),
        '--central-beam-id','center','--canonical-basis-ids','right','down','--basis-transform','1','0','0','-1',
        '--binning','2','--output',str(tmp_path/'result')])
    runpy.run_path(str(root/'examples/a10_scan_geometry_validation.py'),run_name='__main__')
    result=json.loads((tmp_path/'result/scan_geometry.json').read_text())
    assert result['saved_a9_binning']==2 and result['confirmed_conversion_binning']==2
    assert result['a9_array_geometry_matches']
    assert result['maximum_absolute_difference']==pytest.approx(.0886702083897159)
    assert 'passed' not in result
    assert (tmp_path/'result/scan_geometry.md').exists()


def test_a10_requires_physical_selection_and_ignores_diagnostic_size():
    from qcbed.data_preparation import DiskGeometry
    args=dict(central_beam_id='origin',canonical_basis_beam_ids=('upstream_first','upstream_second'),
              zone_axis_origin_px=(10.,10.),basis_transform=((1,0),(0,-1)))
    with pytest.raises(ValueError,match='physical selection'):
        scan_from_catalogue(catalogue(),extraction_config=ExtractionConfig(64,64),**args)
    cfg=ExtractionConfig.from_disk_geometry(DiskGeometry(38,38),final_binning=2,diagnostic_window_shape_px=(64,64))
    a=scan_from_catalogue(catalogue(),extraction_config=cfg,**args)
    b=scan_from_catalogue(catalogue(),extraction_config=replace(cfg,diagnostic_window_shape_px=(20,22)),**args)
    assert a==b
    assert a.points==((-37.,-37.),(37.,-37.),(-37.,37.))
    assert a.scan_geometry.scan_rows==38
