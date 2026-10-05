"""A9 source-equation and ordering regressions; no invented scan geometry."""
from dataclasses import replace
import numpy as np
import pytest
from qcbed.data_preparation import (ArcBeamCatalogue, BootstrapSite, DiskGeometry,
    ReflectionROI, SignalMode, WeightMode, ExtractionConfig,
    build_beam_geometry_mask, construct_signal_fields, prepare_experimental_qcbed)
from qcbed.data_preparation.extraction import _integrate, threshold_support


def site(id='central', center=(12.,12.), index=(0,0), status='reference'):
    return BootstrapSite(index,center,None,None,status,'test',beam_id=id)


def run(pattern=None, sites=None, **kwargs):
    if pattern is None:
        pattern = np.arange(24*24,dtype=float).reshape(24,24)+10
    sites = sites or (site(),)
    return prepare_experimental_qcbed(pattern,ArcBeamCatalogue(tuple(sites),'test',False),
        DiskGeometry(8,6), selected_beam_ids=kwargs.pop('selected_beam_ids',[s.beam_id for s in sites]),
        config=kwargs.pop('config',ExtractionConfig(8,8)),**kwargs)


def test_binary_union_radius_xy_and_boundary():
    geom=DiskGeometry(3.6,1.8)
    r=ReflectionROI('r','a',(4.,3.),geom)
    mask=build_beam_geometry_mask((9,11),(r,r))
    y,x=np.indices(mask.shape)
    np.testing.assert_array_equal(mask,((x-4)/3.6)**2+((y-3)/1.8)**2<1)
    assert mask.dtype==bool and not mask[0,0] and mask[3,7] and not mask[6,4]
    b=run().beams[0]
    assert b.roi.geometry==DiskGeometry(7.2,5.4)


def test_conventional_minimum_scale_uncertainty_and_immutability():
    p=np.arange(576,dtype=float).reshape(24,24)+10
    v=np.full(p.shape,4.)
    before=p.copy()
    r=run(p,variance_map=v,uncertainty_provenance='synthetic variance')
    b=r.beams[0]; local=p[8:16,8:16]
    np.testing.assert_array_equal(r.fields.signal,p)
    minimum=local[b.roi_geometry_mask & (b.support_fraction > 0)].min()
    np.testing.assert_array_equal(b.signal,local-minimum)
    np.testing.assert_array_equal(b.variance,4)
    np.testing.assert_array_equal(p,before); np.testing.assert_array_equal(v,4)
    with pytest.raises(ValueError): b.signal[0,0]=0
    with pytest.raises(ValueError): b.signal.setflags(write=True)
    unchanged=run(p,config=ExtractionConfig(8,8,conventional_minimum_subtraction=False))
    np.testing.assert_array_equal(unchanged.beams[0].signal,local)
    assert unchanged.beams[0].variance is None
    assert unchanged.beams[0].noise_valid_fraction is None


def test_gradient_equations_and_no_wrap():
    y,x=np.indices((9,11)); p=(x*x+2*y*y).astype(float)
    m=np.ones(p.shape,bool); m[3,3]=False
    f=construct_signal_fields(p,m,signal_mode=SignalMode.GRADIENT,
        variance_map=np.full(p.shape,8.),uncertainty_provenance='synthetic')
    offsets=((-1,0),(1,0),(0,-1),(0,1),(-1,-1),(1,1),(-1,1),(1,-1))
    expected=sum(p[4+dy,4+dx] for dx,dy in offsets)/8-p[4,4]
    assert f.signal[4,4]==expected==2.25
    assert f.variance[4,4]==9
    assert f.support_fraction[4,4]==8/9
    assert not f.valid[0,0]
    assert f.signal[0,0]==sum(p[dy,dx] for dx,dy in ((1,0),(0,1),(1,1)))/8


def test_threshold_strict():
    values=np.array([np.nextafter(.875,0),.875,np.nextafter(.875,1)])
    np.testing.assert_array_equal(threshold_support(values),[0,0,1])


def test_noise_mode_requires_inputs_and_does_not_double_inverse_variance():
    with pytest.raises(ValueError,match='requires explicit A4'):
        run(weight_mode=WeightMode.NOISE_CHARACTERIZATION)
    mask=np.ones((24,24),bool); mask[12,12]=False
    r=run(weight_mode=WeightMode.NOISE_CHARACTERIZATION,noise_valid_mask=mask,
          sigma_map=np.full(mask.shape,3.),uncertainty_provenance='synthetic A4')
    b=r.beams[0]
    assert b.final_weight[4,4]==0 and b.final_weight[4,5]==1
    np.testing.assert_array_equal(b.variance,9)
    np.testing.assert_array_equal(b.noise_weight,1)
    assert b.noise_valid_pixels==b.in_frame_pixels-1
    with pytest.raises(ValueError,match='not mixed'): run(noise_valid_mask=mask)


def test_clipping_rejection_and_no_wrap():
    r=run(sites=(site(center=(0.,0.)),))
    b=r.beams[0]
    assert 0<b.coverage_fraction<1 and b.status=='rejected_coverage'
    assert not b.final_weight.any() and not b.valid_mask[0,0]
    assert not r.fields.beam_geometry_mask[-1,-1]


def test_order_duplicates_flatten_and_modes_same_shapes():
    sites=(site('z'),site('a',(16.,16.),(0,0),'observed'))
    r=run(sites=sites,selected_beam_ids=['a','z'])
    assert [b.beam.id for b in r.beams]==['z','a']
    assert r.beams[0].beam.indices is None  # abstract indices are not HKL
    s,v,w=r.blocks()
    np.testing.assert_array_equal(s[:64],r.beams[0].signal.ravel())
    np.testing.assert_array_equal(s[64:],r.beams[1].signal.ravel())
    np.testing.assert_array_equal(r.tilt_beam_arrays()[0][:,1],s[64:])
    g=run(sites=sites,signal_mode=SignalMode.GRADIENT)
    assert [b.signal.shape for b in g.beams]==[b.signal.shape for b in r.beams]
    assert r.beams[0].bounds_xyxy==(8,8,16,16)


@pytest.mark.parametrize('q,e',[(1,1),(2,1),(2,3),(4,3)])
def test_literal_source_expansion_normalization(q,e):
    a=np.arange(64,dtype=float).reshape(8,8)
    config=ExtractionConfig(8,8,final_binning=q,expansion_factor=e)
    expected=np.zeros((8//q,8//q)); variance=expected.copy()
    offset=(e-1)//2
    for y in range(8//q):
        for x in range(8//q):
            values=[a[((y*e+offset)*q+j)//e,((x*e+offset)*q+i)//e]
                    for j in range(q) for i in range(q)]
            expected[y,x]=sum(values)/q**2
            variance[y,x]=sum(values)*e**2/q**4
    np.testing.assert_array_equal(_integrate(a,config),expected)
    np.testing.assert_array_equal(_integrate(a,config,variance=True),variance)


def test_no_silent_fractional_snapping_or_unverified_rounding():
    with pytest.raises(ValueError,match='fractional selection'): run(sites=(site(center=(12.2,12.)),))
    with pytest.raises(ValueError,match='odd expansion'): ExtractionConfig(8,8,expansion_factor=2)
    with pytest.raises(ValueError,match='Rejected'): run(sites=(site(status='rejected'),))


def test_predicted_center_preserved_and_invalid_gradient_stencil_excluded():
    s=replace(site(),observed_center_px=(13.,12.))
    valid=np.ones((24,24),bool); valid[12,12]=False
    r=run(sites=(s,),signal_mode=SignalMode.GRADIENT,detector_valid_mask=valid)
    b=r.beams[0]
    assert b.beam.center_px==(12.,12.) and b.observed_center_px==(13.,12.)
    assert not b.final_weight[3:6,3:6].any()


def test_diagnostics_roundtrip_arrays_and_no_overwrite(tmp_path):
    from qcbed.data_preparation.extraction_diagnostics import write_extraction_diagnostics
    import json
    r=run()
    output=tmp_path/'diagnostics'
    summary=write_extraction_diagnostics(output,r.fields.signal,r,DiskGeometry(8,6),provenance={'test':True})
    assert summary['selected_beam_count']==summary['extracted_beam_count']==1
    assert json.loads((output/'validation_summary.json').read_text())['beams'][0]['variance_mean'] is None
    np.testing.assert_array_equal(np.load(output/'beam_000_signal.npy'),r.beams[0].signal)
    with pytest.raises(FileExistsError):
        write_extraction_diagnostics(output,r.fields.signal,r,DiskGeometry(8,6),provenance={})


def test_literal_gradient_stride_and_variance_variable_map():
    p=np.arange(100,dtype=float).reshape(10,10)**2
    v=p/10
    f=construct_signal_fields(p,np.ones(p.shape,bool),final_binning=2,
        signal_mode=SignalMode.GRADIENT,variance_map=v,uncertainty_provenance='synthetic')
    offsets=((-2,0),(2,0),(0,-2),(0,2),(-2,-2),(2,2),(-2,2),(2,-2))
    assert f.signal[5,5]==sum(p[5+dy,5+dx] for dx,dy in offsets)/8-p[5,5]
    assert f.variance[5,5]==v[5,5]+sum(v[5+dy,5+dx] for dx,dy in offsets)/64


def test_rejected_roi_does_not_support_other_beams():
    r=run(sites=(site(center=(0.,0.)),))
    assert not r.fields.beam_geometry_mask.any()


def test_a4_generated_products_are_accepted_without_refitting():
    from qcbed.data_preparation import characterize_noise, variance_map_from_noise_model, NoiseFitConfig
    rng=np.random.default_rng(42)
    y,x=np.indices((128,128)); pattern=100+10*x+3*y+rng.normal(0,8,(128,128))
    fit=characterize_noise(pattern,config=NoiseFitConfig(sampling_stride=2))
    variance=variance_map_from_noise_model(pattern,fit.model)
    r=run(pattern,sites=(site(center=(64.,64.)),),weight_mode=WeightMode.NOISE_CHARACTERIZATION,
          variance_map=variance,noise_valid_mask=fit.noise_mask,
          uncertainty_provenance='synthetic A4 characterize_noise with explicit config')
    np.testing.assert_array_equal(r.beams[0].variance,variance[60:68,60:68])
    assert r.fields.uncertainty_provenance.startswith('variance supplied directly:')


def test_gradient_does_not_subtract_minimum_and_noise_invalid_is_not_used():
    p=np.ones((24,24)); p[12,12]=10
    mask=np.ones(p.shape,bool); mask[10,10]=False
    r=run(p,signal_mode=SignalMode.GRADIENT,weight_mode=WeightMode.NOISE_CHARACTERIZATION,
        variance_map=np.ones(p.shape),noise_valid_mask=mask,uncertainty_provenance='synthetic')
    assert r.beams[0].signal[4,4]==-9
    assert not r.beams[0].final_weight[1:4,1:4].any()


@pytest.mark.parametrize('radii,center,shape,expected,retained',[
    ((5.5,4.),(40.,40.),(16,16),1.,True),
    ((5.5,4.),(4.,40.),(16,16),62/65,True),
    ((3.,2.),(1.5,40.5),(9,9),.90,True),
    ((3.5,4.),(2.,40.),(16,16),40/45,False),
])
def test_ellipse_coverage_inclusive_90_percent(radii,center,shape,expected,retained):
    s=site(center=center,status='edge_clipped')
    r=prepare_experimental_qcbed(np.full((80,80),7.),ArcBeamCatalogue((s,),'test',False),
        DiskGeometry(*radii),selected_beam_ids=[s.beam_id],
        config=ExtractionConfig(*shape,beam_radius_factor=1.))
    b=r.beams[0]
    assert b.beam_roi_frame_coverage==expected
    assert (b.status=='extracted')==retained
    assert b.signal.shape==shape[::-1]
    assert not b.final_weight[~b.detector_in_frame_mask].any()
    assert not b.signal[~b.detector_in_frame_mask].any()
    if expected<1: assert b.window_frame_coverage<expected
    if retained: assert b.final_weight.any()


def test_clipped_minimum_uses_only_measured_beam_support_and_safe_variance():
    y,x=np.indices((80,80)); p=(100+y+x).astype(float)
    # Very low measured data outside the beam must not set its minimum either.
    p[32,0]=-999
    s=site(center=(4.,40.),status='missing')
    r=prepare_experimental_qcbed(p,ArcBeamCatalogue((s,),'test',False),DiskGeometry(5.5,4.),
        selected_beam_ids=[s.beam_id],config=ExtractionConfig(16,16,beam_radius_factor=1.),
        variance_map=np.full(p.shape,4.),uncertainty_provenance='synthetic variance')
    b=r.beams[0]; support=b.final_weight.astype(bool)
    assert b.status=='extracted' and b.signal[support].min()==0
    assert not b.signal[~b.detector_in_frame_mask].any()
    assert np.all(b.variance[~b.detector_in_frame_mask]==1)
    # The unchanged objective point formula remains finite for missing samples.
    contribution=b.final_weight*(b.signal-5)**2/b.variance
    assert np.isfinite(contribution).all()
    assert not contribution[~b.detector_in_frame_mask].any()


def test_edge_gradient_stencil_and_deterministic_predicted_selection():
    p=np.arange(80*80,dtype=float).reshape(80,80)
    c=ArcBeamCatalogue((site('edge',(4.,40.),status='missing'),site('inside',(40.,40.))), 'test',False)
    kwargs=dict(selected_beam_ids=['inside','edge'],config=ExtractionConfig(16,16,beam_radius_factor=1.),
                signal_mode=SignalMode.GRADIENT)
    a=prepare_experimental_qcbed(p,c,DiskGeometry(5.5,4.),**kwargs)
    b=prepare_experimental_qcbed(p,c,DiskGeometry(5.5,4.),**kwargs)
    assert [v.beam.id for v in a.beams]==['edge','inside']
    assert all(v.status=='extracted' for v in a.beams)
    assert a.beams[0].signal.shape==a.beams[1].signal.shape
    edge=a.beams[0]
    assert not edge.final_weight[:, :5].any() # x<0 and x=0 incomplete stencil
    assert edge.final_weight[:,5:].any()
    for aa,bb in zip(a.blocks(),b.blocks()):
        if aa is not None: np.testing.assert_array_equal(aa,bb)


def test_clipped_binned_placeholders_and_edge_diagnostics(tmp_path):
    from qcbed.data_preparation.extraction_diagnostics import write_extraction_diagnostics
    s=site(center=(4.,40.),status='edge_clipped')
    p=np.full((80,80),20.)
    r=prepare_experimental_qcbed(p,ArcBeamCatalogue((s,),'test',False),DiskGeometry(5.5,4.),
        selected_beam_ids=[s.beam_id],config=ExtractionConfig(16,16,beam_radius_factor=1.,final_binning=2,expansion_factor=3),
        variance_map=np.full(p.shape,4.),uncertainty_provenance='synthetic')
    b=r.beams[0]
    assert b.signal.shape==(8,8) and b.final_weight.any()
    assert not b.final_weight[:,:2].any()
    assert np.all(b.variance[:,:2]==1) and not b.signal[:,:2].any()
    summary=write_extraction_diagnostics(tmp_path/'edge',p,r,DiskGeometry(5.5,4.),provenance={'synthetic':True})
    row=summary['beams'][0]
    assert row['retained_by_a9'] and row['beam_roi_frame_coverage']==62/65
    assert row['a7_status']=='edge_clipped'
    assert (tmp_path/'edge'/'beam_000_edge_support.png').exists()


def test_real_catalogue_geometry_counts_without_reclassifying_observations():
    coordinates=[(0,0),(-1,0),(1,0),(0,-1),(0,1),(-1,-1),(1,1),(-1,1),(1,-1),
                 (-2,0),(2,0),(0,-2),(0,2),(-2,-1),(2,1),(-1,2),(1,-2)]
    sites=tuple(site(str(i),(127.+81*n-10*m,126.+9*n+81*m),(n,m),'missing')
                for i,(n,m) in enumerate(coordinates))
    result=prepare_experimental_qcbed(np.ones((256,256)),ArcBeamCatalogue(sites,'synthetic intensity / real geometry',False),
        DiskGeometry(38,38),selected_beam_ids=[s.beam_id for s in sites],config=ExtractionConfig(64,64))
    assert [b.in_frame_pixels for b in result.beams]==[3681]*9+[0,7,0,24,320,369,369,228]
    assert [b.integer_indices for b in result.beams if b.status=='extracted']==coordinates[:9]
    assert all(b.catalogue_status=='missing' for b in result.beams)


def test_legacy_q2_gradient_on_original_grid_then_q2_semiintegration():
    y,x=np.indices((40,40)); p=(x**4+3*y**3).astype(float)
    variance=1+x.astype(float)+y
    cfg=ExtractionConfig(12,12,final_binning=2,physical_selection_shape_px=(12,12))
    s=site(center=(20.,20.))
    result=prepare_experimental_qcbed(p,ArcBeamCatalogue((s,),'test',False),DiskGeometry(14,14),
        selected_beam_ids=[s.beam_id],config=cfg,signal_mode=SignalMode.GRADIENT,
        variance_map=variance,uncertainty_provenance='synthetic')
    offsets=((-2,0),(2,0),(0,-2),(0,2),(-2,-2),(2,2),(-2,2),(2,-2))
    expected=np.empty((12,12)); v=np.empty_like(expected)
    for j in range(12):
        for i in range(12):
            yy,xx=j+14,i+14
            expected[j,i]=sum(p[yy+dy,xx+dx] for dx,dy in offsets)/8-p[yy,xx]
            v[j,i]=variance[yy,xx]+sum(variance[yy+dy,xx+dx] for dx,dy in offsets)/64
    expected_binned=expected.reshape(6,2,6,2).sum(axis=(1,3))/4
    expected_variance=v.reshape(6,2,6,2).sum(axis=(1,3))/16
    np.testing.assert_array_equal(result.beams[0].signal,expected_binned)
    np.testing.assert_array_equal(result.beams[0].variance,expected_variance)
    one=construct_signal_fields(p,np.ones_like(p,bool),signal_mode=SignalMode.GRADIENT,final_binning=1)
    assert one.signal[20,20]!=result.fields.signal[20,20]
    pre_binned=p.reshape(20,2,20,2).mean(axis=(1,3))
    wrong=construct_signal_fields(pre_binned,np.ones_like(pre_binned,bool),signal_mode=SignalMode.GRADIENT,final_binning=2)
    assert not np.array_equal(wrong.signal[7:13,7:13],expected_binned)
    mask=result.fields.beam_geometry_mask
    yy,xx=20,31
    assert result.fields.support_fraction[yy,xx]==(mask[yy,xx]+sum(float(mask[yy+dy,xx+dx]) for dx,dy in offsets))/9


def test_physical_diameter_selection_and_diagnostic_independence():
    cfg=ExtractionConfig.from_disk_geometry(DiskGeometry(38,38),final_binning=2,diagnostic_window_shape_px=(64,64))
    assert cfg.physical_selection_shape_px==(76,76)
    assert cfg.final_pattern_matching_shape_px==(38,38)
    assert cfg.gradient_step_px==cfg.final_binning==2
    with pytest.raises(ValueError,match='integral'):
        ExtractionConfig.from_disk_geometry(DiskGeometry(38.1,38),final_binning=2)
    with pytest.raises(ValueError,match='must match'):
        ExtractionConfig(64,64,physical_selection_shape_px=(76,76))


def test_conventional_physical_selection_q2_integrates_after_minimum():
    y,x=np.indices((40,40)); p=(100+x+2*y).astype(float)
    cfg=ExtractionConfig.from_disk_geometry(DiskGeometry(6,6),final_binning=2)
    s=site(center=(20.,20.))
    r=prepare_experimental_qcbed(p,ArcBeamCatalogue((s,),'test',False),DiskGeometry(6,6),
        selected_beam_ids=[s.beam_id],config=cfg)
    np.testing.assert_array_equal(r.fields.signal,p)
    b=r.beams[0]; local=p[14:26,14:26]
    minimum=local[b.roi_geometry_mask].min()
    np.testing.assert_array_equal(b.signal,(local-minimum).reshape(6,2,6,2).mean(axis=(1,3)))
    assert b.signal.shape==cfg.final_pattern_matching_shape_px==(6,6)


def test_explicit_nearest_grid_preserves_physical_center_and_samples():
    cfg=ExtractionConfig(8,8,sampling_policy='nearest_detector_grid',conventional_minimum_subtraction=False)
    p=np.arange(576,dtype=float).reshape(24,24)+10
    r=run(p,sites=(site(center=(12.5,12.25)),),config=cfg,
          variance_map=np.full_like(p,4.),uncertainty_provenance='test')
    b=r.beams[0]
    assert b.roi.center_px==b.predicted_center_px==(12.5,12.25)
    assert b.bounds_xyxy==(8,8,16,16)
    np.testing.assert_array_equal(b.signal,p[8:16,8:16])
    np.testing.assert_array_equal(b.variance,4.*np.ones((8,8)))
    from qcbed.data_preparation.extraction import ellipse_mask
    np.testing.assert_array_equal(b.roi_geometry_mask,ellipse_mask((8,8),(4.5,4.25),b.roi.geometry))
    assert any('sampling_offset_px=(-0.5, -0.25)' in v for v in b.provenance)
    strict=run(p,config=ExtractionConfig(8,8))
    nearest=run(p,config=ExtractionConfig(8,8,sampling_policy='nearest_detector_grid'))
    for name in ('signal','variance','final_weight'):
        np.testing.assert_array_equal(getattr(strict.beams[0],name),getattr(nearest.beams[0],name))
