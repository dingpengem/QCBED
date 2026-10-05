"""New A6 inference: synthetic pixel lattices, not legacy numerical equivalence."""
from dataclasses import FrozenInstanceError
import numpy as np
import pytest
from qcbed.data_preparation import (DiskGeometry, LatticeDetectionConfig,
    LatticeCandidate, detect_reciprocal_lattice_seed, detect_lattice_candidates,
    disk_response_map, infer_lattice_seed)


def scene(g1=(30,0), g2=(0,30), difficult=False):
    shape=(200,240); center=np.array([117.,97.]); geometry=DiskGeometry(5,4)
    yy,xx=np.indices(shape); image=3.+.005*xx+.007*yy
    matrix=np.column_stack((g1,g2)); truth=[]
    for i in range(-2,3):
        for j in range(-2,3):
            if difficult and (i,j) in ((1,1),(-1,2),(2,-2)):
                continue
            x,y=center+matrix@[i,j]
            if not (15<x<225 and 15<y<185):
                continue
            r=np.hypot((xx-x)/5,(yy-y)/4)
            image+=50/(1+np.exp(np.clip((r-1)*20,-700,700)))
            if (i,j)!=(0,0): truth.append((x,y))
    if difficult:
        image+=np.random.default_rng(42).normal(0,.15,shape)
    return image,tuple(center),geometry,matrix,np.array(truth)


def assert_same_lattice(result, expected, atol=.08):
    fitted = np.column_stack((result.basis_vector_1_px, result.basis_vector_2_px))
    transform = np.linalg.solve(expected, fitted)
    np.testing.assert_allclose(transform, np.rint(transform), atol=atol)
    assert abs(round(np.linalg.det(transform))) == 1  # primitive observed lattice, not supercell
    assert np.linalg.det(fitted) > 0
    return fitted


@pytest.mark.parametrize('g1,g2', [((30,0),(0,30)), ((32,4),(9,27)), ((38,0),(0,24))])
@pytest.mark.parametrize('difficult', [False, True])
def test_full_pipeline_disk_edges(g1,g2,difficult):
    image, center, geometry, matrix, truth = scene(g1,g2,difficult)
    before = image.copy()
    result = detect_reciprocal_lattice_seed(image, geometry, center)
    assert_same_lattice(result, matrix)
    assert result.support_fraction >= .6
    assert result.fit_rms_residual_px < 1
    assert not hasattr(result, 'zone_axis_origin_px')
    np.testing.assert_array_equal(image, before)
    for a in result.assignments:
        assert np.linalg.norm(np.array(a.candidate.center_px)-center) > 8
        assert a.residual_px >= 0


def test_determinism_candidate_order_ids_and_immutable_result():
    image, center, geometry, matrix, _ = scene((32,4),(9,27),True)
    candidates = detect_lattice_candidates(image, geometry, center)
    one = infer_lattice_seed(candidates, center, geometry)
    two = infer_lattice_seed(tuple(reversed(candidates)), center, geometry)
    assert one == two == detect_reciprocal_lattice_seed(image, geometry, center)
    assert len({c.id for c in candidates}) == len(candidates)
    with pytest.raises(FrozenInstanceError):
        one.reference_center_px = (0,0)


def test_response_linear_background_and_no_intensity_modification():
    y,x = np.indices((70,90)); image = 50.+x*.4+y*.7
    response = disk_response_map(image, DiskGeometry(5,3))
    np.testing.assert_allclose(response[8:-8,12:-12], 0, atol=2e-13)
    assert detect_lattice_candidates(image, DiskGeometry(5,3), (45,35)) == ()


def test_residuals_indices_and_retained_outliers():
    geometry = DiskGeometry(3,2); center = (53.,61.)
    matrix = np.array([[25.,7.],[3.,22.]])
    points = [LatticeCandidate(f'p{i}_{j}', tuple(np.array(center)+matrix@[i,j]), 1.)
              for i in range(-2,3) for j in range(-2,3) if (i,j)!=(0,0)]
    points.append(LatticeCandidate('false', (62.,72.),1000))
    result = infer_lattice_seed(points, center, geometry)
    fitted = assert_same_lattice(result,matrix,atol=1e-10)
    for a in result.assignments:
        residual = np.linalg.norm(np.array(a.candidate.center_px)-center-fitted@a.integer_indices)
        assert a.residual_px == pytest.approx(residual)
        assert a.inlier == (a.candidate.id != 'false')
    assert result.fit_rms_residual_px < 1e-12


def test_collinear_too_few_and_unstable_failures():
    geometry = DiskGeometry(3,2)
    with pytest.raises(ValueError,match='Too few'):
        infer_lattice_seed([], (0,0), geometry)
    line = [LatticeCandidate(str(i),(i*20., i*.001),1) for i in range(1,10)]
    with pytest.raises(ValueError,match='non-collinear'):
        infer_lattice_seed(line,(0,0),geometry)
    random = np.random.default_rng(8).uniform(10,200,(10,2))
    points = [LatticeCandidate(str(i),tuple(p),1) for i,p in enumerate(random)]
    with pytest.raises(ValueError,match='No stable'):
        infer_lattice_seed(points,(0,0),geometry)


@pytest.mark.parametrize('image', [np.ones(5), np.ones((0,5)), np.full((50,60),np.nan), np.full((50,60),np.inf)])
def test_invalid_images(image):
    with pytest.raises(ValueError):
        detect_reciprocal_lattice_seed(image,DiskGeometry(4,3),(10,10))


@pytest.mark.parametrize('kwargs', [dict(max_candidates=True),dict(max_families=1),dict(min_inliers=2),
    dict(displacement_tolerance_px=0),dict(min_basis_sine=1.1),dict(min_inlier_fraction=np.nan)])
def test_invalid_config(kwargs):
    with pytest.raises((TypeError,ValueError)):
        LatticeDetectionConfig(**kwargs)


def test_geometry_and_center_validation():
    with pytest.raises(ValueError,match='fit'):
        disk_response_map(np.ones((20,20)),DiskGeometry(10,10))
    with pytest.raises(ValueError,match='inside'):
        detect_lattice_candidates(np.ones((50,60)),DiskGeometry(4,3),(80,20))
    candidates = [LatticeCandidate('same',(i*10.,20.),1) for i in range(8)]
    with pytest.raises(ValueError,match='Duplicate'):
        infer_lattice_seed(candidates,(0,0),DiskGeometry(3,2))


def test_a7_template_response_against_direct_sum():
    image = np.zeros((60,80)); yy,xx = np.indices(image.shape)
    image += 30*np.exp(-2*(((xx-38.3)/5)**2+((yy-29.2)/4)**2))
    geometry = DiskGeometry(5,4)
    response = disk_response_map(image,geometry)
    y,x = np.mgrid[-8:9,-10:11]; r=np.sqrt((x/5)**2+(y/4)**2)
    inner=np.where(r<1,(1-r*r)**2,0.)
    outer=np.where((r>1)&(r<2),np.sin(np.pi*(r-1))**2,0.)
    expected=np.sum(image[21:38,28:49]*(inner/inner.sum()-outer/outer.sum()))
    assert response[29,38] == pytest.approx(expected,abs=1e-13)


def test_close_false_candidates_do_not_define_shortest_basis():
    center=(60.,60.); geometry=DiskGeometry(3,2)
    matrix=np.array([[26.,8.],[2.,24.]])
    points=[LatticeCandidate(f'{i}_{j}',tuple(np.array(center)+matrix@[i,j]),1)
            for i in range(-2,3) for j in range(-2,3) if (i,j)!=(0,0)]
    points += [LatticeCandidate('close1',(67.,67.),100),LatticeCandidate('close2',(68.,69.),90)]
    result=infer_lattice_seed(points,center,geometry)
    assert_same_lattice(result,matrix,1e-10)
    assert all(not a.inlier for a in result.assignments if a.candidate.id.startswith('close'))



def arc_scene(centers=((90,100),), fractions=(1.,), radius=(20.,20.), interior=False, noise=0.):
    y,x=np.indices((230,310)); image=10.+.002*x+.003*y
    for (cx,cy),fraction in zip(centers,fractions):
        dx,dy=(x-cx)/radius[0],(y-cy)/radius[1]
        r=np.hypot(dx,dy); angle=np.arctan2(dy,dx)
        # Only this fraction of the boundary is bright; no closed visible contour.
        angular=np.exp(-((angle/(np.pi*fraction))**12)) if fraction<1 else 1.
        image+=50*angular/(1+np.exp(np.clip((r-1)*50,-700,700)))
        if interior:
            image+=200*np.exp(-((x-cx+7)**2+(y-cy-3)**2)/8)
    image+=np.random.default_rng(14).normal(0,noise,image.shape)
    return image


@pytest.mark.parametrize('fraction',[1.,.2,.3,.4])
@pytest.mark.parametrize('interior',[False,True])
def test_partial_arcs_and_dynamical_interior(fraction,interior):
    image=arc_scene(fractions=(fraction,),interior=interior,noise=.1); before=image.copy()
    g=DiskGeometry(20,20)
    candidates=detect_lattice_candidates(image,g,(250,180))
    assert candidates
    closest=min(candidates,key=lambda c:np.linalg.norm(np.array(c.center_px)-(90,100)))
    assert np.linalg.norm(np.array(closest.center_px)-(90,100))<=3
    assert closest.supporting_edge_pixels>=12
    assert closest.angular_coverage>=.15
    assert closest.radius_residual_rms_px<=4
    assert g==DiskGeometry(20,20)
    np.testing.assert_array_equal(image,before)
    assert candidates==detect_lattice_candidates(image,g,(250,180))


def test_multiple_arcs_ellipse_radius_error_scaling_and_separation():
    centers=((70,80),(155,80),(240,80))
    image=arc_scene(centers,(.2,.3,.4),radius=(21,16),noise=.1)
    geometry=DiskGeometry(20,15)
    c=detect_lattice_candidates(image,geometry,(155,190))
    for center in centers:
        assert min(np.linalg.norm(np.array(a.center_px)-center) for a in c)<3
    scaled=detect_lattice_candidates(image*1000,geometry,(155,190))
    assert [a.center_px for a in c]==[a.center_px for a in scaled]
    np.testing.assert_allclose([a.vote_strength for a in c],[a.vote_strength for a in scaled],atol=1e-10)


@pytest.mark.parametrize('kwargs',[dict(edge_percentile=100),dict(edge_percentile=True),
    dict(edge_smoothing_sigma_px=0),dict(radius_tolerance_px=-1),dict(radius_samples=2),
    dict(candidate_merge_radius_px=float('nan')),dict(minimum_angular_coverage=1.1)])
def test_invalid_arc_config(kwargs):
    with pytest.raises(ValueError): LatticeDetectionConfig(**kwargs)


def test_vote_boundaries_central_exclusion_and_no_old_configuration():
    from qcbed.data_preparation import arc_center_votes
    from dataclasses import asdict
    image=arc_scene(((70,80),(155,80)),(1.,1.))
    votes,arrays,d=arc_center_votes(image,DiskGeometry(20,20),(70,80))
    assert not arrays['selected_edges'][80,70]
    assert np.isfinite(votes).all() and (votes>=0).all()
    assert d['radius_offsets_px']==[-2,-1,0,1,2]
    assert d['physical_disk_radii']==(20,20)
    assert 'candidate_detection_scale' not in asdict(LatticeDetectionConfig())
    c=detect_lattice_candidates(image,DiskGeometry(20,20),(70,80))
    assert all(np.linalg.norm(np.array(a.center_px)-(70,80))>24 for a in c)


def test_partial_arc_lattice_missing_disks_false_edges_and_false_disk():
    centers=[(155+65*i,115+60*j) for i in (-1,0,1) for j in (-1,0,1)
             if (i,j) not in ((0,0),(1,1))]
    centers.append((45,40))  # a strong non-lattice arc must be rejected globally
    image=arc_scene(centers,[.3+.1*(i%2) for i in range(len(centers))],noise=.15,interior=True)
    y,x=np.indices(image.shape)
    image+=15*(y>195)+20*(x>285)  # unrelated straight boundaries
    geometry=DiskGeometry(20,20)
    result=detect_reciprocal_lattice_seed(image,geometry,(155,115))
    assert_same_lattice(result,np.diag([65.,60.]),atol=.03)
    assert result.fit_rms_residual_px<1.5
    assert any(not a.inlier for a in result.assignments)
    assert result==detect_reciprocal_lattice_seed(image,geometry,(155,115))


def test_real_scale_geometry_38_radius_85_spacing():
    image=arc_scene(((70,90),(155,90),(240,90)),(.3,.4,.3),radius=(38,38),noise=.1)
    geometry=DiskGeometry(38,38)
    candidates=detect_lattice_candidates(image,geometry,(155,190))
    nearest=[min(candidates,key=lambda c:np.linalg.norm(np.array(c.center_px)-target))
             for target in ((70,90),(155,90),(240,90))]
    assert len({c.id for c in nearest})==3
    np.testing.assert_allclose([c.center_px for c in nearest],[(70,90),(155,90),(240,90)],atol=3)


def test_ellipse_normal_vote_formula_and_no_wrapped_votes():
    from qcbed.data_preparation import arc_center_votes
    image=arc_scene(((55,100),),(1.,),radius=(30,18))
    votes,arrays,_=arc_center_votes(image,DiskGeometry(30,18),(250,180),
                                  config=LatticeDetectionConfig(radius_tolerance_px=0))
    y,x=np.unravel_index(votes.argmax(),votes.shape)
    assert abs(x-55)<=1 and abs(y-100)<=1
    # No edge can cast a geometric vote on the far side of this frame.
    assert votes[:,200:].max()<1e-10


@pytest.mark.parametrize('feature',['straight','parallel','curved'])
def test_internal_feature_alone_cannot_pass_arc_validation(feature):
    y,x=np.indices((230,310));cx,cy=90,100
    if feature=='curved':
        # A closed but wrong-radius fringe must not be interpreted as a physical R=20 disk.
        image=200/(1+np.exp(np.clip((np.hypot(x-cx,y-cy)-9)*3,-700,700)))
    else:
        positions=[cy] if feature=='straight' else [cy-8,cy,cy+8]
        image=sum(200*np.exp(-(y-v)**2/2) for v in positions)
    assert detect_lattice_candidates(image,DiskGeometry(20,20),(250,180))==()


@pytest.mark.parametrize('fraction',[.2,.3,.4])
@pytest.mark.parametrize('feature',['straight','parallel','curved'])
def test_partial_outer_arc_survives_stronger_internal_features(fraction,feature):
    image=arc_scene(fractions=(fraction,),noise=.05)
    y,x=np.indices(image.shape);cx,cy=90,100
    if feature=='curved':
        image+=200/(1+np.exp(np.clip((np.hypot(x-cx,y-cy)-9)*3,-700,700)))
    else:
        positions=[cy] if feature=='straight' else [cy-6,cy,cy+6]
        for v in positions:
            image+=200*np.exp(-(y-v)**2/2)*np.exp(-((x-cx)/12)**8)
    cs=detect_lattice_candidates(image,DiskGeometry(20,20),(250,180))
    closest=min(cs,key=lambda c:np.linalg.norm(np.array(c.center_px)-(cx,cy)))
    assert np.linalg.norm(np.array(closest.center_px)-(cx,cy))<=3
    points=np.array(closest.supporting_pixels)
    radii=np.linalg.norm(points-(cx,cy),axis=1)
    assert radii.min()>16  # supporting pixels belong to the outer arc, not internal fringes
    assert closest.center_vote_rms_px<=2.5
    assert closest.mean_normal_residual_deg<=12
    assert closest.arc_length_px>=8
    assert closest.angular_coverage>=.15
    assert closest.confidence>0


def test_arc_polarity_and_support_are_geometric():
    image=arc_scene(fractions=(.3,),interior=True,noise=.05)
    g=DiskGeometry(20,20)
    a=detect_lattice_candidates(image,g,(250,180))
    b=detect_lattice_candidates(-image,g,(250,180))
    assert [c.center_px for c in a]==[c.center_px for c in b]
    assert [c.supporting_pixels for c in a]==[c.supporting_pixels for c in b]
    np.testing.assert_allclose([c.confidence for c in a],[c.confidence for c in b],atol=1e-10)


@pytest.mark.parametrize('kwargs',[dict(center_vote_rms_scale_px=0),dict(minimum_arc_length_px=-1),
    dict(maximum_normal_residual_deg=90),dict(maximum_neighbor_normal_change_deg=True)])
def test_arc_validation_config_errors(kwargs):
    with pytest.raises(ValueError): LatticeDetectionConfig(**kwargs)


def test_tiny_arc_cannot_replace_minimum_geometric_support():
    image=arc_scene(fractions=(.03,))
    assert detect_lattice_candidates(image,DiskGeometry(20,20),(250,180))==()


def disconnected_radius_support(radius=38., center=(110.,110.)):
    # Eight short, disjoint fragments; none meets the former per-fragment count.
    phase=(np.arange(8)[:,None]*np.pi/4+np.arange(-3,4)[None,:]/radius).ravel()
    points=np.unique(np.rint(np.array(center)+radius*np.column_stack((np.cos(phase),np.sin(phase)))),axis=0)
    normals=points-center
    normals/=np.linalg.norm(normals,axis=1)[:,None]
    return dict(edge_x=points[:,0].astype(int),edge_y=points[:,1].astype(int),
                normal_x=normals[:,0],normal_y=normals[:,1])


def test_disconnected_short_fragments_jointly_validate_known_radius_center():
    from qcbed.data_preparation.arc_geometry import _validated_center,_coherent_fragments
    arrays=disconnected_radius_support();cfg=LatticeDetectionConfig()
    points=np.column_stack((arrays['edge_x'],arrays['edge_y']))
    normals=np.column_stack((arrays['normal_x'],arrays['normal_y']))
    fragments=_coherent_fragments(points,normals,cfg)
    assert len(fragments)==8
    assert all(len(f)<cfg.minimum_support_pixels for f in fragments)
    result=_validated_center(np.array([110.,110.]),arrays,DiskGeometry(38,38),cfg)
    assert result is not None and result['fragment_count']==8
    assert result['supporting_edge_pixels']==len(points)
    assert .15<result['angular_coverage']<.3  # absent gaps are not counted as supported arcs
    assert result['radius_residual_rms_px']<1
    assert result['radius_offset_px']==0
    assert result['arc_length_px']==pytest.approx(2*np.pi*38*result['angular_coverage'])
    opposite=dict(arrays,normal_x=-arrays['normal_x'],normal_y=-arrays['normal_y'])
    assert _validated_center(np.array([110.,110.]),opposite,DiskGeometry(38,38),cfg)==result


def test_disconnected_evidence_must_share_nominal_radius_and_normals():
    from qcbed.data_preparation.lattice_detection import _validated_center
    cfg=LatticeDetectionConfig()
    for radius in (33.,35.,41.,43.):
        assert _validated_center(np.array([110.,110.]),disconnected_radius_support(radius),DiskGeometry(38,38),cfg) is None
    arrays=disconnected_radius_support()
    wrong=dict(arrays,normal_x=-arrays['normal_y'],normal_y=arrays['normal_x'])
    assert _validated_center(np.array([110.,110.]),wrong,DiskGeometry(38,38),cfg) is None
    assert _validated_center(np.array([120.,120.]),arrays,DiskGeometry(38,38),cfg) is None


def test_image_with_only_disconnected_short_outer_arcs():
    y,x=np.indices((256,300));dx,dy=x-100.,y-110.
    angle=np.arctan2(dy,dx);radius=np.hypot(dx,dy)
    local_phase=(angle+np.pi/8)%(np.pi/4)-np.pi/8
    visible=1/(1+np.exp(np.clip((abs(local_phase)-.14)*200,-700,700)))
    image=50*visible/(1+np.exp(np.clip((radius-38)*3,-700,700)))
    candidates=detect_lattice_candidates(image,DiskGeometry(38,38),(250,210))
    closest=min(candidates,key=lambda c:np.linalg.norm(np.array(c.center_px)-(100,110)))
    assert np.linalg.norm(np.array(closest.center_px)-(100,110))<=2
    assert closest.fragment_count>1
    assert closest.angular_coverage>=.15


def two_validated_centers():
    image=arc_scene(((70,80),(155,80)),(1.,1.))
    geometry=DiskGeometry(20,20);center=(100.,160.)
    candidates=detect_lattice_candidates(image,geometry,center)
    assert len(candidates)==2
    return image,geometry,center,candidates


def move_arc(candidate,center):
    from dataclasses import replace
    shift=np.array(center)-candidate.center_px
    return replace(candidate,center_px=center,supporting_pixels=tuple(
        tuple(int(v) for v in np.array(point)+shift) for point in candidate.supporting_pixels))


def test_two_validated_arc_centers_minimal_seed_and_provenance(monkeypatch):
    import qcbed.data_preparation.lattice_detection as module
    image,g,center,candidates=two_validated_centers()
    def forbidden(*args,**kwargs): raise AssertionError('Minimal path must not vote on displacement families')
    monkeypatch.setattr(module,'_families',forbidden)
    seed=infer_lattice_seed(candidates,center,g)
    assert seed==infer_lattice_seed(candidates[::-1],center,g)
    assert seed.seed_method=='two_validated_arc_centers'
    assert seed.source_candidates==candidates
    assert seed.initial_basis_vectors_px==tuple(tuple(np.array(c.center_px)-center) for c in candidates)
    assert seed.displacement_families==()
    assert not seed.primitive_basis_certain
    assert len(seed.assignments)==2 and all(a.inlier for a in seed.assignments)
    assert seed.fit_rms_residual_px<1e-10
    assert seed.basis_condition_number<1e6
    fitted=np.column_stack((seed.basis_vector_1_px,seed.basis_vector_2_px))
    np.testing.assert_allclose(fitted,module._canonical(np.array(seed.initial_basis_vectors_px).T))


@pytest.mark.parametrize('bad', ['unvalidated','low_confidence','bad_normal','bad_radius','short'])
def test_minimal_seed_requires_arc_confidence_and_plausible_geometry(bad):
    from dataclasses import replace
    _,g,center,c=two_validated_centers()
    first=c[0]
    if bad=='unvalidated': first=LatticeCandidate(first.id,first.center_px,first.response_score)
    if bad=='low_confidence': first=replace(first,confidence=0)
    if bad=='bad_normal': first=replace(first,mean_normal_residual_deg=30)
    if bad=='bad_radius': first=replace(first,center_px=(first.center_px[0]+10,first.center_px[1]))
    if bad=='short': first=move_arc(first,(105.,160.))
    with pytest.raises(ValueError,match='Minimal'):
        infer_lattice_seed((first,c[1]),center,g)


def test_minimal_seed_collinear_and_ill_conditioned_rejected():
    _,g,center,c=two_validated_centers()
    with pytest.raises(ValueError,match='non-collinear'):
        infer_lattice_seed((move_arc(c[0],(200.,160.)),move_arc(c[1],(300.,160.))),center,g)
    with pytest.raises(ValueError,match='ill-conditioned'):
        infer_lattice_seed((move_arc(c[0],(200000100.,160.)),move_arc(c[1],(200000100.,161.))),center,g,
                           config=LatticeDetectionConfig(min_basis_sine=1e-10))


def test_multi_candidate_seed_keeps_displacement_path():
    image,center,g,_,_=scene()
    seed=detect_reciprocal_lattice_seed(image,g,center)
    assert seed.seed_method=='multi_candidate_displacement_voting'
    assert seed.displacement_families and len(seed.assignments)>2


def test_minimal_seed_is_not_enough_for_unchanged_a7_global_fit():
    from qcbed.data_preparation import refine_lattice_and_generate_catalogue,LatticeRefinementConfig
    image,g,center,c=two_validated_centers()
    seed=infer_lattice_seed(c,center,g)
    # Two actual observations cannot determine A7's free-origin affine fit.
    with pytest.raises(ValueError,match='Too few candidates'):
        refine_lattice_and_generate_catalogue(image,seed,g)
    with pytest.raises(ValueError,match='min_inliers must be >=3'):
        LatticeRefinementConfig(min_inliers=2)
