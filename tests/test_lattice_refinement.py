"""Synthetic A7 affine lattice/catalogue tests, separate from physical HKL."""
from dataclasses import replace, FrozenInstanceError
import numpy as np
import pytest
from qcbed.data_preparation import (DiskGeometry, LatticeCandidate, LatticeAssignment,
    LatticeDetectionConfig, ReciprocalLatticeSeed, LatticeRefinementConfig, BeamStatus,
    refine_reciprocal_lattice, enumerate_lattice_sites, generate_beam_catalogue,
    refine_lattice_and_generate_catalogue)

ORIGIN = np.array([119., 91.])
BASIS = np.array([[30., 8.], [3., 25.]])
GEOMETRY = DiskGeometry(5,4)


def seed(noise=0., outlier=False, missing=False, duplicate=False, wrong=False):
    rng=np.random.default_rng(5); assignments=[]
    for i in range(-2,3):
        for j in range(-2,3):
            if (i,j)==(0,0) or (missing and (i+j)%2==0):
                continue
            point=ORIGIN+BASIS@[i,j]+rng.normal(0,noise,2)
            n=(i,j)
            if wrong and (i,j)==(1,1):
                n=(2,1)
            c=LatticeCandidate(f'p{i}_{j}',tuple(point),10.)
            assignments.append(LatticeAssignment(c,n,0.,not(wrong and (i,j)==(1,1))))
    if outlier:
        point=ORIGIN+BASIS@[1,1]+[4,-3]
        assignments.append(LatticeAssignment(LatticeCandidate('bad',tuple(point),100.),(1,1),0.,True))
    if duplicate:
        point=ORIGIN+BASIS@[1,0]+[.2,.2]
        assignments.append(LatticeAssignment(LatticeCandidate('duplicate',tuple(point),20.),(1,0),0.,True))
    return ReciprocalLatticeSeed(tuple(ORIGIN+[.6,-.4]),tuple(BASIS[:,0]+[.1,-.1]),
        tuple(BASIS[:,1]+[-.1,.1]),tuple(assignments),(),1.,0.,LatticeDetectionConfig())


def matrix(lattice):
    return np.column_stack((lattice.basis_vector_1_px,lattice.basis_vector_2_px))


def assert_equivalent(lattice,tolerance):
    np.testing.assert_allclose(lattice.origin_px,ORIGIN,atol=tolerance)
    transformation=np.linalg.solve(BASIS,matrix(lattice))
    np.testing.assert_allclose(transformation,np.rint(transformation),atol=tolerance)
    assert round(np.linalg.det(transformation))==1


def test_exact_global_fit_and_canonical_labels():
    s=seed(); lattice=refine_reciprocal_lattice(s)
    assert_equivalent(lattice,1e-12)
    assert lattice.fit_rms_residual_px < 1e-12
    assert lattice.initial_rms_residual_px > .5
    for a in lattice.assignments:
        np.testing.assert_allclose(np.array(lattice.origin_px)+matrix(lattice)@a.integer_indices,
                                   a.candidate.center_px,atol=1e-12)
    assert lattice.initial_seed is s
    assert lattice.origin_px != s.reference_center_px
    assert not hasattr(lattice,'zone_axis_origin_px')


def test_noisy_centers_outlier_and_duplicate_preserved():
    s=seed(noise=.12,duplicate=True,outlier=True)
    result=refine_reciprocal_lattice(s)
    assert_equivalent(result,.1)
    assert result.fit_rms_residual_px < result.initial_rms_residual_px
    assert len(result.assignments)==len(s.assignments)
    assert not next(a for a in result.assignments if a.candidate.id=='bad').inlier
    sites=[a.integer_indices for a in result.assignments if a.inlier]
    assert len(sites)==len(set(sites))


def test_reassignment_of_initial_outlier():
    result=refine_reciprocal_lattice(seed(wrong=True))
    a=next(a for a in result.assignments if a.candidate.id=='p1_1')
    assert a.inlier and a.reassigned and a.residual_px < 1e-12


def test_missing_reflections_do_not_remove_predicted_sites():
    full=refine_reciprocal_lattice(seed())
    sparse=refine_reciprocal_lattice(seed(missing=True))
    a=enumerate_lattice_sites(full,(190,250),GEOMETRY)
    b=enumerate_lattice_sites(sparse,(190,250),GEOMETRY)
    assert {n for n,p in a}=={n for n,p in b}
    assert len(b)>len(sparse.assignments)
    assert a[0][0]==(0,0)
    assert len({n for n,p in a})==len(a)
    for n,(x,y) in a:
        assert 5-1e-10 <= x <= 244+1e-10 and 4-1e-10 <= y <= 185+1e-10
    # Independent brute enumeration verifies the oblique corner-box enumeration.
    expected=set()
    for i in range(-20,21):
        for j in range(-20,21):
            x,y=np.array(full.origin_px)+matrix(full)@[i,j]
            if 5-1e-10<=x<=244+1e-10 and 4-1e-10<=y<=185+1e-10:
                expected.add((i,j))
    assert {n for n,p in a}==expected


def image_and_catalogue(lattice,mask=None):
    y,x=np.indices((190,250)); image=np.zeros((190,250))
    for n,p in enumerate_lattice_sites(lattice,image.shape,GEOMETRY):
        amplitude=20. if n==(0,0) else .6 if n==(1,0) else 0.
        image += amplitude*np.exp(-2*(((x-p[0])/5)**2+((y-p[1])/4)**2))
    before=image.copy()
    catalogue=generate_beam_catalogue(image,lattice,GEOMETRY,valid_mask=mask)
    np.testing.assert_array_equal(image,before)
    return image,catalogue


def test_catalogue_classification_and_predicted_observed_separation():
    lattice=refine_reciprocal_lattice(seed())
    image,catalogue=image_and_catalogue(lattice)
    beams={b.lattice_indices:b for b in catalogue.beams}
    assert beams[(0,0)].status==BeamStatus.OBSERVED
    assert beams[(1,0)].status==BeamStatus.WEAK
    assert beams[(-1,0)].status==BeamStatus.MISSING
    assert beams[(1,0)].observed_center_px is None
    assert sum(b.is_central for b in catalogue.beams)==1
    assert sum(catalogue.status_counts.values())==len(catalogue.beams)
    for b in catalogue.beams:
        assert b.observation.center_px==b.predicted_center_px
        assert b.observation.indices is None
    assert catalogue.observations==tuple(b.observation for b in catalogue.beams)
    with pytest.raises(FrozenInstanceError):
        catalogue.beams=()


def test_biased_intensity_cannot_jump_and_mask_rejection():
    lattice=refine_reciprocal_lattice(seed()); y,x=np.indices((190,250))
    origin=np.array(lattice.origin_px)
    image=100*np.exp(-2*(((x-origin[0]-4)/5)**2+((y-origin[1])/4)**2))
    catalogue=generate_beam_catalogue(image,lattice,GEOMETRY)
    central=catalogue.beams[0]
    assert central.status==BeamStatus.REJECTED
    assert central.observed_center_px is None
    assert central.predicted_center_px==lattice.origin_px
    assert central.local_search_radius_px==2.
    mask=np.ones(image.shape,bool); mask[75:105,100:140]=False
    _,masked=image_and_catalogue(lattice,mask)
    assert masked.beams[0].status==BeamStatus.REJECTED
    assert masked.beams[0].mask_coverage < .95


def test_deterministic_order_and_wrapper():
    s=seed(noise=.1)
    result=refine_reciprocal_lattice(s)
    reordered=refine_reciprocal_lattice(replace(s,assignments=tuple(reversed(s.assignments))))
    assert result.origin_px==reordered.origin_px
    assert result.basis_vector_1_px==reordered.basis_vector_1_px
    assert result.assignments==reordered.assignments
    image,_=image_and_catalogue(result)
    a=refine_lattice_and_generate_catalogue(image,s,GEOMETRY)
    b=generate_beam_catalogue(image,reordered,GEOMETRY)
    assert a.beams==b.beams


@pytest.mark.parametrize('kwargs',[dict(max_iterations=0),dict(min_inliers=2),dict(local_step_px=True),
    dict(local_radius_scale=2),dict(minimum_mask_coverage=2),dict(weak_response_threshold=2),
    dict(residual_threshold_px=np.nan)])
def test_invalid_config(kwargs):
    with pytest.raises((ValueError,TypeError)):
        LatticeRefinementConfig(**kwargs)


def test_rank_budget_and_shape_failures():
    s=seed()
    with pytest.raises(ValueError,match='Too few'):
        refine_reciprocal_lattice(replace(s,assignments=s.assignments[:2]))
    with pytest.raises(ValueError,match='singular'):
        refine_reciprocal_lattice(replace(s,basis_vector_2_px=s.basis_vector_1_px))
    with pytest.raises(ValueError,match='converge'):
        refine_reciprocal_lattice(seed(wrong=True),config=LatticeRefinementConfig(max_iterations=1))
    lattice=refine_reciprocal_lattice(s,config=LatticeRefinementConfig(max_enumeration_sites=4))
    with pytest.raises(ValueError,match='truncated'):
        enumerate_lattice_sites(lattice,(190,250),GEOMETRY)
    with pytest.raises(ValueError,match='Central'):
        enumerate_lattice_sites(lattice,(20,20),GEOMETRY)


def test_geometric_outlier_trim_without_duplicate_site():
    s=seed(noise=.1)
    point=ORIGIN+BASIS@[3,1]+[5.,-4.]
    bad=LatticeAssignment(LatticeCandidate('unique_bad',tuple(point),10.),(3,1),0.,True)
    result=refine_reciprocal_lattice(replace(s,assignments=s.assignments+(bad,)))
    assert not next(a for a in result.assignments if a.candidate.id=='unique_bad').inlier
    assert result.iterations>=2
    assert_equivalent(result,.1)


def test_reliable_small_local_shift_kept_separate():
    lattice=refine_reciprocal_lattice(seed())
    y,x=np.indices((190,250)); ox,oy=lattice.origin_px
    image=30*np.exp(-2*(((x-ox-.75)/5)**2+((y-oy+.5)/4)**2))
    beam=generate_beam_catalogue(image,lattice,GEOMETRY).beams[0]
    assert beam.status==BeamStatus.OBSERVED
    assert beam.predicted_center_px==lattice.origin_px
    np.testing.assert_allclose(beam.observed_center_px,(ox+.75,oy-.5),atol=.26)
    assert np.linalg.norm(np.array(beam.observed_center_px)-beam.predicted_center_px)<beam.local_search_radius_px


def test_real_a6_seed_to_a7_catalogue():
    from qcbed.data_preparation import detect_reciprocal_lattice_seed
    y,x=np.indices((190,250)); image=np.zeros((190,250))
    for i in range(-2,3):
        for j in range(-2,3):
            if (i,j) in ((-1,1),(2,-2)):
                continue
            px,py=ORIGIN+BASIS@[i,j]
            # A6 now requires an actual boundary at the physical ellipse radius.
            radial=np.hypot((x-px)/5,(y-py)/4)
            image+=30/(1+np.exp(np.clip((radial-1)*20,-700,700)))
    initial=detect_reciprocal_lattice_seed(image,GEOMETRY,tuple(ORIGIN+[.4,-.3]))
    result=refine_lattice_and_generate_catalogue(image,initial,GEOMETRY)
    assert_equivalent(result.lattice,.02)
    assert len(result.beams)>len(initial.assignments)
    assert result.status_counts['missing']>0
    assert result.beams[0].is_central
