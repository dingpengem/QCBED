"""Independent local boundary evidence, never prediction-as-observation."""
import csv
import numpy as np
import pytest
from qcbed.data_preparation import (DiskGeometry, detect_reciprocal_lattice_seed, BootstrapConfig,
    bootstrap_lattice, LatticeRefinementConfig)
from qcbed.data_preparation.bootstrap_diagnostics import write_bootstrap_diagnostics


def image(centers):
    y,x=np.indices((240,240));result=np.zeros(x.shape)
    for cx,cy in centers:
        result+=50/(1+np.exp(np.clip((np.hypot(x-cx,y-cy)-20)*3,-700,700)))
    return result


@pytest.fixture
def problem():
    g=DiskGeometry(20,20)
    source=image(((180,120),(120,180)))
    seed=detect_reciprocal_lattice_seed(source,g,(120,120))
    assert seed.seed_method=='two_validated_arc_centers'
    return g,source,seed


def test_predictions_without_global_fit_or_whole_image_redetection(problem,monkeypatch):
    import qcbed.data_preparation.lattice_bootstrap as module
    import qcbed.data_preparation.lattice_detection as detection
    g,source,seed=problem
    def forbidden(*a,**kw): raise AssertionError('Must not fit or redetect')
    monkeypatch.setattr(module,'refine_reciprocal_lattice',forbidden)
    monkeypatch.setattr(detection,'detect_lattice_candidates',forbidden)
    monkeypatch.setattr(detection,'arc_center_votes',forbidden)
    before=source.copy()
    result=bootstrap_lattice(source,seed,g)
    assert result.status=='insufficient_observations'
    assert result.newly_confirmed_count==0
    assert len(result.expanded_seed.assignments)==2
    assert len(result.sites)==25
    assert result.sites[0].integer_indices==(0,0)
    assert {s.predicted_center_px for s in result.sites}=={(x,y) for x in (0.,60.,120.,180.,240.) for y in (0.,60.,120.,180.,240.)}
    assert sum(s.status=='defining' for s in result.sites)==2
    assert all(s.observed_center_px is None for s in result.sites if s.status in ('missing','reference'))
    assert result==bootstrap_lattice(source,seed,g)
    np.testing.assert_array_equal(source,before)


def test_independent_third_fourth_and_fit_trigger(problem,monkeypatch):
    import qcbed.data_preparation.lattice_bootstrap as module
    g,_,seed=problem
    original=module.refine_reciprocal_lattice;calls=[]
    def record(expanded,**kwargs):
        calls.append(expanded)
        assert len(expanded.assignments)>=4
        return original(expanded,**kwargs)
    monkeypatch.setattr(module,'refine_reciprocal_lattice',record)
    three=bootstrap_lattice(image(((180,120),(120,180),(62,120))),seed,g)
    assert three.newly_confirmed_count==1 and not calls
    result=bootstrap_lattice(image(((180,120),(120,180),(62,120),(120,62))),seed,g)
    assert result.newly_confirmed_count==2 and len(calls)==1
    assert result.refined_lattice is not None
    assert result.refined_lattice.config.min_inliers==LatticeRefinementConfig().min_inliers==4
    confirmed=[s for s in result.sites if s.status=='confirmed']
    assert all(s.observed_center_px!=s.predicted_center_px and 0<s.shift_px<4 for s in confirmed)
    assert all(s.source=='a7_bootstrap_confirmed' for s in confirmed)
    assert len({a.candidate.id for a in result.expanded_seed.assignments})==4
    assert len({a.integer_indices for a in result.expanded_seed.assignments})==4


def test_bright_internal_features_do_not_confirm_missing_site(problem):
    g,source,seed=problem;y,x=np.indices(source.shape)
    p=source+1000*np.exp(-((x-60)**2+(y-120)**2)/8)
    result=bootstrap_lattice(p,seed,g)
    site=next(s for s in result.sites if s.predicted_center_px==(60.,120.))
    assert site.status=='missing' and site.observed_center_px is None
    assert result.newly_confirmed_count==0


def test_boundary_evidence_is_not_independent_confirmation(problem):
    g,_,seed=problem
    result=bootstrap_lattice(image(((180,120),(120,180),(64,120))),seed,g)
    site=next(s for s in result.sites if s.predicted_center_px==(60.,120.))
    assert site.status=='weak_evidence' and site.observed_center_px is None
    assert site.candidate is not None
    assert result.newly_confirmed_count==0


def test_bootstrap_artifacts_for_incomplete_result(problem,tmp_path):
    g,source,seed=problem
    result=bootstrap_lattice(source,seed,g)
    before=source.copy()
    files=write_bootstrap_diagnostics(tmp_path,source,result,g)
    assert (tmp_path/'bootstrap_overview.png').is_file()
    assert len(list((tmp_path/'bootstrap_beams').glob('*.png')))==len(result.sites)
    with open(files['bootstrap_table']) as f: rows=list(csv.DictReader(f))
    assert len(rows)==len(result.sites)
    assert all(r['observed_x']=='' for r in rows if r['status']=='missing')
    assert sum(r['source']=='a6_defining' for r in rows)==2
    np.testing.assert_array_equal(source,before)


@pytest.mark.parametrize('kwargs',[dict(search_radius_px=0),dict(search_step_px=True),dict(search_step_px=5),dict(search_step_px=.0001)])
def test_invalid_bootstrap_config(kwargs):
    with pytest.raises(ValueError): BootstrapConfig(**kwargs)


def test_search_cannot_jump_to_neighbor(problem):
    g,source,seed=problem
    with pytest.raises(ValueError,match='one quarter'):
        bootstrap_lattice(source,seed,g,config=BootstrapConfig(search_radius_px=30))


def test_shared_enumeration_order_and_coverage(problem):
    from dataclasses import replace
    from qcbed.data_preparation import enumerate_lattice_sites,enumerate_seed_sites
    g,_,seed=problem
    result=bootstrap_lattice(image(((180,120),(120,180),(60,120),(120,60))),seed,g)
    refined=result.refined_lattice
    assert refined is not None
    same_geometry=replace(refined,origin_px=seed.reference_center_px,
        basis_vector_1_px=seed.basis_vector_1_px,basis_vector_2_px=seed.basis_vector_2_px)
    assert enumerate_lattice_sites(same_geometry,(240,240),g)==enumerate_seed_sites(seed,(240,240),g)


def test_diagnostics_preserves_predictions_when_bootstrap_incomplete(problem,tmp_path,monkeypatch):
    from types import SimpleNamespace
    import qcbed.data_preparation.diagnostics as module
    g,source,seed=problem
    path=tmp_path/'pattern.dat';np.savetxt(path,source)
    monkeypatch.setattr(module,'locate_central_disk',lambda *args:SimpleNamespace(central_disk_center_px=seed.reference_center_px))
    result=module.validate_pattern(path,tmp_path/'output',geometry=g,geometry_provenance='synthetic physical radius')
    assert result['status']=='bootstrap_only'
    assert result['a7'] is None
    assert result['bootstrap']['newly_confirmed_count']==0
    assert result['bootstrap']['required_observations']==4
    assert (tmp_path/'output'/'bootstrap_sites.csv').is_file()
    assert len(list((tmp_path/'output'/'bootstrap_beams').glob('*.png')))==25


def test_higher_order_partially_clipped_fourth_observation(problem):
    from qcbed.data_preparation.arc_catalogue import generate_arc_catalogue
    g,_,seed=problem
    source=image(((180,120),(120,180),(60,120),(0,120)))
    result=bootstrap_lattice(source,seed,g)
    edge=next(s for s in result.sites if s.predicted_center_px==(0.,120.))
    assert max(abs(n) for n in edge.integer_indices)==2
    assert edge.status=='confirmed'
    assert .3<edge.usable_boundary_fraction<.5
    assert edge.boundary_fraction==pytest.approx(.5)
    assert not edge.boundary_limited
    assert result.newly_confirmed_count==2
    assert result.refined_lattice is not None
    assert result.refined_lattice.fit_rms_residual_px<1e-10
    assert len(result.refined_lattice.assignments)==4
    catalogue=generate_arc_catalogue(result,source.shape,g)
    assert catalogue.refined and catalogue.status_counts['observed']==4
    assert catalogue==generate_arc_catalogue(result,source.shape,g)
    assert len({s.beam_id for s in catalogue.sites})==len(catalogue.sites)
    for site in catalogue.sites:
        if site.status=='observed':
            assert site.candidate is not None and site.observed_center_px==site.candidate.center_px
            assert site.initial_integer_indices is not None
            assert site.lattice_residual_px<1e-10
        else: assert site.observed_center_px is None


def test_visible_circumference_rejection_and_boundary_flag(problem):
    g,_,seed=problem
    result=bootstrap_lattice(image(((180,120),(120,180),(0,120),(64,120))),seed,g,
                            config=BootstrapConfig(minimum_boundary_fraction=.75))
    clipped=next(s for s in result.sites if s.predicted_center_px==(0.,120.))
    assert clipped.status=='edge_clipped' and clipped.observed_center_px is None
    assert 'circumference' in clipped.reason
    limited=next(s for s in result.sites if s.predicted_center_px==(60.,120.))
    assert limited.boundary_limited and limited.status=='weak_evidence'
    assert limited.observed_center_px is None
    assert result.refined_lattice is None


def test_partial_arc_catalogue_does_not_invent_observations(problem):
    from qcbed.data_preparation.arc_catalogue import generate_arc_catalogue
    g,source,seed=problem
    result=bootstrap_lattice(source,seed,g)
    catalogue=generate_arc_catalogue(result,source.shape,g)
    assert not catalogue.refined and catalogue.status_counts['observed']==2
    assert all(s.lattice_residual_px is None for s in catalogue.sites)
    assert sum(s.observed_center_px is not None for s in catalogue.sites)==2


def test_duplicate_seed_observations_rejected(problem):
    from dataclasses import replace
    g,source,seed=problem
    duplicate=replace(seed,assignments=seed.assignments+(seed.assignments[0],))
    with pytest.raises(ValueError,match='Duplicate'):
        bootstrap_lattice(source,duplicate,g)


def test_visible_ellipse_fraction_and_enumeration_budget():
    from qcbed.data_preparation.boundary_sites import boundary_fraction,enumerate_boundary_sites
    g=DiskGeometry(20,10)
    assert boundary_fraction((0,60),(121,121),g)==pytest.approx(.5)
    assert boundary_fraction((60,60),(121,121),g)==pytest.approx(1)
    assert boundary_fraction((-21,60),(121,121),g)==0
    with pytest.raises(ValueError,match='Enumeration'):
        enumerate_boundary_sites((60,60),np.eye(2),(121,121),g,max_sites=2)


@pytest.mark.parametrize('value',[0,1.1,float('nan'),True])
def test_invalid_boundary_fraction(value):
    with pytest.raises(ValueError): BootstrapConfig(minimum_boundary_fraction=value)


def test_refinement_comparison_and_catalogue_artifacts(problem,tmp_path):
    from qcbed.data_preparation.arc_catalogue import generate_arc_catalogue,refinement_comparison
    g,_,seed=problem
    source=image(((180,120),(120,180),(60,120),(0,120)))
    result=bootstrap_lattice(source,seed,g)
    comparison=refinement_comparison(result)
    assert comparison['rms_after']<1e-10
    np.testing.assert_allclose(comparison['origin_change'],(0,0),atol=1e-10)
    catalogue=generate_arc_catalogue(result,source.shape,g)
    paths=write_bootstrap_diagnostics(tmp_path,source,catalogue,g,catalogue=True)
    with open(paths['beam_table']) as f: rows=list(csv.DictReader(f))
    assert len(rows)==len(catalogue.sites)
    assert sum(r['status']=='observed' for r in rows)==4
    assert all(r['observed_x']=='' for r in rows if r['status']!='observed')
    assert (tmp_path/'beam_catalogue_overview.png').exists()
    assert len(list((tmp_path/'beams').glob('*.png')))==len(rows)
    absent=refinement_comparison(bootstrap_lattice(image(((180,120),(120,180))),seed,g))
    assert all(v is None for v in absent.values())
