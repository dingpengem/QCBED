import numpy as np
import pytest
from qcbed.data_preparation import DiskGeometry
from qcbed.data_preparation.subpixel_arc import refine_arc_center, SubpixelArcConfig


def evidence(center=(50.35,40.70),axes=(20.,20.),partial=False):
    phase=np.concatenate((np.linspace(.1,.7,30),np.linspace(2.,2.7,30),np.linspace(4.,4.5,30))) if partial else np.linspace(0,2*np.pi,100,endpoint=False)
    d=np.column_stack((axes[0]*np.cos(phase),axes[1]*np.sin(phase)))
    normals=d/np.array(axes)**2
    normals/=np.linalg.norm(normals,axis=1)[:,None]
    normals[::2]*=-1
    return d+center,normals


@pytest.mark.parametrize('axes,partial',[((20.,20.),False),((20.,20.),True),((25.,14.),True)])
def test_known_subpixel_fixed_radius_recovery(axes,partial):
    points,normals=evidence(axes=axes,partial=partial)
    result=refine_arc_center((50.,41.),points,normals,DiskGeometry(*axes))
    assert result.accepted
    np.testing.assert_allclose(result.refined_center_px,(50.35,40.70),rtol=0,atol=1e-8)
    assert result.objective_after<result.objective_before
    assert result==refine_arc_center((50.,41.),points,normals,DiskGeometry(*axes))


def test_bound_and_input_immutability():
    points,normals=evidence(center=(52.,43.)); before=points.copy()
    r=refine_arc_center((50.,41.),points,normals,DiskGeometry(20,20),config=SubpixelArcConfig(bound_px=.2))
    assert max(abs(np.array(r.delta_px)))<=.2+1e-14
    assert r.bound_active
    np.testing.assert_array_equal(points,before)


def test_optimizer_failure_and_worse_proposal_keep_initial(monkeypatch):
    import qcbed.data_preparation.subpixel_arc as m
    from types import SimpleNamespace
    points,normals=evidence()
    monkeypatch.setattr(m,'least_squares',lambda *args,**kwargs:SimpleNamespace(
        x=np.array([50.9,41.9]),success=True,jac=np.eye(2),active_mask=np.zeros(2)))
    r=refine_arc_center((50.,41.),points,normals,DiskGeometry(20,20))
    assert not r.accepted and r.refined_center_px==(50.,41.)
    def fail(*args,**kwargs): raise ValueError('failed')
    monkeypatch.setattr(m,'least_squares',fail)
    assert not refine_arc_center((50.,41.),points,normals,DiskGeometry(20,20)).accepted


def test_fixed_support_ignores_internal_features():
    from qcbed.data_preparation.arc_geometry import LatticeCandidate,LatticeDetectionConfig,_validated_center
    from qcbed.data_preparation.subpixel_arc import refine_validated_candidate
    # Integer detector evidence of a subpixel disk; unrelated internal ridges
    # cannot enter the fixed support even if their implied intensity is enormous.
    p,n=evidence(center=(50.35,40.70),axes=(20,20),partial=True)
    p=np.unique(np.rint(p).astype(int),axis=0)
    normals=(p-[50.35,40.70]); normals/=np.linalg.norm(normals,axis=1)[:,None]
    a=dict(edge_x=p[:,0],edge_y=p[:,1],normal_x=normals[:,0],normal_y=normals[:,1])
    cfg=LatticeDetectionConfig()
    metrics=_validated_center(np.array([50.,41.]),a,DiskGeometry(20,20),cfg)
    c=LatticeCandidate('test',(50.,41.),1.,**metrics)
    first=refine_validated_candidate(c,a,DiskGeometry(20,20),cfg)
    augmented={k:np.concatenate((v,np.array([45,46]) if k=='edge_x' else np.array([40,41]) if k=='edge_y' else np.ones(2))) for k,v in a.items()}
    second=refine_validated_candidate(c,augmented,DiskGeometry(20,20),cfg)
    assert first==second
    assert first.subpixel_refinement.accepted


def test_central_opt_in_keeps_initial_evidence_with_strong_internal_fringe():
    from qcbed.data_preparation import locate_central_disk, CentralArcConfig
    y,x=np.indices((128,128)); center=(64.35,63.70)
    radius=np.hypot(x-center[0],y-center[1])
    image=50/(1+np.exp(np.clip((radius-22)*2,-700,700)))
    image+=2000*np.exp(-(y-center[1])**2/2)*np.exp(-((x-center[0])/10)**8)
    options=dict(initial_center_px=(64.,64.),search_config=CentralArcConfig(search_radius_px=2.))
    original=locate_central_disk(image,DiskGeometry(22,22),**options)
    refined=locate_central_disk(image,DiskGeometry(22,22),subpixel_config=SubpixelArcConfig(),**options)
    record=refined.arc_support.subpixel_refinement
    assert record.initial_center_px==original.central_disk_center_px
    assert tuple(record.support_points_px)==tuple(original.arc_support.supporting_pixels)
    assert record.objective_after<=record.objective_before
    np.testing.assert_allclose(refined.central_disk_center_px,center,atol=.3)
    assert original.arc_support.subpixel_refinement is None


def test_unchanged_validator_failure_preserves_original(monkeypatch):
    from qcbed.data_preparation.arc_geometry import LatticeCandidate, LatticeDetectionConfig
    import qcbed.data_preparation.arc_geometry as arc
    from qcbed.data_preparation.subpixel_arc import refine_validated_candidate
    p,n=evidence();p=np.unique(np.rint(p).astype(int),axis=0)
    n=p-[50.35,40.70];n/=np.linalg.norm(n,axis=1)[:,None]
    arrays=dict(edge_x=p[:,0],edge_y=p[:,1],normal_x=n[:,0],normal_y=n[:,1])
    c=LatticeCandidate('c',(50.,41.),1.,supporting_pixels=tuple(map(tuple,p)))
    monkeypatch.setattr(arc,'_validated_center',lambda *args:None)
    refined=refine_validated_candidate(c,arrays,DiskGeometry(20,20),LatticeDetectionConfig())
    assert refined.center_px==c.center_px and not refined.subpixel_refinement.accepted
    assert refined.supporting_pixels==c.supporting_pixels
