from dataclasses import FrozenInstanceError
from itertools import product
import numpy as np
import pytest
from qcbed.crystallography import (DirectLattice,ZoneAxis,OriginalReflectionBasis,
    build_zone_axis_cell,index_experimental_beam,index_experimental_beams,
    original_reflection_basis_vectors)

INDICES=((0,0),(-1,0),(1,0),(0,-1),(0,1),(-1,-1),(1,1),(-1,1),(1,-1))


def context(cell=None,zone=(1,1,0),basis=((1,-1,0),(0,0,2))):
    d=DirectLattice(np.eye(3)*4 if cell is None else cell)
    z=ZoneAxis(zone)
    return OriginalReflectionBasis(*basis,z),build_zone_axis_cell(d,z)


def test_nine_order_original_internal_and_central(monkeypatch):
    from qcbed.crystallography import ReciprocalLattice
    monkeypatch.setattr(ReciprocalLattice,'d_spacing',lambda *a:pytest.fail('indexing must not request d-spacing'))
    basis,cell=context()
    inputs=[(f'beam{i}',n) for i,n in enumerate(INDICES)]
    result=index_experimental_beams(inputs,basis,cell)
    assert [r.stable_id for r in result]==[r[0] for r in inputs]
    assert [r.experimental_index for r in result]==list(INDICES)
    assert [r.original_hkl for r in result]==[(0,0,0),(-1,1,0),(1,-1,0),(0,0,-2),(0,0,2),(-1,1,-2),(1,-1,2),(-1,1,2),(1,-1,-2)]
    assert [r.internal_hk.as_HKL() for r in result]==[(0,0,0),(-1,0,0),(1,0,0),(0,2,0),(0,-2,0),(-1,2,0),(1,-2,0),(-1,-2,0),(1,2,0)]
    with pytest.raises(FrozenInstanceError):result[0].original_hkl=(1,1,1)
    assert result==index_experimental_beams(inputs,basis,cell)


@pytest.mark.parametrize('cell,zone,refs',[
    (np.eye(3)*4,(1,1,0),((1,-1,0),(0,0,2))),
    (np.diag([3.,4.,5.]),(1,1,0),((1,-1,0),(0,0,2))),
    (np.array([[3.,0,0],[-1.5,3*np.sqrt(3)/2,0],[0,0,5.]]),(0,0,1),((1,0,0),(0,1,0))),
    (np.array([[3.,0,0],[1.,4.,0],[.5,-.7,5.]]),(1,1,1),((1,-1,0),(0,1,-1))),
])
def test_metric_independent_indexing_physical_consistency(cell,zone,refs):
    basis,q=context(cell,zone,refs);b=q.original_lattice.reciprocal();g=original_reflection_basis_vectors(b,basis)
    for n in product(range(-3,4),repeat=2):
        r=index_experimental_beam(stable_id='id',experimental_index=n,reflection_basis=basis,zone_cell=q)
        expected=tuple(n[0]*a+n[1]*b for a,b in zip(*refs))
        assert r.original_hkl==expected
        np.testing.assert_allclose(b.g_vector(r.original_hkl),np.array(n)@g,atol=5e-16)
        np.testing.assert_allclose(b.g_vector(r.original_hkl),q.reciprocal_lattice.g_vector(r.internal_hk.as_HKL()),atol=1e-15)


def test_no_automatic_sign_or_order_correction():
    basis,q=context()
    for refs,expected in (((basis.hkl_2,basis.hkl_1),(0,0,2)),
                          ((tuple(-x for x in basis.hkl_1),basis.hkl_2),(-1,1,0))):
        changed=OriginalReflectionBasis(*refs,basis.zone_axis)
        r=index_experimental_beam(stable_id='id',experimental_index=(1,0),reflection_basis=changed,zone_cell=q)
        assert r.original_hkl==expected


@pytest.mark.parametrize('index',[(1.,0),(True,0),(1,2,3),(1,),None])
def test_invalid_experimental_index(index):
    b,q=context()
    with pytest.raises(ValueError):index_experimental_beam(stable_id='id',experimental_index=index,reflection_basis=b,zone_cell=q)


def test_context_ids_duplicates_and_empty():
    b,q=context()
    wrong=build_zone_axis_cell(q.original_lattice,ZoneAxis((0,0,1)))
    with pytest.raises(ValueError,match='zone_axis'):index_experimental_beams([],b,wrong)
    with pytest.raises(ValueError,match='unique'):index_experimental_beams([('id',(1,0)),('id',(0,1))],b,q)
    with pytest.raises(ValueError,match='nonempty'):index_experimental_beams([('',(1,0))],b,q)
    result=index_experimental_beams([('a',(1,0)),('b',(1,0))],b,q)
    assert len(result)==2 and result[0].original_hkl==result[1].original_hkl
    assert index_experimental_beams([],b,q)==()


@pytest.mark.parametrize('refs',[((1,0,0),(0,0,2)),((0,0,0),(0,0,2)),((1,-1,0),(3,-3,0)),((1.,-1,0),(0,0,2)),((1,-1,0,0),(0,0,2))])
def test_invalid_original_basis_reused_validation(refs):
    with pytest.raises(ValueError):context(basis=refs)
