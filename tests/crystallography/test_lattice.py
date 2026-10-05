import numpy as np
import pytest
from qcbed.crystallography import DirectLattice, ReciprocalLattice


def test_cubic():
    direct=DirectLattice(np.eye(3)*4)
    reciprocal=direct.reciprocal()
    np.testing.assert_array_equal(reciprocal.vectors,np.eye(3)*.25)
    assert direct.volume==pytest.approx(64)
    np.testing.assert_array_equal(direct.metric_tensor,np.eye(3)*16)
    np.testing.assert_array_equal(direct.lengths,[4,4,4])
    np.testing.assert_array_equal(direct.angles,[90,90,90])
    for hkl,d in (((1,0,0),4),((1,1,0),4/np.sqrt(2)),((1,1,1),4/np.sqrt(3))):
        assert reciprocal.d_spacing(hkl)==pytest.approx(d)
    np.testing.assert_array_equal(reciprocal.g_vector((-1,2,0)),[-.25,.5,0])
    np.testing.assert_array_equal(reciprocal.g_vector((0,0,0)),[0,0,0])


def test_orthorhombic_and_triclinic():
    np.testing.assert_array_equal(DirectLattice(np.diag([3,4,5])).reciprocal().vectors,np.diag([1/3,1/4,1/5]))
    a=np.array([[3.,0,0],[1.,4.,0],[.5,-.7,5.]])
    d=DirectLattice(a);b=d.reciprocal()
    np.testing.assert_allclose(b.vectors@d.vectors.T,np.eye(3),atol=2e-16)
    # Independent cross-product oracle verifies inverse-transpose convention.
    expected=np.array([np.cross(a[1],a[2]),np.cross(a[2],a[0]),np.cross(a[0],a[1])])/60
    np.testing.assert_allclose(b.vectors,expected,rtol=1e-15,atol=1e-16)
    np.testing.assert_allclose(d.metric_tensor,a@a.T)
    np.testing.assert_allclose(d.angles,[np.degrees(np.arccos(np.dot(a[i],a[j])/np.linalg.norm(a[i])/np.linalg.norm(a[j]))) for i,j in ((1,2),(0,2),(0,1))])


@pytest.mark.parametrize('bad',[np.eye(2),np.ones((3,3)),np.zeros((3,3)),np.diag([-1,1,1]),np.full((3,3),np.nan),np.full((3,3),np.inf), [['1']*3]*3,np.eye(3,dtype=complex)])
@pytest.mark.parametrize('cls',[DirectLattice,ReciprocalLattice])
def test_invalid(bad,cls):
    with pytest.raises((ValueError,TypeError)): cls(bad)


def test_owned_immutable_arrays():
    source=np.eye(3);d=DirectLattice(source);source[0,0]=7
    assert d.vectors[0,0]==1
    for a in (d.vectors,d.lengths,d.angles,d.metric_tensor,d.reciprocal().vectors,d.reciprocal().g_vector((1,0,0))):
        with pytest.raises(ValueError): a.setflags(write=True)


@pytest.mark.parametrize('hkl',[(0,0,0),(1,0),(1.5,0,0),(True,0,0)])
def test_invalid_reflection(hkl):
    with pytest.raises(ValueError): DirectLattice(np.eye(3)).reciprocal().d_spacing(hkl)
