from dataclasses import FrozenInstanceError
import numpy as np
import pytest
from qcbed.crystallography import DirectLattice, ZoneAxis, CrystalOrientation, zone_axis_cartesian, zone_law_value


@pytest.mark.parametrize('uvw,expected',[((2,2,0),(1,1,0)),((-2,-2,0),(1,1,0)),((0,-2,2),(0,1,-1)),((0,0,-4),(0,0,1)),((2,1,1),(2,1,1)),((1,-1,0),(1,-1,0))])
def test_canonical_axis(uvw,expected):
    z=ZoneAxis(uvw);assert z.uvw==expected
    with pytest.raises(FrozenInstanceError): z.uvw=(0,0,1)


@pytest.mark.parametrize('uvw',[(0,0,0),(1,2),(1.,0,0),(1.5,0,0),(True,0,0),('1',0,0),None,(np.nan,0,1)])
def test_invalid_axis(uvw):
    with pytest.raises(ValueError): ZoneAxis(uvw)


@pytest.mark.parametrize('uvw',[(0,0,1),(1,0,0),(1,1,0),(1,1,1),(1,-1,0),(2,1,1)])
@pytest.mark.parametrize('cell',[np.eye(3)*4,np.diag([3.,4.,5.]),np.array([[3.,0,0],[1.,4.,0],[.5,-.7,5.]])])
def test_beam_frame_and_round_trip(uvw,cell):
    d=DirectLattice(cell);z=ZoneAxis(uvw);o=CrystalOrientation(d,z)
    expected=np.asarray(z.uvw)@cell;expected/=np.linalg.norm(expected)
    np.testing.assert_allclose(zone_axis_cartesian(d,z),expected,atol=2e-16)
    r=o.beam_basis
    np.testing.assert_allclose(r@r.T,np.eye(3),atol=5e-16)
    assert np.linalg.det(r)==pytest.approx(1,abs=5e-16)
    np.testing.assert_allclose(np.cross(r[0],r[1]),r[2],atol=3e-16)
    vectors=np.array([[1.2,-3.4,.5],[.2,2.,-1.],[0,0,0]])
    np.testing.assert_allclose(o.beam_to_cartesian(o.cartesian_to_beam(vectors)),vectors,atol=2e-15)
    np.testing.assert_array_equal(r,CrystalOrientation(d,z).beam_basis)
    for a in (r,o.beam_direction_cartesian,o.cartesian_to_beam(vectors)):
        with pytest.raises(ValueError): a.setflags(write=True)


def test_cubic_001_exact_convention():
    o=CrystalOrientation(DirectLattice(np.eye(3)),ZoneAxis((0,0,1)))
    np.testing.assert_array_equal(o.beam_basis,[[0,-1,0],[1,0,0],[0,0,1]])
    np.testing.assert_array_equal(o.cartesian_to_beam([2,3,4]),[-3,2,4])


def test_nonorthogonal_is_not_normalized_indices():
    o=CrystalOrientation(DirectLattice([[3,0,0],[1,4,0],[.5,-.7,5]]),ZoneAxis((1,1,1)))
    assert not np.allclose(o.beam_direction_cartesian,np.ones(3)/np.sqrt(3))


def test_symmetric_reference_axis_roundoff_ties():
    for perturbation in (-1e-15,0,1e-15):
        o=CrystalOrientation(DirectLattice(np.diag([1+perturbation,1,1])),ZoneAxis((1,1,1)))
        assert o.reference_axis_index==0


@pytest.mark.parametrize('values',[[1,2],[1,np.inf,0],['1','2','3'],[1j,0,0],True])
def test_invalid_transform_vectors(values):
    o=CrystalOrientation(DirectLattice(np.eye(3)),ZoneAxis((1,1,0)))
    with pytest.raises(ValueError):o.cartesian_to_beam(values)


def test_zone_law_exact_integer_validation():
    z=ZoneAxis((-2,-2,0))
    assert zone_law_value((1,-1,7),z)==0
    assert zone_law_value((2,1,0),z)==3
    with pytest.raises(ValueError):zone_law_value((1.,0,0),z)
