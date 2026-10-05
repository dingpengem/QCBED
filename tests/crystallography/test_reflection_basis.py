import numpy as np
import pytest
from qcbed.crystallography import (DirectLattice,ZoneAxis,OriginalReflectionBasis,
    reflection_from_basis_coordinates,original_reflection_basis_vectors,
    build_zone_axis_cell,original_hkl_to_internal_HK,internal_HK_to_original_hkl)


def test_user_basis_algebra_and_physical_vectors():
    z=ZoneAxis((1,1,0));basis=OriginalReflectionBasis((1,-1,0),(0,0,2),z)
    d=DirectLattice([[3.,0,0],[1.,4.,0],[.5,-.7,5.]])
    g=original_reflection_basis_vectors(d.reciprocal(),basis)
    np.testing.assert_array_equal(g,np.array([basis.hkl_1,basis.hkl_2])@d.reciprocal().vectors)
    q=build_zone_axis_cell(d,z)
    for n,expected in (((1,0),(1,-1,0)),((0,1),(0,0,2)),((1,1),(1,-1,2)),((-1,1),(-1,1,2)),((0,0),(0,0,0))):
        hkl=reflection_from_basis_coordinates(basis,*n)
        assert hkl==expected
        np.testing.assert_allclose(np.array(hkl)@d.reciprocal().vectors,np.array(n)@g,atol=1e-16)
        assert internal_HK_to_original_hkl(original_hkl_to_internal_HK(hkl,q),q)==hkl
    assert basis.hkl_2==(0,0,2) # nonprimitive basis explicitly preserved
    with pytest.raises(ValueError):g.setflags(write=True)
    with pytest.raises(ValueError):reflection_from_basis_coordinates(basis,1.,0)


@pytest.mark.parametrize('a,b',[((0,0,0),(0,0,2)),((1,0,0),(0,0,2)),((1,-1,0),(2,-2,0)),((2,-2,0),(3,-3,0)),((1.,-1,0),(0,0,2)),((True,-1,0),(0,0,2))])
def test_invalid_basis(a,b):
    with pytest.raises(ValueError):OriginalReflectionBasis(a,b,ZoneAxis((1,1,0)))
