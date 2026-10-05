from itertools import product
import numpy as np
import pytest
from qcbed.crystallography import (DirectLattice,ZoneAxis,build_zone_axis_cell,
    original_hkl_to_zone_HKL,zone_HKL_to_original_hkl,is_zone_reflection_original,
    is_zolz_zone_index,original_hkl_to_internal_HK,internal_HK_to_original_hkl,ZoneReflectionHK)

AXES=((0,0,1),(1,0,0),(1,1,0),(1,1,1),(1,-1,0),(2,1,1),(3,2,1))


@pytest.mark.parametrize('uvw',AXES)
@pytest.mark.parametrize('cell',[np.eye(3)*4,np.diag([3.,4.,5.]),np.array([[3.,0,0],[1.,4.,0],[.5,-.7,5.]])])
def test_exhaustive_index_invariance_and_zone_equivalence(uvw,cell):
    d=DirectLattice(cell);q=build_zone_axis_cell(d,ZoneAxis(uvw))
    for hkl in product(range(-3,4),repeat=3):
        h=original_hkl_to_zone_HKL(hkl,q)
        assert all(type(v) is int for v in h)
        assert zone_HKL_to_original_hkl(h,q)==hkl
        np.testing.assert_allclose(np.array(hkl)@d.reciprocal().vectors,np.array(h)@q.reciprocal_lattice.vectors,atol=4e-15)
        assert is_zone_reflection_original(hkl,q.zone_axis)==is_zolz_zone_index(h)
        assert h[2]==sum(a*b for a,b in zip(hkl,q.zone_axis.uvw))
        if h[2]==0:
            internal=original_hkl_to_internal_HK(hkl,q)
            assert internal.as_HKL()==h
            assert internal_HK_to_original_hkl(internal,q)==hkl
        else:
            with pytest.raises(ValueError,match='outside the zone'):original_hkl_to_internal_HK(hkl,q)


def test_huge_indices_use_python_integers_without_overflow():
    q=build_zone_axis_cell(DirectLattice(np.eye(3)),ZoneAxis((3,2,1)))
    hkl=(10**70,-10**65,17)
    assert zone_HKL_to_original_hkl(original_hkl_to_zone_HKL(hkl,q),q)==hkl


@pytest.mark.parametrize('bad',[(1.,0,0),(True,0,0),(1,2),(np.nan,0,0)])
def test_index_validation(bad):
    q=build_zone_axis_cell(DirectLattice(np.eye(3)),ZoneAxis((1,1,0)))
    for fn in (original_hkl_to_zone_HKL,zone_HKL_to_original_hkl):
        with pytest.raises(ValueError):fn(bad,q)
    with pytest.raises(ValueError):is_zolz_zone_index(bad)


def test_internal_type_and_110_example():
    q=build_zone_axis_cell(DirectLattice(np.eye(3)),ZoneAxis((1,1,0)))
    h=original_hkl_to_internal_HK((1,-1,2),q)
    assert h==ZoneReflectionHK(1,-2)
    assert internal_HK_to_original_hkl(h,q)==(1,-1,2)
    with pytest.raises(TypeError):internal_HK_to_original_hkl((1,-2),q)
    with pytest.raises(ValueError):ZoneReflectionHK(1.,2)
