from dataclasses import FrozenInstanceError
from itertools import product
import numpy as np
import pytest
from qcbed.crystallography import (DirectLattice,ZoneAxis,CrystalOrientation,
    build_zone_axis_cell,original_fractional_to_zone,zone_fractional_to_original,
    AtomicSite,AtomicStructureModel,recast_structure_to_zone_cell)

AXES=((0,0,1),(1,0,0),(1,1,0),(1,1,1),(1,-1,0),(2,1,1),(3,2,1))
CELLS=(np.eye(3)*4,np.diag([3.,4.,5.]),np.array([[3.,0,0],[1.,4.,0],[.5,-.7,5.]]))


def determinant(t):
    a,b,c=map(lambda r:tuple(map(int,r)),t)
    return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])


@pytest.mark.parametrize('uvw',AXES)
@pytest.mark.parametrize('cell',CELLS)
def test_zone_cell_invariants(uvw,cell):
    d=DirectLattice(cell);z=ZoneAxis(uvw);q=build_zone_axis_cell(d,z)
    assert q.transform.dtype.kind=='i' and determinant(q.transform)==1
    assert tuple(q.transform[2])==z.uvw
    np.testing.assert_array_equal(q.transform.astype(object)@q.inverse_transform.astype(object),np.eye(3,dtype=object))
    np.testing.assert_array_equal(q.direct_lattice.vectors,q.transform@cell)
    assert q.direct_lattice.volume==pytest.approx(d.volume,rel=2e-15)
    b=q.reciprocal_lattice.vectors
    np.testing.assert_allclose(b@q.direct_lattice.vectors.T,np.eye(3),atol=2e-15)
    np.testing.assert_allclose(b,q.inverse_transform.T@d.reciprocal().vectors,atol=1e-15)
    o=CrystalOrientation(d,z);c=q.direct_lattice.vectors[2]
    np.testing.assert_allclose(c/np.linalg.norm(c),o.beam_direction_cartesian,atol=3e-16)
    np.testing.assert_allclose(b[:2]@o.beam_direction_cartesian,0,atol=1e-15)
    np.testing.assert_array_equal(q.transform,build_zone_axis_cell(d,z).transform)
    for a in (q.transform,q.inverse_transform,q.direct_lattice.vectors,b):
        with pytest.raises(ValueError):a.setflags(write=True)
    with pytest.raises(FrozenInstanceError):q.zone_axis=ZoneAxis((1,0,0))


def test_completion_all_small_signed_axes():
    # Includes zero components, negative signs, nonprimitive input and w=0.
    for axis in product(range(-3,4),repeat=3):
        if axis==(0,0,0):continue
        z=ZoneAxis(axis);q=build_zone_axis_cell(DirectLattice(np.eye(3)),z)
        assert determinant(q.transform)==1
        assert tuple(q.transform[2])==z.uvw
        np.testing.assert_array_equal(q.transform.astype(object)@q.inverse_transform.astype(object),np.eye(3,dtype=object))


@pytest.mark.parametrize('uvw',AXES)
def test_fractional_and_full_structure_recast(uvw):
    d=DirectLattice(CELLS[2]);q=build_zone_axis_cell(d,ZoneAxis(uvw))
    f=np.array([[-.2,1.7,.35],[.4,-.7,2.1],[-.2,1.7,.35]])
    s=AtomicStructureModel(d,(AtomicSite('Al',tuple(f[0]),.5,.8),AtomicSite('O',tuple(f[1]),1,None),AtomicSite('Al',tuple(f[2]),.5,.8)))
    recast=recast_structure_to_zone_cell(s,q)
    g=original_fractional_to_zone(f,q)
    np.testing.assert_allclose(g@q.direct_lattice.vectors,f@d.vectors,atol=4e-15)
    np.testing.assert_allclose(zone_fractional_to_original(g,q),f,atol=2e-15)
    np.testing.assert_array_equal([a.fractional_position for a in recast.sites],g)
    assert len(recast.sites)==len(s.sites)==3
    assert [(a.symbol,a.occupancy,a.debye_waller_b) for a in recast.sites]==[(a.symbol,a.occupancy,a.debye_waller_b) for a in s.sites]
    np.testing.assert_array_equal([a.fractional_position for a in s.sites],f)
    assert np.any((g<0)|(g>1)) # deliberately unwrapped
    with pytest.raises(ValueError):g.setflags(write=True)


def test_storage_overflow_and_lattice_mismatch():
    with pytest.raises(OverflowError,match='no truncation'):
        build_zone_axis_cell(DirectLattice(np.eye(3)),ZoneAxis((2**64,1,0)))
    d=DirectLattice(np.eye(3));q=build_zone_axis_cell(d,ZoneAxis((1,1,0)))
    s=AtomicStructureModel(DirectLattice(np.eye(3)*2),(AtomicSite('Al',(0,0,0)),))
    with pytest.raises(ValueError,match='original lattice'):recast_structure_to_zone_cell(s,q)
