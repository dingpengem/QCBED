"""Strict geometry proof; no approximate full-period fallback."""
import numpy as np
import pytest
from qcbed.crystallography.lattice import DirectLattice
from qcbed.crystallography.orientation import ZoneAxis
from qcbed.crystallography.zone_cell import ZoneAxisCell
from qcbed.crystallography.structure import AtomicSite, AtomicStructureModel
from qcbed.crystallography.simulation_cell import (
    MultisliceSimulationCell, simulation_cell_compatibility, build_depth_resolved_sample,
)


def cell(rows, axis=(0, 0, 1)):
    return ZoneAxisCell(DirectLattice(rows), ZoneAxis(axis))


@pytest.mark.parametrize('axis,accepted', [((0,0,1),True), ((1,1,0),False), ((1,1,1),False)])
def test_cubic_strict(axis, accepted):
    c = cell(np.eye(3)*4, axis)
    assert simulation_cell_compatibility(c).compatible == accepted
    if not accepted:
        with pytest.raises(ValueError, match='Unsupported depth-resolved'):
            MultisliceSimulationCell(c)


@pytest.mark.parametrize('rows', [np.diag([3,4,5]),
    [[3,0,0],[-1.5,3*np.sqrt(3)/2,0],[0,0,5]],
    [[3,0,5],[1,4,-5],[0,0,5]]])
def test_transverse_metric_and_hk_preserved(rows):
    c = cell(rows)
    sim = MultisliceSimulationCell(c)
    p,q,z = sim.direct_lattice.vectors
    assert np.dot(p,z) == 0
    assert np.dot(q,z) == 0
    np.testing.assert_allclose(sim.direct_lattice.reciprocal().vectors[:2],
                               c.reciprocal_lattice.vectors[:2], atol=1e-15)
    np.testing.assert_allclose(sim.frame @ sim.frame.T, np.eye(3), atol=1e-15)
    assert np.linalg.det(sim.frame) > 0
    with pytest.raises(ValueError):
        sim.transform_from_zone[0,0] = 2


@pytest.mark.parametrize('rows', [ [[3,0,0],[1,4,0],[.5,-.7,5]],
    [[3,0,5+1e-12],[1,4,-5],[0,0,5]] ])
def test_skew_or_nearly_integer_not_approximated(rows):
    assert not simulation_cell_compatibility(cell(rows)).compatible


def test_depth_order_metadata_periodic_ownership_and_cartesian_identity():
    c = cell([[3,0,5],[1,4,-5],[0,0,5]])
    sites = (AtomicSite('Si',(0,0,.5),.7,.3), AtomicSite('O',(0,0,1),1,.2),
             AtomicSite('Si',(2,0,-1.75),.8,.4))
    structure = AtomicStructureModel(c.original_lattice, sites)
    result = build_depth_resolved_sample(structure,c,slice_boundaries_fractional=(0,.5,1))
    assert [s.source_indices for s in result.slices] == [(1,2),(0,)]
    assert result.c_period_translations == (0,1,0)
    assert [s.thickness_angstrom for s in result.slices] == [2.5,2.5]
    for slab in result.slices:
        for i, site in zip(slab.source_indices, slab.sites):
            original = sites[i]
            assert (site.symbol,site.occupancy,site.debye_waller_b) == (
                original.symbol,original.occupancy,original.debye_waller_b)
            actual = np.array(site.fractional_position) @ result.simulation_cell.direct_lattice.vectors
            actual += result.c_period_translations[i]*result.simulation_cell.direct_lattice.vectors[2]
            np.testing.assert_array_equal(actual,np.array(original.fractional_position) @ c.original_lattice.vectors)
    assert result.slices[0].sites[1].fractional_position[0] == 2 # no x wrapping
    assert result.original_structure is structure


@pytest.mark.parametrize('edges', [(),(0,),(0,.5,.5,1),(0,float('nan'),1),(0,.5,2),(False,.5,1)])
def test_explicit_depth_boundaries(edges):
    c = cell(np.eye(3))
    structure = AtomicStructureModel(c.original_lattice,(AtomicSite('Si',(0,0,0)),))
    with pytest.raises(ValueError,match='boundaries'):
        build_depth_resolved_sample(structure,c,slice_boundaries_fractional=edges)
