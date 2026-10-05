from dataclasses import FrozenInstanceError
import numpy as np
import pytest
from qcbed.crystallography import (DirectLattice,AtomicSite,AtomicStructureModel,
    ZoneAxis,CrystalOrientation,project_atomic_structure,project_hkl,zone_law_value,
    plot_projected_atomic_structure)


@pytest.mark.parametrize('cell',[np.eye(3)*4,np.diag([3.,4.,5.]),np.array([[3.,0,0],[1.,4.,0],[.5,-.7,5.]])])
@pytest.mark.parametrize('uvw,hkls',[((0,0,1),[(1,0,0),(0,1,0),(1,1,0),(0,0,1)]),((1,1,0),[(1,-1,0),(0,0,2),(1,1,0)]),((1,1,1),[(1,-1,0),(0,1,-1),(1,1,-2),(1,0,0)]),((2,1,1),[(1,-2,0),(0,1,-1),(1,0,0)])])
def test_zone_law_projection(cell,uvw,hkls):
    d=DirectLattice(cell);o=CrystalOrientation(d,ZoneAxis(uvw));b=d.reciprocal()
    denominator=np.linalg.norm(np.array(o.zone_axis.uvw)@cell)
    for hkl in hkls:
        p=project_hkl(b,o,hkl);law=zone_law_value(hkl,o.zone_axis)
        assert p.beam_parallel_component==pytest.approx(law/denominator,abs=2e-16)
        assert (abs(p.beam_parallel_component)<1e-14)==(law==0)
        np.testing.assert_allclose(np.linalg.norm(p.g_cartesian),np.linalg.norm(p.g_beam),atol=2e-16)
        np.testing.assert_array_equal(p.g_projected_xy,p.g_beam[:2])
        for a in (p.g_cartesian,p.g_beam,p.g_projected_xy):
            with pytest.raises(ValueError): a.setflags(write=True)


@pytest.mark.parametrize('repeats',[(1,1,1),(2,2,2),(1,2,3)])
def test_atom_projection_provenance(repeats):
    d=DirectLattice(np.eye(3)*4)
    sites=(AtomicSite('Al',(-.2,.1,0),.5,.8),AtomicSite('O',(-.2,.1,.5)))
    s=AtomicStructureModel(d,sites);o=CrystalOrientation(d,ZoneAxis((0,0,1)))
    before=d.vectors.copy();r=o.beam_basis.copy()
    p=project_atomic_structure(s,o,repeats)
    assert len(p)==2*np.prod(repeats)
    assert [a.source_site_index for a in p]==[0,1]*int(np.prod(repeats))
    assert p[0].x==p[1].x and p[0].y==p[1].y # overlapping columns NOT merged
    assert p[0].depth!=p[1].depth
    for a in p:
        assert a.source_site is sites[a.source_site_index]
        assert a.symbol==a.source_site.symbol
        expected=(np.array(a.source_site.fractional_position)+a.cell_offset)@d.vectors@r.T
        np.testing.assert_allclose([a.x,a.y,a.depth],expected,atol=1e-15)
        assert np.isfinite(a.depth)
    with pytest.raises(FrozenInstanceError):p[0].x=7
    np.testing.assert_array_equal(d.vectors,before);np.testing.assert_array_equal(o.beam_basis,r)


def test_invalid_projection_geometry_and_repeats():
    s=AtomicStructureModel(DirectLattice(np.eye(3)),(AtomicSite('Al',(0,0,0)),))
    o=CrystalOrientation(s.direct_lattice,ZoneAxis((0,0,1)))
    for repeats in ((1,0,1),(1.,1,1),(True,1,1),None):
        with pytest.raises(ValueError):project_atomic_structure(s,o,repeats)
    other=CrystalOrientation(DirectLattice(np.eye(3)*2),ZoneAxis((0,0,1)))
    with pytest.raises(ValueError,match='same direct lattice'):project_atomic_structure(s,other)


def test_projected_plot_headless(monkeypatch,tmp_path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    s=AtomicStructureModel(DirectLattice(np.eye(3)*4),(AtomicSite('Al',(0,0,0)),AtomicSite('O',(.5,.5,.5))))
    o=CrystalOrientation(s.direct_lattice,ZoneAxis((1,1,1)))
    before=o.beam_basis.copy()
    monkeypatch.setattr(plt,'show',lambda:pytest.fail('unexpected show'))
    fig,ax=plot_projected_atomic_structure(s,o,(2,2,2))
    assert ax.figure is fig and ax.get_aspect()==1
    assert ax.get_title()=='Zone axis [1 1 1]'
    assert sum(len(c.get_offsets()) for c in ax.collections)==16
    assert ax.get_legend_handles_labels()[1]==['Al','O']
    np.testing.assert_array_equal(o.beam_basis,before)
    fig.savefig(tmp_path/'projection.png');plt.close(fig)
