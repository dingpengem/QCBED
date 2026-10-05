import sys
import subprocess
import numpy as np
import pytest
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from qcbed.crystallography import AtomicSite,AtomicStructureModel,DirectLattice,plot_atomic_structure


def structure():
    return AtomicStructureModel(DirectLattice([[3,0,0],[1,4,0],[0,0,5]]),
        (AtomicSite('Al',(-.1,0,0)),AtomicSite('O',(.5,.5,.5))))


@pytest.mark.parametrize('repeats',[(1,1,1),(2,2,2)])
def test_plot(repeats,monkeypatch,tmp_path):
    model=structure(); before=model.direct_lattice.vectors.copy()
    monkeypatch.setattr(plt,'show',lambda: pytest.fail('unexpected interactive display'))
    fig,ax=plot_atomic_structure(model,repeats)
    assert ax.figure is fig and ax.name=='3d'
    assert len(ax.lines)==12*np.prod(repeats)
    assert sum(len(c._offsets3d[0]) for c in ax.collections)==2*np.prod(repeats)
    xyz=np.column_stack(ax.collections[0]._offsets3d)
    np.testing.assert_array_equal(xyz[0],np.array(model.sites[0].fractional_position)@before)
    assert ax.get_legend_handles_labels()[1]==['Al','O']
    np.testing.assert_array_equal(model.direct_lattice.vectors,before)
    fig.savefig(tmp_path/'structure.png');plt.close(fig)


@pytest.mark.parametrize('repeats',[(0,1,1),(1,2),(1.,1,1),(True,1,1)])
def test_invalid_repeats(repeats):
    with pytest.raises(ValueError): plot_atomic_structure(structure(),repeats)


def test_backend_import_does_not_load_frontend_dependencies():
    subprocess.run([sys.executable,'-c',
        'import sys; import qcbed.crystallography; '
        'assert "matplotlib" not in sys.modules; assert "ase" not in sys.modules'],check=True)
