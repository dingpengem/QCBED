import numpy as np
import pytest
ase=pytest.importorskip('ase')
from qcbed.crystallography import atomic_structure_from_ase


def test_ase_order_cell_unwrapped_and_independence():
    a=np.array([[3.,0,0],[1.,4.,0],[.5,-.7,5.]])
    f=np.array([[-.25,1.2,.3],[.7,.4,.1],[-.25,1.2,.3]])
    atoms=ase.Atoms(['Al','O','Al'],scaled_positions=f,cell=a,pbc=True)
    model=atomic_structure_from_ase(atoms)
    np.testing.assert_array_equal(model.direct_lattice.vectors,a)
    np.testing.assert_allclose([s.fractional_position for s in model.sites],f,atol=3e-16,rtol=0)
    assert [s.symbol for s in model.sites]==['Al','O','Al']
    assert all(s.occupancy==1 and s.debye_waller_b is None for s in model.sites)
    atoms.positions[:]=0;atoms.cell[:]=np.eye(3)
    np.testing.assert_array_equal(model.direct_lattice.vectors,a)
    assert model.sites[0].fractional_position[0]<0


@pytest.mark.parametrize('cell',[np.zeros((3,3)),np.diag([1,1,-1])])
def test_invalid_ase_cell(cell):
    with pytest.raises(ValueError): atomic_structure_from_ase(ase.Atoms('Al',cell=cell))


def test_invalid_input_and_empty():
    with pytest.raises(TypeError): atomic_structure_from_ase(None)
    with pytest.raises(ValueError): atomic_structure_from_ase(ase.Atoms(cell=np.eye(3)))
