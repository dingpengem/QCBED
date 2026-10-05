from dataclasses import replace, FrozenInstanceError
import numpy as np
import pytest
from qcbed.crystallography import AtomicSite,AtomicStructureModel,DirectLattice


def test_unwrapped_order_duplicates_immutable():
    coords=[-0.25,1.2,0.3];a=AtomicSite('Al',coords,.5,.8);coords[0]=8
    b=AtomicSite('O',(0,0,0));sites=[a,b,a]
    model=AtomicStructureModel(DirectLattice(np.eye(3)),sites);sites.reverse()
    assert model.sites==(a,b,a)
    assert a.fractional_position==(-.25,1.2,.3)
    with pytest.raises(FrozenInstanceError): a.occupancy=1


@pytest.mark.parametrize('updates',[dict(symbol='X'),dict(symbol='al'),dict(occupancy=0),dict(occupancy=-1),dict(occupancy=1.1),dict(occupancy=np.nan),dict(occupancy=True),dict(debye_waller_b=-1),dict(debye_waller_b=np.inf),dict(fractional_position=(1,2)),dict(fractional_position=(0,np.nan,0))])
def test_invalid_site(updates):
    with pytest.raises(ValueError): replace(AtomicSite('Si',(0,0,0)),**updates)


def test_empty_or_invalid_sites():
    for sites in ((),('Si',)):
        with pytest.raises(ValueError): AtomicStructureModel(DirectLattice(np.eye(3)),sites)
