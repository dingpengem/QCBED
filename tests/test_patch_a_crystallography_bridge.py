from dataclasses import replace
import numpy as np
import pytest
from qcbed.crystallography import (DirectLattice,ZoneAxis,OriginalReflectionBasis,
    build_zone_axis_cell,index_prepared_dataset,index_beam_catalogue,
    internal_scan_basis,adapt_indexed_prepared_dataset)
from qcbed.data_preparation import SignalMode,ArcBeamCatalogue,BootstrapSite
from qcbed.scan import generate_scan_tilts
from test_prepared_dataset_integration import assembled


def indexed_data(mode=SignalMode.CONVENTIONAL,variance=True):
    _,conversion,data=assembled(mode,variance)
    z=ZoneAxis((1,1,0));b=OriginalReflectionBasis((1,-1,0),(0,0,2),z)
    q=build_zone_axis_cell(DirectLattice(np.diag([3.,4.,5.])),z)
    return index_prepared_dataset(data,b,q),conversion


@pytest.mark.parametrize('mode',list(SignalMode))
def test_prepared_bridge_native_adapter_values_and_scan(mode):
    indexed,conversion=indexed_data(mode)
    data=indexed.dataset
    snapshots=[a.copy() for a in (data.intensities,data.variances,data.weights)]
    adapted=adapt_indexed_prepared_dataset(indexed,scan_basis_transform=((1,0),(0,-1)),objective_mode='masked')
    expected=tuple((r.internal_hk.H,r.internal_hk.K) for r in indexed.indexed_beams)
    assert adapted.experimental.beam_indices_hk==expected
    assert adapted.experimental.gradient_option==(1 if mode==SignalMode.GRADIENT else 0)
    assert adapted.scan.basis_1==(1.,0.,0.) and adapted.scan.basis_2==(0.,2.,0.)
    assert adapted.scan.corner_a==conversion.scan_geometry.corner_a
    assert adapted.scan.corner_b==conversion.scan_geometry.corner_b
    assert adapted.scan.corner_c==conversion.scan_geometry.corner_c
    assert (adapted.scan.scan_rows,adapted.scan.scan_columns)==(38,38)
    tilts=generate_scan_tilts(adapted.scan)
    a=np.array(conversion.scan_geometry.corner_a)
    np.testing.assert_array_equal(tilts[0],a@np.array([[1,0],[0,2]]))
    for name,before in zip(('intensities','variances','weights'),snapshots):
        np.testing.assert_array_equal(getattr(data,name),before)
        np.testing.assert_array_equal(getattr(adapted.experimental,name),before)
    assert data.beams[1].indices is None # original Patch A record not relabelled
    assert [r.experimental_index for r in indexed.indexed_beams]==list(data.lattice_coordinates)
    assert [r.stable_id for r in indexed.indexed_beams]==[b.id for b in data.beams]
    assert adapted.experimental.intensities.shape==(1444,9)


def test_missing_variance_does_not_block_indexing_or_infer_uncertainty():
    indexed,_=indexed_data(variance=False)
    assert len(indexed.indexed_beams)==9 and indexed.dataset.variances is None
    with pytest.raises(ValueError,match='variance is absent'):
        adapt_indexed_prepared_dataset(indexed,scan_basis_transform=((1,0),(0,-1)),objective_mode='masked')


def test_catalogue_order_and_no_status_filtering():
    indexed,_=indexed_data()
    data=indexed.dataset
    sites=tuple(BootstrapSite(n,b.center_px,None,None,'missing','test',beam_id=b.id)
                for b,n in zip(data.beams,data.lattice_coordinates))
    cat=ArcBeamCatalogue(tuple(reversed(sites)),'test',False)
    result=index_beam_catalogue(cat,indexed.reflection_basis,indexed.zone_cell)
    assert result==tuple(reversed(indexed.indexed_beams))


def test_missing_coordinates_not_taken_from_beam_physical_indices():
    indexed,_=indexed_data()
    data=replace(indexed.dataset,lattice_coordinates=())
    with pytest.raises(ValueError,match='explicit'):
        index_prepared_dataset(data,indexed.reflection_basis,indexed.zone_cell)


@pytest.mark.parametrize('transform',[((1,0),(0,0)),((2,0),(0,1)),((1.,0),(0,1)),((True,0),(0,1)),((1,0,0),(0,1,0))])
def test_invalid_explicit_scan_transform(transform):
    indexed,_=indexed_data()
    with pytest.raises(ValueError):internal_scan_basis(indexed,scan_basis_transform=transform)


def test_scan_transform_is_explicit_and_not_beam_reindexing():
    indexed,_=indexed_data()
    first=internal_scan_basis(indexed,scan_basis_transform=((1,0),(0,1)))
    second=internal_scan_basis(indexed,scan_basis_transform=((0,1),(-1,0)))
    assert first==((1,0,0),(0,-2,0))
    assert second==((0,2,0),(1,0,0))
    assert indexed.indexed_beams[4].original_hkl==(0,0,2)
    with pytest.raises(TypeError):adapt_indexed_prepared_dataset(indexed,objective_mode='masked')
