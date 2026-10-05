"""Patch A assembly/adapter tests: no inferred real noise or physical indexing."""
from dataclasses import replace
import json
import numpy as np
import pytest

from qcbed.data_preparation import (
    ArcBeamCatalogue, BootstrapSite, DiskGeometry, ExtractionConfig, SignalMode,
    WeightMode, prepare_experimental_qcbed, scan_from_catalogue,
    assemble_prepared_dataset, adapt_prepared_dataset, prepared_dataset_manifest)
from qcbed.experimental import parse_experimental_intensities
from qcbed.refinement import pack_refinement_parameters
from qcbed.refinement_trial import evaluate_refinement_trial
from qcbed.execution import ForwardExecutionConfig
from examples.patch_a_closure import synthetic_input


INDICES=((0,0),(-1,0),(1,0),(0,-1),(0,1),(-1,-1),(1,1),(-1,1),(1,-1))


def assembled(mode, variance):
    # Deliberate synthetic geometry matching reference dimensions, no DM-input inference.
    origin=np.array([127.,126.]); g1=np.array([81.,9.]); g2=np.array([-10.,81.])
    sites=tuple(BootstrapSite(n,tuple(origin+n[0]*g1+n[1]*g2),None,None,
        'reference' if n==(0,0) else 'missing','synthetic geometry',beam_id=f'b{i}')
        for i,n in enumerate(INDICES))
    catalogue=ArcBeamCatalogue(sites,'synthetic',False)
    geometry=DiskGeometry(38,38); config=ExtractionConfig.from_disk_geometry(geometry,final_binning=2)
    y,x=np.indices((256,256),dtype=float); pattern=10+x*x+2*y*y
    before=pattern.copy()
    extra={} if not variance else dict(variance_map=np.full(pattern.shape,4.),uncertainty_provenance='synthetic')
    result=prepare_experimental_qcbed(pattern,catalogue,geometry,
        selected_beam_ids=[s.beam_id for s in reversed(sites)],config=config,signal_mode=mode,**extra)
    conversion=scan_from_catalogue(catalogue,central_beam_id='b0',canonical_basis_beam_ids=('b2','b4'),
        zone_axis_origin_px=tuple(origin),extraction_config=config,basis_transform=((1,0),(0,-1)))
    dataset=assemble_prepared_dataset(result,conversion,intensity_units='synthetic counts',physical_geometry=geometry)
    np.testing.assert_array_equal(before,pattern)
    return result,conversion,dataset


@pytest.mark.parametrize('mode',list(SignalMode))
@pytest.mark.parametrize('variance',[False,True])
def test_assembly_shape_values_order_and_immutability(mode,variance):
    extraction,conversion,data=assembled(mode,variance)
    assert data.intensities.shape==(1444,9)
    assert data.intensities.size==12996
    assert data.lattice_coordinates==INDICES
    assert [b.id for b in data.beams]==[f'b{i}' for i in range(9)]
    assert all(b.indices is None for b in data.beams)
    assert data.scan_geometry==conversion.scan_geometry
    assert data.extraction_config.physical_selection_shape_px==(76,76)
    assert data.physical_disk_geometry==DiskGeometry(38,38)
    assert data.weight_mode is WeightMode.BEAM_MASK
    for actual,expected in zip((data.intensities,data.variances,data.weights),extraction.tilt_beam_arrays()):
        if expected is None:
            assert actual is None
        else:
            np.testing.assert_array_equal(actual,expected)
            with pytest.raises(ValueError): actual.setflags(write=True)
    for i,b in enumerate(extraction.beams):
        np.testing.assert_array_equal(data.intensities[:,i],b.signal.ravel(order='C'))
        np.testing.assert_array_equal(data.weights[:,i],b.final_weight.ravel(order='C'))
        assert data.center_sources[i]==b.center_source
    manifest=prepared_dataset_manifest(data)
    encoded=json.loads(json.dumps(manifest))
    assert encoded['beam_mapping'][8]['lattice_coordinates']==[1,-1]
    assert encoded['raw_pixels_per_beam']==5776
    assert encoded['flattened_signal_size']==12996
    assert encoded['scan_geometry']['corner_a']==list(conversion.scan_geometry.corner_a)
    assert encoded['arrays']['variances'] is None if not variance else encoded['arrays']['variances']['shape']==[1444,9]


def test_absent_variance_is_dataset_valid_but_not_engine_ready():
    _,_,data=assembled(SignalMode.CONVENTIONAL,False)
    with pytest.raises(ValueError,match='variance is absent'):
        adapt_prepared_dataset(data,beam_hk_assignments=(),basis_1=(1,0,0),basis_2=(0,1,0),objective_mode='masked')


@pytest.mark.parametrize('mode',list(SignalMode))
def test_direct_vs_legacy_loader_full_trial_exact(mode,tmp_path):
    data,adapter,state=synthetic_input(mode)
    e=adapter.experimental
    path=tmp_path/'synthetic.dat'
    lines=[f'{e.scan_columns} {e.scan_rows} {e.gradient_option}']
    for name in ('intensities','variances','weights'):
        values=getattr(e,name)
        for beam in range(values.shape[1]):
            for row in values[:,beam].reshape(e.scan_rows,e.scan_columns):
                lines.append(' '.join(format(float(v),'.17g') for v in row))
    path.write_text('\n'.join(lines)+'\n')
    loaded=parse_experimental_intensities(path,beam_indices_hk=e.beam_indices_hk,scan=adapter.scan)
    assert loaded.beam_indices_hk==e.beam_indices_hk
    for name in ('intensities','variances','weights'):
        np.testing.assert_array_equal(getattr(e,name),getattr(loaded,name))
        np.testing.assert_array_equal(getattr(e,name),getattr(data,name))
    assert adapter.scan.corner_a==data.scan_geometry.corner_a
    vector=pack_refinement_parameters(state)
    kwargs=dict(objective_mode=adapter.objective_mode,execution=ForwardExecutionConfig('batched',64))
    a=evaluate_refinement_trial(state,vector,e,**kwargs)
    b=evaluate_refinement_trial(state,vector,loaded,**kwargs)
    assert np.isfinite(a.objective) and np.isfinite(a.scaling)
    assert a.objective==b.objective and a.scaling==b.scaling
    np.testing.assert_array_equal(a.amplitudes,b.amplitudes)
    np.testing.assert_array_equal(a.calculated_intensities,b.calculated_intensities)
    np.testing.assert_array_equal(a.objective_intensities,b.objective_intensities)
    # Numerical state and input arrays survive evaluation unchanged.
    np.testing.assert_array_equal(vector,pack_refinement_parameters(state))
    assert adapter.scan==a.physical_state.scan


def test_adapter_mapping_duplicates_readiness_and_masked_zero_variance():
    data,_,_=synthetic_input(SignalMode.CONVENTIONAL)
    assignments=tuple((b.id,(0,0)) for b in reversed(data.beams))
    kwargs=dict(beam_hk_assignments=assignments,basis_1=(1,0,0),basis_2=(0,1,0),objective_mode='masked')
    adapted=adapt_prepared_dataset(data,**kwargs)
    assert adapted.experimental.beam_indices_hk==((0,0),)*9
    assert [r[0] for r in adapted.beam_mapping]==[b.id for b in data.beams]
    assert [r[1] for r in adapted.beam_mapping]==list(data.lattice_coordinates)
    assert [r[3] for r in adapted.beam_mapping]==list(range(9))
    for assignments in ((),kwargs['beam_hk_assignments'][:-1],kwargs['beam_hk_assignments']+(kwargs['beam_hk_assignments'][0],)):
        with pytest.raises(ValueError,match='exactly once'):
            adapt_prepared_dataset(data,**{**kwargs,'beam_hk_assignments':assignments})
    v=data.variances.copy();w=data.weights.copy();v[0,0]=0;w[0,0]=0
    altered=replace(data,variances=v,weights=w)
    accepted=adapt_prepared_dataset(altered,**kwargs)
    assert accepted.experimental.variances[0,0]==0 # no placeholder or scientific repair
    with pytest.raises(ValueError,match='positive variance'):
        adapt_prepared_dataset(altered,**{**kwargs,'objective_mode':'legacy'})
    w[0,0]=1
    with pytest.raises(ValueError,match='positive variance'):
        adapt_prepared_dataset(replace(data,variances=v,weights=w),**kwargs)
    np.testing.assert_array_equal(altered.weights,accepted.experimental.weights)


def test_metadata_validation_and_no_mutable_aliases():
    _,_,data=assembled(SignalMode.CONVENTIONAL,False)
    source=data.intensities.copy(); snapshot=replace(data,intensities=source)
    source[:]=0
    np.testing.assert_array_equal(snapshot.intensities,data.intensities)
    for updates in (dict(center_sources=('one',)),dict(lattice_coordinates=((0,0),)),
                    dict(weight_mode='beam_mask'),dict(physical_disk_geometry=(38,38)),
                    dict(extraction_config=ExtractionConfig(2,2))):
        with pytest.raises((ValueError,TypeError)): replace(data,**updates)


def test_rejected_slots_keep_original_block_provenance():
    extraction,conversion,_=assembled(SignalMode.CONVENTIONAL,False)
    # Removing an already rejected block from the final dataset must not relabel others.
    rejected=replace(extraction.beams[0],status='rejected_coverage')
    extraction=replace(extraction,beams=(rejected,)+extraction.beams[1:])
    dataset=assemble_prepared_dataset(extraction,conversion,intensity_units='counts')
    assert dataset.extraction_block_indices==tuple(range(1,9))
    assert dataset.lattice_coordinates==INDICES[1:]
    manifest=prepared_dataset_manifest(dataset)
    assert manifest['beam_mapping'][0]['original_extraction_block']==1
    assert manifest['beam_mapping'][0]['refinement_column']==0
