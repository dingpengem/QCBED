from dataclasses import replace
import numpy as np
import pytest
from qcbed.crystallography import (DirectLattice, ZoneAxis, ZoneAxisCell,
    AtomicSite, AtomicStructureModel, IndexedBeamReflection, ZoneReflectionHK,
    build_zone_axis_simulation_sample)
from qcbed.models import MeshConfig, SliceGeometry
from qcbed.legacy_io import LegacyInput, Species, AtomSite as LegacyAtom
from qcbed.simulation import (SamplePreparationParameters, SliceAppliedParameters,
    prepare_multislice_sample, simulate_one_tilt)


def inputs():
    lattice = DirectLattice(np.diag([4.,4.,4.]))
    cell = ZoneAxisCell(lattice,ZoneAxis((0,0,1)))
    structure = AtomicStructureModel(lattice,(AtomicSite('Si',(.1,.2,.1),.8,.3),
                                              AtomicSite('Si',(.3,.4,.6),.7,.3)))
    beams = tuple(IndexedBeamReflection(str(i),(h,k),(h,k,0),ZoneReflectionHK(h,k))
                  for i,(h,k) in enumerate(((0,0),(1,0),(1,0))))
    options = dict(slice_boundaries_fractional=(0,.5,1),depth_angstrom=5.,
        mesh=MeshConfig(8,8),indexed_beams=beams,energy_kv=200.,dsap_inverse_angstrom=.5,
        absorption_coefficients=None,reciprocal_ordering='unshifted',phase_grating_method='direct')
    return structure,cell,options


def test_bridge_manual_legacy_complex_wave_equivalence():
    structure,cell,options=inputs()
    result=build_zone_axis_simulation_sample(structure,cell,**options)
    model=LegacyInput('manual',200.,False,(0,0,1),MeshConfig(8,8),
        (SliceGeometry(2.,cell_a_angstrom=4.,cell_b_angstrom=4.,lattice_gamma_degrees=90.),)*2,
        (Species('Si',14),),((LegacyAtom(1,.1,.2,.8,'S'),),(LegacyAtom(1,.3,.4,.7,'S'),)),
        '', '',absorption_raw='NO')
    parameters=SamplePreparationParameters(200.,(SliceAppliedParameters((.3,),0.,1.),)*2,
        .5,None,(0,1,0),((0,0),(1,0),(1,0)),'unshifted','direct','legacy','unity')
    assert result.parameters == parameters
    assert result.model.atoms_by_slice == model.atoms_by_slice
    assert result.effective_depth_angstrom == 6.
    assert result.requested_depth_angstrom == 5.
    reference=prepare_multislice_sample(model,parameters)
    for h,k in ((0.,0.),(.17,-.11)):
        kwargs=dict(tilt_h_undilated_fractional=h,tilt_k_undilated_fractional=k)
        a,b=simulate_one_tilt(result.prepared,**kwargs),simulate_one_tilt(reference,**kwargs)
        np.testing.assert_array_equal(a.final_reciprocal_wave,b.final_reciprocal_wave)
        np.testing.assert_array_equal(a.amplitudes,b.amplitudes)
        np.testing.assert_array_equal(a.intensities,b.intensities)
    for beam in result.indexed_beams:
        g=np.array(beam.original_hkl)@cell.original_lattice.reciprocal().vectors
        gi=np.array(beam.internal_hk.as_HKL())@result.depth_sample.simulation_cell.direct_lattice.reciprocal().vectors
        np.testing.assert_array_equal(g,gi)


@pytest.mark.parametrize('b',[None,.7])
def test_unsupported_thermal_not_fabricated(b):
    structure,cell,options=inputs()
    structure=replace(structure,sites=(structure.sites[0],replace(structure.sites[1],debye_waller_b=b)))
    with pytest.raises(ValueError,match='per-species B'):
        build_zone_axis_simulation_sample(structure,cell,**options)


def test_mismatched_hk_rejected():
    structure,cell,options=inputs()
    options['indexed_beams']=(replace(options['indexed_beams'][0],internal_hk=ZoneReflectionHK(1,0)),)
    with pytest.raises(ValueError,match='HK does not match'):
        build_zone_axis_simulation_sample(structure,cell,**options)
