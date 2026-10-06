"""One-period slicing through the unchanged numerical engine and source depth loop."""
from dataclasses import replace
from pathlib import Path
import shutil
import subprocess
import numpy as np
import pytest
from qcbed.crystallography import (DirectLattice,ZoneAxis,ZoneAxisCell,AtomicSite,
    AtomicStructureModel,build_zone_axis_simulation_sample)
from qcbed.models import MeshConfig,SliceGeometry
from qcbed.legacy_io import LegacyInput,Species,AtomSite
from qcbed.simulation import SamplePreparationParameters,SliceAppliedParameters,prepare_multislice_sample,build_tilt_propagators
from qcbed.multislice import simulate_slice_sequence
from qcbed.depth_geometry import requested_depth_to_sample_geometry


def build(structure,cell,edges,depth=5.,absorption=None,mesh=MeshConfig(8,8)):
    return build_zone_axis_simulation_sample(structure,cell,slice_boundaries_fractional=edges,
        depth_angstrom=depth,mesh=mesh,indexed_beams=(),energy_kv=200.,dsap_inverse_angstrom=.5,
        absorption_coefficients=absorption,reciprocal_ordering='unshifted',phase_grating_method='direct')


@pytest.mark.parametrize('edges,count', [((0,1),1),((0,.5,1),2),((0,.25,.5,.75,1),4)])
def test_stored_slices_assignment_and_repetition(edges,count):
    lattice=DirectLattice(np.eye(3)*4);cell=ZoneAxisCell(lattice,ZoneAxis((0,0,1)))
    structure=AtomicStructureModel(lattice,tuple(AtomicSite('Si',(.1,.2,z),1.,.3) for z in (0.,.25,.5,.75,1.,-.25)))
    sample=build(structure,cell,edges)
    assert len(sample.model.slices)==len(sample.prepared.slice_types)==count
    assert sorted(i for s in sample.depth_sample.slices for i in s.source_indices)==list(range(6))
    assert sample.depth_sample.c_period_translations==(0,0,0,0,1,-1)
    assert all(s.lower_fraction<=site.fractional_position[2]<s.upper_fraction for s in sample.depth_sample.slices for site in s.sites)
    assert sum(s.thickness_angstrom for s in sample.model.slices)==4.
    assert sample.parameters.slice_sequence_zero_based==tuple(i%count for i in range(count+1))
    assert sample.effective_depth_angstrom==4.+4./count


@pytest.mark.parametrize('absorbing',[False,True])
def test_one_slice_manual_presliced_complex_equivalence(monkeypatch,absorbing):
    from qcbed.absorption import AbsorptionCoefficients
    import qcbed.phase_grating as pg
    potentials=[];original=pg.projected_potential_direct
    def capture(*args):
        value=original(*args);potentials.append(value.copy());return value
    monkeypatch.setattr(pg,'projected_potential_direct',capture)
    lattice=DirectLattice(np.eye(3)*4);cell=ZoneAxisCell(lattice,ZoneAxis((0,0,1)))
    structure=AtomicStructureModel(lattice,(AtomicSite('Si',(.1,.2,.1),.8,.3),AtomicSite('Si',(.3,.4,.6),.7,.3)))
    absorption=(AbsorptionCoefficients(.1,1.,1.),) if absorbing else None
    actual=build(structure,cell,(0,1),depth=9.,absorption=absorption)
    manual=LegacyInput('explicit one-slice reference',200.,False,(0,0,1),MeshConfig(8,8),
        (SliceGeometry(4.,cell_a_angstrom=4.,cell_b_angstrom=4.,lattice_gamma_degrees=90.),),
        (Species('Si',14),),((AtomSite(1,.1,.2,.8,'S'),AtomSite(1,.3,.4,.7,'S')),),'','',absorption_raw='YE' if absorbing else 'NO')
    parameters=SamplePreparationParameters(200.,(SliceAppliedParameters((.3,),0.,1.),),
        .5,absorption,(0,0,0),(),'unshifted','direct','legacy','unity')
    reference=prepare_multislice_sample(manual,parameters)
    assert actual.parameters==parameters
    assert actual.effective_depth_angstrom==12.
    np.testing.assert_array_equal(potentials[0],potentials[1])
    assert potentials[0].dtype==np.complex128 and np.isfinite(potentials[0]).all()
    np.testing.assert_array_equal(actual.prepared.slice_types[0].phase_grating,reference.slice_types[0].phase_grating)
    maximum=0.
    for h,k in ((0.,0.),(.17,-.11)):
        options=dict(tilt_h_undilated_fractional=h,tilt_k_undilated_fractional=k)
        kernels=[build_tilt_propagators(p,**options) for p in (actual.prepared,reference)]
        for n in (1,2,3): # Compare the complex wave after each slice occurrence.
            waves=[simulate_slice_sequence((p.slice_types[0].phase_grating,),kernel,(0,)*n,
                beam_indices_hk=((0,0),(1,0),(0,1)),reciprocal_ordering='unshifted',propagator_ordering='unshifted')
                for p,kernel in zip((actual.prepared,reference),kernels)]
            error=float(np.max(np.abs(waves[0].final_reciprocal_wave-waves[1].final_reciprocal_wave)))
            maximum=max(maximum,error)
            np.testing.assert_array_equal(waves[0].final_reciprocal_wave,waves[1].final_reciprocal_wave)
            np.testing.assert_array_equal(waves[0].intensities,waves[1].intensities)
    print(f'one-slice bridge/manual complex-wave max absolute error (absorption={absorbing}): {maximum}')


def test_non_001_mesh_is_transverse_recast_cell():
    lattice=DirectLattice(np.diag([3.,4.,5.]));cell=ZoneAxisCell(lattice,ZoneAxis((1,0,0)))
    structure=AtomicStructureModel(lattice,(AtomicSite('Si',(.1,.2,.3),1.,.3),))
    sample=build(structure,cell,(0,1),mesh=MeshConfig(8,6))
    p,q,c=sample.depth_sample.simulation_cell.direct_lattice.vectors
    np.testing.assert_array_equal(p,(0,-4,0));np.testing.assert_array_equal(q,(0,0,-5));np.testing.assert_array_equal(c,(3,0,0))
    g=sample.model.slices[0];assert (g.cell_a_angstrom,g.cell_b_angstrom,g.thickness_angstrom)==(4.,5.,3.)
    assert sample.model.atoms_by_slice[0][0].x_fractional==-.2
    assert sample.model.atoms_by_slice[0][0].y_fractional==-.3
    prepared=sample.prepared.slice_types[0]
    assert prepared.phase_grating.shape==(6,8) # rows Q, columns P; detector shape is independent.
    assert (g.cell_a_angstrom/8,g.cell_b_angstrom/6)==(.5,5./6)
    # Changing transverse sampling must not change the period or sequence.
    finer=build(structure,cell,(0,1),mesh=MeshConfig(16,12))
    assert finer.parameters.slice_sequence_zero_based==sample.parameters.slice_sequence_zero_based==(0,0)
    assert finer.prepared.slice_types[0].phase_grating.shape==(12,16)


def test_one_slice_does_not_bypass_geometry_proof():
    lattice=DirectLattice(np.eye(3)*4);cell=ZoneAxisCell(lattice,ZoneAxis((1,1,0)))
    structure=AtomicStructureModel(lattice,(AtomicSite('Si',(0,0,0),1.,.3),))
    with pytest.raises(ValueError,match='Unsupported depth-resolved geometry'):
        build(structure,cell,(0,1))


def test_compiled_original_fortran_sequence_loop(tmp_path):
    compiler=shutil.which('gfortran')
    if compiler is None:pytest.skip('gfortran unavailable for executable source-loop comparison')
    from private_reference import private_reference_root
    source=(private_reference_root()/'legacy_fortran/oneit.f').read_text()
    source=source.split('C DETERMINE THE SLICE SEQUENCE FOR THE WHOLE SPECIMEN',1)[1]
    body=source[source.index('          K=0'):source.index('          MS=K')]+ '          MS=K\n'
    harness='''program sequence_reference
implicit none
integer :: k,sln,j,i,l,nregions,ms
integer,allocatable :: rnsl(:),mseq(:),rslices(:)
real(kind=8),allocatable :: depth(:),deltaz(:),applieddilation(:),slices(:)
real(kind=8) :: thick
character(len=80) :: arg
nregions=1
allocate(rnsl(1),rslices(1),depth(1),deltaz(1),applieddilation(1),slices(1))
call get_command_argument(1,arg)
read(arg,*) depth(1)
call get_command_argument(2,arg)
read(arg,*) deltaz(1)
allocate(mseq(ceiling(depth(1)/deltaz(1))+2))
rnsl=1
applieddilation=1d0
'''+body+'''write(*,*) ms,thick
write(*,*) mseq(:ms)
end program
'''
    path=tmp_path/'sequence.f90';path.write_text(harness);exe=tmp_path/'sequence'
    subprocess.run([compiler,'-O0',str(path),'-o',str(exe)],check=True,capture_output=True,text=True)
    for thickness in (4.,.1):
        for depth in (.01,.3,4.,8.,float(np.nextafter(8.,np.inf)),9.):
            lines=subprocess.check_output([str(exe),repr(depth),repr(thickness)],text=True).splitlines()
            n,effective=lines[0].split();sequence=tuple(int(v)-1 for v in lines[1].split())
            actual=requested_depth_to_sample_geometry(depth,(thickness,))
            assert int(n)==actual.occurrence_count and sequence==actual.sequence
            assert float(effective)==actual.effective_depth
    print('Compiled original ONEIT MSEQ loop: exact sequence/effective-depth agreement for 12 one-slice cases')


def test_modern_depth_refinement_rebuild_roundtrip_and_resume(tmp_path):
    from test_modern_grouped_run import modern_config
    from qcbed import prepare_run_snapshot,run_snapshot,save_yaml_config,load_yaml_config
    from qcbed.crystallography_refinement import apply_crystallographic_parameters
    from qcbed.grouped_run import evaluate_grouped,grouped_identity,load_grouped_checkpoint,check_grouped
    c=modern_config(tmp_path,mode=0)
    c=replace(c,crystallography=replace(c.crystallography,slice_boundaries_fractional=(0.,1.)))
    snapshot=prepare_run_snapshot(c,base_directory=tmp_path)
    sample=snapshot.compiled.crystallography.simulation_sample
    assert len(sample.model.slices)==snapshot.compiled.state.regions[0].slice_count==1
    assert len(check_grouped(snapshot.compiled,snapshot.experimental))==2
    vector=np.asarray(snapshot.compiled.packed_starting_vector,dtype=np.float64)
    index=next(p.vector_index for p in snapshot.compiled.layout.parameters if p.category=='depth')
    changed=vector.copy();changed[index]=9.
    rebuilt=apply_crystallographic_parameters(snapshot.compiled,changed)
    later=rebuilt.crystallography.simulation_sample
    assert len(later.prepared.slice_types)==1
    assert later.parameters.slice_sequence_zero_based==(0,0,0)
    np.testing.assert_array_equal(sample.prepared.slice_types[0].phase_grating,later.prepared.slice_types[0].phase_grating)
    assert sample.parameters.slice_sequence_zero_based==(0,)
    extent=requested_depth_to_sample_geometry(9.,(sample.model.slices[0].thickness_angstrom,))
    assert later.effective_depth_angstrom==extent.effective_depth
    first=evaluate_grouped(snapshot,vector)[1]
    moved=evaluate_grouped(snapshot,changed)[1]
    repeated=evaluate_grouped(snapshot,vector)[1]
    assert np.isfinite(moved['amplitudes']).all()
    assert not np.array_equal(first['amplitudes'],moved['amplitudes'])
    np.testing.assert_array_equal(first['amplitudes'],repeated['amplitudes'])
    save_yaml_config(c,tmp_path/'one.yml',base_directory=tmp_path)
    loaded=load_yaml_config(tmp_path/'one.yml')
    assert loaded.crystallography.slice_boundaries_fractional==(0.,1.)
    restored=prepare_run_snapshot(loaded,base_directory=tmp_path)
    assert grouped_identity(restored)==grouped_identity(snapshot)
    # Real initialization/checkpoint/resume, cancelled before optimizer evaluation.
    assert run_snapshot(snapshot,cancel_requested=lambda:True).status=='cancelled'
    assert run_snapshot(restored,resume=True,cancel_requested=lambda:True).status=='cancelled'
    cp=load_grouped_checkpoint(tmp_path/'run',snapshot)
    assert cp['restart_count']==1
    two=replace(c,crystallography=replace(c.crystallography,slice_boundaries_fractional=(0.,.5,1.)))
    alternative=prepare_run_snapshot(two,base_directory=tmp_path)
    assert grouped_identity(alternative)!=grouped_identity(snapshot)
    with pytest.raises(ValueError,match='identity mismatch'):
        load_grouped_checkpoint(tmp_path/'run',alternative)
