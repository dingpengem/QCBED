"""Profiling only: three unchanged production trial calls, no cache/optimizer."""
from collections import defaultdict
from contextlib import ExitStack
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import platform
import sys
from time import perf_counter
from unittest.mock import patch

import numpy as np

import qcbed.refinement_trial as trial
import qcbed.simulation as simulation
from qcbed.experimental import parse_experimental_intensities
from qcbed.refinement import build_refinement_parameter_layout, pack_refinement_parameters, apply_refinement_parameters
from qcbed.refinement_driver import parse_legacy_refinement
from qcbed.refinement_io import _CaptureRecords
from qcbed.legacy_io import _parse_for007_records


def thickness_only(state):
    """Copy flags, not physical values; preserve the source depth range pair."""
    layout=build_refinement_parameter_layout(state)
    depths=[p for p in layout.parameters if p.category=="depth"]
    if len(state.regions)!=1 or len(depths)!=1 or depths[0].target!=(0,):
        raise ValueError("this profiling case requires one active region depth")
    model=replace(state.model,energy_refined=False,atoms_by_slice=tuple(
        tuple(replace(a,x_mask=0,y_mask=0) for a in atoms) for atoms in state.model.atoms_by_slice))
    result=replace(state,model=model,
        regions=tuple(replace(r,thickness_mask=1,vacancy_masks=tuple(0 for _ in r.vacancy_masks),
            debye_waller_masks=tuple(0 for _ in r.debye_waller_masks)) for r in state.regions),
        dilation_mask=0,absorption_masks=tuple((0,0,0) for _ in state.absorption_masks),
        scan_masks=("NN",)*3,
        reflections=tuple(tuple(replace(r,mask="NN") for r in row) for row in state.reflections),
        ranges=(depths[0].range_factors,) if depths[0].range_factors is not None else ())
    only=build_refinement_parameter_layout(result).parameters
    assert len(only)==1 and only[0].category=="depth" and only[0].target==(0,)
    return result,depths[0].vector_index


def profile_case(state, vector, experimental):
    times=defaultdict(float); counts=defaultdict(int)
    evidence={}; kernels_hash=hashlib.sha256()

    def sample_evidence(sample):
        evidence['wavelength_angstrom']=sample.wavelength_angstrom
        evidence['grating_hashes']=[hashlib.sha256(s.phase_grating.tobytes()).hexdigest() for s in sample.slice_types]
        evidence['grid_hashes']=[hashlib.sha256(s.reciprocal_grid.stl_inverse_angstrom.tobytes()).hexdigest()
                                 for s in sample.slice_types]
        evidence['beam_positions']=sample.beam_positions_row_col

    def sequence_evidence(value):
        evidence['occurrences']=len(value[1])
        evidence['sequence']=value[1]
        evidence['realized_depth_angstrom']=sum(state.model.slices[i].thickness_angstrom for i in value[1])

    def kernels_evidence(kernels):
        for kernel in kernels: kernels_hash.update(kernel.tobytes())

    def wrapper(name, fn, after=None):
        def measured(*args,**kwargs):
            start=perf_counter()
            try:
                result=fn(*args,**kwargs)
            finally:
                times[name]+=perf_counter()-start; counts[name]+=1
            if after is not None: after(result)
            return result
        return measured

    probes=[('mapping',trial.apply_refinement_parameters,None),
            ('sequence_and_applied_state',trial._prepare_trial,sequence_evidence),
            ('sample_preparation',trial.prepare_multislice_sample,sample_evidence),
            ('propagator_construction',trial.build_tilt_propagators,kernels_evidence),
            ('tilt_multislice',trial.simulate_slice_sequence,None),
            ('gradient',trial._gradient,None),
            ('scan_generation',trial.generate_scan_tilts,
             lambda a:evidence.update(scan_hash=hashlib.sha256(a.tobytes()).hexdigest())),
            ('ynorm_objective',trial.compute_selected_objective,None),
            ('elastic_structure_factors',simulation.elastic_structure_factors,None),
            ('absorptive_structure_factors',simulation.absorptive_structure_factors,None),
            ('phase_gratings',simulation.phase_grating_from_structure_factors,None),
            ('wavelength',simulation.electron_wavelength_angstrom,None),
            ('beam_mapping',simulation.beam_indices_to_row_col,None)]
    with ExitStack() as stack:
        # Patch imported aliases of the SAME function; do not substitute kernels.
        for name,fn,after in probes:
            measured=wrapper(name,fn,after)
            for module_name,module in tuple(sys.modules.items()):
                if module_name.startswith('qcbed.') and module is not None:
                    for key,value in tuple(vars(module).items()):
                        if value is fn: stack.enter_context(patch.object(module,key,measured))
        for name in ('fft2','ifft2'):
            stack.enter_context(patch.object(np.fft,name,wrapper('numpy_'+name,getattr(np.fft,name))))
        start=perf_counter()
        result=trial.evaluate_refinement_trial(state,vector,experimental,objective_mode='masked')
        wall=perf_counter()-start
    evidence['propagator_hash']=kernels_hash.hexdigest()
    evidence['amplitudes_hash']=hashlib.sha256(result.amplitudes.tobytes()).hexdigest()
    top=('mapping','sequence_and_applied_state','sample_preparation','propagator_construction',
         'tilt_multislice','gradient','scan_generation','ynorm_objective')
    evidence.update(wall_seconds=wall,times=dict(times),counts=dict(counts),
        unattributed_wall_seconds=wall-sum(times[n] for n in top),
        objective=float(result.objective),ynorm=float(result.scaling),depth=float(vector[0]))
    return result,evidence


def main():
    root=Path(__file__).resolve().parents[1]
    source=root/'tests/fixtures/legacy_reference'/'FOR007.DAT'
    baseline,_=parse_legacy_refinement(source,seed=2026)
    original=repr(baseline)
    state,index=thickness_only(baseline)
    state_before=repr(state)
    records=_CaptureRecords(source); _parse_for007_records(records)
    row=records.values['output beam count and indices']
    beams=tuple(zip(row[1::2],row[2::2]))
    experimental=parse_experimental_intensities(source.parent/'in.dat',beam_indices_hk=beams,scan=state.scan)
    vector=pack_refinement_parameters(state)
    step=state.model.slices[0].thickness_angstrom
    changed=vector.copy(); changed[0]+=step
    mapped=apply_refinement_parameters(state,changed)
    assert replace(mapped,regions=(replace(mapped.regions[0],depth_angstrom=state.regions[0].depth_angstrom),))==state
    assert all(not a.x_mask and not a.y_mask for row in state.model.atoms_by_slice for a in row)
    destination=root/'tmp'/'thickness_profile'; destination.mkdir(parents=True,exist_ok=True)
    report=dict(interpreter=sys.executable,python=sys.version,numpy=np.__version__,platform=platform.platform(),
        original_parameter_index_zero_based=index,original_layout_count=len(build_refinement_parameter_layout(baseline).parameters),
        thickness_only_count=1,starting_depth=float(vector[0]),perturbation=step,
        mesh=[state.model.mesh.mesh_y,state.model.mesh.mesh_x],tilts=1521,beams=len(beams),cases={})
    first=None
    for label,v in (('initial',vector),('repeat',vector),('changed',changed)):
        print(f'Start {label}: depth={v[0]}',flush=True)
        result,metrics=profile_case(state,v,experimental)
        if first is None: first=result
        elif label=='repeat':
            np.testing.assert_array_equal(result.amplitudes,first.amplitudes)
            np.testing.assert_array_equal(result.calculated_intensities,first.calculated_intensities)
            assert result.objective==first.objective and result.scaling==first.scaling
        metrics['max_amplitude_difference_from_initial']=float(np.max(np.abs(result.amplitudes-first.amplitudes)))
        assert repr(baseline)==original and repr(state)==state_before
        assert vector[0]==state.regions[0].depth_angstrom
        report['cases'][label]=metrics
        report['baseline_unchanged']=True
        (destination/'profile.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
        print(label,json.dumps({k:v for k,v in metrics.items() if k not in ('sequence',)}),flush=True)
    assert report['cases']['changed']['occurrences']==report['cases']['initial']['occurrences']+1
    for key in ('wavelength_angstrom','grating_hashes','grid_hashes','beam_positions','scan_hash','propagator_hash'):
        assert report['cases']['initial'][key]==report['cases']['repeat'][key]==report['cases']['changed'][key]
    print(destination/'profile.json',flush=True)


if __name__=='__main__':
    main()
