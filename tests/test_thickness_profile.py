"""Profiling setup checks only; no production fast path or optimizer."""
from dataclasses import replace
import importlib.util
from pathlib import Path

import numpy as np

from qcbed.refinement_driver import parse_legacy_refinement
from qcbed.refinement import build_refinement_parameter_layout, apply_refinement_parameters
from qcbed.refinement_trial import _prepare_trial


def test_reference_thickness_only_layout_and_fixed_physics():
    root=Path(__file__).parents[1]
    spec=importlib.util.spec_from_file_location('thickness_profile',root/'benchmarks'/'profile_thickness_only.py')
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    baseline,_=parse_legacy_refinement(root/'tests/fixtures/legacy_reference'/'FOR007.DAT',seed=2026)
    before=repr(baseline)
    state,index=module.thickness_only(baseline)
    layout=build_refinement_parameter_layout(state).parameters
    assert index==0 and len(layout)==1
    assert (layout[0].category,layout[0].target,layout[0].initial_value)==('depth',(0,),840.)
    vector=np.array([842.0245]); saved=vector.copy()
    changed=apply_refinement_parameters(state,vector)
    assert replace(changed,regions=(replace(changed.regions[0],depth_angstrom=840.),))==state
    assert state.model.slices==baseline.model.slices
    assert state.absorption==baseline.absorption and state.scan==baseline.scan
    assert state.regions[0].debye_waller==baseline.regions[0].debye_waller
    initial_sequence=_prepare_trial(state,original_voltage=200.)[1]
    changed_sequence=_prepare_trial(changed,original_voltage=200.)[1]
    assert initial_sequence==(0,)*415 and changed_sequence==(0,)*416
    plateau=apply_refinement_parameters(state,np.array([840.1]))
    assert _prepare_trial(plateau,original_voltage=200.)[1]==initial_sequence
    np.testing.assert_array_equal(vector,saved)
    assert repr(baseline)==before
