"""Audit-only B7e regression: do not reinterpret vacancy dilation as strain."""
from dataclasses import replace
import numpy as np
from examples.b7e_gradient_audit import run_audit
from examples.b7a_refinement_layout import real_config
from qcbed import compile_config_from_sources
from qcbed.refinement import apply_refinement_parameters, pack_refinement_parameters
from qcbed.refinement_trial import _prepare_trial


def test_legacy_dilation_mapping_has_no_zero_vacancy_effect():
    state = compile_config_from_sources(real_config()).state
    state = replace(state, dilation_mask=1,
        regions=tuple(replace(r, debye_waller_masks=(0,)*len(r.debye_waller)) for r in state.regions))
    for raw, expected in ((0.,0.), (.25,.25), (1.,1.), (-2.,.5)):
        mapped = apply_refinement_parameters(state, np.array([raw], dtype=np.float64))
        assert mapped.dilation == expected
        # Existing supported trial subset explicitly has unit applied dilation.
        applied, _, _, _ = _prepare_trial(mapped, original_voltage=state.model.energy_kv)
        assert all(s.dilation == 1. for s in applied)
    assert pack_refinement_parameters(replace(state, dilation=0.))[0] == 1e-10


def test_gradient_baseline_and_deferred_gate():
    report = run_audit()
    assert report['status'] == 'PASS (audit only)'
    assert report['dilation_status'] == 'DEFERRED'
    assert report['global_lattice_scale_bounds'] is None
    assert not report['optimizer_executed']
    assert report['canonical_configuration_unchanged']
    assert report['baseline']['active_count'] == 7281
    assert report['active_weights'] == [1.0]
    np.testing.assert_allclose(
        [report['baseline']['ynorm'], report['baseline']['objective']],
        [13447206.106343817, 1.1896911365294378], rtol=1e-14, atol=0)
