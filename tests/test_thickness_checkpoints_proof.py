"""Pre-implementation proof: literal checkpoints against independent full runs."""
from dataclasses import replace

import numpy as np

from qcbed.refinement_trial import _prepare_trial, _gradient, evaluate_refinement_trial
from qcbed.refinement import apply_refinement_parameters
from qcbed.simulation import SamplePreparationParameters, prepare_multislice_sample, build_tilt_propagators
from qcbed.scan import generate_scan_tilts
from qcbed.multislice import extract_beam_amplitudes, beam_intensities
from qcbed.wave_propagation import (apply_phase_grating, real_to_reciprocal_wave,
    apply_propagator, reciprocal_to_real_wave)
from qcbed.masked_objective import compute_selected_objective
from test_refinement_driver import problem
from test_refinement import data


def proof_problem():
    s = problem()
    return replace(s, iref=-1., scan=replace(s.scan, corner_b=(.3, .2), corner_c=(-.2, .4)))


def literal_checkpoints(s, e, counts):
    mapped = apply_refinement_parameters(s, np.array([2.*max(counts)]))
    applied, sequence, _, coefficients = _prepare_trial(mapped, original_voltage=s.model.energy_kv)
    params = SamplePreparationParameters(mapped.model.energy_kv, applied, mapped.dsap,
        coefficients, sequence, e.beam_indices_hk, "centered", "fft", "legacy", "unity")
    base = prepare_multislice_sample(mapped.model, params)
    tilts = generate_scan_tilts(s.scan)
    out = {n: np.empty((len(tilts), len(e.beam_indices_hk)), dtype=np.complex128) for n in counts}
    for j, (h, k) in enumerate(tilts):
        kernels = build_tilt_propagators(base, tilt_h_undilated_fractional=h, tilt_k_undilated_fractional=k)
        wave = None
        for n, index in enumerate(sequence, 1):
            grating = base.slice_types[index].phase_grating
            transmitted = grating if wave is None else apply_phase_grating(wave, grating)
            reciprocal = real_to_reciprocal_wave(transmitted, ordering="centered")
            sampled = apply_propagator(reciprocal, kernels[index], wave_ordering="centered", propagator_ordering="centered")
            if n in out:
                out[n][j] = extract_beam_amplitudes(sampled, e.beam_indices_hk, ordering="centered")
            if n < len(sequence):
                wave = reciprocal_to_real_wave(sampled, ordering="centered")
    return out


def test_same_occurrence_different_depth_full_equality():
    s, e = proof_problem(), data(gradient=1)
    results = [evaluate_refinement_trial(s, np.array([d]), e, objective_mode="masked")
               for d in (829., 829.5, 830.)]
    for result in results:
        assert len(_prepare_trial(result.physical_state, original_voltage=200.)[1]) == 415
        for name in ("amplitudes", "calculated_intensities", "objective_intensities"):
            np.testing.assert_array_equal(getattr(result, name), getattr(results[0], name))
        assert result.scaling == results[0].scaling
        assert result.objective == results[0].objective


def test_literal_checkpoint_vs_independent_runs():
    s, e = proof_problem(), data(gradient=1)
    counts = (413, 414, 415, 416, 417)
    checkpoints = literal_checkpoints(s, e, counts)
    for n in counts:
        full = evaluate_refinement_trial(s, np.array([2.*n]), e, objective_mode="masked")
        np.testing.assert_array_equal(checkpoints[n], full.amplitudes)
        intensities = beam_intensities(checkpoints[n].ravel()).reshape(full.amplitudes.shape)
        np.testing.assert_array_equal(intensities, full.calculated_intensities)
        gradient = _gradient(intensities, e.scan_rows, e.scan_columns)
        np.testing.assert_array_equal(gradient, full.objective_intensities)
        objective = compute_selected_objective(gradient, e, iref=s.iref, objective_mode="masked")
        assert objective.scaling == full.scaling
        assert objective.objective == full.objective
