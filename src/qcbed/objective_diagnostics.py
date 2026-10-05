"""B6b fixed-input diagnostics: existing B5 forward engine and objective only.

No optimizer, parameter mapping, sample reconstruction, or model fitting. The
attached prepared sample is used by identity, after the B6a provenance gate.
"""
from dataclasses import asdict
import numpy as np
from .input_preflight import preflight_inputs
from .scan import generate_scan_tilts
from .tilt_series import simulate_tilt_series_batched
from .refinement_trial import _gradient
from .masked_objective import compute_selected_objective


def statistics(values):
    a = np.asarray(values)
    if not a.size:
        return dict(count=0)
    return dict(shape=list(a.shape), min=float(a.min()), max=float(a.max()),
                median=float(np.median(a)), nan_count=int(np.isnan(a).sum()),
                inf_count=int(np.isinf(a).sum()), negative_count=int((a < 0).sum()))


def residual_diagnostics(calculated, experimental, scaling, objective):
    """Decompose existing masked root error, without fitting or normalization.

    Beam contributions add to objective squared. Their error fractions times the
    root error provide an explicitly labelled additive allocation of that error.
    Zero-weight entries are never divided by variance.
    """
    x, v, w = experimental.intensities, experimental.variances, experimental.weights
    active = w != 0
    scaled = scaling*calculated
    residual = scaled-x
    errors, norms = np.zeros_like(x), np.zeros_like(x)
    errors[active] = w[active]*residual[active]*residual[active]/v[active]
    norms[active] = w[active]*x[active]*x[active]/v[active]
    error, norm = np.float64(0), np.float64(0)
    for e,n in zip(errors.T.ravel(), norms.T.ravel()):
        error += e; norm += n
    if not np.isfinite(error) or not np.isfinite(norm) or norm <= 0:
        raise ValueError('Invalid diagnostic residual decomposition')
    reconstructed = float(np.sqrt(error/norm))
    if not np.isclose(reconstructed, objective, rtol=1e-14, atol=1e-14):
        raise ValueError('Residual decomposition differs from engine objective')
    beams = []
    for b in range(x.shape[1]):
        a = active[:, b]
        observed = x[a,b]; theory = calculated[a,b]
        e = float(sum(errors[:,b], np.float64(0)))
        fraction = e/float(error) if error != 0 else 0.
        observed_sum, calculated_sum = float(observed.sum()), float(theory.sum())
        beams.append(dict(active_count=int(a.sum()), experimental=statistics(observed),
            objective_signal=statistics(theory), experimental_sum=observed_sum,
            experimental_mean=float(observed.mean()), calculated_sum=calculated_sum,
            calculated_mean=float(theory.mean()),
            observed_to_calculated_sum_ratio=observed_sum/calculated_sum if calculated_sum else None,
            observed_to_scaled_sum_ratio=observed_sum/(scaling*calculated_sum) if scaling*calculated_sum else None,
            weighted_squared_error=e, squared_objective_contribution=e/float(norm),
            fraction_of_total_error=fraction, allocated_root_objective=fraction*float(objective)))
    return dict(weighted_squared_error=float(error), weighted_experimental_norm=float(norm),
        reconstructed_objective=reconstructed, zero_weight_contribution_exact=bool(np.all(errors[~active]==0)),
        beams=beams), dict(scaled_calculated=scaled, residual=residual,
                          weighted_squared_residual=errors, weights=w, variances=v, experimental=x)


def evaluate_preflight_objectives(compiled, datasets, *, baseline_manifests,
                                 input_kind, physical_input_confirmation=None, verify_repeat=True):
    """One fixed evaluation per mode; optional verification repeat (default enabled).

    Currently restricted to the existing B4b attached sample and batched backend.
    Reject other execution policies instead of silently overriding configuration.
    Preflight checks the same sample instance that is passed to propagation.
    """
    if type(verify_repeat) is not bool:
        raise TypeError('verify_repeat must be bool')
    preflight = preflight_inputs(compiled, datasets, baseline_manifests=baseline_manifests,
        input_kind=input_kind, physical_input_confirmation=physical_input_confirmation)
    if compiled.execution.execution_backend != 'batched':
        raise ValueError('B6b diagnostic requires configured batched execution')
    if compiled.runner.objective_mode != 'masked':
        raise ValueError('B6b diagnostic currently requires explicit masked objective')
    crystal = compiled.crystallography
    sample = crystal.simulation_sample.prepared
    source = compiled.source_configuration
    report = dict(status='PASS', preflight=preflight,
        source_configuration=asdict(source), resolved_configuration=asdict(compiled.configuration),
        prepared_parameters=asdict(sample.parameters),
        active_parameter_count=compiled.parameter_count,
        parameter_layout=[repr(p) for p in compiled.layout.parameters],
        packed_starting_vector=compiled.packed_starting_vector.tolist(),
        parameter_policy='All current values fixed; no optimizer or parameter updates',
        forward_policy='Exact attached B5 sample; its direct gratings and unshifted ordering retained',
        modes={}, legacy_comparison='No matched historical firstout.dat comparison claimed')
    artifacts = {}
    for mode, data in datasets.items():
        adapted = crystal.adapt_dataset(data, objective_mode=compiled.runner.objective_mode)
        exp = adapted.experimental
        tilts = generate_scan_tilts(adapted.scan)
        def once():
            series = simulate_tilt_series_batched(sample, tilts, batch_size=compiled.execution.batch_size)
            if series.beam_indices_hk != exp.beam_indices_hk or not np.array_equal(series.tilts,tilts):
                raise ValueError('Forward beam/tilt ordering changed')
            c = series.intensities
            if c.shape != exp.intensities.shape or not np.isfinite(c).all() or np.any(c < 0):
                raise ValueError('Invalid calculated intensities before objective')
            if not np.isfinite(series.amplitudes).all():
                raise ValueError('Nonfinite complex amplitudes')
            values = _gradient(c,exp.scan_rows,exp.scan_columns) if exp.gradient_option==1 else c
            result = compute_selected_objective(values,exp,iref=compiled.state.iref,
                                                objective_mode=compiled.runner.objective_mode)
            if not np.isfinite(result.objective) or result.objective < 0 or not np.isfinite(result.scaling):
                raise ValueError('Invalid objective/YNORM; no numerical repair attempted')
            return series, values, result
        first, values, result = once()
        second, repeated_values, repeated_result = once() if verify_repeat else (first,values,result)
        pairs = dict(amplitudes=(first.amplitudes,second.amplitudes),
                     calculated_intensities=(first.intensities,second.intensities),
                     objective_intensities=(values,repeated_values),
                     ynorm=(result.scaling,repeated_result.scaling),
                     objective=(result.objective,repeated_result.objective))
        repeat_errors = {k:float(np.max(np.abs(np.asarray(a)-np.asarray(b)))) for k,(a,b) in pairs.items()}
        if any(not np.array_equal(a,b) for a,b in pairs.values()):
            raise ValueError(f'Nonidentical deterministic repeat: {repeat_errors}')
        if not verify_repeat:
            repeat_errors = {}  # Not measured, never report a fabricated zero repeat.
        decomposition, arrays = residual_diagnostics(values,exp,result.scaling,result.objective)
        for index,row in enumerate(decomposition['beams']):
            row.update(preflight['beam_table'][index])
            row['calculated_raw_intensity'] = statistics(first.intensities[exp.weights[:,index]!=0,index])
        arrays.update(complex_amplitudes=first.amplitudes, calculated_intensities=first.intensities,
                      objective_intensities=values, tilts=tilts)
        artifacts[mode] = arrays
        report['modes'][mode] = dict(gradient_option=exp.gradient_option,
            primary_evaluations=1, determinism_repeat_evaluations=int(verify_repeat),
            ynorm=float(result.scaling), objective=float(result.objective),
            calculated_intensities=statistics(first.intensities), objective_signal=statistics(values),
            amplitudes=dict(shape=list(first.amplitudes.shape),
                finite_real_count=int(np.isfinite(first.amplitudes.real).sum()),
                finite_imaginary_count=int(np.isfinite(first.amplitudes.imag).sum()),
                max_absolute=float(np.abs(first.amplitudes).max())),
            active_count=int((exp.weights!=0).sum()), repeat_max_absolute_errors=repeat_errors,
            decomposition=decomposition)
    return report, artifacts
