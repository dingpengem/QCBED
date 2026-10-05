"""Opt-in current-path profiling and exact reuse census; never caches results.

Run from the repository root with qcbed-python. Only this benchmark installs
temporary wrappers/line tracing. Production imports and execution are untouched.
"""
from collections import defaultdict
from contextlib import ExitStack
from dataclasses import fields, is_dataclass, asdict
from datetime import datetime
import argparse
import hashlib
import importlib.util
import inspect
import json
import os
from pathlib import Path
import platform
import struct
import sys
from time import perf_counter
from unittest.mock import patch

import numpy as np
import psutil

import qcbed.refinement_trial as trial
import qcbed.tilt_series as series
from qcbed.refinement import pack_refinement_parameters, build_refinement_parameter_layout
from qcbed.refinement_driver import parse_legacy_refinement, run_legacy_refinement
from qcbed.refinement_subset import select_refinement_categories
from qcbed.refinement_bounds import RefinementBoundRules, resolve_refinement_bounds
from qcbed.experimental import parse_experimental_intensities
from qcbed.refinement_io import _CaptureRecords
from qcbed.legacy_io import _parse_for007_records


def exact(value):
    """Hashable typed value identity, including exact IEEE bytes (no rounding).

    Sets compare these complete keys, not digests. Signed zeros stay distinct.
    Dataclass metadata is conservatively included where present.
    """
    if isinstance(value, np.ndarray):
        return ('array', value.dtype.str, value.shape, value.tobytes(order='C'))
    if is_dataclass(value):
        return (type(value).__module__, type(value).__name__,
                tuple((f.name, exact(getattr(value, f.name))) for f in fields(value)))
    if isinstance(value, (float, np.floating)):
        return ('float64', struct.pack('!d', float(value)))
    if isinstance(value, (tuple, list)):
        return (type(value).__name__, tuple(exact(v) for v in value))
    if isinstance(value, dict):
        return ('dict', tuple((k, exact(v)) for k, v in sorted(value.items())))
    if value is None or isinstance(value, (str, int, bool)):
        return (type(value).__name__, value)
    raise TypeError(f'unsupported identity value: {type(value)}')


def signatures(mapped, applied_state, prepared, tilts, execution):
    """Conservative dependencies for the supported deterministic trial subset."""
    applied, sequence, region_ids, coefficients = applied_state
    p = prepared.parameters
    preparation = exact((mapped.model, applied, coefficients, mapped.dsap,
                         mapped.reflections, mapped.symmetry, region_ids,
                         p.phase_grating_method, p.boundary_policy, p.reciprocal_ordering))
    propagators = exact((mapped.model.mesh, mapped.model.slices,
                         p.applied_voltage_kv, tuple(s.dilation for s in applied),
                         tilts, p.boundary_policy, p.reciprocal_ordering))
    return dict(mapped=exact(mapped), sequence=exact(sequence), scan=exact(mapped.scan),
                prepared=preparation, propagator=propagators,
                forward=(preparation, propagators, exact(sequence),
                         exact(p.beam_indices_hk), exact(execution)))


class Census:
    """All-history opportunity counts; no computed arrays/results retained."""
    def __init__(self):
        self.seen = defaultdict(set)
        self.repeats = defaultdict(int)
        self.eligible = defaultdict(int)

    def add(self, keys):
        repeated = {}
        for name, key in keys.items():
            repeated[name] = key in self.seen[name]
            self.repeats[name] += repeated[name]
            self.eligible[name] += 1
            self.seen[name].add(key)
        return repeated

    def report(self, total):
        return {name: dict(eligible=self.eligible[name], unique=len(values),
                          repeats=self.repeats[name],
                          potential_hit_rate=self.repeats[name]/self.eligible[name],
                          repeats_per_all_trials=self.repeats[name]/total)
                for name, values in self.seen.items()}


def _line_regions():
    """Locate inline phases in original code, fail loudly if source changes."""
    def lines(fn, markers):
        source, start = inspect.getsourcelines(fn)
        found = {}
        for name, marker in markers.items():
            matches = [start+i for i, line in enumerate(source) if line.strip().startswith(marker)]
            if len(matches) != 1:
                raise RuntimeError(f'profiling marker changed: {marker}')
            found[name] = matches[0]
        return found
    chunk = lines(series._chunk, dict(propagation='sequence =', extraction='return reciprocal['))
    setup = lines(trial._prepare_trial, dict(start='extent = requested_depth_to_sample_geometry(', end='offset +='))
    return {series._chunk.__code__: lambda line: (
        'beam_extraction' if line == chunk['extraction'] else
        'propagator_construction' if line < chunk['propagation'] else 'propagation'),
        trial._prepare_trial.__code__: lambda line: (
            'sequence' if setup['start'] <= line <= setup['end'] else 'applied_state_setup')}


class TrialProfiler:
    """Single-thread benchmark context. Restores all hooks even on failure.

    Inline phases require line tracing; times are instrumented, not estimates
    of uninstrumented wall time. Nested inclusive timers are reported separately.
    """
    def __init__(self, on_trial=None):
        self.rows = []
        self.census = Census()
        self.active = None
        self.on_trial = on_trial

    def __enter__(self):
        if sys.gettrace() is not None:
            raise RuntimeError('run profiling outside a debugger/coverage tracer')
        self.stack = ExitStack()
        regions = _line_regions()

        def trace(frame, event, arg):
            if event != 'call' or frame.f_code not in regions or self.active is None:
                return None
            classify = regions[frame.f_code]
            previous = [None, perf_counter()]

            def local(frame, event, arg):
                now = perf_counter()
                if previous[0] is not None:
                    self.active['seconds'][previous[0]] += now-previous[1]
                previous[:] = [classify(frame.f_lineno) if event == 'line' else None, now]
                return local
            return local

        def wrap(fn, name, capture=None):
            def measured(*args, **kwargs):
                if self.active is None:
                    return fn(*args, **kwargs)
                start = perf_counter()
                try:
                    value = fn(*args, **kwargs)
                finally:
                    self.active['seconds'][name] += perf_counter()-start
                if capture:
                    self.active[capture] = value
                return value
            return measured

        original = trial._evaluate_with_invalid_policy

        def evaluate(*args, **kwargs):
            if self.active is not None:
                raise RuntimeError('nested trial profiling is unsupported')
            active = dict(seconds=defaultdict(float))
            self.active = active
            start = perf_counter()
            try:
                result = original(*args, **kwargs)
            finally:
                elapsed = perf_counter()-start
                self.active = None
            signature_start = perf_counter()
            mapped = getattr(result, 'physical_state', None)
            keys = {} if mapped is None else dict(mapped=exact(mapped), scan=exact(mapped.scan))
            if 'prepared' in active and 'tilts' in active and isinstance(result, trial.RefinementTrialResult):
                keys = signatures(mapped, active['applied'], active['prepared'],
                                  active['tilts'], kwargs.get('execution'))
            repeated = self.census.add(keys)
            times = dict(active['seconds'])
            # Inclusive timers below are NOT additional disjoint components.
            components = ('mapping', 'sequence', 'applied_state_setup', 'sample_preparation',
                          'tilt_generation', 'propagator_construction', 'propagation',
                          'beam_extraction', 'intensity', 'gradient', 'objective')
            accounted = sum(times.get(k, 0.) for k in components)
            row = dict(index=len(self.rows)+1, total_seconds=elapsed, seconds=times,
                       other_seconds=elapsed-accounted, repeated=repeated,
                       classification=type(result).__name__, objective=float(result.objective),
                       occurrences=len(active['applied'][1]) if 'applied' in active else None,
                       signature_seconds=perf_counter()-signature_start,
                       rss_mib=psutil.Process().memory_info().rss/2**20)
            self.rows.append(row)
            if self.on_trial is not None:
                self.on_trial(row)
            return result

        try:
            for module, name, label, capture in (
                (trial, 'apply_refinement_parameters', 'mapping', 'mapped'),
                (trial, '_prepare_trial', 'sequence_setup_inclusive', 'applied'),
                (trial, 'prepare_multislice_sample', 'sample_preparation', 'prepared'),
                (trial, 'generate_scan_tilts', 'tilt_generation', 'tilts'),
                (trial, 'simulate_tilt_series_batched', 'forward_inclusive', None),
                (series, 'beam_intensities', 'intensity', None),
                (trial, '_gradient', 'gradient', None),
                (trial, 'compute_selected_objective', 'objective', None)):
                self.stack.enter_context(patch.object(module, name, wrap(getattr(module, name), label, capture)))
            self.stack.enter_context(patch.object(trial, '_evaluate_with_invalid_policy', evaluate))
            sys.settrace(trace)
        except BaseException:
            self.stack.close()
            raise
        return self

    def __exit__(self, *args):
        sys.settrace(None)
        self.stack.close()

    def report(self):
        return dict(total_trials=len(self.rows), census=self.census.report(len(self.rows)), trials=self.rows)


def load_workload():
    """Use the editable runner's current categories, bounds and ANSI controls."""
    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location('refinement_example', root/'examples/run_full_refinement.py')
    example = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(example)
    state, original = parse_legacy_refinement(example.INPUT_FILE, seed=example.SEED,
                                             experimental_header_path=example.EXPERIMENTAL_FILE)
    state, original = select_refinement_categories(state, original, example.REFINEMENT_CATEGORIES)
    rules = RefinementBoundRules(example.SPECIFIC_BOUNDS, example.CATEGORY_BOUNDS,
                                example.CATEGORY_RELATIVE_BOUNDS, example.CATEGORY_RELATIVE_FALLBACK)
    config = example.make_config(original, refinement_bounds=resolve_refinement_bounds(state, rules),
                                 bound_rules=rules)
    records = _CaptureRecords(example.INPUT_FILE)
    _parse_for007_records(records)
    row = records.values['output beam count and indices']
    experimental = parse_experimental_intensities(example.EXPERIMENTAL_FILE,
        beam_indices_hk=tuple(zip(row[1::2], row[2::2])), scan=state.scan)
    return state, experimental, config


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    destination = args.output or root/'tmp/refinement_reuse'/datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    destination = destination.resolve()
    if any(destination.is_relative_to(root/p) for p in ('legacy_fortran', 'reference_data', 'src', 'tests')):
        raise ValueError('benchmark output must be separate from source/reference data')
    destination.mkdir(parents=True, exist_ok=False)
    state, experimental, config = load_workload()
    kwargs = dict(objective_mode=config.objective_mode, execution=config.refinement.execution,
                  bounds_policy=config.bounds_policy, refinement_bounds=config.refinement_bounds,
                  invalid_trial_policy=config.invalid_trial_policy)
    vector = pack_refinement_parameters(state)
    report = dict(platform=platform.platform(), machine=platform.machine(), interpreter=sys.executable,
                  numpy=np.__version__, psutil=psutil.__version__, configuration=asdict(config),
                  layout=asdict(build_refinement_parameter_layout(state)),
                  thread_environment={k:os.environ.get(k) for k in
                    ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS')},
                  source_sha256={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in [*sorted((root/'src/qcbed').glob('*.py')), Path(__file__),
                              root/'examples/run_full_refinement.py', root/'reference_data/FOR007.DAT',
                              root/'reference_data/in.dat']})
    def save():
        (destination/'results.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(destination, flush=True)
    controls = []
    for i in range(2):
        start = perf_counter()
        baseline = trial.evaluate_refinement_trial(state, vector, experimental, **kwargs)
        controls.append(perf_counter()-start)
    with TrialProfiler() as probe:
        measured = trial.evaluate_refinement_trial(state, vector, experimental, **kwargs)
    for name in ('amplitudes', 'calculated_intensities', 'objective_intensities'):
        np.testing.assert_array_equal(getattr(baseline, name), getattr(measured, name))
    assert baseline.objective == measured.objective and baseline.scaling == measured.scaling
    report.update(uninstrumented_control_seconds=controls, detailed_control=probe.report(),
                  control_max_amplitude_error=float(np.max(abs(baseline.amplitudes-measured.amplitudes))),
                  control_max_intensity_error=float(np.max(abs(baseline.calculated_intensities-measured.calculated_intensities))))
    save()
    # Run the original driver, including its stage-boundary evaluations. No
    # runner preview/verification calls are mixed into the optimizer census.
    def record(row):
        print(f"trial {row['index']}: {row['objective']}", flush=True)
        report['refinement'] = profile.report()
        save()
    with TrialProfiler(on_trial=record) as profile:
        start = perf_counter()
        result = run_legacy_refinement(state, experimental, config.refinement,
            objective_mode=config.objective_mode, invalid_trial_policy=config.invalid_trial_policy,
            bounds_policy=config.bounds_policy, refinement_bounds=config.refinement_bounds)
        report['driver_wall_seconds'] = perf_counter()-start
    report.update(refinement=profile.report(), driver_evaluations=result.evaluations,
                  initial_objective=float(result.initial_objective), final_objective=float(result.final_objective))
    assert result.evaluations == len(profile.rows)
    save()
    print(json.dumps(report['refinement']['census'], indent=2), flush=True)


if __name__ == '__main__':
    main()
