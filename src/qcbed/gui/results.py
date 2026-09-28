"""Read-only durable artifact inspection; never evaluates a numerical trial."""
from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np
from qcbed.full_runner import load_refinement_checkpoint

COMPLETED_EVENTS = frozenset(('trial', 'initial_retained', 'retained_boundary',
                             'snapshot_verification', 'verification'))


def completed_events(history):
    """Finite completed evaluations only; stage acceptance is not an evaluation."""
    result = {}
    for row in history:
        if row.get('classification') not in COMPLETED_EVENTS:
            continue
        try:
            if not np.isfinite(float(row['objective'])):
                continue
            vector = row.get('physical_active_vector', row.get('active_vector'))
            if vector is None or not np.all(np.isfinite(np.asarray(vector, dtype=float))):
                continue
        except (ValueError, TypeError, KeyError):
            continue
        result[row['global_evaluation_index']] = row
    return tuple(result.values())


@dataclass(frozen=True)
class RunResults:
    directory: Path
    summary: dict
    parameters: tuple
    history: tuple
    comparison: dict | None
    warnings: tuple


def load_results(directory, *, evaluation=None):
    directory = Path(directory).resolve()
    def read(name):
        with (directory / name).open(encoding='utf-8') as stream:
            value = json.load(stream)
        if not isinstance(value, dict):
            raise ValueError(f'{name}: expected an object')
        return value
    try:
        cp = load_refinement_checkpoint(directory)
        if cp.get('checkpoint_identity', {}).get('format') != 'qcbed-canonical-run':
            raise ValueError('Not a recognizable canonical run')
        config = read('run_config.json')
        summary = read('summary.json') if (directory / 'summary.json').exists() else {}
        warnings, history = [], []
        if (directory / 'progress.jsonl').exists():
            # Only newline-terminated records are committed by the writer.
            lines = (directory / 'progress.jsonl').read_text(encoding='utf-8').splitlines(keepends=True)
            for number, line in enumerate(lines, 1):
                try:
                    if not line.endswith('\n'):
                        raise ValueError('Incomplete record')
                    row = json.loads(line)
                    if not isinstance(row, dict) or 'global_evaluation_index' not in row or 'classification' not in row:
                        raise ValueError('Invalid event')
                    history.append(row)
                except ValueError:
                    warnings.append(f'Unreadable progress record at line {number}; subsequent history omitted.')
                    break
        valid = completed_events(history)
        selected = next((r for r in valid if r['global_evaluation_index'] == evaluation), None) if evaluation is not None else (valid[-1] if valid else None)
        if evaluation is not None and selected is None:
            raise ValueError(f'Evaluation {evaluation} has no valid completed record')
        selected_index = selected['global_evaluation_index'] if selected else cp.get('retained_at_evaluation')
        vector = selected.get('physical_active_vector', selected.get('active_vector')) if selected else cp['retained_vector']
        objective = selected['objective'] if selected else cp.get('retained_objective')
        ynorm = selected.get('ynorm') if selected else cp.get('retained_ynorm')
        partial = summary.get('status') != 'completed' or not cp['completed'] or (directory / '.run.lock').exists()
        summary = dict(summary, partial=partial, output_directory=str(directory),
            checkpoint_status='Present (optimizer retained restart state; independent of displayed evaluation)',
            final_state_status='Present' if (directory / 'final_parameters.json').exists() else 'Unavailable',
            retained_objective=objective, retained_ynorm=ynorm,
            displayed_evaluation=selected_index, latest_valid_evaluation=valid[-1]['global_evaluation_index'] if valid else None,
            historical=bool(valid and selected_index != valid[-1]['global_evaluation_index']),
            recorded_evaluations=max([cp['evaluation_count']] + [r['global_evaluation_index'] for r in history]),
            recorded_elapsed_seconds=max([cp['elapsed_seconds'], summary.get('runtime_seconds', 0)] +
                                         [r.get('elapsed_seconds', 0) for r in history]),
            recorded_stage=selected.get('stage') if selected else cp.get('stage_name'),
            restart_count=cp['restart_count'])
        summary.setdefault('status', 'incomplete')
        if (directory / '.run.lock').exists():
            summary['status'] = 'Running (run lock present)'
        specs = cp['layout']['parameters']
        if len(specs) != len(vector):
            raise ValueError('Checkpoint parameter/vector length mismatch')
        targets = (config.get('results_metadata') or {}).get('targets', [])
        if targets and len(targets) != len(specs):
            raise ValueError('Stored target metadata length mismatch')
        bounds = config.get('refinement_bounds', {}).get('entries', [])
        parameters = []
        for index, (spec, value) in enumerate(zip(specs, vector)):
            info = targets[index] if targets else {}
            bound = next((b for b in bounds if b['category'] == spec['category'] and b['target'] == spec['target']), None)
            parameters.append(dict(label=info.get('label', f"{spec['category']} — legacy target {spec['target']}"),
                category=spec['category'], initial=spec['initial_value'],
                final=value if objective is not None else None,
                delta=value-spec['initial_value'] if objective is not None else None,
                units=info.get('units', 'unavailable'), bounds=bound))
        if specs and not targets:
            warnings.append('Stable target IDs were not stored by this older run; legacy physical target coordinates are shown.')
        if selected and 'physical_active_vector' not in selected and selected.get('vector_representation') != 'mapped physical active values':
            warnings.append('Older event stores an optimizer vector; mapped physical values were not separately recorded.')
        comparison = None
        # Atomic current-pattern file is usable only after its matching complete event.
        candidates = [(directory / 'latest_evaluation.npz', None)]
        for metadata_path in sorted(directory.glob('*_snapshot.json'), key=lambda p: p.stat().st_mtime_ns, reverse=True):
            try:
                metadata = read(metadata_path.name)
            except (ValueError, OSError) as exc:
                warnings.append(f'{metadata_path.name}: {exc}')
                continue
            if selected_index in (metadata.get('evaluation_index'), metadata.get('retained_at_evaluation')):
                candidates.append((metadata_path.with_suffix('.npz'), metadata))
        for path, metadata in candidates:
            if not path.exists():
                continue
            try:
                with np.load(path, allow_pickle=False) as arrays:
                    if metadata is None and int(arrays['evaluation_index']) != selected_index:
                        continue
                    if metadata is not None and 'evaluation_index' in arrays and int(arrays['evaluation_index']) != metadata.get('evaluation_index'):
                        continue  # Snapshot files are individually atomic; ignore a mismatched pair.
                    stored_vector = arrays['physical_vector'] if 'physical_vector' in arrays else arrays['vector']
                    if not np.array_equal(stored_vector, vector):
                        continue
                    if float(arrays['objective']) != objective or float(arrays['ynorm']) != ynorm:
                        continue
                    if not {'experimental', 'scaled', 'beams'} <= set(arrays.files):
                        continue
                    comparison = {k: arrays[k].copy() for k in ('experimental', 'scaled', 'beams')}
                    if 'scan_shape' in arrays:
                        comparison['scan_shape'] = arrays['scan_shape'].copy()
                    if 'scan_json' in arrays:
                        comparison['scan'] = json.loads(str(arrays['scan_json']))
                if comparison['experimental'].shape != comparison['scaled'].shape or comparison['scaled'].ndim != 2 or comparison['beams'].shape != (comparison['scaled'].shape[1], 2):
                    raise ValueError('Comparison array dimensions disagree')
                for array in comparison.values():
                    if isinstance(array, np.ndarray):
                        array.setflags(write=False)
                break
            except (OSError, ValueError, KeyError) as exc:
                comparison = None
                warnings.append(f'{path.name}: {exc}')
        if comparison is None:
            warnings.append('Pattern unavailable for this evaluation. No values were recomputed.')
        return RunResults(directory, summary, tuple(parameters), tuple(history), comparison, tuple(warnings))
    except (OSError, ValueError, KeyError, TypeError, AttributeError) as exc:
        raise ValueError(f'Cannot inspect canonical run {directory}: {exc}') from exc
