"""Frontend-independent immutable run snapshots and canonical execution."""
from dataclasses import asdict, dataclass
import hashlib
import json

from .config_paths import resolve_config_paths
from .config_compiler import CompiledConfig, compile_config
from .experimental_config import load_experimental_data
from .experimental import ExperimentalIntensityData
from .full_runner import run_full_refinement


def canonical_checkpoint_identity(compiled):
    # Compiled layout is authoritative: declaration order of bindings is irrelevant.
    return dict(format='qcbed-canonical-run', version=1,
        digest=hashlib.sha256(json.dumps(dict(state=asdict(compiled.state),
            runner=asdict(compiled.runner), layout=asdict(compiled.layout),
            experiment=asdict(compiled.experiment), run=asdict(compiled.configuration.run)),
            sort_keys=True,allow_nan=False).encode()).hexdigest())


@dataclass(frozen=True)
class RunSnapshot:
    """Compiled immutable config plus an immutable copy of parsed experimental data.

    Later configuration edits or edits to the source file cannot affect this run.
    Resume remains subject to the shared runner's configuration/data fingerprints.
    """
    compiled: CompiledConfig
    experimental: ExperimentalIntensityData

    def __post_init__(self):
        if not isinstance(self.compiled,CompiledConfig) or not isinstance(self.experimental,ExperimentalIntensityData):
            raise TypeError('RunSnapshot requires CompiledConfig and ExperimentalIntensityData')


def prepare_run_snapshot(config, *, base_directory):
    compiled = compile_config(resolve_config_paths(config,base_directory=base_directory))
    return RunSnapshot(compiled,load_experimental_data(compiled,base_directory=base_directory))


def run_snapshot(snapshot: RunSnapshot, *, resume=False, progress_hook=None):
    if not isinstance(snapshot,RunSnapshot):
        raise TypeError('snapshot must be RunSnapshot')
    compiled = snapshot.compiled
    from .config_targets import available_refinement_targets, target_label
    targets = [dict(label=target_label(compiled.configuration, i), category=i.category.value,
                    identifiers=i.target.identifiers, units=i.units)
               for i in available_refinement_targets(compiled) if i.active]
    return run_full_refinement(compiled.state,snapshot.experimental,compiled.runner,
        output_directory=compiled.configuration.run.output_directory,
        final_output_path=compiled.configuration.run.output_filename,
        checkpoint_identity=canonical_checkpoint_identity(compiled),resume=resume,progress_hook=progress_hook,
        results_metadata=dict(targets=targets))


def run_config(config, *, base_directory, resume=False, progress_hook=None):
    """Run an immutable canonical snapshot, with no YAML or temporary file."""
    return run_snapshot(prepare_run_snapshot(config,base_directory=base_directory),
                        resume=resume,progress_hook=progress_hook)
