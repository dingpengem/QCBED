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
    identity = dict(format='qcbed-canonical-run', version=1,
        digest=hashlib.sha256(json.dumps(dict(state=asdict(compiled.state),
            runner=asdict(compiled.runner), layout=asdict(compiled.layout),
            experiment=asdict(compiled.experiment), run={k:v for k,v in asdict(compiled.configuration.run).items() if k!='max_stage_iterations' or v is not None}),
            sort_keys=True,allow_nan=False).encode()).hexdigest())
    if compiled.crystallography is not None and compiled.crystallography.site_groups:
        # Opt-in ownership is part of identity even when numerical arrays match.
        ownership = dict(original_lattice=compiled.crystallography.original_structure.direct_lattice.vectors.tolist(),
            original_sites=[asdict(s) for s in compiled.crystallography.original_structure.sites],
            zone_axis=compiled.crystallography.zone_axis.uvw,
            original_basis=(compiled.crystallography.reflection_basis.hkl_1,compiled.crystallography.reflection_basis.hkl_2),
            groups=[asdict(g) for g in compiled.crystallography.site_groups],
            target_aliases=[(asdict(a),asdict(b)) for a,b in compiled.crystallography.target_aliases
                if a.category.value == "debye_waller" or any(
                    p.category == a.category.value for p in compiled.layout.parameters)])
        if compiled.source_configuration.crystallography.projected_elastic:
            ownership['projected_elastic']=[{k:v for k,v in asdict(r).items() if k!='equality_group' or v is not None} for r in compiled.source_configuration.crystallography.projected_elastic]
            ownership['projected_elastic_adapter']='slice-eg-v1'
        return dict(format='qcbed-canonical-run',version=2,
            digest=hashlib.sha256(json.dumps(dict(base=identity,ownership=ownership),
                sort_keys=True,allow_nan=False).encode()).hexdigest())
    return identity


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
    from .crystallography_config import compile_config_from_sources
    compiled = compile_config_from_sources(resolve_config_paths(config,base_directory=base_directory))
    return RunSnapshot(compiled,load_experimental_data(compiled,base_directory=base_directory))


def run_snapshot(snapshot: RunSnapshot, *, resume=False, progress_hook=None, cancel_requested=None):
    if not isinstance(snapshot,RunSnapshot):
        raise TypeError('snapshot must be RunSnapshot')
    compiled = snapshot.compiled
    if compiled.crystallography is not None and compiled.crystallography.site_groups:
        from .grouped_run import run_grouped
        return run_grouped(snapshot,resume=resume,progress_hook=progress_hook,cancel_requested=cancel_requested)
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
