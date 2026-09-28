"""Explicit document-relative paths, shared by all configuration frontends."""
from dataclasses import replace
from pathlib import Path

from .config_compiler import compile_config


def resolve_config_paths(config, *, base_directory):
    """Return a snapshot with absolute paths; the base must itself be absolute.

    No existence checks, directory creation, environment or tilde expansion.
    Raw filename metadata keeps its trailing whitespace and consistency checks.
    """
    base = Path(base_directory)
    if not base.is_absolute():
        raise ValueError('base_directory must be absolute')
    compile_config(config)
    def resolved(value):
        return None if value is None else str((base / value).resolve())
    def raw_path(raw, value):
        return None if raw is None else resolved(value)+raw[len(raw.rstrip()):]
    result = replace(config,
        experiment=replace(config.experiment,path=resolved(config.experiment.path)),
        run=replace(config.run,output_directory=resolved(config.run.output_directory),
                    output_filename=resolved(config.run.output_filename)),
        legacy_export=replace(config.legacy_export,filename=resolved(config.legacy_export.filename),
            experimental_filename_raw=raw_path(config.legacy_export.experimental_filename_raw,config.experiment.path),
            output_filename_raw=raw_path(config.legacy_export.output_filename_raw,config.run.output_filename)))
    compile_config(result)
    return result
