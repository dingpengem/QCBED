"""Canonical YAML saving; format conversion only, never numerical execution."""
from dataclasses import asdict
from enum import Enum
import os
from pathlib import Path
import tempfile

from .config_paths import resolve_config_paths


def _plain(value):
    if isinstance(value,Enum): return value.value
    if isinstance(value,dict): return {k:_plain(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)): return [_plain(v) for v in value]
    return value


def save_yaml_config(config, path, *, base_directory, overwrite=False):
    """Save As exclusively, or atomically replace on explicit overwrite=True.

    Both source-relative fields and a relative destination use the required
    absolute base_directory. YAML stores resolved paths. Loading reproduces
    resolve_config_paths(config, base_directory=...), not its relative spelling.
    No comments/formatting preservation. Returns the resolved configuration.
    """
    import yaml
    resolved = resolve_config_paths(config,base_directory=base_directory)
    destination = (Path(base_directory)/path).resolve()
    root = Path(__file__).resolve().parents[2]
    if any(destination.is_relative_to(root/name) for name in ('legacy_fortran','reference_data')):
        raise ValueError('cannot write protected legacy/reference files')
    if destination == Path(resolved.experiment.path):
        raise ValueError('configuration destination is the experimental data file')
    text = yaml.safe_dump(_plain(asdict(resolved)),sort_keys=False)
    # Stage beside the destination: replacement/link is on the same filesystem.
    descriptor,temporary = tempfile.mkstemp(prefix='.'+destination.name+'.',dir=destination.parent)
    try:
        with os.fdopen(descriptor,'w',encoding='utf-8') as stream:
            stream.write(text)
            stream.flush()
            os.fsync(stream.fileno())
        if overwrite:
            os.replace(temporary,destination)
        else:
            os.link(temporary,destination)  # Atomic exclusive publication.
    finally:
        Path(temporary).unlink(missing_ok=True)
    return resolved
