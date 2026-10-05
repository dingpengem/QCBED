"""Reusable modern output locations and explicit legacy output backups."""
from pathlib import Path
from uuid import uuid4
from qcbed import resolve_config_paths


def existing_outputs(config,base):
    c=resolve_config_paths(config,base_directory=base)
    directory=Path(c.run.output_directory);file=Path(c.run.output_filename)
    return tuple(p for p in (directory,file) if p.exists() and not (p==file and file!=directory and directory.exists() and file.is_relative_to(directory)))


def check_replacement(config,base):
    c=resolve_config_paths(config,base_directory=base)
    directory=Path(c.run.output_directory)
    if directory.exists() and not directory.is_dir():raise ValueError('Run output directory is not a directory')
    if (directory/'.run.lock').exists():raise ValueError('Run lock exists; stop the active writer before replacing outputs')
    protected=[Path(base).resolve(),Path(__file__).resolve().parents[3],Path.home()]
    inputs=[Path(c.experiment.path).resolve()]
    if c.crystallography:inputs.append(Path(c.crystallography.structure.path).resolve())
    from .data_model import provenance
    meta=provenance(c)
    if meta.get('source'):inputs.append(Path(meta['source']['path']).resolve())
    if meta.get('overview'):inputs.extend(Path(meta['overview'][k]).resolve() for k in ('pattern','manifest'))
    for target in existing_outputs(c,base):
        resolved=target.resolve()
        if target.is_symlink():raise ValueError('Refusing to replace a symbolic-link output path')
        if any(p==resolved or p.is_relative_to(resolved) for p in protected+inputs):
            raise ValueError('Output replacement would move project/source/input files; select an output-only path')
    return c


def backup_outputs(config,base):
    c=check_replacement(config,base);moved=[]
    try:
        for path in existing_outputs(c,base):
            backup=path.with_name(path.name+'.backup-'+uuid4().hex)
            path.rename(backup);moved.append((path,backup))
    except BaseException:
        for path,backup in reversed(moved):backup.rename(path)
        raise
    return tuple(moved)


def validate_output_directory(path, base, *, create=False):
    """Probe writability without interpreting folder contents as run metadata."""
    import tempfile
    if not str(path).strip():raise ValueError('Select an output directory')
    directory=Path(path).expanduser()
    if not directory.is_absolute():directory=Path(base)/directory
    directory=directory.resolve()
    if directory.exists() and not directory.is_dir():raise ValueError('Output path is a file; select a directory')
    parent=directory
    while not parent.exists():parent=parent.parent
    if not parent.is_dir():raise ValueError('Cannot create output directory: an ancestor is a file')
    try:
        if create:directory.mkdir(parents=True,exist_ok=True)
        with tempfile.TemporaryFile(dir=directory if directory.exists() else parent) as probe:
            probe.write(b'QCBED write probe');probe.flush()
    except OSError as exc:
        raise ValueError(f'Cannot write/create output directory: {exc.strerror or exc}') from exc
    return directory


def modern_execution_config(config,base,*,resume=False):
    """Each Fresh Run owns a new child; never move or replace selected-folder data."""
    from dataclasses import replace
    from datetime import datetime,timezone
    import re
    c=resolve_config_paths(config,base_directory=base)
    root=Path(c.run.output_directory)
    if resume:
        candidates=sorted((p for p in root.iterdir() if p.is_dir() and not p.is_symlink()
            and re.fullmatch(r'qcbed-run-\d{8}T\d{12}Z-[0-9a-f]{32}',p.name)
            and (p/'checkpoint_latest.json').is_file()),reverse=True) if root.is_dir() else []
        directory=candidates[0] if candidates else root
    else:
        stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        directory=root/f'qcbed-run-{stamp}-{uuid4().hex}'
    return replace(c,run=replace(c.run,output_directory=str(directory),output_filename=str(directory/'final_out.dat')))
