"""GUI readiness wording; canonical validation remains the compatibility authority."""
from dataclasses import replace
from pathlib import Path

from qcbed import (validate_config, resolve_config_paths, ValidationPhase,
                   ValidationIssue, ValidationLocation)

NO_CHECKPOINT = (
    'Resume: Not available — no previous run found'
)


def readiness(config, base_directory, *, resume=False):
    if config.crystallography is not None:
        from .output_paths import validate_output_directory,modern_execution_config
        from qcbed.modern_frontend import validate_modern
        try:
            validate_output_directory(config.run.output_directory,base_directory)
            execution=modern_execution_config(config,base_directory,resume=resume)
        except (ValueError,OSError) as exc:
            return (ValidationIssue('output_not_writable','error',ValidationPhase.FILES,
                ValidationLocation('run',('run',),'output_directory'),str(exc)),)
        if resume and not (Path(execution.run.output_directory)/'checkpoint_latest.json').is_file():
            return (ValidationIssue('resume_checkpoint_missing','error',ValidationPhase.EXECUTION,
                ValidationLocation('run',('run',),'output_directory'),NO_CHECKPOINT),)
        return validate_modern(execution,base_directory,check_files=True,resume=resume)
    phases = (ValidationPhase.FIELD, ValidationPhase.COMPILE,
              ValidationPhase.EXECUTION, ValidationPhase.FILES)
    validate = lambda **kw: validate_config(config,base_directory=base_directory,phases=phases,**kw)
    fresh = validate()
    if any(i.phase in (ValidationPhase.FIELD, ValidationPhase.COMPILE) and
           i.severity == 'error' for i in fresh):
        return fresh
    if base_directory is None or not Path(base_directory).is_absolute():
        return fresh
    resolved = resolve_config_paths(config, base_directory=base_directory)
    directory = Path(resolved.run.output_directory)
    checkpoint = directory / 'checkpoint_latest.json'
    # Missing checkpoints are a normal readiness state, not an attempted load.
    if not checkpoint.exists():
        message = (NO_CHECKPOINT if directory.exists() else
                   'Resume is unavailable: the output directory does not exist.')
        resume_issues = tuple(i for i in fresh if i.code != 'output_exists') + (
            ValidationIssue('resume_checkpoint_missing', 'error', ValidationPhase.EXECUTION,
                ValidationLocation('run', ('run',), 'output_directory'), message),)
    else:
        resume_issues = validate(resume=True)
        resume_issues = tuple(
            replace(i, message='Resume is unavailable: ' + i.message)
            if i.code == 'resume_incompatible' else i for i in resume_issues)
    if resume:
        return resume_issues
    fresh=tuple(i for i in fresh if i.code!='output_exists')
    try:
        from .output_paths import check_replacement
        check_replacement(config,base_directory)
    except (ValueError,OSError) as exc:
        fresh+=(ValidationIssue('output_replacement_unsafe','error',ValidationPhase.FILES,
            ValidationLocation('run',('run',),'output_directory'),str(exc)),)
    return fresh
