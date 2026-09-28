"""GUI readiness wording; canonical validation remains the compatibility authority."""
from dataclasses import replace
from pathlib import Path

from qcbed import (validate_config, resolve_config_paths, ValidationPhase,
                   ValidationIssue, ValidationLocation)

NO_CHECKPOINT = (
    'The output directory already exists but does not contain a resumable checkpoint. '
    'Choose a new output directory or remove/rename the existing directory.'
)


def readiness(config, base_directory, *, resume=False):
    phases = (ValidationPhase.FIELD, ValidationPhase.COMPILE,
              ValidationPhase.EXECUTION, ValidationPhase.FILES)
    fresh = validate_config(config, base_directory=base_directory, phases=phases)
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
        resume_issues = validate_config(config, base_directory=base_directory,
                                       resume=True, phases=phases)
        resume_issues = tuple(
            replace(i, message='Resume is unavailable: ' + i.message)
            if i.code == 'resume_incompatible' else i for i in resume_issues)
        reasons = [i.message for i in resume_issues if i.severity == 'error']
        message = (
            'Fresh run cannot use an existing run output directory. ' +
            ('Resume is unavailable. ' + '; '.join(reasons) +
             ' Choose a new output directory or remove/rename the existing directory.'
             if reasons else 'Choose a new output directory or use Resume.')
        )
    if resume:
        return resume_issues
    return tuple(replace(i, message=message if i.field == 'output_directory' else
                         'Fresh run cannot use an existing output file. Choose a new output filename.')
                 if i.code == 'output_exists' else i for i in fresh)
