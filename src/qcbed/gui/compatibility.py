"""Read-only compatibility planning and explicit Stage 5 frontend adapters."""
from pathlib import Path
from qcbed import compile_config, resolve_config_paths, validate_config, ValidationPhase
from qcbed.legacy_config import import_legacy_config, export_for007, legacy_export_report
from .document import EditorDocument


def import_document(source, experimental, policy):
    source = Path(source).resolve()
    config = import_legacy_config(source, experimental_path=source.parent/experimental, policy=policy)
    config = resolve_config_paths(config,base_directory=source.parent)
    document = EditorDocument(config,source.parent/'migrated.yml')
    document.saved_config = None
    document.needs_save_as = True
    return document


def compatibility_report(config, *, base_directory, bundle=True, filename_limits=None):
    resolved = resolve_config_paths(config,base_directory=base_directory)
    if resolved.crystallography is not None:
        from qcbed import ValidationIssue, ValidationLocation
        issue = ValidationIssue('crystallography_export_unsupported','error',ValidationPhase.EXPORT,
            ValidationLocation('crystallography',('crystallography',),None),
            'Original structure/hkl/site-group provenance is not representable; modern FOR007 export is deferred')
        return dict(represented=(),omissions={},diagnostics=(issue,),errors=(issue,))
    diagnostics = validate_config(resolved,phases=(ValidationPhase.EXPORT,),
                                  export_bundle=bundle,filename_limits=filename_limits)
    omissions = legacy_export_report(compile_config(resolved))
    errors = tuple(i for i in diagnostics if i.severity == 'error')
    return dict(represented=(
        'Microscope, ordered species/regions/slices/atoms, depths and scan geometry',
        'Ordered beams/reflections including duplicates and supported physical multipliers',
        'Compiled masks, linked absorption, D1/D2, PLAMDA and legacy optimizer controls',
        'DSMAX and retained legacy metadata',
        'NUMP/NUML/GRADIENTOPT and beam-major data in separate in.dat (bundle only)'),
        omissions=omissions, diagnostics=diagnostics, errors=errors)


def export_document(document, destination, *, bundle=False, filename_limits=None):
    """Exclusive publication: Stage 5 owns protection and representation rules."""
    config = resolve_config_paths(document.config,base_directory=document.base_directory)
    report = compatibility_report(config,base_directory=document.base_directory,
                                  bundle=bundle,filename_limits=filename_limits)
    if report['errors']:
        raise ValueError('\n'.join(i.message for i in report['errors']))
    path = Path(destination)
    if not path.is_absolute(): path = document.base_directory/path
    if path.resolve() == Path(config.experiment.path):
        raise ValueError('Cannot export over the canonical experimental data path')
    return export_for007(config,path,bundle=bundle,filename_limits=filename_limits)
