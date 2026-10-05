"""Frontend-independent modern capability/readiness and stable-ID edits."""
from dataclasses import replace
from pathlib import Path
import os
from .config_validation import ValidationIssue,ValidationLocation,ValidationPhase,validate_config
from .crystallography_config import compile_config_from_sources,load_config_structure
from .config_paths import resolve_config_paths
from .config_editing import Id,update_object

CAPABILITIES=(
    ('Shared Debye–Waller B','READY','site-group','Å²','Explicit group B','Explicit target bounds; reference .5–1.5','Original/recast/all slice occurrences; Dynamic active-vector execution'),
    ('Requested depth','READY','sample','Å','Crystallography depth','Explicit target bounds; reference600–1200','Strict cyclic whole slices; effective depth derived'),
    ('Chemical occupancy / vacancy','READ-ONLY / DISPLAY ONLY','site-group','fraction','Structure/group occupancy','No validated refinement bounds','Chemical occupancy is not projected half-weight or legacy VAC; vacancy explicit independent coefficient refinement'),
    ('Absorption a/b/c','READY','site-group/species','legacy coefficients','Explicit group coefficients or disabled','Explicit nonnegative target bounds','One group per element; linked C not supported in modern bridge; explicit independent coefficient refinement'),
    ('Accelerating voltage','READY WITHOUT ABSORPTION','global','kV','Microscope','Explicit positive target bounds','Legacy absorption voltage update is unsupported; voltage refinement requires absorption disabled'),
    ('Scan A/B/C','READY','scan','basis coefficients','Prepared A10 geometry','Explicit signed target bounds','Explicit opt-in refinement of six A10 basis coefficients; preparation arrays remain fixed'),
    ('Reflection multipliers','NOT APPROPRIATE FOR CURRENT MODERN PATH','reflection','dimensionless','Unity in B5','Legacy only','No original-HKL multiplier binding in strict bridge'),
    ('Atom x/y','NOT APPROPRIATE FOR CURRENT MODERN PATH','original sites','fractional','Original structure','None','No symmetry-aware physical coordinate mode'),
    ('Legacy DILATION','DEFERRED','vacancy coupling','dimensionless','Legacy coefficient0','Coefficient domain0–1 only','Vacancy coupling, not lattice scale/strain; unit applied scale at zero vacancy'),
)


def validate_modern(config,base_directory,*,check_files=False,resume=False,execution=True):
    def issue(code,message,field='structure',phase=ValidationPhase.COMPILE,section='crystallography'):
        return ValidationIssue(code,'error',phase,ValidationLocation(section,(section,),field),str(message))
    if check_files and not config.run.output_directory.strip():
        return (issue('output_not_selected','Select and validate an output directory in D',
                      'output_directory',ValidationPhase.FILES,'run'),)
    try:
        resolved=resolve_config_paths(config,base_directory=base_directory)
        structure=load_config_structure(resolved)
        issues=validate_config(resolved,phases=(ValidationPhase.FIELD,ValidationPhase.COMPILE),structure=structure)
        if issues:return issues
        c=compile_config_from_sources(resolved)
    except (OSError,ValueError,TypeError) as exc:
        return (issue('crystallography_error',exc,getattr(exc,'field','structure')),)
    found=[]
    from .grouped_run import check_grouped,load_grouped_checkpoint
    if execution:
        try:check_grouped(c)
        except (ValueError,KeyError) as exc:found.append(issue('grouped_capability',exc,'parameters',ValidationPhase.EXECUTION,'refinement'))
    if not check_files:return tuple(found)
    from .canonical_run import RunSnapshot
    from .experimental_config import load_experimental_data
    data=None
    try:
        data=load_experimental_data(c,base_directory=base_directory)
    except (OSError,ValueError,TypeError) as exc:
        found.append(issue('experimental_data_error',exc,'path',ValidationPhase.FILES,'experiment'))
    if data is not None and not any(i.code=='grouped_capability' for i in found):
        try:check_grouped(c,data)
        except (ValueError,KeyError) as exc:
            found.append(issue('grouped_capability',exc,'parameters',ValidationPhase.EXECUTION,'refinement'))
    root=Path(__file__).resolve().parents[2]
    for field in ('output_directory','output_filename'):
        path=Path(getattr(resolved.run,field))
        if any(path.is_relative_to(root/n) for n in ('legacy_fortran','reference_data')):
            found.append(issue('protected_output','Output is in protected source/reference data',field,ValidationPhase.FILES,'run'))
        if not resume and path.exists():found.append(issue('output_exists','Fresh run requires a new output path',field,ValidationPhase.FILES,'run'))
        parent=path.parent
        while not parent.exists():parent=parent.parent
        if not parent.is_dir() or not os.access(parent,os.W_OK):found.append(issue('output_not_writable','Output parent is not writable',field,ValidationPhase.FILES,'run'))
    if resume and data is not None:
        try:
            if (Path(resolved.run.output_directory)/'.run.lock').exists():raise ValueError('Run lock exists; another writer or an interrupted process requires investigation')
            load_grouped_checkpoint(resolved.run.output_directory,RunSnapshot(c,data))
        except ValueError as exc:found.append(issue('resume_incompatible',exc,'output_directory',ValidationPhase.FILES,'run'))
    return tuple(found)


def edit_group(config,group_id,*,debye_waller_b):
    """Only the validated field; identity, membership, occupancy and absorption stay fixed."""
    return update_object(config,('crystallography','site_groups',Id(group_id)),debye_waller_b=debye_waller_b)


def set_target_bounds(config,target,lower,upper):
    from .configuration import TargetBound
    old=config.refinement.bounds
    replacement=TargetBound(target,lower,upper)
    items=tuple(replacement if b.target==target else b for b in old.specific) if any(b.target==target for b in old.specific) else old.specific+(replacement,)
    return update_object(config,('refinement','bounds'),specific=items)
