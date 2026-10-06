"""Structured frontend diagnostics, without changing existing exception APIs.

Default validation is filesystem-independent. File and export checks are explicit
phases. Readiness is a preflight, not a guarantee that every numerical trial will
succeed. Existing compiler/evaluator/runner checks remain authoritative.
"""
from dataclasses import MISSING, dataclass, fields, is_dataclass
from enum import Enum
from pathlib import Path
from typing import get_type_hints
import os

from .configuration import QcbedConfig, _typed
from .config_editing import Id
from .config_compiler import compile_config, _target_map


class ValidationPhase(str,Enum):
    FIELD = 'field'
    COMPILE = 'compile'
    EXECUTION = 'execution'
    FILES = 'files'
    EXPORT = 'export'


@dataclass(frozen=True)
class ValidationLocation:
    section: str
    object_path: tuple = ()
    field: str | None = None
    occurrence: int | None = None  # disambiguates duplicate IDs in one collection

    def __post_init__(self):
        object.__setattr__(self,'object_path',tuple(self.object_path))


@dataclass(frozen=True)
class ValidationIssue:
    code: str
    severity: str
    phase: ValidationPhase
    location: ValidationLocation
    message: str
    related_locations: tuple[ValidationLocation, ...] = ()

    def __post_init__(self):
        if self.severity not in ('error','warning','info'):
            raise ValueError('severity must be error, warning or info')
        object.__setattr__(self,'phase',ValidationPhase(self.phase))
        object.__setattr__(self,'related_locations',tuple(self.related_locations))

    @property
    def section(self): return self.location.section
    @property
    def object_path(self): return self.location.object_path
    @property
    def field(self): return self.location.field


def validate_fields(model_type, values, *, section='config', object_path=()):
    """Validate constructor keyword types before constructing an immutable object.

    Nested values are existing typed objects; no second dict-based config schema.
    Local domain checks reuse the constructor; cross-object checks use validate_config.
    """
    declarations = {f.name:f for f in fields(model_type)}
    hints = get_type_hints(model_type)
    issues = []
    for name,value in values.items():
        location = ValidationLocation(section,tuple(object_path),name)
        if name not in declarations:
            issues.append(ValidationIssue('unknown_field','error',ValidationPhase.FIELD,location,'Unknown field'))
            continue
        try: _typed(value,hints[name],name)
        except (TypeError,ValueError) as exc:
            issues.append(ValidationIssue('invalid_field_type','error',ValidationPhase.FIELD,location,str(exc)))
    for name,f in declarations.items():
        if name not in values and f.default is MISSING and f.default_factory is MISSING:
            issues.append(ValidationIssue('missing_field','error',ValidationPhase.FIELD,
                ValidationLocation(section,tuple(object_path),name),'Required field is missing'))
    if not issues:
        try:
            model_type(**values)  # Reuse constructor domain checks; discard this temporary object.
        except (TypeError,ValueError) as exc:
            issues.append(ValidationIssue('invalid_field_value','error',ValidationPhase.FIELD,
                ValidationLocation(section,tuple(object_path)),str(exc)))
    return tuple(issues)


def validate_config(config, *, phases=(ValidationPhase.FIELD,ValidationPhase.COMPILE,ValidationPhase.EXECUTION),
                    base_directory=None, resume=False, export_bundle=True, filename_limits=None, structure=None):
    """Return immutable diagnostic records; never write files or run propagation.

    The FILES phase requires an explicit absolute base_directory. EXPORT validates
    a prospective bundle by default; standalone widths must be supplied explicitly.
    Generic compiler fallbacks are retained for domain errors lacking a location.
    """
    phases = {ValidationPhase(p) for p in phases}
    issues = []
    def add(code,phase,section,path,field,message,*,related=(),occurrence=None,severity='error'):
        issues.append(ValidationIssue(code,severity,phase,ValidationLocation(section,tuple(path),field,occurrence),message,tuple(related)))
    if not isinstance(config,QcbedConfig):
        add('invalid_config_type',ValidationPhase.FIELD,'config',(),None,'Expected QcbedConfig; use validate_fields for uncommitted field input')
        return tuple(issues)
    if not phases: return ()
    def field_walk(value,section,path):
        values = {f.name:getattr(value,f.name) for f in fields(value)}
        found = validate_fields(type(value),values,section=section,object_path=path)
        issues.extend(found)
        if found: return
        for name,child in values.items():
            child_section = name if not path else section
            if is_dataclass(child): field_walk(child,child_section,path+(name,))
            elif isinstance(child,tuple):
                for i,item in enumerate(child):
                    if is_dataclass(item):
                        selector = Id(item.id) if hasattr(item,'id') else (Id(item.species,'species') if hasattr(item,'species') else i)
                        field_walk(item,child_section,path+(name,selector))
    field_walk(config,'config',())
    if issues or phases == {ValidationPhase.FIELD}: return tuple(issues)
    if config.crystallography is not None:
        from .crystallography_config import compile_config_from_sources, CrystallographyConfigError
        try:
            from .config_paths import resolve_config_paths
            resolved = resolve_config_paths(config,base_directory=base_directory) if base_directory is not None else config
            compiled = (compile_config_from_sources(resolved) if ValidationPhase.FILES in phases
                        else compile_config(resolved,structure=structure))
        except (ValueError,TypeError,OSError) as exc:
            field = exc.field if isinstance(exc,CrystallographyConfigError) else 'structure'
            phase = ValidationPhase.FILES if field == 'structure.path' else ValidationPhase.COMPILE
            add('crystallography_error',phase,'crystallography',('crystallography',),field,str(exc))
            return tuple(issues)
        if ValidationPhase.EXPORT in phases:
            add('crystallography_export_unsupported',ValidationPhase.EXPORT,'crystallography',
                ('crystallography',),None,'Original structure/hkl provenance is not representable; crystallography export is deferred')
        if compiled.crystallography.site_groups and phases & {ValidationPhase.EXECUTION, ValidationPhase.FILES}:
            # Physical ownership must not be discarded by validating a derived
            # legacy sample against the legacy checkpoint format.
            if ValidationPhase.FILES in phases:
                from .modern_frontend import validate_modern
                return tuple(issues) + validate_modern(resolved,base_directory,
                    check_files=True,resume=resume,execution=ValidationPhase.EXECUTION in phases)
            from .grouped_run import check_grouped
            try:
                check_grouped(compiled)
            except ValueError as exc:
                add('grouped_capability',ValidationPhase.EXECUTION,'refinement',
                    ('refinement',),'parameters',str(exc))
            return tuple(issues)
        # Pure modern configuration validation does not require a run destination.
        preview=compiled.configuration
        if not preview.run.output_directory.strip() and not phases & {ValidationPhase.EXECUTION,ValidationPhase.FILES}:
            from dataclasses import replace
            preview=replace(preview,run=replace(preview.run,output_directory='.'))
        return tuple(issues) + validate_config(preview,phases=phases-{ValidationPhase.EXPORT},base_directory=base_directory,
                               resume=resume,export_bundle=export_bundle,filename_limits=filename_limits)
    # Every dependent phase requires a valid compilation; expose its blockers too.
    C = ValidationPhase.COMPILE
    def unique(items,path,key='id',seen=None):
        seen = {} if seen is None else seen
        for i,item in enumerate(items):
            value = getattr(item,key)
            loc = ValidationLocation('sample',path+(Id(value,key),),key,i)
            if not value:
                add('empty_id',C,'sample',loc.object_path,key,'ID must not be empty',occurrence=i)
            elif value in seen:
                add('duplicate_id',C,'sample',loc.object_path,key,f'Duplicate {key}: {value}',related=(seen[value],),occurrence=i)
            else: seen[value] = loc
        return seen
    sample = config.sample
    if sample is None:
        add("missing_sample",C,"sample",("sample",),None,"sample is required without crystallography")
        return tuple(issues)
    if config.schema_version != 1:
        add('unsupported_schema_version',C,'config',(),'schema_version','Only canonical schema version 1 is supported')
    if not sample.species:
        add('empty_species',C,'sample',('sample',),'species','At least one species is required')
    if not sample.regions:
        add('empty_regions',C,'sample',('sample',),'regions','At least one region is required')
    for sp in sample.species:
        if sp.atomic_number <= 0:
            add('invalid_atomic_number',C,'sample',('sample','species',Id(sp.id)),'atomic_number','Atomic number must be positive')
    legacy = config.legacy_export
    for name,value in (('title_raw',config.title),('output_filename_raw',config.run.output_filename),
                       ('experimental_filename_raw',config.experiment.path)):
        raw = getattr(legacy,name)
        if raw is not None and raw.rstrip() != value:
            add('legacy_metadata_mismatch',C,'legacy_export',('legacy_export',),name,'Raw legacy metadata disagrees with the canonical value')
    unique(sample.species,('sample','species'))
    unique(sample.regions,('sample','regions'))
    species = {sp.id for sp in sample.species}
    slices_seen = {}
    for region in sample.regions:
        rp = ('sample','regions',Id(region.id))
        if not region.slices:
            add('empty_slices',C,'sample',rp,'slices','Each region requires at least one stored slice')
        unique(region.slices,rp+('slices',),seen=slices_seen)
        unique(region.reflections,rp+('reflections',))
        unique(region.debye_waller,rp+('debye_waller',),key='species')
        thermal = {v.species for v in region.debye_waller}
        if thermal != species:
            add('thermal_species_mismatch',C,'sample',rp,'debye_waller',f'Expected one entry per species; missing={sorted(species-thermal)}, unknown={sorted(thermal-species)}')
        for stored in region.slices:
            sp = rp+('slices',Id(stored.id))
            unique(stored.atoms,sp+('atoms',))
            for atom in stored.atoms:
                if atom.species not in species:
                    add('unknown_species',C,'sample',sp+('atoms',Id(atom.id)),'species',f'Unknown species {atom.species}')
    unique(sample.absorption.coefficients,('sample','absorption','coefficients'),key='species')
    if sample.absorption.enabled and {v.species for v in sample.absorption.coefficients} != species:
        add('absorption_species_mismatch',C,'sample',('sample','absorption'),'coefficients','Absorption requires exactly one entry per species')
    for name in ('scan_columns','scan_rows'):
        if getattr(config.scan,name) <= 0:
            add('invalid_dimension',C,'scan',('scan',),name,'Scan dimension must be positive')
    if not issues:
        try:
            coefficients = tuple(next(c for c in sample.absorption.coefficients if c.species == sp.id)
                for sp in sample.species) if sample.absorption.enabled else ()
            targets = _target_map(config,coefficients)
            for name,items in (('parameters',config.refinement.parameters),('specific',config.refinement.bounds.specific)):
                seen = {}
                base = ('refinement','parameters') if name == 'parameters' else ('refinement','bounds','specific')
                for i,item in enumerate(items):
                    loc = ValidationLocation('refinement',base+(item.target,),'target',i)
                    if item.target not in targets:
                        add('unknown_target',C,'refinement',loc.object_path,'target','Unresolved parameter target',occurrence=i)
                    if item.target in seen:
                        add('duplicate_target',C,'refinement',loc.object_path,'target','Duplicate parameter target',related=(seen[item.target],),occurrence=i)
                    seen[item.target] = loc
                    if name == 'parameters' and item.target in targets:
                        category,t = targets[item.target]
                        if category == 'absorption' and t[1] == 2 and coefficients[t[0]].c_mode == 'linked_to_b':
                            add('linked_target',C,'refinement',loc.object_path,'target','Linked absorption C cannot be independently refined')
                        if category == 'atom_coordinate':
                            stored = tuple(s for r in sample.regions for s in r.slices)[t[0]]
                            if stored.record_type == 'G':
                                add('gas_target',C,'refinement',loc.object_path,'target','G records cannot have coordinate refinement masks')
        except (ValueError,TypeError) as exc:
            add('compile_error',C,'config',(),None,str(exc))
    if issues: return tuple(issues)
    try: compiled = compile_config(config)
    except (ValueError,TypeError) as exc:
        add('compile_error',C,'config',(),None,str(exc))
        return tuple(issues)
    if ValidationPhase.EXECUTION in phases:
        from .refinement import apply_refinement_parameters
        from .refinement_categories import LEGACY_OPTIMIZER_PARAMETER_LIMIT
        E = ValidationPhase.EXECUTION
        mapped = apply_refinement_parameters(compiled.state,compiled.packed_starting_vector)
        if mapped.model.energy_kv <= 0:
            add('nonpositive_voltage',E,'microscope',('microscope',),'energy_kv','Mapped electron energy must be positive')
        for i,r in enumerate(mapped.regions):
            path = ('sample','regions',Id(sample.regions[i].id))
            if r.vacancy_flag != 0 or r.vacancy_values != (0.,) or r.free_surface[0] != 0:
                add('unsupported_region_state',E,'sample',path,'free_surface' if r.free_surface[0] != 0 else 'vacancy','Current evaluator does not support this vacancy/free-surface state')
            if r.depth_angstrom <= 0:
                add('nonpositive_depth',E,'sample',path,'depth_angstrom','Mapped depth must be positive')
        for region in sample.regions:
            for stored in region.slices:
                if stored.record_type == 'G' and stored.atoms:
                    add('unsupported_gas_atoms',E,'sample',('sample','regions',Id(region.id),'slices',Id(stored.id)),'record_type','G atom randomization is unsupported')
        if sample.absorption.enabled and compiled.state.model.energy_refined and compiled.optimizer.loops > 0:
            add('voltage_absorption_refinement',E,'microscope',('microscope',),'energy_kv','Changing voltage with absorption is unsupported')
        for name in ('scan_columns','scan_rows'):
            if getattr(config.scan,name) == 1:
                add('singleton_scan',E,'scan',('scan',),name,'Singleton scan dimensions are undefined by the existing scan generator')
        n = compiled.parameter_count
        if n > LEGACY_OPTIMIZER_PARAMETER_LIMIT or (compiled.optimizer.loops > 0 and n == 0):
            add('optimizer_parameter_count',E,'refinement',('refinement',),'parameters',f'Current driver requires at most {LEGACY_OPTIMIZER_PARAMETER_LIMIT} parameters, and a nonempty layout for optimization')
        if 'AMOEBA' in compiled.optimizer.stages:
            for info,p in zip(compiled.layout.parameters,compiled.range_factors):
                if p[0] == p[1]:
                    add('singular_range_factors',E,'refinement',('refinement',),'parameters',f'Singular D1/D2 for {info.category} {info.target}')
    if ValidationPhase.FILES in phases:
        F = ValidationPhase.FILES
        if base_directory is None or not Path(base_directory).is_absolute():
            add('base_directory_required',F,'project',(),'base_directory','An explicit absolute base_directory is required')
        else:
            from .config_paths import resolve_config_paths
            from .experimental_config import load_experimental_data
            resolved = compile_config(resolve_config_paths(config,base_directory=base_directory))
            data = None
            try: data = load_experimental_data(resolved,base_directory=base_directory)
            except (OSError,ValueError,TypeError) as exc:
                add('experimental_data_error',F,'experiment',('experiment',),'path',str(exc))
            root = Path(__file__).resolve().parents[2]
            for field,value in (('output_directory',resolved.configuration.run.output_directory),('output_filename',resolved.configuration.run.output_filename)):
                path = Path(value)
                if any(path.is_relative_to(root/n) for n in ('legacy_fortran','reference_data','tests/fixtures')):
                    add('protected_output',F,'run',('run',),field,'Output is inside a protected source/reference directory')
                if not resume and path.exists():
                    add('output_exists',F,'run',('run',),field,'Fresh run requires a new output path')
                parent = path.parent
                while not parent.exists(): parent = parent.parent
                if not parent.is_dir() or not os.access(parent,os.W_OK):
                    add('output_not_writable',F,'run',('run',),field,'Output parent is not a writable directory')
            if resume and data is not None:
                from .canonical_run import canonical_checkpoint_identity
                from .full_runner import load_refinement_checkpoint, _digest
                try:
                    cp = load_refinement_checkpoint(resolved.configuration.run.output_directory)
                    if cp.get('checkpoint_identity') != canonical_checkpoint_identity(resolved) or cp['problem_digest'] != _digest(resolved.state,data) or cp['input_files'] != {}:
                        raise ValueError('Canonical checkpoint configuration/data identity mismatch')
                except (OSError,ValueError,KeyError,TypeError) as exc:
                    add('resume_incompatible',F,'run',('run',),'output_directory',str(exc))
    if ValidationPhase.EXPORT in phases:
        from .legacy_config import _for007_text, legacy_export_report
        X = ValidationPhase.EXPORT
        if config.legacy_export.dsmax_inverse_angstrom is None:
            add('missing_dsmax',X,'legacy_export',('legacy_export',),'dsmax_inverse_angstrom','DSMAX must be explicitly supplied for FOR007 export')
        else:
            try:
                if not export_bundle:
                    if config.experiment.header != 'validate': raise ValueError('Standalone export requires header-validation mode')
                    if filename_limits is None or len(filename_limits) != 2 or any(type(v) is not int or v < 1 for v in filename_limits):
                        raise ValueError('Standalone export requires verified filename limits')
                    if len(config.experiment.path) > filename_limits[0] or len(config.run.output_filename) > filename_limits[1]:
                        raise ValueError('Filenames exceed target widths')
                _for007_text(compiled,experimental_filename='in.dat' if export_bundle else config.experiment.path,
                             output_filename='out.dat' if export_bundle else config.run.output_filename)
            except (ValueError,TypeError) as exc:
                add('legacy_not_representable',X,'legacy_export',('legacy_export',),None,str(exc))
        for field,value in legacy_export_report(compiled).items():
            add('legacy_omission',X,'legacy_export',('legacy_export',),field,str(value),severity='warning')
    return tuple(issues)
