"""Toolkit-independent document adapter; scientific rules belong to Stage 6A."""
from dataclasses import dataclass, fields, replace
from pathlib import Path
import re

from qcbed import (Id, get_object, update_object, replace_object, set_refinement, validate_fields,
    validate_config, compile_config, available_refinement_targets, load_yaml_config,
    save_yaml_config, resolve_config_paths, ValidationIssue, ValidationLocation, ValidationPhase)


@dataclass(frozen=True)
class FieldSpec:
    section: str
    path: tuple
    name: str
    label: str
    kind: str

    @property
    def key(self): return self.path+(self.name,)


def editable_fields(config):
    result = [FieldSpec('Project',(),'title','Project title','text'),
        FieldSpec('Microscope',('microscope',),'energy_kv','Energy (kV)','real'),
        FieldSpec('Microscope',('microscope',),'zone_axis','Zone axis (3 integers)','int3')]
    for r in config.sample.regions:
        result.append(FieldSpec('Regions',('sample','regions',Id(r.id)),'depth_angstrom',f'Region {r.id} depth (Å)','real'))
    for name,label,kind in (
        ('scan_columns','Scan columns / NUMP','int'),('scan_rows','Scan rows / NUML','int'),
        ('basis_1','Basis 1 (3 components)','real3'),('basis_2','Basis 2 (3 components)','real3'),
        ('corner_a','Corner A (2 components)','real2'),('corner_b','Corner B (2 components)','real2'),
        ('corner_c','Corner C (2 components)','real2')):
        result.append(FieldSpec('Scan',('scan',),name,label,kind))
    result.extend((FieldSpec('Experimental data',('experiment',),'path','Experimental data file','file'),
        FieldSpec('Execution',('execution',),'execution_backend','Execution backend','backend'),
        FieldSpec('Execution',('execution',),'batch_size','Batch size (blank for serial)','optional_int'),
        FieldSpec('Outputs',('run',),'output_directory','Output directory','directory'),
        FieldSpec('Outputs',('run',),'output_filename','Output filename','save_file')))
    return tuple(result)


def display_value(value):
    if value is None: return ''
    if isinstance(value,tuple): return ', '.join(display_value(v) for v in value)
    return repr(value) if isinstance(value,float) else str(value)


def parse_text(text,kind):
    """Parse widget syntax only. Existing validators enforce scientific rules."""
    if kind in ('text','file','directory','save_file','backend'): return text
    if kind == 'optional_int' and not text.strip(): return None
    if kind in ('int','optional_int'):
        if not re.fullmatch(r'[+-]?[0-9]+',text.strip()): raise ValueError('Enter an integer, without decimals')
        return int(text)
    if kind == 'real': return float(text)
    if kind in ('real2','real3','int3'):
        parts = re.split(r'[\s,]+',text.strip())
        if len(parts) != int(kind[-1]): raise ValueError(f'Enter exactly {kind[-1]} components')
        return tuple(parse_text(p,'int' if kind == 'int3' else 'real') for p in parts)
    raise ValueError(f'Unknown widget input type {kind}')


def input_issue(message,section,path=(),field=None,code='invalid_editor_input',phase=ValidationPhase.FIELD):
    return ValidationIssue(code,'error',phase,ValidationLocation(section,tuple(path),field),str(message))


from qcbed.config_targets import target_label

def diagnostic_presentation(issue):
    """Editor wording only; retain the original diagnostic and runner severity."""
    phase = issue.phase
    message = issue.message
    if issue.code == 'output_exists':
        phase = ValidationPhase.EXECUTION
    labels = {ValidationPhase.FIELD: 'Configuration validity',
              ValidationPhase.COMPILE: 'Compile validity',
              ValidationPhase.EXECUTION: 'Execution readiness',
              ValidationPhase.FILES: 'External file/output readiness',
              ValidationPhase.EXPORT: 'Legacy export readiness'}
    severity = ('Run blocker' if issue.severity == 'error' and phase in
                (ValidationPhase.EXECUTION, ValidationPhase.FILES) else issue.severity)
    return severity, labels[phase], message


class EditorDocument:
    def __init__(self,config,path):
        self.path = Path(path).resolve()
        self.base_directory = self.path.parent
        self.config = config
        self.saved_config = config
        self.needs_save_as = False
        self.undo_stack = []
        self.redo_stack = []
        self.files_checked = False
        self.issues = ()
        self.validate()

    @classmethod
    def open(cls,path):
        import yaml
        try:
            return cls(load_yaml_config(path),path)
        except yaml.YAMLError as exc:
            raise ValueError(f'Invalid YAML in {path}: {exc}') from exc

    @property
    def dirty(self): return self.config != self.saved_config

    def apply(self,texts,refinement=()):
        """Atomically commit all valid drafts; keep scientific snapshot on failure.

        texts maps FieldSpec.key to text; refinement contains
        (ParameterTarget, selected, D1 text, D2 text, PLAMDA text).
        """
        candidate = self.config
        issues = []
        grouped = {}
        for spec in editable_fields(self.config):
            if spec.key not in texts: continue
            try: value = parse_text(texts[spec.key],spec.kind)
            except (ValueError,OverflowError) as exc:
                issues.append(input_issue(exc,spec.path[0] if spec.path else 'project',spec.path,spec.name))
                continue
            grouped.setdefault(spec.path,{})[spec.name] = value
        for path,changes in grouped.items():
            old = get_object(candidate,path)
            values = {f.name:getattr(old,f.name) for f in fields(old)}
            values.update(changes)
            if path == ('execution',) and values['execution_backend'] == 'batched' and values['batch_size'] is None:
                issues.append(input_issue('Enter an explicit batch size for batched execution', 'execution', path, 'batch_size'))
                continue
            local = validate_fields(type(old),values,section=path[0] if path else 'project',object_path=path)
            if local: issues.extend(local); continue
            try: candidate = replace_object(candidate,path,replace(old,**changes))
            except (ValueError,TypeError) as exc: issues.append(input_issue(exc,path[0] if path else 'project',path))
        # Raw export spellings describe the old associated value. Clear only when
        # that user-edited value changes; never silently alter unrelated metadata.
        raw_changes = {}
        for name,before,after in (
            ('title_raw',self.config.title,candidate.title),
            ('experimental_filename_raw',self.config.experiment.path,candidate.experiment.path),
            ('output_filename_raw',self.config.run.output_filename,candidate.run.output_filename)):
            if before != after: raw_changes[name] = None
        if raw_changes: candidate = update_object(candidate,('legacy_export',),**raw_changes)
        if not issues:
            for target,selected,d1,d2,step in refinement:
                previous = next((b for b in candidate.refinement.parameters if b.target == target),None)
                try:
                    if not selected:
                        if previous: candidate = set_refinement(candidate,target,enabled=False)
                        elif any(v.strip() for v in (d1,d2,step)):
                            raise ValueError('Select the binding to apply search settings, or clear these draft values')
                        continue
                    factors = (parse_text(d1,'real'),parse_text(d2,'real'))
                    plamda = parse_text(step,'real')
                    if previous and previous.range_factors == factors and previous.plamda == plamda: continue
                    candidate = set_refinement(candidate,target,enabled=True,range_factors=factors,plamda=plamda)
                except (ValueError,TypeError,OverflowError) as exc:
                    issues.append(input_issue(exc,'refinement',('refinement','parameters',target),'parameters'))
        if not issues:
            issues = list(validate_config(candidate,phases=(ValidationPhase.FIELD,ValidationPhase.COMPILE)))
        if any(i.severity == 'error' for i in issues):
            self.issues = tuple(issues) + self.committed_issues
            return False
        try:
            resolved = resolve_config_paths(candidate,base_directory=self.base_directory)
        except (ValueError, OSError, RuntimeError) as exc:
            self.issues = (input_issue(exc, 'project', field='path', phase=ValidationPhase.FILES),) + self.committed_issues
            return False
        self._record(resolved)
        self.files_checked = False
        self.validate()
        return True

    def _record(self, config):
        if config != self.config:
            self.undo_stack.append(self.config)
            self.redo_stack.clear()
            self.config = config

    def commit_structure(self, config):
        issues = validate_config(config, phases=(ValidationPhase.FIELD, ValidationPhase.COMPILE))
        if any(i.severity == 'error' for i in issues):
            self.issues = issues
            return False
        self._record(config)
        self.validate()
        return True

    def undo(self):
        if not self.undo_stack: return False
        self.redo_stack.append(self.config)
        self.config = self.undo_stack.pop()
        self.validate()
        return True

    def redo(self):
        if not self.redo_stack: return False
        self.undo_stack.append(self.config)
        self.config = self.redo_stack.pop()
        self.validate()
        return True

    def validate(self,*,check_files=False):
        phases = (ValidationPhase.FIELD,ValidationPhase.COMPILE,ValidationPhase.EXECUTION)
        if check_files: phases += (ValidationPhase.FILES,)
        if check_files:
            from .readiness import readiness
            self.issues = readiness(self.config, self.base_directory)
        else:
            self.issues = validate_config(self.config,phases=phases,base_directory=self.base_directory)
        self.committed_issues = self.issues
        self.files_checked = check_files
        return self.issues

    def summary(self):
        compiled = compile_config(self.config)
        execution_errors = any(i.severity == 'error' and (i.phase == ValidationPhase.EXECUTION or i.code == 'output_exists') for i in self.committed_issues)
        file_errors = any(i.severity == 'error' and i.phase == ValidationPhase.FILES for i in self.committed_issues)
        return dict(regions=len(self.config.sample.regions),stored_slices=len(compiled.state.model.slices),
            tilts=self.config.scan.scan_rows*self.config.scan.scan_columns,
            beams=len(self.config.experiment.beam_indices_hk),active_parameters=compiled.parameter_count,
            configuration='Valid (committed snapshot)',
            compile='Valid (committed snapshot)',
            execution='Blocked' if execution_errors else 'Preflight passed; numerical trials not evaluated',
            files=('Blocked' if file_errors else 'Preflight passed') if self.files_checked else 'Not checked — use Validate',
            bounds=compiled.bounds)

    def save(self,path=None,*,overwrite=False):
        destination = self.path if path is None else Path(path)
        resolved = save_yaml_config(self.config,destination,base_directory=self.base_directory,overwrite=overwrite)
        self.path = (self.base_directory/destination).resolve()
        self.base_directory = self.path.parent
        self.config = self.saved_config = resolved
        self.needs_save_as = False
        return resolved
