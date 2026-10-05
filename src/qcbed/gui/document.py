"""Toolkit-independent document adapter over canonical scientific validation."""
from dataclasses import dataclass, fields, replace
from pathlib import Path
import re

from qcbed import (Id, get_object, update_object, replace_object, set_refinement, validate_fields,
    validate_config, load_yaml_config,
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
    for r in (() if config.sample is None else config.sample.regions):
        result.append(FieldSpec('Regions',('sample','regions',Id(r.id)),'depth_angstrom',f'Region {r.id} depth (Å)','real'))
    for name,label,kind in (
        ('scan_columns','Scan columns / NUMP','int'),('scan_rows','Scan rows / NUML','int'),
        ('basis_1','Basis 1 (3 components)','real3'),('basis_2','Basis 2 (3 components)','real3'),
        ('corner_a','Corner A (2 components)','real2'),('corner_b','Corner B (2 components)','real2'),
        ('corner_c','Corner C (2 components)','real2')):
        if hasattr(config.scan,name):result.append(FieldSpec('Scan',('scan',),name,label,kind))
    result.extend((FieldSpec('Experimental data',('experiment',),'path','Experimental data file','file'),
        FieldSpec('Execution',('execution',),'execution_backend','Execution mode','backend'),
        FieldSpec('Execution',('execution',),'batch_size','Batch size (blank for serial)','optional_int'),
        FieldSpec('Outputs',('run',),'output_directory','Output directory','directory'),
        FieldSpec('Outputs',('run',),'output_filename','Output filename','save_file')))
    if config.crystallography is not None:
        result=[f for f in result if f.key not in (('microscope','zone_axis'),('experiment','path'),('run','output_filename'))]
        result.extend((FieldSpec('Structure',('crystallography','structure'),'path','Structure file','file'),
            FieldSpec('Structure',('crystallography','structure'),'format','Structure format (blank: infer from file)','optional_text'),
            FieldSpec('Structure',('crystallography',),'zone_axis','Zone axis [u,v,w]','int3'),
            FieldSpec('Structure',('crystallography','cbed_basis'),'g1_hkl','g1 hkl (original/conventional cell)','int3'),
            FieldSpec('Structure',('crystallography','cbed_basis'),'g2_hkl','g2 hkl (original/conventional cell)','int3'),
            FieldSpec('Numerical / Simulation cell',('crystallography','mesh'),'mesh_x','Transverse simulation mesh X (points)','int'),
            FieldSpec('Numerical / Simulation cell',('crystallography','mesh'),'mesh_y','Transverse simulation mesh Y (points)','int'),
            FieldSpec('Numerical / Simulation cell',('crystallography',),'dsap_inverse_angstrom','Reciprocal cutoff (Å⁻¹)','real'),
            FieldSpec('Numerical / Simulation cell',('crystallography',),'slice_boundaries_fractional','Slice boundaries','boundaries'),
            FieldSpec('Numerical / Simulation cell',('crystallography',),'depth_angstrom','Requested depth (Å)','real'),
            FieldSpec('Microscope',('experiment',),'gradient_option','Signal mode','signal')))
    for name,label,kind in (('method','Workflow','workflow'),('loops','Loops','int'),
        ('maxfsimp','Stage iteration / soft evaluation budget','int'),('maxfpara','Coordinate-search budget','int'),('seed','Random seed','int'),
        ('maxtemp','ANNEAL maximum temperature','real'),('mintemp','ANNEAL minimum temperature','real'),
        ('coolrate','Cooling rate','real'),('ftol','Objective-spread tolerance','real'),('chimin','Minimum objective','real')):
        result.append(FieldSpec('Optimizer / Run',('refinement','optimizer'),name,label,kind))
    if config.crystallography is not None:
        result.extend((
            FieldSpec('Execution',('run',),'max_total_evaluations','Maximum total evaluations (blank: unlimited)','optional_int'),
            FieldSpec('Execution',('run',),'max_wall_time_seconds','Maximum wall time (seconds; blank: unlimited)','optional_real')))
    return tuple(result)


def display_value(value):
    if value is None: return ''
    if isinstance(value,tuple): return ', '.join(display_value(v) for v in value)
    return repr(value) if isinstance(value,float) else str(value)


def parse_text(text,kind):
    """Parse widget syntax only. Existing validators enforce scientific rules."""
    if kind=='boundaries':
        import math
        parts=re.split(r'[\s,]+',text.strip().strip('[]').strip())
        try: values=tuple(float(v) for v in parts)
        except ValueError: raise ValueError('Slice boundaries must start at 0, end at 1, and contain at least 2 strictly increasing values.') from None
        if (len(values)<2 or not all(math.isfinite(v) for v in values) or
                values[0]!=0 or values[-1]!=1 or
                any(b<=a for a,b in zip(values,values[1:]))):
            raise ValueError('Slice boundaries must start at 0, end at 1, and contain at least 2 strictly increasing values.')
        return values
    if kind=='optional_text':return text.strip() or None
    if kind=='signal':return {'Conventional':0,'Gradient':1}[text]
    if kind=='workflow':return {'ANNEAL → AMOEBA':'ANSI','ANNEAL':'ANNE','AMOEBA':'SIMP','PARAMIN':'PARA','SIMPARA':'SPAR'}.get(text,text)
    if kind in ('text','file','directory','save_file','backend'): return text
    if kind == 'optional_int' and not text.strip(): return None
    if kind in ('int','optional_int'):
        if not re.fullmatch(r'[+-]?[0-9]+',text.strip()): raise ValueError('Enter an integer, without decimals')
        return int(text)
    if kind == 'optional_real':return float(text) if text.strip() else None
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
        if config.crystallography is not None:
            config=replace(config,run=replace(config.run,
                output_filename=str(Path(config.run.output_directory)/'final_out.dat')))
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

    def apply(self,texts,refinement=(),bounds=(),*,stage=None):
        """Atomically commit all valid drafts; keep scientific snapshot on failure.

        texts maps FieldSpec.key to text; refinement contains
        (ParameterTarget, selected, D1 text, D2 text, PLAMDA text).
        """
        if stage in (0,1,"simulation","execution"):refinement,bounds=(),()
        candidate = self.config
        issues = []
        grouped = {}
        for spec in editable_fields(self.config):
            if stage=='simulation' and spec.section not in ('Microscope','Numerical / Simulation cell'):continue
            if stage=='refinement':continue
            if stage=='execution' and spec.section != 'Optimizer / Run':continue
            if stage==0 and spec.section not in ('Data Prep','Scan','Experimental data'):continue
            if stage==1 and spec.section not in ('Structure','Site groups'):continue
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
        if candidate.crystallography is not None:
            candidate=update_object(candidate,('microscope',),zone_axis=candidate.crystallography.zone_axis)
            candidate=replace(candidate,run=replace(candidate.run,
                output_filename=str(Path(candidate.run.output_directory)/'final_out.dat')))
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
                        elif candidate.crystallography is None and any(v.strip() for v in (d1,d2,step)):
                            raise ValueError('Select the binding to apply search settings, or clear these draft values')
                        continue
                    missing=[name for name,value in (('D1',d1),('D2',d2),('PLAMDA',step)) if not value.strip()]
                    if missing:
                        label='Requested depth' if target.category.value=='depth' else (
                            'Site group '+target.identifiers[-1]+' Debye–Waller B' if target.category.value=='debye_waller' else target.category.value)
                        raise ValueError(label+' refinement requires '+', '.join(missing)+'.')
                    factors = (parse_text(d1,'real'),parse_text(d2,'real'))
                    plamda = parse_text(step,'real')
                    if previous and previous.range_factors == factors and previous.plamda == plamda: continue
                    candidate = set_refinement(candidate,target,enabled=True,range_factors=factors,plamda=plamda)
                except (ValueError,TypeError,OverflowError) as exc:
                    issues.append(input_issue(exc,'refinement',('refinement','parameters',target),'parameters'))
        if not issues:
            for target,lo,hi in bounds:
                try:
                    from qcbed.modern_frontend import set_target_bounds
                    candidate=set_target_bounds(candidate,target,parse_text(lo,'real') if lo.strip() else None,parse_text(hi,'real') if hi.strip() else None)
                except (ValueError,TypeError) as exc:issues.append(input_issue(exc,'refinement',('refinement','bounds',target),'bounds'))
            if stage==1 and candidate.crystallography:
                from .workflow import validate_atomic_model
                try:validate_atomic_model(candidate)
                except (ValueError,TypeError,OSError) as exc:
                    issues.append(input_issue(exc,'crystallography',('crystallography',)))
            elif stage=='simulation':
                # Validate the forward model without requiring optimizer bindings.
                from qcbed.crystallography_config import compile_config_from_sources
                try:
                    compile_config_from_sources(replace(candidate,refinement=replace(candidate.refinement,
                        parameters=(),bounds=type(candidate.refinement.bounds)())))
                except (ValueError,TypeError,OSError) as exc:
                    issues.append(input_issue(exc,'simulation'))
            elif stage=='refinement':
                from qcbed.crystallography_config import compile_config_from_sources
                from qcbed.crystallography_refinement import physical_target_layout
                import math
                try:
                    layout=physical_target_layout(compile_config_from_sources(candidate))
                    if not layout:raise ValueError('Select at least one complete supported refinement binding.')
                    for row in layout:
                        lo,hi,value=row['lower'],row['upper'],row['initial_value']
                        if lo is None or hi is None or not all(math.isfinite(v) for v in (lo,hi)) or not lo<hi or not lo<=value<=hi:
                            raise ValueError(row['target_id']+': Hard bounds must be ordered and contain the starting value.')
                except (ValueError,TypeError,OSError) as exc:
                    issues.append(input_issue(exc,'refinement'))
            elif stage not in (0,1):
                issues.extend(self._validate_config(candidate,execution=False))
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

    @staticmethod
    def thermal_message(config):
        for g in config.crystallography.site_groups:
            if g.debye_waller_b is None:
                return f'Debye–Waller B is required for site group {g.id}.'
        return 'Define physical site groups and Debye–Waller B (Å²) before validating B. Atomic Model.'

    @staticmethod
    def thermal_incomplete(config):
        c=config.crystallography
        return c is not None and (not c.site_groups and not c.debye_waller or
            bool(c.site_groups) and any(g.debye_waller_b is None for g in c.site_groups))

    def _validate_config(self,config,*,execution=True):
        if config.crystallography is not None:
            try: parse_text(', '.join(map(str,config.crystallography.slice_boundaries_fractional)),'boundaries')
            except ValueError as exc:
                return (input_issue(exc,'crystallography',('crystallography',),'slice_boundaries_fractional'),)
        if self.thermal_incomplete(config):
            return (input_issue(self.thermal_message(config),
                'crystallography',('crystallography',),'site_groups'),)

        if config.crystallography is not None:
            from qcbed.modern_frontend import validate_modern
            return validate_modern(config,self.base_directory,execution=execution)
        phases=(ValidationPhase.FIELD,ValidationPhase.COMPILE)+( (ValidationPhase.EXECUTION,) if execution else ())
        return validate_config(config,phases=phases,base_directory=self.base_directory)

    def _record(self, config):
        if config != self.config:
            self.undo_stack.append(self.config)
            self.redo_stack.clear()
            self.config = config

    def commit_preparation(self, config):
        """Commit independently validated A inputs while B may be incomplete."""
        from .data_model import read_data
        old=self.config
        crystal=config.crystallography
        if crystal is not None and old.crystallography is not None:
            crystal=replace(crystal,beams=old.crystallography.beams,
                            scan_basis_transform=old.crystallography.scan_basis_transform)
        comparison=replace(config,experiment=old.experiment,scan=old.scan,
                           crystallography=crystal,legacy_export=old.legacy_export)
        if comparison!=old:raise ValueError('Preparation edit may change only A inputs and their legacy filename metadata')
        read_data(config,self.base_directory)
        self._record(config);self.files_checked=False;self.validate()
        return True

    def load_structure_file(self, path):
        """Commit a readable ASE structure independently of downstream model readiness."""
        from qcbed.crystallography_config import load_config_structure
        path = Path(path).expanduser()
        if not path.is_absolute():
            path = self.base_directory / path
        if str(path.resolve()) == self.config.crystallography.structure.path:
            return
        source = replace(self.config.crystallography.structure,
                         path=str(path.resolve()), format=None)
        candidate = replace(self.config, crystallography=replace(
            self.config.crystallography, structure=source))
        structure=load_config_structure(candidate)  # Failure leaves the document/history unchanged.
        from qcbed.configuration import AbsorptionConfig, AtomicSiteGroup
        # Group by explicit source metadata, not inferred crystallographic symmetry.
        keys=tuple(dict.fromkeys((s.symbol,s.occupancy,s.debye_waller_b) for s in structure.sites))
        import numpy as np
        fcc_al=(len(keys)==1 and keys[0][0]=='Al' and len(structure.sites)==4 and
            np.allclose(structure.direct_lattice.lengths,structure.direct_lattice.lengths[0]) and
            np.allclose(structure.direct_lattice.angles,(90,90,90)) and
            sorted(s.fractional_position for s in structure.sites)==
            [(0.,0.,0.),(0.,.5,.5),(.5,0.,.5),(.5,.5,0.)])
        groups=tuple(AtomicSiteGroup(('Al_fcc' if fcc_al else f'{symbol}_group_{i+1}'),symbol,
            tuple(j for j,s in enumerate(structure.sites)
                  if (s.symbol,s.occupancy,s.debye_waller_b)==(symbol,occupancy,b)),
            b,occupancy) for i,(symbol,occupancy,b) in enumerate(keys))
        candidate=replace(candidate,crystallography=replace(candidate.crystallography,
            site_groups=groups,debye_waller=(),absorption=AbsorptionConfig()),
            refinement=replace(candidate.refinement,parameters=(),
                bounds=replace(candidate.refinement.bounds,specific=())))

        self._record(candidate)
        self.files_checked = False
        self.validate()

    def edit_site_group(self, index, group):
        from qcbed.crystallography_config import load_config_structure
        groups=list(self.config.crystallography.site_groups)
        old=groups[index]
        if group.id!=old.id and (self.config.refinement.parameters or self.config.refinement.bounds.specific):
            raise ValueError('Remove refinement bindings/bounds before renaming a group.')
        groups[index]=group
        if len({g.id for g in groups})!=len(groups):raise ValueError('Duplicate group ID')
        structure=load_config_structure(self.config)
        members=[i for g in groups for i in g.site_indices]
        if sorted(members)!=list(range(len(structure.sites))):
            raise ValueError('Groups must cover each source site exactly once')
        if any(structure.sites[i].symbol!=g.element for g in groups for i in g.site_indices):
            raise ValueError('Group element must match its member sites')
        self._record(replace(self.config,crystallography=replace(
            self.config.crystallography,site_groups=tuple(groups))))
        self.validate()

    def commit_structure(self, config):
        issues = self._validate_config(config,execution=False)
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

    def validate(self,*,check_files=False,execution=True):
        if check_files and not self.thermal_incomplete(self.config):
            from .readiness import readiness
            self.issues = readiness(self.config, self.base_directory)
        else:
            self.issues = self._validate_config(self.config,execution=execution)
        self.committed_issues = self.issues
        self.files_checked = check_files
        return self.issues

    def summary(self):
        if self.config.crystallography is not None:
            try: parse_text(', '.join(map(str,self.config.crystallography.slice_boundaries_fractional)),'boundaries')
            except ValueError:
                return dict(configuration='Incomplete — define slice boundaries in C / Numerical / Simulation cell',
                    compile='Not attempted — explicit slice boundaries required',bounds=None)
        if self.thermal_incomplete(self.config):
            return dict(configuration=self.thermal_message(self.config),
                compile='Not attempted — Atomic Model incomplete',bounds=None)

        from qcbed.crystallography_config import compile_config_from_sources
        try:compiled = compile_config_from_sources(self.config)
        except (OSError,ValueError,TypeError) as exc:return dict(configuration='Invalid',compile=str(exc),bounds=None)
        effective=compiled.configuration
        execution_errors = any(i.severity == 'error' and (i.phase == ValidationPhase.EXECUTION or i.code == 'output_exists') for i in self.committed_issues)
        file_errors = any(i.severity == 'error' and i.phase == ValidationPhase.FILES for i in self.committed_issues)
        return dict(regions=len(effective.sample.regions),stored_slices=len(compiled.state.model.slices),
            tilts=self.config.scan.scan_rows*self.config.scan.scan_columns,
            beams=len(compiled.experiment.beam_indices_hk),active_parameters=compiled.parameter_count,
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
