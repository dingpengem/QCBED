"""One-time legacy migration and optional export; never used by canonical runs."""
from dataclasses import asdict, dataclass, replace
from enum import Enum
from pathlib import Path

from . import configuration as c
from .config_compiler import CompiledConfig, compile_config
from .experimental import parse_experimental_intensities
from .experimental_config import load_experimental_data
from .legacy_io import _parse_for007_records
from .refinement_io import _CaptureRecords, _parse_refinement_input
from .refinement import build_refinement_parameter_layout


@dataclass(frozen=True)
class MigrationPolicy:
    """Required explicit choices for settings that legacy inputs do not define."""
    seed: int
    execution: c.ForwardExecutionConfig
    run: c.RunConfig
    objective_mode: str
    invalid_trial_policy: str
    bounds: c.BoundsConfig


class _MigrationRecords(_CaptureRecords):
    def __init__(self, path):
        super().__init__(path)
        self.kinds = []

    def read(self, section):
        read = super().read(section)
        if section.startswith('slice ') and section.endswith(' atom header'):
            token = read.token
            def capture():
                value = token()
                self.kinds.append(value.value)
                return value
            read.token = capture
        return read


def import_legacy_config(path, *, experimental_path, policy: MigrationPolicy) -> c.QcbedConfig:
    """Read legacy inputs once; explicit data path and Python policy are required.

    Deterministic ordinal IDs preserve records, including duplicate reflections.
    Supplied run settings select new canonical destinations, not the old BMNAME.
    """
    if not isinstance(policy, MigrationPolicy):
        raise TypeError('policy must be MigrationPolicy')
    path, experimental_path = Path(path), Path(experimental_path).resolve()
    state, controls, records, _ = _parse_refinement_input(path,experimental_header_path=experimental_path)
    if any(line.strip() for line in records.lines[records.position:]):
        raise ValueError('unrecognized trailing legacy records; resolve before migration')
    prefix = _MigrationRecords(path)
    _parse_for007_records(prefix)
    beams = prefix.values['output beam count and indices'][1:]
    beams = tuple(zip(beams[::2],beams[1::2]))
    data = parse_experimental_intensities(experimental_path,beam_indices_hk=beams,scan=state.scan)
    m = state.model
    species = tuple(c.SpeciesConfig(f'species-{i}',s.label,s.atomic_number,s.label_raw) for i,s in enumerate(m.species))
    regions, flat_slices = [], []
    offset = 0
    for i,r in enumerate(state.regions):
        slices = []
        for j in range(offset,offset+r.slice_count):
            stored = c.StoredSliceConfig(f'slice-{j}',m.slices[j],tuple(c.AtomConfig(f'atom-{k}',
                species[a.species_index_one_based-1].id,a.x_fractional,a.y_fractional,a.starting_occupancy)
                for k,a in enumerate(m.atoms_by_slice[j])),prefix.kinds[j],prefix.kinds[j])
            slices.append(stored)
            flat_slices.append(stored)
        regions.append(c.RegionConfig(f'region-{i}',r.depth_angstrom,tuple(slices),
            tuple(c.ThermalFactorConfig(sp.id,v) for sp,v in zip(species,r.debye_waller)),
            tuple(c.ReflectionConfig(f'reflection-{j}',v.h,v.k,v.elastic,v.absorptive)
                  for j,v in enumerate(state.reflections[i])),state.symmetry[i],
            c.VacancyConfig(c.VacancyKind(r.vacancy_flag),r.vacancy_values),
            c.SurfaceConfig(c.SurfaceKind(r.free_surface[0]),r.free_surface[1:])))
        offset += r.slice_count
    def identifiers(p):
        t = p.target
        if p.category in ('voltage','dilation'): return ()
        if p.category == 'depth': return (regions[t[0]].id,)
        if p.category == 'constant_vacancy': return (regions[t[0]].id,)
        if p.category == 'vacancy_coefficient': return (regions[t[0]].id,('a','b')[t[1]])
        if p.category == 'absorption': return (species[t[0]].id,('a','b','c')[t[1]])
        if p.category == 'debye_waller': return (regions[t[0]].id,species[t[1]].id)
        if p.category == 'atom_coordinate': return (flat_slices[t[0]].id,flat_slices[t[0]].atoms[t[1]].id,('x','y')[t[2]])
        if p.category == 'scan_coordinate': return (('A','B','C')[t[0]],('first','second')[t[1]])
        if p.category == 'reflection_multiplier': return (regions[t[0]].id,regions[t[0]].reflections[t[1]].id,('elastic','absorptive')[t[2]])
        raise ValueError(f'unsupported migration target {p.category}')
    layout = build_refinement_parameter_layout(state)
    plamda = controls.pop('plamda')
    raw = c.QcbedConfig(m.title,c.MicroscopeConfig(m.energy_kv,m.zone_axis),
        c.SampleConfig(m.mesh,species,tuple(regions),state.dsap,state.dilation,
            c.AbsorptionConfig(m.absorption_enabled,tuple(c.AbsorptionSpeciesConfig(species[i].id,*row,
                'linked_to_b' if state.absorption_masks[i][2] == -1 else 'independent')
                for i,row in enumerate(state.absorption)),m.absorption_raw)),
        c.ScanConfig(**asdict(state.scan)),
        c.ExperimentalDataConfig(str(experimental_path),beams,data.gradient_option,header='validate'),
        c.RefinementConfig(c.OptimizerConfig(**controls,seed=policy.seed),tuple(
            c.ParameterBinding(c.ParameterTarget(c.ParameterCategory(p.category),identifiers(p)),state.ranges[i],plamda[i])
            for i,p in enumerate(layout.parameters)),state.iref,policy.objective_mode,policy.invalid_trial_policy,
            bounds=policy.bounds),policy.execution,policy.run,
        c.LegacyExportConfig(prefix.values['DSMAX'][0],title_raw=m.title_raw,
            output_filename_raw=policy.run.output_filename+m.output_filename_raw[len(m.output_filename_raw.rstrip()):],
            experimental_filename_raw=str(experimental_path)+m.experimental_intensity_filename_raw[len(m.experimental_intensity_filename_raw.rstrip()):]))
    compile_config(raw)
    return raw


def _new_path(path):
    path = Path(path).resolve()
    root = Path(__file__).resolve().parents[2]
    if any(path.is_relative_to(root/name) for name in ('legacy_fortran','reference_data')):
        raise ValueError('cannot write protected legacy/reference files')
    if path.exists():
        raise FileExistsError(path)
    return path


def _plain(value):
    if isinstance(value,Enum): return value.value
    if isinstance(value,dict): return {k:_plain(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)): return [_plain(v) for v in value]
    return value


def write_yaml_config(config, path):
    """Exclusive YAML write; absolute I/O paths preserve migration destinations."""
    import yaml
    compile_config(config)
    destination = _new_path(path)
    def raw_path(raw, value):
        return None if raw is None else str(Path(value).resolve())+raw[len(raw.rstrip()):]
    config = replace(config,experiment=replace(config.experiment,path=str(Path(config.experiment.path).resolve())),
        run=replace(config.run,output_directory=str(Path(config.run.output_directory).resolve()),
                    output_filename=str(Path(config.run.output_filename).resolve())),
        legacy_export=replace(config.legacy_export,
            filename=None if config.legacy_export.filename is None else str(Path(config.legacy_export.filename).resolve()),
            experimental_filename_raw=raw_path(config.legacy_export.experimental_filename_raw,config.experiment.path),
            output_filename_raw=raw_path(config.legacy_export.output_filename_raw,config.run.output_filename)))
    with destination.open('x',encoding='utf-8') as stream:
        stream.write(yaml.safe_dump(_plain(asdict(config)),sort_keys=False))


def legacy_export_report(compiled):
    """Explicit values omitted from FOR007, returned even for default policies."""
    raw = compiled.configuration
    return dict(execution=asdict(raw.execution),seed=raw.refinement.optimizer.seed,
        bounds=asdict(raw.refinement.bounds),objective_mode=raw.refinement.objective_mode,
        invalid_trial_policy=raw.refinement.invalid_trial_policy,run=asdict(raw.run),
        stable_ids='not representable; reimport assigns ordinal IDs',
        active_categories='compiled masks exported; inactive bindings/category selectors not retained',
        experimental_header='NUMP/NUML/GRADIENTOPT require a separate header-bearing data file',
        executable_limits='legacy binary capacities/string widths are build-dependent; not inferred from MULTIS.FTN')


def _for007_text(compiled, *, experimental_filename, output_filename):
    s,m = compiled.state,compiled.state.model
    metadata = compiled.configuration.legacy_export
    if metadata.dsmax_inverse_angstrom is None:
        raise ValueError('DSMAX is required for export; specify legacy_export.dsmax_inverse_angstrom')
    if any(r.depth_angstrom <= 0 for r in s.regions):
        raise ValueError('legacy input requires positive region depths; no folding on export')
    def line(value):
        if not value.isascii():
            raise ValueError('legacy text export requires ASCII; no implicit encoding conversion')
        if '\n' in value or '\r' in value:
            raise ValueError('legacy text fields cannot contain newlines')
        return value
    if m.energy_kv <= 0: raise ValueError('legacy voltage must be positive')
    if any(v != int(v) for v in (*s.scan.basis_1,*s.scan.basis_2)):
        raise ValueError('legacy scan basis requires integer-valued components; no rounding on export')
    if len(m.title) > 80: raise ValueError('legacy title exceeds 80 characters')
    lines = [line(m.title_raw)]
    def emit(*values):
        lines.append(' '.join(str(v) if isinstance(v,int) else repr(float(v)) for v in values))
    emit(m.energy_kv,int(m.energy_refined),*m.zone_axis)
    emit(len(m.slices),m.mesh.mesh_x,m.mesh.mesh_y)
    for g in m.slices:
        emit(g.cell_a_angstrom,g.cell_b_angstrom,g.lattice_gamma_degrees,g.thickness_angstrom,g.x_shift_fractional,g.y_shift_fractional)
    emit(len(m.species))
    for sp in m.species:
        if len(sp.label) != 2: raise ValueError('legacy species labels require exactly two characters')
        raw = sp.label_raw or sp.label
        if raw[:2].ljust(2) != sp.label: raise ValueError('species raw label disagrees with A2 label')
        lines.append(line(raw)); emit(sp.atomic_number)
    emit(s.dsap); emit(metadata.dsmax_inverse_angstrom); emit(len(s.regions))
    emit(*(v for r in s.regions for v in (r.depth_angstrom,r.thickness_mask)))
    emit(*(r.slice_count for r in s.regions))
    for r in s.regions:
        emit(r.vacancy_flag); emit(*(v for pair in zip(r.vacancy_values,r.vacancy_masks) for v in pair))
    emit(s.dilation,s.dilation_mask)
    lines.append(line(m.absorption_raw))
    if m.absorption_enabled:
        for values,masks in zip(s.absorption,s.absorption_masks):
            emit(*(v for pair in zip(values,masks) for v in pair))
    for r in s.regions:
        emit(*r.debye_waller); emit(*r.debye_waller_masks)
        emit(int(r.free_surface[0]),*r.free_surface[1:])
    stored = [v for r in compiled.configuration.sample.regions for v in r.slices]
    for definition,atoms in zip(stored,m.atoms_by_slice):
        lines.append(f'{len(atoms)} {definition.record_type}')
        for a in atoms:
            if definition.record_type == 'S':
                emit(a.species_index_one_based,a.x_fractional,a.x_mask,a.y_fractional,a.y_mask,a.starting_occupancy)
            else: emit(a.species_index_one_based,a.x_fractional,a.y_fractional,a.starting_occupancy)
    lines.append(line(output_filename))
    beams = compiled.experiment.beam_indices_hk
    emit(len(beams),*(v for beam in beams for v in beam))
    lines.append(line(experimental_filename))
    emit(*(int(v) for v in (*s.scan.basis_1,*s.scan.basis_2)))
    for corner in (s.scan.corner_a,s.scan.corner_b,s.scan.corner_c): emit(*corner)
    lines.extend(s.scan_masks); emit(s.iref)
    for sym,reflections in zip(s.symmetry,s.reflections):
        emit(len(reflections)); lines.append(sym)
        for r in reflections: emit(r.h,r.k,r.elastic,r.absorptive)
        lines.extend(r.mask for r in reflections)
    o = compiled.optimizer
    if len(o.method) != 4 or any(v.isspace() for v in o.method):
        raise ValueError('legacy METHOD export requires an exact four-character selector')
    lines.append(f'{o.method} {o.paraopt} {o.loud}')
    for name in ('loops','maxfsimp','maxfpara','ftol','chimin','maxtemp','mintemp','coolrate'): emit(getattr(o,name))
    emit(*o.plamda)
    for pair in s.ranges: emit(*pair)
    return '\n'.join(lines)+'\n'


def export_for007(config, path, *, bundle=False, experimental_filename=None, output_filename=None, filename_limits=None):
    """Exclusive optional export. Bundle creates a new directory with its own data.

    Returns a report of Python-only policies. No existing file is overwritten.
    Standalone export requires configured header-validation mode. Explicit name
    overrides are relocation choices, not modifications to canonical objects.
    """
    compiled = config if isinstance(config,CompiledConfig) else compile_config(config)
    destination = _new_path(path)
    if bundle and (experimental_filename is not None or output_filename is not None):
        raise ValueError('bundle uses fixed separate names in.dat and out.dat')
    if not bundle and compiled.experiment.header != 'validate':
        raise ValueError('headerless experimental data requires bundle=True')
    exp = 'in.dat' if bundle else (experimental_filename or compiled.experiment.path)
    out = 'out.dat' if bundle else (output_filename or compiled.configuration.run.output_filename)
    if not bundle:
        if filename_limits is None or len(filename_limits) != 2 or any(type(v) is not int or v < 1 for v in filename_limits):
            raise ValueError('standalone export requires verified filename_limits=(EXPINT_width, BMNAME_width); use bundle for reference names')
        if len(exp) > filename_limits[0] or len(out) > filename_limits[1]:
            raise ValueError('legacy filenames exceed explicitly supplied target widths; use explicit relocation names or bundle')
    text = _for007_text(compiled,experimental_filename=exp,output_filename=out)
    report = legacy_export_report(compiled)
    if bundle:
        data = load_experimental_data(compiled,base_directory=Path.cwd())
        body = [f'{data.scan_columns} {data.scan_rows} {data.gradient_option}']
        for block in (data.intensities,data.variances,data.weights):
            for beam in range(block.shape[1]):
                for row in range(data.scan_rows):
                    body.append(' '.join(repr(float(v)) for v in block[row*data.scan_columns:(row+1)*data.scan_columns,beam]))
        destination.mkdir()
        (destination/'in.dat').write_text('\n'.join(body)+'\n',encoding='utf-8')
        (destination/'FOR007.dat').write_text(text,encoding='utf-8')
        import json
        (destination/'compatibility.json').write_text(json.dumps(_plain(report),indent=2)+'\n',encoding='utf-8')
    else:
        with destination.open('x',encoding='utf-8') as stream: stream.write(text)
    return report
