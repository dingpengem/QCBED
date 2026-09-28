"""Pure canonical-to-existing-state adapter. No file reads or physics execution."""
from dataclasses import asdict, dataclass, replace

from .configuration import QcbedConfig, ParameterTarget, ParameterCategory, ExperimentalDataConfig
from .legacy_io import LegacyInput, Species, AtomSite
from .scan import LegacyScanDefinition
from .refinement import (RefinementState, RefinementRegion, RefinementReflection,
                         RefinementParameterLayout, build_refinement_parameter_layout,
                         pack_refinement_parameters)
from .refinement_driver import LegacyRefinementConfig
from .refinement_subset import select_refinement_categories
from .refinement_bounds import (RefinementBound, RefinementBoundRules,
                                RefinementBounds, resolve_refinement_bounds)
from .full_runner import FullRefinementConfig


@dataclass(frozen=True)
class CompiledConfig:
    configuration: QcbedConfig
    state: RefinementState
    optimizer: LegacyRefinementConfig
    layout: RefinementParameterLayout
    bounds: RefinementBounds
    bound_rules: RefinementBoundRules
    runner: FullRefinementConfig
    experiment: ExperimentalDataConfig

    @property
    def parameter_count(self):
        """Active instance count, never the number of categories or records."""
        return len(self.layout.parameters)

    @property
    def packed_starting_vector(self):
        """Existing ARRANGEPARS projection; never written into raw configuration."""
        return pack_refinement_parameters(self.state)

    @property
    def execution(self):
        return self.optimizer.execution

    @property
    def range_factors(self):
        return self.state.ranges

    @property
    def plamda(self):
        return self.optimizer.plamda


def _index(items, field='id'):
    result = {}
    for i, item in enumerate(items):
        key = getattr(item, field)
        if not key or key in result:
            raise ValueError(f'empty or duplicate {field}: {key!r}')
        result[key] = i
    return result


def _target_map(config, coefficients):
    """Single source of stable-target enumeration for compilation and frontends."""
    s = config.sample
    species_ids = _index(s.species)
    slices = tuple(stored for region in s.regions for stored in region.slices)
    slice_ids = _index(slices)
    targets = {}
    def target(category, ids, indices):
        targets[ParameterTarget(ParameterCategory(category), ids)] = (category, indices)
    target('voltage', (), ())
    target('dilation', (), ())
    for r, region in enumerate(s.regions):
        if not region.slices:
            raise ValueError('each region requires at least one stored slice')
        target('depth', (region.id,), (r,))
        if region.vacancy.kind == 0:
            target('constant_vacancy', (region.id,), (r, 0))
        else:
            for j, name in enumerate(('a', 'b')):
                target('vacancy_coefficient', (region.id, name), (r, j))
        thermal_ids = _index(region.debye_waller, 'species')
        if set(thermal_ids) != set(species_ids):
            raise ValueError('Debye-Waller requires exactly one entry per species')
        for sp, j in species_ids.items():
            target('debye_waller', (region.id, sp), (r, j))
        _index(region.reflections)
        for j, reflection in enumerate(region.reflections):
            for k, name in enumerate(('elastic', 'absorptive')):
                target('reflection_multiplier', (region.id, reflection.id, name), (r, j, k))
    for i, coeff in enumerate(coefficients):
        for j, name in enumerate(('a', 'b', 'c')):
            target('absorption', (coeff.species, name), (i, j))
    for stored in slices:
        for j, atom in enumerate(stored.atoms):
            for k, name in enumerate(('x', 'y')):
                target('atom_coordinate', (stored.id, atom.id, name), (slice_ids[stored.id], j, k))
    for i, corner in enumerate(('A', 'B', 'C')):
        for j, name in enumerate(('first', 'second')):
            target('scan_coordinate', (corner, name), (i, j))

    return targets


def compile_config(config: QcbedConfig) -> CompiledConfig:
    """Resolve stable references/masks in source order, retaining input values.

    Bounds use the existing mapped-start resolver on a separate state. The raw
    baseline is never mapped. Unsupported numerical states remain representable
    and continue to be rejected by the existing trial evaluator.
    """
    if not isinstance(config, QcbedConfig):
        raise TypeError('config must be QcbedConfig')
    if config.schema_version != 1:
        raise ValueError('unsupported configuration schema_version')
    s, ref = config.sample, config.refinement
    if not s.species or not s.regions:
        raise ValueError('species and regions must not be empty')
    species_ids = _index(s.species)
    _index(s.regions)
    slices = tuple(stored for region in s.regions for stored in region.slices)
    slice_ids = _index(slices)
    for sp in s.species:
        if sp.atomic_number <= 0:
            raise ValueError('atomic_number must be positive')
    for stored in slices:
        _index(stored.atoms)
        if any(getattr(stored.geometry, field) is None for field in
               ('cell_a_angstrom', 'cell_b_angstrom', 'lattice_gamma_degrees')):
            raise ValueError('canonical slice geometry requires an explicit lateral cell')
        for atom in stored.atoms:
            if atom.species not in species_ids:
                raise ValueError(f'unknown atom species: {atom.species}')
    absorption = s.absorption
    if absorption.enabled:
        coeff_index = _index(absorption.coefficients, 'species')
        if set(coeff_index) != set(species_ids):
            raise ValueError('absorption requires exactly one entry per species')
        coefficients = tuple(absorption.coefficients[coeff_index[sp.id]] for sp in s.species)
    else:
        if absorption.coefficients:
            raise ValueError('disabled absorption requires empty coefficients')
        coefficients = ()
    raw_absorption = absorption.legacy_raw
    if raw_absorption is None:
        raw_absorption = 'YES' if absorption.enabled else 'NO'
    if (raw_absorption[:2] == 'YE') != absorption.enabled:
        raise ValueError('absorption legacy_raw disagrees with enabled')

    targets = _target_map(config, coefficients)

    def resolve(value):
        if value not in targets:
            raise ValueError(f'unknown parameter target: {value}')
        return targets[value]
    bindings = {}
    for binding in ref.parameters:
        key = resolve(binding.target)
        if key in bindings:
            raise ValueError(f'duplicate parameter binding: {binding.target}')
        category, indices = key
        if category == 'absorption' and indices[1] == 2 and coefficients[indices[0]].c_mode == 'linked_to_b':
            raise ValueError('linked absorption C cannot be independently refined')
        if category == 'atom_coordinate' and slices[indices[0]].record_type == 'G':
            raise ValueError('G records cannot have atom refinement masks')
        bindings[key] = binding
    def mask(category, *indices):
        return int((category, indices) in bindings)
    def pair(category, *indices):
        return ''.join('Y' if mask(category, *indices, k) else 'N' for k in range(2))

    legacy = config.legacy_export
    for raw, value, name in ((legacy.title_raw, config.title, 'title'),
                            (legacy.output_filename_raw, config.run.output_filename, 'output filename'),
                            (legacy.experimental_filename_raw, config.experiment.path, 'experimental filename')):
        if raw is not None and raw.rstrip() != value:
            raise ValueError(f'legacy raw {name} disagrees with configured value')
    if not config.experiment.path.strip() or not config.run.output_filename.strip() or not config.run.output_directory.strip():
        raise ValueError('input/output paths must not be empty')
    model = LegacyInput(config.title, config.microscope.energy_kv, bool(mask('voltage')),
        config.microscope.zone_axis, s.mesh, tuple(st.geometry for st in slices),
        tuple(Species(sp.label, sp.atomic_number, sp.label_raw) for sp in s.species),
        tuple(tuple(AtomSite(species_ids[a.species]+1, a.x_fractional, a.y_fractional,
                            a.starting_occupancy, st.record_type,
                            mask('atom_coordinate', i, j, 0), mask('atom_coordinate', i, j, 1),
                            st.record_type_raw or st.record_type)
                    for j, a in enumerate(st.atoms)) for i, st in enumerate(slices)),
        config.run.output_filename, config.experiment.path,
        config.title if legacy.title_raw is None else legacy.title_raw,
        raw_absorption,
        config.run.output_filename if legacy.output_filename_raw is None else legacy.output_filename_raw,
        config.experiment.path if legacy.experimental_filename_raw is None else legacy.experimental_filename_raw)
    regions = []
    for i, region in enumerate(s.regions):
        thermal = {v.species: v.value_angstrom_squared for v in region.debye_waller}
        vacancy_category = 'vacancy_coefficient' if region.vacancy.kind else 'constant_vacancy'
        regions.append(RefinementRegion(region.depth_angstrom, mask('depth', i), len(region.slices),
            int(region.vacancy.kind), region.vacancy.values,
            tuple(mask(vacancy_category, i, j) for j in range(len(region.vacancy.values))),
            tuple(thermal[sp.id] for sp in s.species),
            tuple(mask('debye_waller', i, j) for j in range(len(s.species))),
            (int(region.free_surface.kind), *region.free_surface.coefficients)))
    scan = LegacyScanDefinition(**asdict(config.scan))
    state = RefinementState(model, tuple(regions), s.dilation, mask('dilation'),
        tuple((v.a, v.b, v.c) for v in coefficients),
        tuple((mask('absorption', i, 0), mask('absorption', i, 1),
               -1 if v.c_mode == 'linked_to_b' else mask('absorption', i, 2)) for i, v in enumerate(coefficients)),
        scan, tuple(pair('scan_coordinate', i) for i in range(3)),
        tuple(tuple(RefinementReflection(v.h, v.k, v.elastic, v.absorptive,
                    pair('reflection_multiplier', i, j)) for j, v in enumerate(r.reflections))
              for i, r in enumerate(s.regions)), tuple(r.symmetry for r in s.regions),
        s.dsap_inverse_angstrom, ref.iref)
    order = build_refinement_parameter_layout(state).parameters
    if len(order) != len(bindings):
        raise ValueError('parameter bindings do not match source-defined masks')
    state = replace(state, ranges=tuple(bindings[p.category, p.target].range_factors for p in order))
    optimizer = LegacyRefinementConfig(**asdict(ref.optimizer),
        plamda=tuple(bindings[p.category, p.target].plamda for p in order), execution=config.execution)
    categories = None if ref.active_categories is None else tuple(c.value for c in ref.active_categories)
    if categories is not None:
        if len(set(categories)) != len(categories):
            raise ValueError('duplicate active category')
        state, optimizer = select_refinement_categories(state, optimizer, categories)
    rules = RefinementBoundRules(
        tuple(RefinementBound(*resolve(b.target), b.lower, b.upper) for b in ref.bounds.specific),
        tuple((b.category.value, (b.lower, b.upper)) for b in ref.bounds.absolute),
        tuple((b.category.value, b.fraction) for b in ref.bounds.relative),
        tuple((b.category.value, b.absolute_fallback) for b in ref.bounds.relative))
    bounds = resolve_refinement_bounds(state, rules)
    runner = FullRefinementConfig(optimizer, objective_mode=ref.objective_mode,
        invalid_trial_policy=ref.invalid_trial_policy, refinement_categories=categories,
        snapshot_every=config.run.snapshot_every, checkpoint_every=config.run.checkpoint_every,
        max_total_evaluations=config.run.max_total_evaluations,
        max_wall_time=config.run.max_wall_time_seconds, bounds_policy=ref.bounds.policy,
        refinement_bounds=bounds, bound_rules=rules)
    return CompiledConfig(config, state, optimizer, build_refinement_parameter_layout(state),
                          bounds, rules, runner, config.experiment)
