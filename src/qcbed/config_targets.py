"""Physical target catalogue backed by the compiler's own target enumeration."""
from dataclasses import dataclass, replace

from .configuration import ParameterTarget
from .config_compiler import CompiledConfig, compile_config, _target_map
from .refinement_categories import ParameterCategory, PARAMETER_SEMANTICS


@dataclass(frozen=True)
class RefinementTargetInfo:
    target: ParameterTarget
    label: str
    units: str | None
    starting_value: float
    selected: bool
    active: bool
    available: bool
    reason: str | None = None
    references: tuple = ()  # (identifier position, stable-ID object path)

    @property
    def category(self):
        return self.target.category

    @property
    def stable_ids(self):
        """Referenced object IDs only; excludes axis/component tokens."""
        return tuple(self.target.identifiers[i] for i,_ in self.references)


def _value_and_support(raw, category, t, coefficients, slices):
    s = raw.sample
    units = {'voltage':'kV','depth':'angstrom','debye_waller':'angstrom squared'}
    reason = None
    if category == 'voltage':
        value = raw.microscope.energy_kv
        if s.absorption.enabled: reason = 'Voltage refinement with absorption is unsupported.'
    elif category == 'dilation': value = s.dilation
    elif category == 'depth': value = s.regions[t[0]].depth_angstrom
    elif category in ('constant_vacancy','vacancy_coefficient'):
        value = s.regions[t[0]].vacancy.values[t[1]]
        reason = 'Nonzero vacancy/profile trials are unsupported by the current evaluator.'
    elif category == 'debye_waller':
        value = next(v.value_angstrom_squared for v in s.regions[t[0]].debye_waller if v.species == s.species[t[1]].id)
    elif category == 'absorption':
        coeff = coefficients[t[0]]
        value = (coeff.a,coeff.b,coeff.c)[t[1]]
        if t[1] == 2 and coeff.c_mode == 'linked_to_b': reason = 'C is linked to B; independent refinement is forbidden.'
    elif category == 'atom_coordinate':
        atom = slices[t[0]].atoms[t[1]]
        value = (atom.x_fractional,atom.y_fractional)[t[2]]
        if slices[t[0]].record_type == 'G': reason = 'G atom coordinates cannot have refinement masks.'
    elif category == 'scan_coordinate':
        value = (raw.scan.corner_a,raw.scan.corner_b,raw.scan.corner_c)[t[0]][t[1]]
    elif category == 'reflection_multiplier':
        reflection = s.regions[t[0]].reflections[t[1]]
        value = (reflection.elastic,reflection.absorptive)[t[2]]
    else:
        raise ValueError(f'No catalogue value semantics for {category}')
    return value,units.get(category,'dimensionless'),reason


def available_refinement_targets(config):
    """Enumerate compiler-valid targets, including unavailable targets with reasons.

    Raw starting values are not packed, folded or mapped. Selected denotes a
    binding; active additionally accounts for the category filter. Runtime
    readiness of the complete experiment is reported by validate_config.
    """
    from .config_editing import _target_references
    if isinstance(config,CompiledConfig):
        compiled = config
    elif config.crystallography is not None:
        from .crystallography_config import compile_config_from_sources
        compiled = compile_config_from_sources(config)
    else:
        compiled = compile_config(config)
    raw = compiled.configuration
    coefficients = tuple(next(c for c in raw.sample.absorption.coefficients if c.species == sp.id)
                         for sp in raw.sample.species) if raw.sample.absorption.enabled else ()
    slices = tuple(v for r in raw.sample.regions for v in r.slices)
    targets = _target_map(raw,coefficients)
    selected = {b.target for b in raw.refinement.parameters}
    active = {(p.category,p.target) for p in compiled.layout.parameters}
    result = []
    for target,(category,indices) in sorted(targets.items(),key=lambda item:
            (PARAMETER_SEMANTICS[item[1][0]].source_order_group,item[1][1])):
        value,units,reason = _value_and_support(raw,category,indices,coefficients,slices)
        label = 'Vacancy-coupling coefficient (legacy DILATION)' if category=='dilation' else category.replace('_',' ').capitalize()
        if target.identifiers: label += ': '+' / '.join(target.identifiers)
        result.append(RefinementTargetInfo(target,label,units,value,target in selected,
            (category,indices) in active,reason is None,reason,_target_references(raw,target)))
    if compiled.crystallography is not None and compiled.crystallography.site_groups:
        aliases = {internal:physical for physical,internal in compiled.crystallography.target_aliases}
        unsupported = {
            'atom_coordinate': 'Missing original-site coordinate binding across recasting/slicing; stored-slice indices are not physical site identities.',
            'reflection_multiplier': 'Strict crystallography bridge uses unity multipliers; original-hkl multiplier mapping is not implemented.',
            'constant_vacancy': 'Legacy vacancy is not chemical occupancy; nonzero vacancy propagation is unsupported.',
            'vacancy_coefficient': 'Nonzero vacancy/profile propagation is unsupported.',
            'dilation': 'Vacancy-coupling coefficient, not lattice dilation; nonzero vacancy propagation is unsupported.',
        }
        modern = []
        for info in result:
            target = aliases.get(info.target)
            if target is None:
                modern.append(replace(info, available=False, reason=unsupported.get(info.category.value, 'No modern physical target mapping.')))
                continue
            label = info.label
            refs = ()
            if info.category.value == 'depth': label = 'Sample → requested depth'
            elif info.category.value in ('debye_waller', 'absorption'):
                label = f'Site group {target.identifiers[1]} → ' + ('Debye–Waller B' if info.category.value == 'debye_waller' else 'absorption '+target.identifiers[2])
                refs = ((1, ('crystallography', 'site_groups', target.identifiers[1])),)
            elif info.category.value == 'voltage': label = 'Microscope → accelerating voltage'
            elif info.category.value == 'scan_coordinate': label = f'Scan {target.identifiers[0]} → {target.identifiers[1]}'
            modern.append(replace(info, target=target, label=label, references=refs))
        result = modern
    return tuple(result)


def target_label(config,info):
    """Presentation of catalogue identities only; no target enumeration here."""
    ids = info.target.identifiers
    category = info.category.value
    if info.label == 'Sample → requested depth': return info.label
    if category == 'depth': return f'Region {ids[0]} → depth'
    if category == 'debye_waller' and ids[0] == 'site_group': return info.label
    if category == 'debye_waller': return f'Region {ids[0]} → species {ids[1]} → Debye–Waller'
    if category == 'reflection_multiplier':
        region = next(r for r in config.sample.regions if r.id == ids[0])
        reflection = next(r for r in region.reflections if r.id == ids[1])
        return f'Region {ids[0]} → reflection ({reflection.h}, {reflection.k}) [{reflection.id}] → {ids[2]} multiplier'
    return info.label

