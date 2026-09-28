"""Physical target catalogue backed by the compiler's own target enumeration."""
from dataclasses import dataclass

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
    compiled = config if isinstance(config,CompiledConfig) else compile_config(config)
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
        label = category.replace('_',' ').capitalize()
        if target.identifiers: label += ': '+' / '.join(target.identifiers)
        result.append(RefinementTargetInfo(target,label,units,value,target in selected,
            (category,indices) in active,reason is None,reason,_target_references(raw,target)))
    return tuple(result)


def target_label(config,info):
    """Presentation of catalogue identities only; no target enumeration here."""
    ids = info.target.identifiers
    category = info.category.value
    if category == 'depth': return f'Region {ids[0]} → depth'
    if category == 'debye_waller': return f'Region {ids[0]} → species {ids[1]} → Debye–Waller'
    if category == 'reflection_multiplier':
        region = next(r for r in config.sample.regions if r.id == ids[0])
        reflection = next(r for r in region.reflections if r.id == ids[1])
        return f'Region {ids[0]} → reflection ({reflection.h}, {reflection.k}) [{reflection.id}] → {ids[2]} multiplier'
    return info.label

