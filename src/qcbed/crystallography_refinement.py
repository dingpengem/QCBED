"""Physical group/scalar application through the existing dynamic parameter layout.

Supported source values are rebuilt through the original crystallographic bridge.
Bounds reject the original finite vector before legacy scalar mapping can fold it.
"""
from dataclasses import replace, asdict
import numpy as np
from .refinement import apply_refinement_parameters
from .config_compiler import compile_config
from .config_targets import available_refinement_targets


def _binding_keys(compiled):
    from .config_compiler import _target_map
    raw = compiled.configuration
    coefficients = tuple(next(c for c in raw.sample.absorption.coefficients if c.species == sp.id)
                         for sp in raw.sample.species) if raw.sample.absorption.enabled else ()
    targets = _target_map(raw, coefficients)
    return {targets[internal]: physical for physical, internal in compiled.crystallography.target_aliases}


def apply_crystallographic_parameters(compiled, vector):
    """Return a new compiled physical snapshot; baseline remains untouched."""
    crystal = compiled.crystallography
    if crystal is None or not crystal.site_groups:
        raise ValueError('explicit physical site groups required')
    values = np.asarray(vector)
    if values.dtype != np.float64 or values.shape != (compiled.parameter_count,) or not np.isfinite(values).all():
        raise ValueError('require finite float64 vector matching compiled layout')
    for info in available_refinement_targets(compiled):
        if info.active and not info.available: raise ValueError(info.reason)
    species = [s.id for s in compiled.configuration.sample.species]
    bindings = _binding_keys(compiled)
    bounds = {(b.category,b.target):b for b in compiled.bounds.entries}
    for parameter,value in zip(compiled.layout.parameters,values):
        key=(parameter.category,parameter.target)
        if parameter.category == 'projected_elastic':
            row = next(r for r in physical_target_layout(compiled) if r['internal_target']==list(parameter.target) and r['category']=='projected_elastic')
            if row['lower'] is None or row['upper'] is None or not row['lower'] <= value <= row['upper']:
                raise ValueError(row['target_id']+': value outside bounds')
            continue
        if key not in bindings:
            raise ValueError('Target has no implemented modern physical application')
        bound=bounds.get(key)
        if compiled.runner.bounds_policy != 'reject' or bound is None or bound.lower is None or bound.upper is None:
            raise ValueError('Active targets require explicit finite lower/upper rejection bounds')
        if value < bound.lower or value > bound.upper:
            raise ValueError(f'{bindings[key].identifiers}: value {value} outside bounds [{bound.lower},{bound.upper}]')
        if parameter.category != 'scan_coordinate' and value < 0:
            raise ValueError('Physical scalar must be nonnegative; never abs-fold a negative input')
    ordinary = np.array([v for p,v in zip(compiled.layout.parameters,values) if p.category!='projected_elastic'],dtype=np.float64)
    mapped = apply_refinement_parameters(compiled.state,ordinary,coordinate_mode='physical')
    if compiled.bounds.violations(mapped):
        raise ValueError('mapped physical values violate compiled bounds')
    groups={g.id:g for g in crystal.site_groups}
    for parameter in compiled.layout.parameters:
        if parameter.category == 'projected_elastic': continue
        physical=bindings[parameter.category,parameter.target]
        if parameter.category == 'debye_waller':
            group_id = physical.identifiers[1]
            groups[group_id]=replace(groups[group_id],debye_waller_b=mapped.regions[0].debye_waller[parameter.target[1]])
        elif parameter.category == 'absorption':
            group_id = physical.identifiers[1]
            groups[group_id]=replace(groups[group_id],absorption=mapped.absorption[parameter.target[0]])
    source=compiled.source_configuration
    categories = {p.category for p in compiled.layout.parameters}
    config=replace(source,
        microscope=replace(source.microscope,energy_kv=mapped.model.energy_kv) if 'voltage' in categories else source.microscope,
        scan=replace(source.scan,corner_a=mapped.scan.corner_a,corner_b=mapped.scan.corner_b,corner_c=mapped.scan.corner_c) if 'scan_coordinate' in categories else source.scan,
        crystallography=replace(source.crystallography,
        site_groups=tuple(groups[g.id] for g in crystal.site_groups),
        depth_angstrom=mapped.regions[0].depth_angstrom))
    from .projected_elastic import trial_records
    config=replace(config,crystallography=replace(config.crystallography,projected_elastic=trial_records(compiled,values)))
    return compile_config(config,structure=crystal.original_structure)


def physical_target_layout(compiled):
    """Human/GUI target metadata in existing numerical layout order."""
    infos={i.target:i for i in available_refinement_targets(compiled)}
    crystal=compiled.crystallography
    species=[s.id for s in compiled.configuration.sample.species]
    lookup=_binding_keys(compiled)
    bounds={(b.category,b.target):b for b in compiled.bounds.entries}
    result=[]
    for p in compiled.layout.parameters:
        if p.category == 'projected_elastic': continue
        target=lookup[p.category,p.target]; info=infos[target]; b=bounds[p.category,p.target]
        depth = p.category == 'depth'
        from .modern_yaml import target_name
        result.append(dict(target_id=target_name(target),
            scope={'depth':'sample','voltage':'microscope','scan_coordinate':'scan'}.get(p.category,'site-group'),site_group=target.identifiers[1] if target.identifiers and target.identifiers[0]=='site_group' else None,category=p.category,
            initial_value=info.starting_value,lower=b.lower,upper=b.upper,units=info.units,
            binding=asdict(target), internal_target=list(p.target), selected=True))
    from .projected_elastic import layout_rows
    return tuple(result) + layout_rows(compiled)


def site_occurrence_provenance(compiled):
    """Compact original ownership plus stored slices/sequence; copies reconstruct exactly."""
    sample=compiled.crystallography.simulation_sample
    return dict(groups=[asdict(g) for g in compiled.crystallography.site_groups],
        stored_slices=[[(s.original_site_index,s.physical_group_id) for s in slab.sites]
                       for slab in sample.depth_sample.slices],
        sequence=list(sample.parameters.slice_sequence_zero_based))
