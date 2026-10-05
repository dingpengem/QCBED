"""Physical site ownership and alias resolution into the existing target layout."""
from dataclasses import replace
from ..configuration import (AtomicSiteGroup, ParameterTarget, ParameterCategory,
    ThermalFactorConfig, AbsorptionConfig, AbsorptionSpeciesConfig)


def group_target(group_id):
    return ParameterTarget(ParameterCategory.DEBYE_WALLER, ('site_group', group_id))


def requested_depth_target():
    return ParameterTarget(ParameterCategory.DEPTH, ('sample','requested_depth'))


def resolve_site_groups(structure, groups):
    """Apply explicit QCBED metadata to original sites; retain source order.

    The present scattering bridge supports one group per element. Independent
    same-element groups need a future species/scattering mapping, not silent merging.
    """
    groups = tuple(groups)
    if any(not isinstance(g, AtomicSiteGroup) for g in groups):
        raise TypeError('site_groups requires AtomicSiteGroup records')
    if len({g.id for g in groups}) != len(groups):
        raise ValueError('duplicate physical site-group ID')
    owners = {}
    for group in groups:
        for index in group.site_indices:
            if index >= len(structure.sites):
                raise ValueError(f'{group.id}: original site index {index} does not exist')
            if index in owners:
                raise ValueError(f'original site {index} belongs to conflicting groups')
            if structure.sites[index].symbol != group.element:
                raise ValueError(f'{group.id}: element differs from member site {index}')
            owners[index] = group
    if set(owners) != set(range(len(structure.sites))):
        raise ValueError('explicit site groups must cover every original site exactly once')
    if len({g.element for g in groups}) != len(groups):
        raise ValueError('independent same-element groups are deferred: backend requires one shared B per element')
    if any(g.debye_waller_b is None for g in groups):
        raise ValueError('every site group requires explicit B for scattering')
    return replace(structure, sites=tuple(replace(site,
        occupancy=owners[i].occupancy, debye_waller_b=owners[i].debye_waller_b,
        physical_group_id=owners[i].id, original_site_index=i)
        for i,site in enumerate(structure.sites)))


def group_compiler_inputs(config, structure):
    """Return effective species metadata and canonical→internal target aliases."""
    c = config.crystallography
    applied = resolve_site_groups(structure, c.site_groups)
    symbols = tuple(dict.fromkeys(s.symbol for s in structure.sites))
    by_element = {g.element:g for g in c.site_groups}
    if c.debye_waller or c.absorption.coefficients or c.absorption.enabled:
        raise ValueError('site_groups own B/absorption; leave old species metadata empty to avoid conflicting authorities')
    supplied = [by_element[s].absorption is not None for s in symbols]
    if any(supplied) and not all(supplied):
        raise ValueError('absorption metadata must cover all site groups or be absent for all')
    absorption = AbsorptionConfig(bool(all(supplied)), tuple(
        AbsorptionSpeciesConfig(s,*by_element[s].absorption) for s in symbols) if all(supplied) else ())
    thermal = tuple(ThermalFactorConfig(s,by_element[s].debye_waller_b) for s in symbols)
    aliases = tuple((group_target(by_element[s].id),
                    ParameterTarget(ParameterCategory.DEBYE_WALLER,('crystal',s))) for s in symbols)
    aliases += ((requested_depth_target(),ParameterTarget(ParameterCategory.DEPTH,('crystal',))),)
    aliases += ((ParameterTarget(ParameterCategory.VOLTAGE, ()), ParameterTarget(ParameterCategory.VOLTAGE, ())),)
    aliases += tuple((ParameterTarget(ParameterCategory.SCAN_COORDINATE, (corner, axis)),
                      ParameterTarget(ParameterCategory.SCAN_COORDINATE, (corner, axis)))
                     for corner in ('A', 'B', 'C') for axis in ('first', 'second'))
    if absorption.enabled:
        aliases += tuple((ParameterTarget(ParameterCategory.ABSORPTION, ('site_group', by_element[s].id, component)),
                          ParameterTarget(ParameterCategory.ABSORPTION, (s, component)))
                         for s in symbols for component in ('a', 'b', 'c'))
    return applied, thermal, absorption, aliases


def remap_refinement(refinement, aliases):
    """Use existing masks/layout/ranges/bounds; only resolve physical stable IDs."""
    mapping = dict(aliases)
    internal = set(mapping.values())
    if any(b.target not in mapping for b in refinement.parameters):
        raise ValueError("Target has no implemented modern physical mapping")
    def target(t):
        if t in mapping: return mapping[t]
        if t in internal:
            raise ValueError('grouped Debye-Waller must bind by physical site-group identity')
        if t.category == ParameterCategory.DEBYE_WALLER and t not in mapping:
            raise ValueError(f'unknown physical site-group target: {t.identifiers}')
        return mapping.get(t,t)
    return replace(refinement,
        parameters=tuple(replace(b,target=target(b.target)) for b in refinement.parameters),
        bounds=replace(refinement.bounds,specific=tuple(replace(b,target=target(b.target))
            for b in refinement.bounds.specific)))
