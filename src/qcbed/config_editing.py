"""Pure stable-ID editing of immutable canonical snapshots, without GUI state.

Edits preserve field types; cross-object completeness is checked separately by
validate_config/compile_config. This permits explicit multi-step species edits.
No operation silently changes physical values or establishes resume compatibility.
"""
from dataclasses import dataclass, fields, is_dataclass, replace

from .configuration import ParameterBinding, ParameterTarget


@dataclass(frozen=True)
class Id:
    value: str
    field: str = 'id'


@dataclass(frozen=True)
class Dependency:
    object_path: tuple
    field: str
    message: str


class DependencyError(ValueError):
    def __init__(self, dependencies):
        self.dependencies = tuple(dependencies)
        super().__init__('Deletion requires explicit resolution of dependent references')


def get_object(config, path):
    value = config
    for step in path:
        if isinstance(step,Id):
            matches = [v for v in value if getattr(v,step.field,None) == step.value]
            if len(matches) != 1: raise ValueError(f'Expected one object for {step}')
            value = matches[0]
        elif isinstance(step,str): value = getattr(value,step)
        else: raise TypeError('paths contain field names and Id selectors, never row indices')
    return value


def _replace_at(value,path,replacement):
    if not path: return replacement
    step,*tail = path
    if isinstance(step,Id):
        old = get_object(value,(step,))
        return tuple(_replace_at(v,tail,replacement) if v is old else v for v in value)
    if not isinstance(step,str): raise TypeError('paths contain field names and Id selectors')
    return replace(value,**{step:_replace_at(getattr(value,step),tail,replacement)})


def replace_object(config,path,value):
    """Replace a value/object; object IDs cannot be renamed by this operation."""
    path = tuple(path)
    old = get_object(config,path)
    if hasattr(old,'id') and (not hasattr(value,'id') or old.id != value.id):
        raise ValueError('replace preserves stable IDs; use explicit reference-aware migration to rename')
    if type(old) is not type(value): raise TypeError('replacement must have the same type')
    return _replace_at(config,path,value)


def update_object(config,path,**changes):
    return replace_object(config,path,replace(get_object(config,path),**changes))


def add_object(config,collection_path,value,*,before_id=None):
    """Insert in scientific order. Default appends; before_id is explicit."""
    collection_path = tuple(collection_path)
    items = get_object(config,collection_path)
    if not isinstance(items,tuple) or not hasattr(value,'id'):
        raise TypeError('add_object requires an ID-bearing ordered collection')
    if any(item.id == value.id for item in items): raise ValueError(f'duplicate ID {value.id}')
    from .configuration import RegionConfig, StoredSliceConfig
    new_slices = value.slices if isinstance(value,RegionConfig) else ((value,) if isinstance(value,StoredSliceConfig) else ())
    old_ids = {s.id for r in config.sample.regions for s in r.slices}
    new_ids = [s.id for s in new_slices]
    if len(new_ids) != len(set(new_ids)) or old_ids.intersection(new_ids):
        raise ValueError('stored slice IDs must be globally unique')
    if before_id is None: position = len(items)
    else:
        positions = [i for i,item in enumerate(items) if item.id == before_id]
        if len(positions) != 1: raise ValueError('before_id must identify one existing object')
        position = positions[0]
    return _replace_at(config,collection_path,items[:position]+(value,)+items[position:])


def reorder_objects(config,collection_path,ordered_ids):
    """Explicit scientific reordering, never a presentation-only sort."""
    items = get_object(config,collection_path)
    mapping = {v.id:v for v in items}
    ordered_ids = tuple(ordered_ids)
    if len(mapping) != len(items) or len(set(ordered_ids)) != len(items) or set(ordered_ids) != set(mapping):
        raise ValueError('ordered_ids must be an exact permutation of unique existing IDs')
    return _replace_at(config,tuple(collection_path),tuple(mapping[i] for i in ordered_ids))


def _target_references(config,target):
    """Identifier positions and scoped object addresses, not vector indices."""
    ids,c = target.identifiers,target.category.value
    if config.crystallography is not None:
        if c in ('debye_waller','absorption') and ids and ids[0]=='site_group':
            return ((1,('crystallography','site_groups',Id(ids[1]))),)
        if c=='depth' and ids==('sample','requested_depth'): return ()
    from .refinement_categories import PARAMETER_SEMANTICS
    arity = 1 if c == 'constant_vacancy' else PARAMETER_SEMANTICS[c].target_arity
    if len(ids) != arity:
        raise ValueError('Resolve malformed refinement targets before dependency-aware editing')
    region = lambda i: ('sample','regions',Id(i))
    species = lambda i: ('sample','species',Id(i))
    if c in ('depth','constant_vacancy','vacancy_coefficient'): return ((0,region(ids[0])),)
    if c == 'debye_waller': return ((0,region(ids[0])),(1,species(ids[1])))
    if c == 'absorption': return ((0,species(ids[0])),)
    if c == 'reflection_multiplier': return ((0,region(ids[0])),(1,region(ids[0])+('reflections',Id(ids[1]))))
    if c == 'atom_coordinate':
        matches = [(r,s) for r in config.sample.regions for s in r.slices if s.id == ids[0]]
        if len(matches) != 1: return ()
        r,s = matches[0]
        path = region(r.id)+('slices',Id(s.id))
        return ((0,path),(1,path+('atoms',Id(ids[1]))))
    return ()


def deletion_dependencies(config,path):
    """Report references outside the subtree being removed, including inactive bindings."""
    path = tuple(path)
    get_object(config,path)
    result = []
    def referenced(reference): return reference[:len(path)] == path
    for name,items in (('parameters',config.refinement.parameters),('specific',config.refinement.bounds.specific)):
        for binding in items:
            if any(referenced(p) for _,p in _target_references(config,binding.target)):
                base = ('refinement','parameters') if name == 'parameters' else ('refinement','bounds','specific')
                result.append(Dependency(base+(binding.target,), 'target','Refinement binding or bound refers to deleted object'))
    if len(path) == 3 and path[:2] == ('sample','species') and isinstance(path[-1],Id):
        species = path[-1].value
        for r in config.sample.regions:
            base = ('sample','regions',Id(r.id))
            for s in r.slices:
                for a in s.atoms:
                    if a.species == species:
                        result.append(Dependency(base+('slices',Id(s.id),'atoms',Id(a.id)),'species','Atom uses this species'))
            for v in r.debye_waller:
                if v.species == species: result.append(Dependency(base+('debye_waller',Id(species,'species')),'species','Thermal entry uses this species'))
        for v in config.sample.absorption.coefficients:
            if v.species == species:
                result.append(Dependency(('sample','absorption','coefficients',Id(species,'species')),'species','Absorption entry uses this species'))
    return tuple(result)


def delete_object(config,path):
    """Refuse dependent deletion; callers must explicitly remove/reassign references."""
    path = tuple(path)
    if not path or not isinstance(path[-1],Id): raise ValueError('delete requires an ID-selected object')
    dependencies = deletion_dependencies(config,path)
    if dependencies: raise DependencyError(dependencies)
    items = get_object(config,path[:-1])
    old = get_object(config,path)
    return _replace_at(config,path[:-1],tuple(v for v in items if v is not old))


def duplicate_object(config,path,*,copy_refinement=False):
    """Duplicate a subtree with deterministic fresh IDs; originals never change.

    By default no bindings/bounds are copied. copy_refinement explicitly remaps
    associated bindings and specific bounds. Shared species references remain
    shared when copying regions/slices/atoms. Returns a new snapshot.
    """
    path = tuple(path)
    source = get_object(config,path)
    if not path or not isinstance(path[-1],Id) or not hasattr(source,'id'):
        raise ValueError('duplicate requires an ID-bearing object')
    used = set()
    def collect(value):
        if is_dataclass(value):
            if hasattr(value,'id'): used.add(value.id)
            for f in fields(value): collect(getattr(value,f.name))
        elif isinstance(value,tuple):
            for v in value: collect(v)
    collect(config)
    renamed = {}
    def clone(value,address):
        if not is_dataclass(value): return value
        changes = {}
        if hasattr(value,'id'):
            n = 1
            while f'{value.id}-copy-{n}' in used: n += 1
            new_id = f'{value.id}-copy-{n}'
            used.add(new_id); changes['id'] = new_id; renamed[address] = new_id
        for f in fields(value):
            members = getattr(value,f.name)
            if isinstance(members,tuple) and any(hasattr(v,'id') for v in members):
                changes[f.name] = tuple(clone(v,address+(f.name,Id(v.id))) for v in members)
        return replace(value,**changes)
    duplicate = clone(source,path)
    result = add_object(config,path[:-1],duplicate)
    if copy_refinement:
        def copied(items):
            extras = []
            for item in items:
                ids = list(item.target.identifiers)
                for i,ref in _target_references(config,item.target):
                    if ref in renamed: ids[i] = renamed[ref]
                if tuple(ids) != item.target.identifiers:
                    extras.append(replace(item,target=replace(item.target,identifiers=tuple(ids))))
            return items+tuple(extras)
        result = replace(result,refinement=replace(result.refinement,
            parameters=copied(config.refinement.parameters),
            bounds=replace(result.refinement.bounds,specific=copied(config.refinement.bounds.specific))))
    return result


def set_refinement(config,target,*,enabled,range_factors=None,plamda=None):
    """Manage a binding by physical identity. Never alter the category filter.

    New bindings require explicit factors and PLAMDA; existing settings are kept
    if omitted. Disabling is allowed even for dangling targets in a draft.
    """
    if type(enabled) is not bool: raise TypeError('enabled must be bool')
    if not isinstance(target,ParameterTarget): raise TypeError('target must be ParameterTarget')
    items = config.refinement.parameters
    matches = [b for b in items if b.target == target]
    if len(matches) > 1: raise ValueError('duplicate bindings must be resolved first')
    if not enabled:
        return replace(config,refinement=replace(config.refinement,parameters=tuple(b for b in items if b.target != target)))
    # Discover targets using the compiler catalogue with bindings/filter removed;
    # selection errors must not prevent correcting an otherwise valid draft.
    from .config_targets import available_refinement_targets
    from .configuration import BoundsConfig
    base = replace(config,refinement=replace(config.refinement,parameters=(),active_categories=None,bounds=BoundsConfig()))
    info = next((v for v in available_refinement_targets(base) if v.target == target),None)
    if info is None: raise ValueError('unknown refinement target')
    if not info.available: raise ValueError(info.reason)
    previous = matches[0] if matches else None
    factors = previous.range_factors if range_factors is None and previous else range_factors
    step = previous.plamda if plamda is None and previous else plamda
    if factors is None or step is None: raise ValueError('new bindings require explicit range_factors and plamda')
    binding = ParameterBinding(target,factors,step)
    new = tuple(binding if b.target == target else b for b in items) if previous else items+(binding,)
    return replace(config,refinement=replace(config.refinement,parameters=new))
