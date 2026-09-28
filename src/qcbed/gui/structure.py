"""Transactional structural editing. Canonical helpers own identity/dependency rules."""
from dataclasses import replace
from qcbed import (Id, get_object, replace_object, add_object, delete_object,
    duplicate_object, reorder_objects, deletion_dependencies, set_refinement)


def selector(record):
    for field in ('id', 'species', 'target', 'category'):
        if hasattr(record, field): return Id(getattr(record, field), field)
    raise TypeError('Record has no canonical identity')


def reference_path(path):
    from qcbed.configuration import ParameterTarget
    return tuple(Id(p, 'target') if isinstance(p, ParameterTarget) else p for p in path)


def removal_plan(config, path):
    """Transitive external references, including bindings on dependent atoms."""
    result = []
    seen = set()
    def visit(current):
        for dependency in deletion_dependencies(config, current):
            address = reference_path(dependency.object_path)
            if address not in seen:
                seen.add(address)
                visit(address)
                result.append(dependency)
    visit(tuple(path))
    return tuple(result)


class StructuralDraft:
    def __init__(self, config): self.config = config

    def update(self, path, **changes):
        old = get_object(self.config, path)
        # Raw spellings belong to the associated presentation value only.
        if 'label' in changes and changes['label'] != old.label: changes['label_raw'] = ''
        if 'record_type' in changes and changes['record_type'] != old.record_type: changes['record_type_raw'] = ''
        if 'enabled' in changes and hasattr(old, 'legacy_raw'): changes['legacy_raw'] = None
        self.config = replace_object(self.config, path, replace(old, **changes))

    def add(self, path, record):
        if hasattr(record, 'id'):
            self.config = add_object(self.config, path, record)
        else:
            items = get_object(self.config, path)
            key = selector(record)
            if any(selector(item) == key for item in items): raise ValueError('Duplicate record identity')
            self.config = replace_object(self.config, path, items + (record,))

    def duplicate(self, path, *, copy_refinement=False):
        self.config = duplicate_object(self.config, path, copy_refinement=copy_refinement)

    def move(self, path, offset):
        items = get_object(self.config, path[:-1])
        ids = [v.id for v in items]
        index = ids.index(path[-1].value)
        destination = index + offset
        if not 0 <= destination < len(ids): return
        ids[index], ids[destination] = ids[destination], ids[index]
        self.config = reorder_objects(self.config, path[:-1], ids)

    def delete(self, path, *, remove_dependencies=False):
        plan = removal_plan(self.config, path)
        if plan and not remove_dependencies:
            from qcbed.config_editing import DependencyError
            raise DependencyError(plan)
        candidate = self.config
        for dependency in plan:
            address = reference_path(dependency.object_path)
            candidate = delete_object(candidate, address)
        self.config = delete_object(candidate, path)


def resolved_bound_rows(config):
    """Present compiler-resolved limits and the configured rule responsible."""
    from qcbed import compile_config, available_refinement_targets
    from .document import target_label
    compiled = compile_config(config)
    bounds = {(b.category,b.target):b for b in compiled.bounds.entries}
    rules = config.refinement.bounds
    rows = []
    active = (info for info in available_refinement_targets(compiled) if info.active)
    for info, parameter in zip(active,compiled.layout.parameters):
        bound = bounds.get((parameter.category,parameter.target))
        if any(b.target == info.target for b in rules.specific): source = 'specific target'
        elif any(b.category == info.category for b in rules.absolute): source = 'category absolute'
        elif any(b.category == info.category for b in rules.relative): source = 'category relative'
        else: source = 'unbounded'
        rows.append((target_label(config,info),bound.lower if bound else None,
                     bound.upper if bound else None,source+f'; policy={rules.policy}'))
    return tuple(rows)
