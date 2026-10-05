"""Present existing project hard-bound rules, without inventing physical limits."""
from qcbed.refinement_bounds import RefinementBoundRules, resolve_refinement_bounds
from qcbed.refinement import RefinementParameterSpec, RefinementParameterLayout
from qcbed.refinement_categories import PARAMETER_SEMANTICS


def configured_bounds(config, info):
    target=info.target
    for bound in config.refinement.bounds.specific:
        if bound.target==target:
            return bound.lower,bound.upper,'Configured target bounds'
    rules=config.refinement.bounds
    category=target.category.value
    absolute=tuple((b.category.value,(b.lower,b.upper)) for b in rules.absolute)
    relative=tuple((b.category.value,b.fraction) for b in rules.relative)
    fallback=tuple((b.category.value,b.absolute_fallback) for b in rules.relative)
    spec=RefinementParameterSpec(0,category,(0,)*PARAMETER_SEMANTICS[category].target_arity,info.starting_value,None)
    resolved=resolve_refinement_bounds(RefinementParameterLayout((spec,)),
        RefinementBoundRules(category_bounds=absolute,category_relative_bounds=relative,
                            relative_bound_absolute_fallback=fallback))
    if resolved.entries:
        b=resolved.entries[0]
        return b.lower,b.upper,'Configured category bounds' if category in dict(absolute) else 'Configured relative bounds'
    return None,None,'No canonical hard-bound default — enter specimen-appropriate bounds'
