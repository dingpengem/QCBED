"""Version-stable frontend search metadata; not physical bounds or new algorithms.

Reference profiles are documented in docs/modern_refinement_catalogue.md.
Explicit bindings always take precedence. Keep these values stable for replay.
"""
from .configuration import ParameterBinding


def recommended_binding(target):
    category=target.category.value
    steps={'depth':.5,'debye_waller':.0002,'absorption':.005,'scan_coordinate':.0005}
    if category in steps:
        return ParameterBinding(target,(1.05,.95),steps[category])
    if category=='voltage':
        # Existing dynamic voltage regression profile, not a fitted optimum.
        return ParameterBinding(target,(.9,1.1),.01)
    raise ValueError(f'No audited recommended search metadata for {category}; specify advanced settings.')


def binding_settings(binding):
    return (*binding.range_factors,binding.plamda)
