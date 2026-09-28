"""User-selected active categories, derived without changing input semantics."""

from dataclasses import replace
from .refinement import build_refinement_parameter_layout
from .refinement_categories import PARAMETER_SEMANTICS


def select_refinement_categories(state, config, categories):
    """Return immutable state/config copies with excluded masks disabled.

    Original physical values, source parameter order and mapping stay unchanged.
    Range pairs and PLAMDA are selected using original active-vector positions.
    """
    if isinstance(categories, str):
        raise TypeError("categories must be a collection of category names")
    selected = frozenset(categories)
    if not selected or not selected <= PARAMETER_SEMANTICS.keys():
        raise ValueError("categories must select known refinement categories")
    layout = build_refinement_parameter_layout(state)
    indices = [p.vector_index for p in layout.parameters if p.category in selected]
    if not indices or len(config.plamda) != len(layout.parameters):
        raise ValueError("subset requires active parameters and matching original PLAMDA")
    model = replace(state.model, energy_refined=state.model.energy_refined and "voltage" in selected,
        atoms_by_slice=tuple(tuple(replace(a,
            x_mask=a.x_mask if "atom_coordinate" in selected else 0,
            y_mask=a.y_mask if "atom_coordinate" in selected else 0) for a in row) for row in state.model.atoms_by_slice))
    regions = tuple(replace(r, thickness_mask=r.thickness_mask if "depth" in selected else 0,
        vacancy_masks=r.vacancy_masks if ("vacancy_coefficient" if r.vacancy_flag else "constant_vacancy") in selected else (0,)*len(r.vacancy_masks),
        debye_waller_masks=r.debye_waller_masks if "debye_waller" in selected else (0,)*len(r.debye_waller_masks)) for r in state.regions)
    subset = replace(state, model=model, regions=regions,
        dilation_mask=state.dilation_mask if "dilation" in selected else 0,
        absorption_masks=state.absorption_masks if "absorption" in selected else tuple((0,0,0) for _ in state.absorption_masks),
        scan_masks=state.scan_masks if "scan_coordinate" in selected else ("NN",)*3,
        reflections=state.reflections if "reflection_multiplier" in selected else tuple(tuple(replace(r, mask="NN") for r in row) for row in state.reflections),
        ranges=tuple(state.ranges[i] for i in indices) if state.ranges else ())
    return subset, replace(config, plamda=tuple(config.plamda[i] for i in indices))
