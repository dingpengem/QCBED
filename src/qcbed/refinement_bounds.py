"""User-selected production constraints, separate from legacy parameter mapping."""
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, replace
from collections.abc import Mapping

from .refinement import (RefinementParameterSpec, RefinementParameterLayout, RefinementState,
    build_refinement_parameter_layout, pack_refinement_parameters, apply_refinement_parameters, _entries)
from .structure_factors import _finite
from .refinement_categories import PARAMETER_SEMANTICS


@dataclass(frozen=True)
class RefinementBound:
    category: str
    target: tuple[int, ...]
    lower: float | None = None
    upper: float | None = None

    def __post_init__(self):
        spec = RefinementParameterSpec(0, self.category, self.target, 0., None)
        object.__setattr__(self, "target", spec.target)
        for name in ("lower", "upper"):
            if getattr(self, name) is not None:
                object.__setattr__(self, name, _finite(getattr(self, name), name))
        if self.lower is None and self.upper is None:
            raise ValueError("bound requires lower or upper; omit unbounded parameters")
        if self.lower is not None and self.upper is not None and self.lower > self.upper:
            raise ValueError("bound requires lower <= upper")


@dataclass(frozen=True)
class BoundViolation:
    category: str
    target: tuple[int, ...]
    value: float
    lower: float | None
    upper: float | None


@dataclass(frozen=True)
class RefinementBounds:
    entries: tuple[RefinementBound, ...] = ()

    def __post_init__(self):
        object.__setattr__(self, "entries", tuple(self.entries))
        if any(not isinstance(b, RefinementBound) for b in self.entries):
            raise TypeError("entries must contain RefinementBound values")
        if len({(b.category, b.target) for b in self.entries}) != len(self.entries):
            raise ValueError("duplicate bound category/target")

    @classmethod
    def from_dict(cls, value):
        """Restore the dataclass JSON representation saved in run configuration."""
        return cls(tuple(RefinementBound(**b) for b in value["entries"]))

    def validate(self, state):
        """Validate actual targets even if inactive; inactive bounds are metadata."""
        known = {("voltage", ()), ("dilation", ())}
        for r, region in enumerate(state.regions):
            known.add(("depth", (r,)))
            category = "vacancy_coefficient" if region.vacancy_flag else "constant_vacancy"
            known.update((category, (r, j)) for j in range(len(region.vacancy_values)))
            known.update(("debye_waller", (r, j)) for j in range(len(region.debye_waller)))
        known.update(("absorption", (i, j)) for i, row in enumerate(state.absorption) for j in range(len(row)))
        known.update(("atom_coordinate", (i, j, k)) for i, row in enumerate(state.model.atoms_by_slice)
                     for j in range(len(row)) for k in range(2))
        known.update(("scan_coordinate", (i, j)) for i in range(3) for j in range(2))
        known.update(("reflection_multiplier", (r, i, j)) for r, row in enumerate(state.reflections)
                     for i in range(len(row)) for j in range(2))
        for b in self.entries:
            if (b.category, b.target) not in known:
                raise ValueError(f"unknown bound parameter target: {b.category} {b.target}")

    def violations(self, mapped):
        active = {(c, t): v for c, t, v in _entries(mapped)}
        return tuple(BoundViolation(b.category, b.target, active[b.category, b.target], b.lower, b.upper)
            for b in self.entries if (b.category, b.target) in active
            and ((b.lower is not None and active[b.category, b.target] < b.lower)
                 or (b.upper is not None and active[b.category, b.target] > b.upper)))


@dataclass(frozen=True)
class RefinementBoundRules:
    """Mapping inputs are copied into immutable tuples; no live dict references."""
    specific_bounds: tuple = ()
    category_bounds: tuple = ()
    category_relative_bounds: tuple = ()
    relative_bound_absolute_fallback: tuple = ()

    def __post_init__(self):
        specific = self.specific_bounds
        if isinstance(specific, Mapping):
            specific = tuple(RefinementBound(c, t, *pair) for (c, t), pair in specific.items())
        specific = RefinementBounds(specific).entries
        object.__setattr__(self, "specific_bounds", specific)
        for field in ("category_bounds", "category_relative_bounds", "relative_bound_absolute_fallback"):
            value = getattr(self, field)
            items = tuple(value.items()) if isinstance(value, Mapping) else tuple(value)
            normalized = []
            for category, setting in items:
                if category not in PARAMETER_SEMANTICS:
                    raise ValueError(f"unknown bound category: {category}")
                if field == "category_bounds":
                    bound = RefinementBound(category, (0,)*PARAMETER_SEMANTICS[category].target_arity, *setting)
                    setting = (bound.lower, bound.upper)
                else:
                    setting = _finite(setting, field)
                    if setting <= 0:
                        raise ValueError(f"{field} must be positive")
                normalized.append((category, setting))
            if len({c for c, _ in normalized}) != len(normalized):
                raise ValueError(f"duplicate category in {field}")
            object.__setattr__(self, field, tuple(normalized))
        relative = dict(self.category_relative_bounds)
        fallback = dict(self.relative_bound_absolute_fallback)
        if set(relative) != set(fallback):
            raise ValueError("each relative category requires exactly one explicit absolute fallback")

    @classmethod
    def from_dict(cls, value):
        value = dict(value)
        value["specific_bounds"] = tuple(RefinementBound(**b) for b in value["specific_bounds"])
        return cls(**value)


def resolve_refinement_bounds(layout, rules):
    """Resolve specific > category absolute > category relative > unbounded.

    Accept a RefinementState to derive starting physical values through existing
    source mapping. Alternatively accept a layout whose initial_value fields
    already contain mapped physical starting values (not simplex coordinates).
    Inactive specific targets remain in rules, not the active resolved table.
    """
    if not isinstance(rules, RefinementBoundRules):
        raise TypeError("rules must be RefinementBoundRules")
    if isinstance(layout, RefinementState):
        state = layout
        RefinementBounds(rules.specific_bounds).validate(state)
        original = build_refinement_parameter_layout(state)
        mapped = apply_refinement_parameters(state, pack_refinement_parameters(state))
        values = {(c, t): v for c, t, v in _entries(mapped)}
        layout = RefinementParameterLayout(tuple(replace(p, initial_value=values[p.category, p.target])
            for p in original.parameters))
    if not isinstance(layout, RefinementParameterLayout):
        raise TypeError("layout must be RefinementParameterLayout or RefinementState")
    specific = {(b.category, b.target): b for b in rules.specific_bounds}
    absolute = dict(rules.category_bounds)
    relative = dict(rules.category_relative_bounds)
    fallback = dict(rules.relative_bound_absolute_fallback)
    resolved = []
    for p in layout.parameters:
        key = (p.category, p.target)
        if key in specific:
            resolved.append(specific[key])
        elif p.category in absolute:
            resolved.append(RefinementBound(*key, *absolute[p.category]))
        elif p.category in relative:
            width = max(abs(p.initial_value)*relative[p.category], fallback[p.category])
            resolved.append(RefinementBound(*key, p.initial_value-width, p.initial_value+width))
    return RefinementBounds(tuple(resolved))


def scan_relative_bounds(state, *, relative_fraction, absolute_fallback):
    """Active scan coordinates: width=max(abs(start)*fraction, absolute_fallback).

    Both controls are explicit user choices. Coordinates/fallback are undilated
    fractional reciprocal-cell coefficients (dimensionless); widths are positive.
    """
    fraction = _finite(relative_fraction, "relative_fraction")
    fallback = _finite(absolute_fallback, "absolute_fallback")
    if fraction <= 0 or fallback <= 0:
        raise ValueError("relative_fraction and absolute_fallback must be positive")
    bounds = []
    for c, t, v in _entries(state):
        if c == "scan_coordinate":
            width = max(abs(v)*fraction, fallback)
            bounds.append(RefinementBound(c, t, v-width, v+width))
    return RefinementBounds(tuple(bounds))


_scope = ContextVar("qcbed_refinement_bounds", default=("none", RefinementBounds(), None))


def _validate_policy(policy, bounds, state):
    if policy not in ("none", "reject"):
        raise ValueError("bounds_policy must be 'none' or 'reject'")
    if bounds is None:
        bounds = RefinementBounds()
    if not isinstance(bounds, RefinementBounds):
        raise TypeError("refinement_bounds must be RefinementBounds")
    bounds.validate(state)
    return bounds


@contextmanager
def refinement_bound_policy(state, *, bounds_policy="none", refinement_bounds=None,
                            invalid_trial_policy=None):
    """Scoped trial policy for unchanged optimizer wrappers; reset even on error.

    Standalone optimizers need invalid_trial_policy='reject' to accept inf scores.
    This scope's optional numerical policy keeps numerical rejection independent.
    The driver configures both automatically. No global or worker mutation.
    """
    bounds = _validate_policy(bounds_policy, refinement_bounds, state)
    if invalid_trial_policy not in (None, "raise", "reject"):
        raise ValueError("invalid_trial_policy must be 'raise' or 'reject'")
    token = _scope.set((bounds_policy, bounds, invalid_trial_policy))
    try:
        yield
    finally:
        _scope.reset(token)
