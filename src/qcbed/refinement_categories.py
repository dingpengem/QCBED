"""Physical taxonomy, distinct from concrete targets and active vector size.

Order groups describe _entries, not a sorting replacement. Vacancy variants
share a group and interleave by region. Arrays are always sized from the layout.
"""
from dataclasses import dataclass
from enum import Enum
from types import MappingProxyType


class ParameterCategory(str, Enum):
    VOLTAGE = 'voltage'
    DEPTH = 'depth'
    CONSTANT_VACANCY = 'constant_vacancy'
    VACANCY_COEFFICIENT = 'vacancy_coefficient'
    DILATION = 'dilation'
    ABSORPTION = 'absorption'
    DEBYE_WALLER = 'debye_waller'
    ATOM_COORDINATE = 'atom_coordinate'
    SCAN_COORDINATE = 'scan_coordinate'
    REFLECTION_MULTIPLIER = 'reflection_multiplier'


@dataclass(frozen=True)
class ParameterSemantics:
    target_arity: int
    source_order_group: int
    write_back_mapped_value: bool = False


PARAMETER_SEMANTICS = MappingProxyType({
    'voltage': ParameterSemantics(0, 0),
    'depth': ParameterSemantics(1, 1),
    'constant_vacancy': ParameterSemantics(2, 2, True),
    'vacancy_coefficient': ParameterSemantics(2, 2, True),
    'dilation': ParameterSemantics(0, 3),
    'absorption': ParameterSemantics(2, 4, True),
    'debye_waller': ParameterSemantics(2, 5, True),
    'atom_coordinate': ParameterSemantics(3, 6),
    'scan_coordinate': ParameterSemantics(2, 7),
    'reflection_multiplier': ParameterSemantics(3, 8),
})

# Existing adapter guard, NOT a compiler/layout limit or allocation size.
# Removing this historical public-API restriction is a separate change.
LEGACY_OPTIMIZER_PARAMETER_LIMIT = 1000
