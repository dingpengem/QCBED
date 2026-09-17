"""Immutable active-parameter mapping from ARRANGEPARS and UPDPARS.

Mapping covers all packed scalar categories. Trial preparation is a separately
validated subset; it never guesses unresolved vacancy/free-surface state.
"""

from dataclasses import dataclass, replace
from numbers import Integral

import numpy as np

from .legacy_io import LegacyInput
from .scan import LegacyScanDefinition
from .multislice import _snapshot
from .structure_factors import _finite


def _flag(value, name, allowed=(0, 1)):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral) or value not in allowed:
        raise ValueError(f"{name} must be one of {allowed}, excluding bool")
    return int(value)


def _pair_mask(value):
    if value not in ("NN", "YN", "NY", "YY"):
        raise ValueError("mask must be NN, YN, NY or YY")
    return value


@dataclass(frozen=True)
class RefinementRegion:
    depth_angstrom: float
    thickness_mask: int
    slice_count: int
    vacancy_flag: int
    vacancy_values: tuple[float, ...]
    vacancy_masks: tuple[int, ...]
    debye_waller: tuple[float, ...]
    debye_waller_masks: tuple[int, ...]
    free_surface: tuple[float, ...]  # FSTYPE, BINDA, BINDB, BINDC

    def __post_init__(self):
        object.__setattr__(self, "depth_angstrom", _finite(self.depth_angstrom, "depth_angstrom"))
        _flag(self.thickness_mask, "thickness_mask")
        _flag(self.vacancy_flag, "vacancy_flag", (0, 1, 2))
        if isinstance(self.slice_count, (bool, np.bool_)) or not isinstance(self.slice_count, Integral) or self.slice_count < 1:
            raise ValueError("slice_count must be a positive integer")
        for field in ("vacancy_values", "debye_waller", "free_surface"):
            object.__setattr__(self, field, tuple(_finite(v, field) for v in getattr(self, field)))
        for field in ("vacancy_masks", "debye_waller_masks"):
            object.__setattr__(self, field, tuple(_flag(v, field) for v in getattr(self, field)))
        n = 2 if self.vacancy_flag else 1
        if len(self.vacancy_values) != n or len(self.vacancy_masks) != n:
            raise ValueError("vacancy values/masks must match vacancy_flag")
        if len(self.debye_waller) != len(self.debye_waller_masks):
            raise ValueError("Debye-Waller values/masks must match")
        if len(self.free_surface) != 4 or self.free_surface[0] not in (0, 1, 2):
            raise ValueError("free_surface must contain FSTYPE (0/1/2), A, B, C")


@dataclass(frozen=True)
class RefinementReflection:
    h: int
    k: int
    elastic: float
    absorptive: float
    mask: str

    def __post_init__(self):
        for name in ("h", "k"):
            value = getattr(self, name)
            if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral):
                raise TypeError(f"{name} must be an integer")
        _pair_mask(self.mask)
        for name in ("elastic", "absorptive"):
            object.__setattr__(self, name, _finite(getattr(self, name), name))


@dataclass(frozen=True)
class RefinementState:
    model: LegacyInput
    regions: tuple[RefinementRegion, ...]
    dilation: float
    dilation_mask: int
    absorption: tuple[tuple[float, float, float], ...]
    absorption_masks: tuple[tuple[int, int, int], ...]
    scan: LegacyScanDefinition
    scan_masks: tuple[str, str, str]
    reflections: tuple[tuple[RefinementReflection, ...], ...]
    symmetry: tuple[str, ...]
    dsap: float
    iref: float
    ranges: tuple[tuple[float, float], ...] = ()  # D1,D2, not enforced bounds

    def __post_init__(self):
        if not isinstance(self.model, LegacyInput) or not isinstance(self.scan, LegacyScanDefinition):
            raise TypeError("model/scan must be LegacyInput/LegacyScanDefinition")
        object.__setattr__(self, "model", replace(self.model,
            slices=tuple(self.model.slices), species=tuple(self.model.species),
            zone_axis=tuple(self.model.zone_axis),
            atoms_by_slice=tuple(tuple(row) for row in self.model.atoms_by_slice)))
        for name in ("regions", "scan_masks", "symmetry"):
            object.__setattr__(self, name, tuple(getattr(self, name)))
        for name in ("absorption", "absorption_masks", "reflections", "ranges"):
            object.__setattr__(self, name, tuple(tuple(row) for row in getattr(self, name)))
        if not self.regions or any(not isinstance(r, RefinementRegion) for r in self.regions):
            raise TypeError("regions must contain RefinementRegion values")
        if sum(r.slice_count for r in self.regions) != len(self.model.slices):
            raise ValueError("region slice counts must cover model slices")
        if len(self.reflections) != len(self.regions) or len(self.symmetry) != len(self.regions):
            raise ValueError("reflection/symmetry groups must match regions")
        if any(not isinstance(v, RefinementReflection) for row in self.reflections for v in row):
            raise TypeError("reflections must contain RefinementReflection values")
        if any(s not in ("NONE", "INVE", "SQUA", "RECT", "HEXA") for s in self.symmetry):
            raise ValueError("unsupported symmetry")
        if len(self.scan_masks) != 3:
            raise ValueError("three scan masks required")
        for value in self.scan_masks:
            _pair_mask(value)
        t = len(self.model.species)
        if any(len(r.debye_waller) != t for r in self.regions):
            raise ValueError("Debye-Waller counts must match species")
        n = t if self.model.absorption_enabled else 0
        if len(self.absorption) != n or len(self.absorption_masks) != n:
            raise ValueError("absorption values/masks must match raw flag and species")
        for row, masks in zip(self.absorption, self.absorption_masks):
            if len(row) != 3 or len(masks) != 3:
                raise ValueError("absorption requires A,B,C")
            for j, (v, mask) in enumerate(zip(row, masks)):
                _finite(v, "absorption")
                _flag(mask, "absorption mask", (-1, 0, 1) if j == 2 else (0, 1))
        for row in self.ranges:
            if len(row) != 2:
                raise ValueError("ranges require D1,D2 pairs")
            for v in row:
                _finite(v, "range")
        for name in ("dilation", "dsap", "iref"):
            object.__setattr__(self, name, _finite(getattr(self, name), name))
        _flag(self.dilation_mask, "dilation_mask")


@dataclass(frozen=True)
class RefinementParameterSpec:
    """Zero-based vector index and zero-based physical target indices."""
    vector_index: int
    category: str
    target: tuple[int, ...]
    initial_value: float
    range_factors: tuple[float, float] | None

    def __post_init__(self):
        if isinstance(self.vector_index, (bool, np.bool_)) or not isinstance(self.vector_index, Integral) or self.vector_index < 0:
            raise ValueError("vector_index must be a nonnegative integer")
        target = tuple(self.target)
        lengths = dict(voltage=0, depth=1, vacancy_coefficient=2, constant_vacancy=2,
                       dilation=0, absorption=2, debye_waller=2, atom_coordinate=3,
                       scan_coordinate=2, reflection_multiplier=3)
        if self.category not in lengths or len(target) != lengths[self.category]:
            raise ValueError("unknown category or incorrect physical target")
        if any(isinstance(i, (bool, np.bool_)) or not isinstance(i, Integral) or i < 0 for i in target):
            raise ValueError("target indices must be nonnegative integers")
        object.__setattr__(self, "target", target)
        object.__setattr__(self, "initial_value", _finite(self.initial_value, "initial_value"))
        if self.range_factors is not None:
            pair = tuple(_finite(v, "range factor") for v in self.range_factors)
            if len(pair) != 2:
                raise ValueError("range_factors requires D1,D2")
            object.__setattr__(self, "range_factors", pair)


@dataclass(frozen=True)
class RefinementParameterLayout:
    parameters: tuple[RefinementParameterSpec, ...]

    def __post_init__(self):
        object.__setattr__(self, "parameters", tuple(self.parameters))
        if any(not isinstance(p, RefinementParameterSpec) or p.vector_index != i
               for i, p in enumerate(self.parameters)):
            raise ValueError("layout requires consecutive source-ordered parameter specs")


def _entries(s):
    if s.model.energy_refined:
        yield "voltage", (), s.model.energy_kv
    for r, region in enumerate(s.regions):
        if region.thickness_mask:
            yield "depth", (r,), region.depth_angstrom
    for r, region in enumerate(s.regions):
        for j, (v, mask) in enumerate(zip(region.vacancy_values, region.vacancy_masks)):
            if mask:
                yield ("vacancy_coefficient" if region.vacancy_flag else "constant_vacancy"), (r, j), v
    if s.dilation_mask:
        yield "dilation", (), s.dilation
    for i, (row, masks) in enumerate(zip(s.absorption, s.absorption_masks)):
        for j, (v, mask) in enumerate(zip(row, masks)):
            if mask == 1:
                yield "absorption", (i, j), v
    for r, region in enumerate(s.regions):
        for j, (v, mask) in enumerate(zip(region.debye_waller, region.debye_waller_masks)):
            if mask:
                yield "debye_waller", (r, j), v
    for i, atoms in enumerate(s.model.atoms_by_slice):
        for j, atom in enumerate(atoms):
            for k, (v, mask) in enumerate(((atom.x_fractional, atom.x_mask), (atom.y_fractional, atom.y_mask))):
                if mask == 1:
                    yield "atom_coordinate", (i, j, k), v
    for i, (corner, mask) in enumerate(zip((s.scan.corner_a, s.scan.corner_b, s.scan.corner_c), s.scan_masks)):
        for j, v in enumerate(corner):
            if mask[j] == "Y":
                yield "scan_coordinate", (i, j), v
    for r, row in enumerate(s.reflections):
        for i, reflection in enumerate(row):
            for j, v in enumerate((reflection.elastic, reflection.absorptive)):
                if reflection.mask[j] == "Y":
                    yield "reflection_multiplier", (r, i, j), v


def _fold_dilation(v):
    v = abs(v)
    return 1/v if v > 1 else v


def _packed(category, v):
    if category in ("vacancy_coefficient", "absorption"):
        v = abs(v)
    if category == "dilation":
        v = _fold_dilation(v)
    return 1e-10 if abs(v) < 1e-10 else v


def build_refinement_parameter_layout(state: RefinementState) -> RefinementParameterLayout:
    if not isinstance(state, RefinementState):
        raise TypeError("state must be RefinementState")
    entries = tuple(_entries(state))
    if state.ranges and len(state.ranges) != len(entries):
        raise ValueError("range count must match active parameter count")
    return RefinementParameterLayout(tuple(RefinementParameterSpec(i, c, t, _packed(c, v),
        state.ranges[i] if state.ranges else None) for i, (c, t, v) in enumerate(entries)))


def pack_refinement_parameters(state: RefinementState) -> np.ndarray:
    """ARRANGEPARS value projection, returned as immutable float64 storage."""
    return _snapshot(np.array([p.initial_value for p in build_refinement_parameter_layout(state).parameters], dtype=np.float64))


def _vector(vector, n):
    if not isinstance(vector, np.ndarray) or vector.dtype != np.float64:
        raise TypeError("trial_vector must be a float64 NumPy array")
    if vector.shape != (n,) or not np.all(np.isfinite(vector)):
        raise ValueError(f"trial_vector must be finite with shape ({n},)")
    return vector.copy()


def _scal(layout):
    a, b = [], []
    for p in layout.parameters:
        if p.range_factors is None:
            raise ValueError("legacy_simplex requires D1,D2 range factors for every parameter")
        d1, d2 = p.range_factors
        if d1 == d2 or p.initial_value == 0:
            raise ValueError("singular SCAL conversion")
        a.append(1.0/(d2-d1)/p.initial_value)
        b.append(d1*p.initial_value)
    return np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)


def apply_refinement_parameters(state: RefinementState, trial_vector: np.ndarray, *,
                                coordinate_mode="physical") -> RefinementState:
    """UPDPARS scalar assignment on a copy; dependent physics is prepared later.

    legacy_simplex uses SCAL anchored to this baseline's packed starting values.
    No D1/D2 clipping occurs. Input vector is never modified or repacked.
    """
    layout = build_refinement_parameter_layout(state)
    values = _vector(trial_vector, len(layout.parameters))
    if coordinate_mode == "legacy_simplex":
        a, b = _scal(layout)
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            values = values/a+b
    elif coordinate_mode != "physical":
        raise ValueError("coordinate_mode must be physical or legacy_simplex")
    regions = list(state.regions)
    absorption = [list(row) for row in state.absorption]
    atoms = [list(row) for row in state.model.atoms_by_slice]
    corners = [list(row) for row in (state.scan.corner_a, state.scan.corner_b, state.scan.corner_c)]
    reflections = [list(row) for row in state.reflections]
    energy, dilation = state.model.energy_kv, state.dilation
    for p, value in zip(layout.parameters, values):
        c, t = p.category, p.target
        v = float(value)
        if c in ("voltage", "depth", "vacancy_coefficient", "absorption", "debye_waller"):
            v = abs(v)
        if c == "voltage":
            energy = v
        elif c == "depth":
            regions[t[0]] = replace(regions[t[0]], depth_angstrom=v)
        elif c in ("vacancy_coefficient", "constant_vacancy"):
            if c == "constant_vacancy":
                if v >= 1:
                    v = 1.0
                if v <= 0:
                    v = abs(v)
            row = list(regions[t[0]].vacancy_values)
            row[t[1]] = v
            regions[t[0]] = replace(regions[t[0]], vacancy_values=tuple(row))
        elif c == "dilation":
            dilation = _fold_dilation(v)
        elif c == "absorption":
            absorption[t[0]][t[1]] = v
        elif c == "debye_waller":
            row = list(regions[t[0]].debye_waller)
            row[t[1]] = v
            regions[t[0]] = replace(regions[t[0]], debye_waller=tuple(row))
        elif c == "atom_coordinate":
            atoms[t[0]][t[1]] = replace(atoms[t[0]][t[1]], **{("x_fractional", "y_fractional")[t[2]]: v})
        elif c == "scan_coordinate":
            corners[t[0]][t[1]] = v
        elif c == "reflection_multiplier":
            reflections[t[0]][t[1]] = replace(reflections[t[0]][t[1]], **{("elastic", "absorptive")[t[2]]: v})
    for i, masks in enumerate(state.absorption_masks):
        if masks[2] == -1:
            absorption[i][2] = absorption[i][1]
    model = replace(state.model, energy_kv=energy, atoms_by_slice=tuple(tuple(row) for row in atoms))
    scan = replace(state.scan, corner_a=tuple(corners[0]), corner_b=tuple(corners[1]), corner_c=tuple(corners[2]))
    return replace(state, model=model, regions=tuple(regions), dilation=dilation,
                   absorption=tuple(tuple(row) for row in absorption), scan=scan,
                   reflections=tuple(tuple(row) for row in reflections))
