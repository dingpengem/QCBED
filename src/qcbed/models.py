"""Typed configuration models for confirmed QCBED concepts."""

from __future__ import annotations

from dataclasses import dataclass, field
import math



def _require_finite(value: float, name: str) -> float:
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


@dataclass(frozen=True)
class MeshConfig:
    """Real-space mesh dimensions, expressed as x and y point counts."""

    mesh_x: int
    mesh_y: int

    def __post_init__(self) -> None:
        if isinstance(self.mesh_x, bool) or not isinstance(self.mesh_x, int):
            raise TypeError("mesh_x must be an integer")
        if isinstance(self.mesh_y, bool) or not isinstance(self.mesh_y, int):
            raise TypeError("mesh_y must be an integer")
        if self.mesh_x <= 0 or self.mesh_y <= 0:
            raise ValueError("mesh_x and mesh_y must be positive")


@dataclass(frozen=True)
class MicroscopeConfig:
    """Microscope settings with electron energy in kV."""

    energy_kv: float

    def __post_init__(self) -> None:
        energy_kv = _require_finite(float(self.energy_kv), "energy_kv")
        if energy_kv <= 0:
            raise ValueError("energy_kv must be positive")
        object.__setattr__(self, "energy_kv", energy_kv)


@dataclass(frozen=True)
class SliceGeometry:
    """Slice lengths are in angstroms; shifts are fractional direct-lattice translations."""

    thickness_angstrom: float
    x_shift_fractional: float = 0.0
    y_shift_fractional: float = 0.0
    cell_a_angstrom: float | None = None
    cell_b_angstrom: float | None = None
    lattice_gamma_degrees: float | None = None

    def __post_init__(self) -> None:
        thickness = _require_finite(float(self.thickness_angstrom), "thickness_angstrom")
        if thickness <= 0:
            raise ValueError("thickness_angstrom must be positive")
        object.__setattr__(self, "thickness_angstrom", thickness)
        object.__setattr__(self, "x_shift_fractional", _require_finite(float(self.x_shift_fractional), "x_shift_fractional"))
        object.__setattr__(self, "y_shift_fractional", _require_finite(float(self.y_shift_fractional), "y_shift_fractional"))
        for name in ("cell_a_angstrom", "cell_b_angstrom", "lattice_gamma_degrees"):
            value = getattr(self, name)
            if value is not None:
                object.__setattr__(self, name, _require_finite(float(value), name))
        for name in ("cell_a_angstrom", "cell_b_angstrom"):
            value = getattr(self, name)
            if value is not None and value <= 0:
                raise ValueError(f"{name} must be positive")
        if self.lattice_gamma_degrees is not None and not 0 < self.lattice_gamma_degrees < 180:
            raise ValueError("lattice_gamma_degrees must be between 0 and 180 degrees")


@dataclass(frozen=True)
class TiltCoordinate:
    """One tilt in the legacy internal coordinate convention."""

    first: float
    second: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "first", _require_finite(float(self.first), "first"))
        object.__setattr__(self, "second", _require_finite(float(self.second), "second"))


@dataclass(frozen=True)
class SimulationConfig:
    """Validated configuration shared by a future simulation workflow."""

    mesh: MeshConfig
    microscope: MicroscopeConfig
    slices: tuple[SliceGeometry, ...]
    tilts: tuple[TiltCoordinate, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        if not isinstance(self.mesh, MeshConfig):
            raise TypeError("mesh must be a MeshConfig")
        if not isinstance(self.microscope, MicroscopeConfig):
            raise TypeError("microscope must be a MicroscopeConfig")
        if not self.slices:
            raise ValueError("slices must not be empty")
        if any(not isinstance(item, SliceGeometry) for item in self.slices):
            raise TypeError("slices must contain SliceGeometry values")
        if any(not isinstance(item, TiltCoordinate) for item in self.tilts):
            raise TypeError("tilts must contain TiltCoordinate values")
