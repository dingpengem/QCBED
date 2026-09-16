"""Source-backed 2D reciprocal geometry; no propagation or cutoff masks."""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from numbers import Real

import numpy as np
from numpy.typing import NDArray

from .grids import reciprocal_index_grid


FloatGrid = NDArray[np.float64]
_SOURCE_PI = 3.14159265358


def _finite(value: float, name: str, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a real scalar, not bool")
    try:
        result = float(value)
    except OverflowError as exc:
        raise ValueError(f"{name} must be finite") from exc
    if not math.isfinite(result) or (positive and result <= 0):
        raise ValueError(f"{name} must be finite" + (" and positive" if positive else ""))
    return result


def _immutable_grid(values: FloatGrid) -> FloatGrid:
    """Own immutable byte storage; even setflags(write=True) cannot mutate it."""
    values = np.asarray(values, dtype=np.float64, order="C")
    if not np.all(np.isfinite(values)):
        raise ValueError("reciprocal geometry must remain finite in float64")
    return np.frombuffer(values.tobytes(order="C"), dtype=np.float64).reshape(values.shape)


@dataclass(frozen=True)
class DirectLattice2D:
    """Direct a/b lengths in angstroms and included gamma in degrees.

    Cartesian frame: a=(a,0), b=(b*cos(gamma),b*sqrt(1-cos(gamma)^2)).
    The angle uses the readable source's printed pi literal, not math.pi.
    """

    a_angstrom: float
    b_angstrom: float
    gamma_degrees: float
    cos_gamma: float = field(init=False)
    sin_gamma: float = field(init=False)
    area_angstrom_squared: float = field(init=False)

    def __post_init__(self) -> None:
        a = _finite(self.a_angstrom, "a_angstrom", positive=True)
        b = _finite(self.b_angstrom, "b_angstrom", positive=True)
        gamma = _finite(self.gamma_degrees, "gamma_degrees")
        if not 0 < gamma < 180:
            raise ValueError("gamma_degrees must be between 0 and 180")
        cosine = math.cos(_SOURCE_PI * gamma / 180.0)
        sine_squared = 1.0 - cosine * cosine
        if sine_squared <= 0:
            raise ValueError("gamma_degrees is degenerate at float64 precision")
        sine = math.sqrt(sine_squared)
        area = _finite(a * b * sine, "area_angstrom_squared", positive=True)
        for name, value in (("a_angstrom", a), ("b_angstrom", b),
                            ("gamma_degrees", gamma), ("cos_gamma", cosine),
                            ("sin_gamma", sine), ("area_angstrom_squared", area)):
            object.__setattr__(self, name, value)

    @property
    def a_vector_angstrom(self) -> tuple[float, float]:
        return self.a_angstrom, 0.0

    @property
    def b_vector_angstrom(self) -> tuple[float, float]:
        return self.b_angstrom * self.cos_gamma, self.b_angstrom * self.sin_gamma


@dataclass(frozen=True)
class ReciprocalLattice2D:
    """CELL2R lengths (inverse angstroms) and dimensionless cosine; no 2*pi."""

    direct: DirectLattice2D
    a_star_inverse_angstrom: float = field(init=False)
    b_star_inverse_angstrom: float = field(init=False)
    cos_gamma_star: float = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.direct, DirectLattice2D):
            raise TypeError("direct must be DirectLattice2D")
        area = self.direct.area_angstrom_squared
        a_star = _finite(self.direct.b_angstrom / area, "a_star_inverse_angstrom", positive=True)
        b_star = _finite(self.direct.a_angstrom / area, "b_star_inverse_angstrom", positive=True)
        object.__setattr__(self, "a_star_inverse_angstrom", a_star)
        object.__setattr__(self, "b_star_inverse_angstrom", b_star)
        object.__setattr__(self, "cos_gamma_star", -self.direct.cos_gamma)

    @property
    def a_star_vector_inverse_angstrom(self) -> tuple[float, float]:
        return (self.a_star_inverse_angstrom * self.direct.sin_gamma,
                -self.a_star_inverse_angstrom * self.direct.cos_gamma)

    @property
    def b_star_vector_inverse_angstrom(self) -> tuple[float, float]:
        return 0.0, self.b_star_inverse_angstrom


@dataclass(frozen=True, eq=False)
class ReciprocalIndexGrid:
    """Dimensionless integer-valued H/K, stored as immutable float64 grids.

    Axis 0 is row/K, axis 1 column/H. Shapes are (mesh_y, mesh_x).
    Ordering comes directly from the existing generic grid helper.
    """

    mesh_x: int
    mesh_y: int
    centered: bool = True
    h: FloatGrid = field(init=False, repr=False)
    k: FloatGrid = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if not isinstance(self.centered, bool):
            raise TypeError("centered must be bool")
        h, k = reciprocal_index_grid(self.mesh_x, self.mesh_y, centered=self.centered)
        object.__setattr__(self, "h", _immutable_grid(h))
        object.__setattr__(self, "k", _immutable_grid(k))


@dataclass(frozen=True, eq=False)
class ReciprocalVectorGrid:
    """Physical g components, magnitude and STL, all in inverse angstroms.

    Arrays are separate from dimensionless indices and immutable. Dilation
    is APPLIEDDILATION: g scales by 1/dilation and |g|^2 by 1/dilation^2.
    Magnitude uses the source metric and STL=sqrt(0.25*|g|^2).
    """

    lattice: ReciprocalLattice2D
    indices: ReciprocalIndexGrid
    dilation: float = 1.0
    gx_inverse_angstrom: FloatGrid = field(init=False, repr=False)
    gy_inverse_angstrom: FloatGrid = field(init=False, repr=False)
    magnitude_inverse_angstrom: FloatGrid = field(init=False, repr=False)
    stl_inverse_angstrom: FloatGrid = field(init=False, repr=False)

    def __post_init__(self) -> None:
        if not isinstance(self.lattice, ReciprocalLattice2D):
            raise TypeError("lattice must be ReciprocalLattice2D")
        if not isinstance(self.indices, ReciprocalIndexGrid):
            raise TypeError("indices must be ReciprocalIndexGrid")
        dilation = _finite(self.dilation, "dilation", positive=True)
        dilation_squared = _finite(dilation * dilation, "dilation squared", positive=True)
        object.__setattr__(self, "dilation", dilation)
        h, k = self.indices.h, self.indices.k
        a, b = self.lattice.a_star_inverse_angstrom, self.lattice.b_star_inverse_angstrom
        cosine = self.lattice.cos_gamma_star
        ax, ay = self.lattice.a_star_vector_inverse_angstrom
        bx, by = self.lattice.b_star_vector_inverse_angstrom
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            try:
                gx = (h * ax + k * bx) / dilation
                gy = (h * ay + k * by) / dilation
                squared = (h*h*a*a + k*k*b*b + 2.0*h*k*a*b*cosine) / dilation_squared
                magnitude = np.sqrt(squared)
                stl = np.sqrt(0.25 * squared)
            except FloatingPointError as exc:
                raise ValueError("reciprocal metric is not representable in float64") from exc
        for name, values in (("gx_inverse_angstrom", gx), ("gy_inverse_angstrom", gy),
                             ("magnitude_inverse_angstrom", magnitude),
                             ("stl_inverse_angstrom", stl)):
            object.__setattr__(self, name, _immutable_grid(values))
