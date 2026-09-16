"""One-tilt execution of precomputed slices, ending at reciprocal beam sampling."""

from collections.abc import Iterable
from dataclasses import dataclass, field
from numbers import Integral

import numpy as np
from numpy.typing import NDArray

from .fft_backend import ComplexGrid, Ordering, _array, _ordering
from .grids import row_col_to_pmap_index, validate_mesh_dimensions
from .wave_propagation import (
    _matching_orderings, apply_phase_grating, apply_propagator,
    real_to_reciprocal_wave, reciprocal_to_real_wave,
)


BeamIndices = tuple[tuple[int, int], ...]


def _integer(value: int, name: str) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral):
        raise TypeError(f"{name} must be an integer, excluding bool")
    return int(value)


def _items(values: Iterable, name: str) -> tuple:
    if isinstance(values, (str, bytes)):
        raise TypeError(f"{name} must be an iterable of entries, not text")
    try:
        return tuple(values)
    except TypeError as exc:
        raise TypeError(f"{name} must be an iterable of entries") from exc


def _sequence(values: Iterable[int], slice_type_count: int, *, one_based: bool) -> tuple[int, ...]:
    count = _integer(slice_type_count, "slice_type_count")
    if count <= 0:
        raise ValueError("slice_type_count must be positive")
    name = "sequence_one_based" if one_based else "slice_sequence_zero_based"
    entries = _items(values, name)
    if not entries:
        raise ValueError(f"{name} must not be empty")
    lower = 1 if one_based else 0
    result = []
    for position, value in enumerate(entries):
        value = _integer(value, f"{name}[{position}]")
        if not lower <= value < count+lower:
            raise ValueError(f"{name}[{position}] must be in {lower}..{count+lower-1}; got {value}")
        result.append(value-lower)
    return tuple(result)


def legacy_slice_sequence_to_zero_based(
    sequence_one_based: Iterable[int], slice_type_count: int,
) -> tuple[int, ...]:
    """Validate active legacy MSEQ entries, then subtract one; no dummy tail."""
    return _sequence(sequence_one_based, slice_type_count, one_based=True)


def unit_incident_wave(mesh_x: int, mesh_y: int) -> ComplexGrid:
    """Construct a unit-amplitude real-space plane wave on (mesh_y,mesh_x)."""
    nx, ny = validate_mesh_dimensions(mesh_x, mesh_y)
    return np.ones((ny, nx), dtype=np.complex128)


def _beams(values: Iterable[tuple[int, int]]) -> BeamIndices:
    pairs = []
    for position, pair in enumerate(_items(values, "beam_indices_hk")):
        pair = _items(pair, f"beam_indices_hk[{position}]")
        if len(pair) != 2:
            raise ValueError(f"beam_indices_hk[{position}] must contain exactly H,K")
        pairs.append((_integer(pair[0], f"beam H[{position}]"),
                      _integer(pair[1], f"beam K[{position}]")))
    return tuple(pairs)


def beam_indices_to_row_col(
    beam_indices_hk: Iterable[tuple[int, int]], *, mesh_x: int, mesh_y: int,
    ordering: Ordering,
) -> tuple[tuple[int, int], ...]:
    """Map canonical signed H/K to row/column, retaining order and duplicates.

    Each signed index is in -(N//2)..(N-1)//2. Out-of-grid aliases, including
    even positive Nyquist, are rejected before any modulo mapping.
    """
    nx, ny = validate_mesh_dimensions(mesh_x, mesh_y)
    _ordering(ordering)
    positions = []
    for h, k in _beams(beam_indices_hk):
        if not -(nx//2) <= h <= (nx-1)//2 or not -(ny//2) <= k <= (ny-1)//2:
            raise ValueError(f"beam (H,K)=({h},{k}) is outside signed grid bounds "
                             f"H={-(nx//2)}..{(nx-1)//2}, K={-(ny//2)}..{(ny-1)//2}")
        positions.append((k+ny//2, h+nx//2) if ordering == "centered" else (k % ny, h % nx))
    return tuple(positions)


def beam_indices_to_legacy_ib(
    beam_indices_hk: Iterable[tuple[int, int]], *, mesh_x: int, mesh_y: int,
) -> tuple[int, ...]:
    """Return one-based legacy PMAP IB values, always using unshifted bins."""
    positions = beam_indices_to_row_col(beam_indices_hk, mesh_x=mesh_x, mesh_y=mesh_y,
                                       ordering="unshifted")
    return tuple(row_col_to_pmap_index(row, col, mesh_x, mesh_y) for row, col in positions)


def extract_beam_amplitudes(
    reciprocal_wave: ComplexGrid, beam_indices_hk: Iterable[tuple[int, int]], *,
    ordering: Ordering,
) -> NDArray[np.complex128]:
    """Copy complex amplitudes in request order, without intensity scaling."""
    wave = _array(reciprocal_wave, "reciprocal_wave")
    ny, nx = wave.shape
    positions = beam_indices_to_row_col(beam_indices_hk, mesh_x=nx, mesh_y=ny, ordering=ordering)
    return np.array([wave[row, col] for row, col in positions], dtype=np.complex128)


def beam_intensities(amplitudes: NDArray[np.complex128]) -> NDArray[np.float64]:
    """Source real**2+imag**2, without RMT, YNORM, or unit-total rescaling."""
    if not isinstance(amplitudes, np.ndarray) or amplitudes.dtype != np.complex128:
        raise TypeError("amplitudes must be a complex128 NumPy array")
    if amplitudes.ndim != 1:
        raise ValueError("amplitudes must have one-dimensional shape (beam_count,)")
    if not np.all(np.isfinite(amplitudes)):
        raise ValueError("amplitudes must be finite")
    with np.errstate(over="ignore", invalid="ignore"):
        intensities = amplitudes.real**2 + amplitudes.imag**2
    if not np.all(np.isfinite(intensities)):
        raise ValueError("beam intensities must remain finite in float64")
    return intensities


def _snapshot(value: np.ndarray) -> np.ndarray:
    return np.frombuffer(value.tobytes(order="C"), dtype=value.dtype).reshape(value.shape)


@dataclass(frozen=True)
class MultisliceResult:
    """Immutable final sampling state; no unused final real wave or wave history."""

    final_reciprocal_wave: ComplexGrid
    beam_indices_hk: BeamIndices
    slice_sequence_zero_based: tuple[int, ...]
    slice_type_count: int
    reciprocal_ordering: Ordering
    amplitudes: NDArray[np.complex128] = field(init=False)
    intensities: NDArray[np.float64] = field(init=False)

    def __post_init__(self) -> None:
        sequence = _sequence(self.slice_sequence_zero_based, self.slice_type_count, one_based=False)
        beams = _beams(self.beam_indices_hk)
        wave = _array(self.final_reciprocal_wave, "final_reciprocal_wave")
        amplitudes = extract_beam_amplitudes(wave, beams, ordering=self.reciprocal_ordering)
        intensities = beam_intensities(amplitudes)
        object.__setattr__(self, "slice_sequence_zero_based", sequence)
        object.__setattr__(self, "slice_type_count", int(self.slice_type_count))
        object.__setattr__(self, "beam_indices_hk", beams)
        for name, value in (("final_reciprocal_wave", wave), ("amplitudes", amplitudes),
                            ("intensities", intensities)):
            object.__setattr__(self, name, _snapshot(value))


def simulate_slice_sequence(
    phase_gratings: Iterable[ComplexGrid], propagators: Iterable[ComplexGrid],
    slice_sequence_zero_based: Iterable[int], *,
    beam_indices_hk: Iterable[tuple[int, int]], reciprocal_ordering: Ordering,
    propagator_ordering: Ordering, incident_wave: ComplexGrid | None = None,
) -> MultisliceResult:
    """Execute precomputed slices for one tilt, stopping after final PRARRAY.

    Collections share stored slice-type indices. Sequence indices are strictly
    zero-based; use the named legacy converter for MSEQ. Both order labels
    are mandatory and must agree. Inputs are never mutated; a custom incident
    wave is a Python extension to the default unit real-space plane wave.
    No construction, sequence padding, final inverse transform or extra RMT.
    """
    _matching_orderings(reciprocal_ordering, propagator_ordering)
    gratings = _items(phase_gratings, "phase_gratings")
    kernels = _items(propagators, "propagators")
    if not gratings or len(gratings) != len(kernels):
        raise ValueError("phase_gratings and propagators must have equal nonzero collection lengths")
    sequence = _sequence(slice_sequence_zero_based, len(gratings), one_based=False)
    shape = _array(gratings[0], "phase_gratings[0]").shape
    for name, collection in (("phase_gratings", gratings), ("propagators", kernels)):
        for index, value in enumerate(collection):
            if _array(value, f"{name}[{index}]").shape != shape:
                raise ValueError(f"{name}[{index}] shape must match all stored slices {shape}")
    beams = _beams(beam_indices_hk)
    beam_indices_to_row_col(beams, mesh_x=shape[1], mesh_y=shape[0], ordering=reciprocal_ordering)
    real_wave = None if incident_wave is None else _array(incident_wave, "incident_wave")
    if real_wave is not None and real_wave.shape != shape:
        raise ValueError("incident_wave shape must match stored slices")
    for occurrence, slice_id in enumerate(sequence):
        # Unit illumination times the first grating is already that grating.
        transmitted = gratings[slice_id] if real_wave is None else apply_phase_grating(real_wave, gratings[slice_id])
        reciprocal = real_to_reciprocal_wave(transmitted, ordering=reciprocal_ordering)
        sampled = apply_propagator(reciprocal, kernels[slice_id], wave_ordering=reciprocal_ordering,
                                   propagator_ordering=propagator_ordering)
        if occurrence == len(sequence)-1:
            return MultisliceResult(sampled, beams, sequence, len(gratings), reciprocal_ordering)
        real_wave = reciprocal_to_real_wave(sampled, ordering=reciprocal_ordering)
    raise AssertionError("validated nonempty sequence did not produce a result")
