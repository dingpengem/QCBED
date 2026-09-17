"""One grating/propagator step using the readable source's wave Fourier signs."""

from dataclasses import dataclass

import numpy as np

from .fft_backend import (
    ComplexGrid, Ordering, _array, _ordering,
    wave_real_to_reciprocal, wave_reciprocal_to_real,
)


def _product(left: ComplexGrid, right: ComplexGrid, left_name: str, right_name: str) -> ComplexGrid:
    left = _array(left, left_name)
    right = _array(right, right_name)
    if left.shape != right.shape:
        raise ValueError(f"{right_name} shape must match {left_name}; broadcasting is not allowed")
    with np.errstate(over="ignore", invalid="ignore"):
        result = left * right
    return _array(result, f"{left_name} product").copy(order="C")


def _matching_orderings(wave_ordering: Ordering, propagator_ordering: Ordering) -> None:
    _ordering(wave_ordering)
    _ordering(propagator_ordering)
    if wave_ordering != propagator_ordering:
        raise ValueError("wave ordering and propagator ordering must match; reorder explicitly")


def apply_phase_grating(real_wave: ComplexGrid, phase_grating: ComplexGrid) -> ComplexGrid:
    """Multiply real-space samples without shifting, scaling or modifying inputs."""
    return _product(real_wave, phase_grating, "real_wave", "phase_grating")


def real_to_reciprocal_wave(real_wave: ComplexGrid, *, ordering: Ordering) -> ComplexGrid:
    """Source negative-exponent transform / (Nx*Ny); ordering labels bins."""
    return wave_real_to_reciprocal(real_wave, ordering=ordering)


def apply_propagator(
    reciprocal_wave: ComplexGrid, propagator: ComplexGrid, *,
    wave_ordering: Ordering, propagator_ordering: Ordering,
) -> ComplexGrid:
    """Multiply matching reciprocal bins, including supplied zeros; no extra RMT."""
    _matching_orderings(wave_ordering, propagator_ordering)
    return _product(reciprocal_wave, propagator, "reciprocal_wave", "propagator")


def reciprocal_to_real_wave(reciprocal_wave: ComplexGrid, *, ordering: Ordering) -> ComplexGrid:
    """Source unnormalized positive-exponent sum; no real-space shift."""
    return wave_reciprocal_to_real(reciprocal_wave, ordering=ordering)


@dataclass(frozen=True)
class WaveStepResult:
    """Immutable snapshots after propagation, before the next phase grating.

    reciprocal_wave is the post-propagator beam-sampling intermediate in
    ordering; real_wave is its real-space reconstruction. Arrays are copied
    to immutable byte storage rather than exposing reusable work buffers.
    """

    real_wave: ComplexGrid
    reciprocal_wave: ComplexGrid
    ordering: Ordering

    def __post_init__(self) -> None:
        _ordering(self.ordering)
        real = _array(self.real_wave, "real_wave")
        reciprocal = _array(self.reciprocal_wave, "reciprocal_wave")
        if real.shape != reciprocal.shape:
            raise ValueError("real_wave and reciprocal_wave shape must match")
        for name, value in (("real_wave", real), ("reciprocal_wave", reciprocal)):
            snapshot = np.frombuffer(value.tobytes(order="C"), dtype=np.complex128).reshape(value.shape)
            object.__setattr__(self, name, snapshot)


def propagate_one_step(
    real_wave: ComplexGrid, phase_grating: ComplexGrid, propagator: ComplexGrid, *,
    ordering: Ordering, propagator_ordering: Ordering,
) -> WaveStepResult:
    """Grating -> normalized reciprocal transform -> propagator -> real sum.

    ordering selects the reciprocal intermediate and must equal the explicit
    propagator_ordering. Inputs have identical (mesh_y, mesh_x) shapes. The
    propagator supplies any boundary zeros. No extra RMT, next grating,
    sequence iteration or historical CDFT equivalence is implied.
    """
    _matching_orderings(ordering, propagator_ordering)
    transmitted = apply_phase_grating(real_wave, phase_grating)
    reciprocal = real_to_reciprocal_wave(transmitted, ordering=ordering)
    sampled = apply_propagator(reciprocal, propagator, wave_ordering=ordering,
                               propagator_ordering=propagator_ordering)
    propagated = reciprocal_to_real_wave(sampled, ordering=ordering)
    return WaveStepResult(propagated, sampled, ordering)
