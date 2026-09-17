"""Independent plane-wave tilt experiments, serial or chunked NumPy FFT."""

from dataclasses import dataclass, field
from typing import Literal

import numpy as np
from numpy.typing import NDArray

from .fft_backend import (
    _batch_array, batched_wave_real_to_reciprocal,
    batched_wave_reciprocal_to_real,
)
from .multislice import BeamIndices, _beams, _integer, _snapshot, beam_intensities
from .propagator import propagation_boundary_mask, propagation_kernel
from .simulation import PreparedMultisliceSample, _geometry_parameters, simulate_one_tilt


def _tilts(values: NDArray[np.float64]) -> NDArray[np.float64]:
    if not isinstance(values, np.ndarray) or values.dtype != np.float64:
        raise TypeError("tilts must be a float64 NumPy array")
    if values.ndim != 2 or values.shape[1] != 2 or values.shape[0] == 0:
        raise ValueError("tilts must have nonempty shape (n_tilts, 2), columns H,K")
    if not np.all(np.isfinite(values)):
        raise ValueError("tilts must be finite")
    return values


def _batch_size(value: int) -> int:
    value = _integer(value, "batch_size")
    if value <= 0:
        raise ValueError("batch_size must be positive")
    return value


@dataclass(frozen=True, eq=False)
class ParallelExecutionMetadata:
    """Pool limits and distinct worker processes that completed chunk tasks."""

    requested_workers: int
    worker_limit: int
    actual_workers: int
    process_chunk_size: int
    start_method: Literal["spawn"] = "spawn"

    def __post_init__(self) -> None:
        for name in ("requested_workers", "worker_limit", "actual_workers", "process_chunk_size"):
            value = _integer(getattr(self, name), name)
            if value <= 0:
                raise ValueError(f"{name} must be positive")
            object.__setattr__(self, name, value)
        if not self.actual_workers <= self.worker_limit <= self.requested_workers:
            raise ValueError("require actual_workers <= worker_limit <= requested_workers")
        if not isinstance(self.start_method, str) or self.start_method != "spawn":
            raise ValueError("start_method must be 'spawn'")


@dataclass(frozen=True, eq=False)
class TiltSeriesResult:
    """Immutable snapshots in request order, with no implicit disk reshape.

    Tilts are undilated fractional reciprocal-cell coefficients, H then K.
    Amplitudes/intensities have shape (n_tilts, n_beams), including duplicates.
    batch_size records the requested maximum, not the final chunk length.
    """

    tilts: NDArray[np.float64]
    amplitudes: NDArray[np.complex128]
    beam_indices_hk: BeamIndices
    execution_mode: Literal["serial", "batched", "parallel"]
    batch_size: int | None = None
    parallel: ParallelExecutionMetadata | None = None
    intensities: NDArray[np.float64] = field(init=False)

    def __post_init__(self) -> None:
        tilts = _tilts(self.tilts)
        beams = _beams(self.beam_indices_hk)
        values = self.amplitudes
        if not isinstance(values, np.ndarray) or values.dtype != np.complex128:
            raise TypeError("amplitudes must be a complex128 NumPy array")
        if values.shape != (len(tilts), len(beams)):
            raise ValueError("amplitudes must have shape (n_tilts, n_beams)")
        intensities = beam_intensities(values.reshape(-1)).reshape(values.shape)
        if not isinstance(self.execution_mode, str):
            raise ValueError("execution_mode must be 'serial', 'batched' or 'parallel'")
        if self.execution_mode == "serial":
            if self.batch_size is not None:
                raise ValueError("serial execution requires batch_size=None")
        elif self.execution_mode in ("batched", "parallel"):
            object.__setattr__(self, "batch_size", _batch_size(self.batch_size))
        else:
            raise ValueError("execution_mode must be 'serial', 'batched' or 'parallel'")
        if self.execution_mode == "parallel":
            if not isinstance(self.parallel, ParallelExecutionMetadata):
                raise TypeError("parallel execution requires ParallelExecutionMetadata")
        elif self.parallel is not None:
            raise ValueError("parallel metadata requires execution_mode='parallel'")
        object.__setattr__(self, "tilts", _snapshot(tilts))
        object.__setattr__(self, "amplitudes", _snapshot(values))
        object.__setattr__(self, "intensities", _snapshot(intensities))
        object.__setattr__(self, "beam_indices_hk", beams)


def _inputs(prepared: PreparedMultisliceSample, tilts: NDArray[np.float64]) -> NDArray[np.float64]:
    if not isinstance(prepared, PreparedMultisliceSample):
        raise TypeError("prepared_sample must be PreparedMultisliceSample")
    return _tilts(tilts)


def simulate_tilt_series_serial(
    prepared_sample: PreparedMultisliceSample, tilts: NDArray[np.float64],
) -> TiltSeriesResult:
    """Call simulate_one_tilt independently for each explicit float64 H,K row."""
    values = _inputs(prepared_sample, tilts)
    beams = prepared_sample.parameters.beam_indices_hk
    amplitudes = np.empty((len(values), len(beams)), dtype=np.complex128)
    for i, (h, k) in enumerate(values):
        amplitudes[i] = simulate_one_tilt(prepared_sample,
            tilt_h_undilated_fractional=h, tilt_k_undilated_fractional=k).amplitudes
    return TiltSeriesResult(values, amplitudes, beams, "serial")


def _chunk(prepared, tilts, mask, permutation, beam_rows, beam_cols):
    """Working waves and kernels live only for this chunk; bins are unshifted."""
    count = len(tilts)
    shape = (count, prepared.mesh.mesh_y, prepared.mesh.mesh_x)
    kernels = []
    for stored in prepared.slice_types:
        batch = np.empty(shape, dtype=np.complex128)
        for i, (h, k) in enumerate(tilts):
            dilation = stored.reciprocal_grid.dilation
            parameters = _geometry_parameters(stored.geometry, prepared.wavelength_angstrom,
                tilt_h_fractional=h/dilation, tilt_k_fractional=k/dilation)
            kernel = propagation_kernel(stored.reciprocal_grid, parameters,
                ordering=prepared.parameters.reciprocal_ordering, boundary="none")
            batch[i] = kernel if permutation is None else kernel[permutation]
        batch[:, ~mask] = 0j
        kernels.append(batch)
    sequence = prepared.parameters.slice_sequence_zero_based
    # Unit illumination times the first grating, independently broadcast.
    wave = np.broadcast_to(prepared.slice_types[sequence[0]].phase_grating, shape)
    for occurrence, index in enumerate(sequence):
        if occurrence:
            wave = wave * prepared.slice_types[index].phase_grating
        reciprocal = batched_wave_real_to_reciprocal(wave)
        reciprocal *= kernels[index]
        _batch_array(reciprocal, "propagated reciprocal waves")
        if occurrence == len(sequence)-1:
            return reciprocal[:, beam_rows, beam_cols].copy()
        wave = batched_wave_reciprocal_to_real(reciprocal)
    raise AssertionError("prepared sequence is empty")


def simulate_tilt_series_batched(
    prepared_sample: PreparedMultisliceSample, tilts: NDArray[np.float64], *, batch_size: int,
) -> TiltSeriesResult:
    """Chunk independent experiments; FFT axes are strictly (-2,-1).

    Accepts either prepared reciprocal ordering. Kernels and beam locations
    are permuted to unshifted FFT bins without touching the batch axis.
    No final inverse transform, RMT factor, or coherent combination of tilts.
    """
    values = _inputs(prepared_sample, tilts)
    size = _batch_size(batch_size)
    prepared = prepared_sample
    beams = prepared.parameters.beam_indices_hk
    amplitudes = np.empty((len(values), len(beams)), dtype=np.complex128)
    mask = propagation_boundary_mask(prepared.slice_types[0].reciprocal_grid)
    positions = np.asarray(prepared.beam_positions_row_col, dtype=np.intp).reshape(-1, 2)
    permutation = None
    if prepared.parameters.reciprocal_ordering == "centered":
        ny, nx = mask.shape
        permutation = np.ix_((np.arange(ny)+ny//2) % ny, (np.arange(nx)+nx//2) % nx)
        mask = mask[permutation]
        positions = (positions-np.array([ny//2, nx//2])) % np.array([ny, nx])
    for start in range(0, len(values), size):
        stop = min(start+size, len(values))
        amplitudes[start:stop] = _chunk(prepared, values[start:stop], mask, permutation,
                                       positions[:, 0], positions[:, 1])
    return TiltSeriesResult(values, amplitudes, beams, "batched", size)
