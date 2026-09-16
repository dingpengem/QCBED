"""NumPy CPU transforms of the Python Fourier specification, not legacy CDFT."""

from typing import Literal

import numpy as np
from numpy.typing import NDArray

from .grids import reciprocal_index_grid


ComplexGrid = NDArray[np.complex128]
Ordering = Literal["centered", "unshifted"]


def _array(value: ComplexGrid, name: str) -> ComplexGrid:
    if not isinstance(value, np.ndarray) or value.dtype != np.complex128:
        raise TypeError(f"{name} must be a complex128 NumPy array")
    if value.ndim != 2 or 0 in value.shape:
        raise ValueError(f"{name} must have nonempty shape (mesh_y, mesh_x)")
    if not np.all(np.isfinite(value)):
        raise ValueError(f"{name} must be finite")
    return value


def _ordering(value: Ordering) -> None:
    if not isinstance(value, str) or value not in ("centered", "unshifted"):
        raise ValueError("ordering must be 'centered' or 'unshifted'")


def _centered_bins(shape: tuple[int, int]) -> tuple[np.ndarray, np.ndarray]:
    ny, nx = shape
    h, k = reciprocal_index_grid(nx, ny, centered=True)
    return k[:, 0].astype(np.intp) % ny, h[0, :].astype(np.intp) % nx


def centered_to_unshifted(coefficients: ComplexGrid) -> ComplexGrid:
    """Reorder centered reciprocal coefficients into FFT bins, without masking.

    Project-defined signed H/K map to (K mod Ny,H mod Nx). Even negative
    Nyquist values are preserved. Odd/even rectangular shapes are supported.
    Returns a new complex128 array; no Fourier transform or scaling occurs.
    """
    values = _array(coefficients, "coefficients")
    rows, cols = _centered_bins(values.shape)
    result = np.empty(values.shape, dtype=np.complex128)
    result[np.ix_(rows, cols)] = values
    return result


def unshifted_to_centered(coefficients: ComplexGrid) -> ComplexGrid:
    """Reorder unshifted FFT bins into project-centered reciprocal order.

    Exact inverse permutation of centered_to_unshifted, including odd sizes.
    Returns a new complex128 array without masking, transforming or scaling.
    """
    values = _array(coefficients, "coefficients")
    rows, cols = _centered_bins(values.shape)
    return values[np.ix_(rows, cols)].copy(order="C")


def reciprocal_to_real_sum(
    coefficients: ComplexGrid, *, ordering: Ordering = "unshifted",
) -> ComplexGrid:
    """Reciprocal coefficients -> unnormalized negative-exponent real-space sum.

    V[row,col]=sum C[H,K]*exp(-2*pi*i*(H*col/Nx+K*row/Ny)). Ordering
    describes the input only; output starts at fractional (0,0), with no
    spatial shift. Uses fft2(norm='backward'), with no masks or extra factors.
    """
    _ordering(ordering)
    values = _array(coefficients, "coefficients")
    bins = centered_to_unshifted(values) if ordering == "centered" else values
    with np.errstate(over="ignore", invalid="ignore"):
        result = np.fft.fft2(bins, axes=(0, 1), norm="backward")
    return _array(result, "real-space sum").copy(order="C")


def real_to_reciprocal_coefficients(
    real_space: ComplexGrid, *, ordering: Ordering = "unshifted",
) -> ComplexGrid:
    """Inverse Python transform: positive exponent, divided by Nx*Ny.

    Input is a complex real-space field sampled at (col/Nx,row/Ny), without
    spatial centering. Uses ifft2(norm='backward'); ordering selects the
    returned reciprocal coefficient order. This is not yet a validated
    legacy real-to-reciprocal wave-propagation operation.
    """
    _ordering(ordering)
    values = _array(real_space, "real_space")
    with np.errstate(over="ignore", invalid="ignore"):
        bins = np.fft.ifft2(values, axes=(0, 1), norm="backward")
    _array(bins, "reciprocal coefficients")
    return unshifted_to_centered(bins) if ordering == "centered" else bins.copy(order="C")


def _batch_array(value: ComplexGrid, name: str) -> ComplexGrid:
    if not isinstance(value, np.ndarray) or value.dtype != np.complex128:
        raise TypeError(f"{name} must be a complex128 NumPy array")
    if value.ndim != 3 or 0 in value.shape:
        raise ValueError(f"{name} must have nonempty shape (batch, mesh_y, mesh_x)")
    if not np.all(np.isfinite(value)):
        raise ValueError(f"{name} must be finite")
    return value


def batched_reciprocal_to_real_sum(coefficients: ComplexGrid) -> ComplexGrid:
    """Independent unshifted coefficient sums; negative exponent, no scaling.

    Only the final two axes are transformed. The leading axis labels separate
    experiments and is never transformed or reduced. No ordering permutation.
    """
    values = _batch_array(coefficients, "coefficients")
    with np.errstate(over="ignore", invalid="ignore"):
        result = np.fft.fft2(values, axes=(-2, -1), norm="backward")
    return _batch_array(result, "real-space sums")


def batched_real_to_reciprocal_coefficients(real_space: ComplexGrid) -> ComplexGrid:
    """Independent positive-exponent transforms divided by mesh_x*mesh_y.

    Returns unshifted bins with shape (batch, mesh_y, mesh_x). There is no
    division by batch size and no transform along the batch axis.
    """
    values = _batch_array(real_space, "real_space")
    with np.errstate(over="ignore", invalid="ignore"):
        result = np.fft.ifft2(values, axes=(-2, -1), norm="backward")
    return _batch_array(result, "reciprocal coefficients")
