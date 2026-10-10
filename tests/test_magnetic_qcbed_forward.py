"""Research-level adapter tests; no modification of baseline simulation."""
import numpy as np
import pytest

from qcbed.magnetism import add_magnetic_gratings, magnetic_phase_grating, project_az_midpoint
from qcbed.magnetism.research import e, hbar


def test_zero_Az_leaves_existing_complex_gratings_bitwise():
    g = np.full((10, 8), 0.9 + 0.01j, dtype=np.complex128)
    combined = add_magnetic_gratings((g,), (np.zeros(g.shape),))[0]
    assert np.array_equal(g, combined)
    assert not np.shares_memory(g, combined)


def test_simple_magnetic_phase_sign():
    phase = np.full((4, 4), np.pi * hbar / e)
    grating = magnetic_phase_grating(phase)
    assert np.allclose(grating, -1 + 0j, atol=1e-14)


def test_projected_midpoint_is_weber_not_tesla():
    a = np.ones((2, 6, 6), dtype=float) * 2e-12  # T*m
    values = project_az_midpoint(a, [3e-10, 4e-10])
    assert np.allclose(values[0], 6e-22)
    assert np.allclose(values[1], 8e-22)


def test_invalid_nonphysical_shapes_and_types_rejected():
    with pytest.raises(ValueError):
        add_magnetic_gratings((np.ones((5, 5), dtype=np.complex128),), (np.ones((4, 5)),))
    with pytest.raises(ValueError):
        add_magnetic_gratings((np.ones((5, 5)),), (np.zeros((5, 5)),))
