"""Decimal-table references, not executable Fortran equivalence tests."""

import math
from pathlib import Path
import re

import numpy as np
import pytest

from qcbed import elastic_atomic_scattering_factor
from qcbed import scattering


# Independently evaluated with exact rational Lagrange interpolation using
# explicit source stencils. No production cubic coefficient algebra was used.
# Stencils and table values are recorded in docs/scattering.md.
CASES = [
    (1, 0.0, 25.327466761),
    (1, 0.025, 24.824622984184895),
    (1, 0.375, 6.886204263203125),
    (1, 1.5, 0.526658099),
    (13, 0.0, 281.953595001),
    (13, 0.025, 275.340445007875),
    (13, 0.375, 48.5632630038125),
    (13, 6.0, 0.430902081),
    (98, 0.0, 806.313549569),
    (98, 0.025, 802.5855584),
    (98, 0.375, 246.79875153955555),
    (98, 6.0, 2.9426290975555554),
]


@pytest.mark.parametrize("z, stl, expected", CASES)
def test_independent_elastic_values(z, stl, expected):
    value = elastic_atomic_scattering_factor(z, stl)
    assert type(value) is float
    assert value == pytest.approx(expected, rel=2e-12, abs=1e-12)


def test_explicit_option_and_numpy_scalars():
    assert elastic_atomic_scattering_factor(np.int64(13), np.float64(.375), option=np.int64(2)) == pytest.approx(
        48.5632630038125, rel=2e-12, abs=1e-12
    )


@pytest.mark.parametrize("z", [0, -1, 99])
def test_atomic_number_bounds(z):
    with pytest.raises(ValueError, match="atomic_number must be in 1..98"):
        elastic_atomic_scattering_factor(z, 0.0)


@pytest.mark.parametrize("z", [True, False, np.bool_(True), 13.0, math.nan, math.inf, "13"])
def test_atomic_number_type(z):
    with pytest.raises(TypeError, match="atomic_number must be an integer"):
        elastic_atomic_scattering_factor(z, 0.0)


@pytest.mark.parametrize("stl", [-.1, math.nan, math.inf, -math.inf, 6.000001])
def test_scattering_coordinate_bounds(stl):
    with pytest.raises(ValueError, match="stl_inverse_angstrom must be finite and in 0..6.0"):
        elastic_atomic_scattering_factor(13, stl)


def test_hydrogen_coordinate_limit():
    with pytest.raises(ValueError, match="0..1.5"):
        elastic_atomic_scattering_factor(1, 1.500001)


@pytest.mark.parametrize("stl", [True, False, np.bool_(False), ".1", None, .1+0j, [.1]])
def test_scattering_coordinate_type(stl):
    with pytest.raises(TypeError, match="stl_inverse_angstrom must be a real scalar"):
        elastic_atomic_scattering_factor(13, stl)


@pytest.mark.parametrize("option", [-1, 0, 1, 3])
def test_unsupported_options(option):
    with pytest.raises(ValueError, match="option must be 2"):
        elastic_atomic_scattering_factor(13, 0.0, option=option)


@pytest.mark.parametrize("option", [True, False, np.bool_(True), 2.0, math.nan, math.inf, "2"])
def test_option_type(option):
    with pytest.raises(TypeError, match="option must be an integer"):
        elastic_atomic_scattering_factor(13, 0.0, option=option)


def test_all_table_entries_preserve_source_order_and_sentinels():
    # Verify all copied data against readable source, including unused X-ray
    # rows and negative sentinels. This is a provenance check, not execution.
    source = (Path(__file__).parents[1] / "legacy_fortran" / "fftsg.f").read_text()
    for source_name, table in [("EFARRAY", scattering._EF), ("XFARRAY", scattering._XF)]:
        rows = re.findall(rf"{source_name}\((\d+),:\)\s*=\s*\(/(.*?)/\)", source, re.S)
        assert len(rows) == len(table) == 98
        for z, body in rows:
            expected = tuple(float(v) for v in body.replace('*', '').split(','))
            assert len(expected) == 62
            assert table[int(z)-1] == expected
    for source_name, table in [("STLARRAY", scattering._STL), ("EMARRAY", scattering._METHOD)]:
        body = re.search(rf"{source_name}(?:\(:\))?\s*=\s*\(/(.*?)/\)", source, re.S)[1]
        assert table == tuple(float(v) for v in body.replace('*', '').split(','))
