import math

import pytest
import numpy as np

from qcbed import electron_wavelength_angstrom as public_wavelength
from qcbed.electron import electron_wavelength_angstrom


@pytest.mark.parametrize(
    "voltage_kv, expected_angstrom",
    [
        (100.0, 0.0370145195150564),
        (200.0, 0.0250794525129636),
        (300.0, 0.0196875824004297),
    ],
)
def test_wavelength_independently_checked_decimal_values(
    voltage_kv: float, expected_angstrom: float
) -> None:
    # Independently checked using 60-digit Decimal arithmetic on the printed
    # source constants. These are decimal-equation references, not executable
    # Fortran outputs. Keep expectations literal rather than recomputing here.
    assert electron_wavelength_angstrom(voltage_kv) == pytest.approx(
        expected_angstrom, rel=1e-14, abs=0.0
    )


@pytest.mark.parametrize("voltage_kv", [0.0, -0.0, -1.0, math.nan, math.inf, -math.inf])
def test_wavelength_rejects_invalid_voltage(voltage_kv: float) -> None:
    with pytest.raises(ValueError, match="voltage_kv must be finite and positive in kV"):
        electron_wavelength_angstrom(voltage_kv)


@pytest.mark.parametrize("value", [True, False, np.bool_(True), "200", "200 V", None, 200+0j, [200]])
def test_wavelength_rejects_bool_and_non_real_inputs(value) -> None:
    with pytest.raises(TypeError, match="voltage_kv must be a real numeric scalar in kV, not bool"):
        electron_wavelength_angstrom(value)


@pytest.mark.parametrize("value", [200, 200.0, np.float64(200)])
def test_public_api_returns_python_float(value) -> None:
    assert public_wavelength is electron_wavelength_angstrom
    assert type(public_wavelength(voltage_kv=value)) is float


def test_wavelength_decreases_with_voltage() -> None:
    values = [electron_wavelength_angstrom(v) for v in (60, 100, 200, 300, 1000)]
    assert all(left > right > 0 for left, right in zip(values, values[1:]))


def test_voltage_units_are_never_inferred_from_magnitude() -> None:
    # 200000 is 200000 kV here, never silently converted from volts to 200 kV.
    assert electron_wavelength_angstrom(200000.0) == pytest.approx(
        0.000061834814305985503, rel=1e-14, abs=0.0
    )
    assert electron_wavelength_angstrom(voltage_kv=200.0) == pytest.approx(
        0.0250794525129636, rel=1e-14, abs=0.0
    )
