import pytest

from qcbed import (
    MeshConfig,
    MicroscopeConfig,
    SimulationConfig,
    SliceGeometry,
    TiltCoordinate,
)


def test_models_validate_confirmed_conventions() -> None:
    mesh = MeshConfig(mesh_x=4, mesh_y=3)
    microscope = MicroscopeConfig(energy_kv=200.0)
    slices = (SliceGeometry(thickness_angstrom=1.5),)
    tilts = (TiltCoordinate(first=0.0, second=0.25),)

    config = SimulationConfig(mesh, microscope, slices, tilts)

    assert config.mesh == mesh
    assert config.microscope.energy_kv == 200.0
    assert config.slices[0].thickness_angstrom == 1.5
    assert config.tilts[0].second == 0.25


@pytest.mark.parametrize("mesh_x, mesh_y", [(0, 3), (4, 0), (-1, 2)])
def test_mesh_config_rejects_invalid_dimensions(mesh_x: int, mesh_y: int) -> None:
    with pytest.raises(ValueError):
        MeshConfig(mesh_x, mesh_y)


def test_dataclasses_reject_invalid_values() -> None:
    with pytest.raises(ValueError):
        MicroscopeConfig(energy_kv=0.0)
    with pytest.raises(ValueError):
        SliceGeometry(thickness_angstrom=-1.0)
    with pytest.raises(ValueError):
        SimulationConfig(MeshConfig(2, 2), MicroscopeConfig(200.0), ())
