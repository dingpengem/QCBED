from qcbed import MeshConfig


def test_package_import() -> None:
    assert MeshConfig(2, 3).mesh_x == 2
