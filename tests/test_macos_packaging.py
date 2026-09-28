"""Packaging is additive; frozen bootstrap dispatches multiprocessing first."""
from pathlib import Path
import runpy
import sys
import types
import tomllib

ROOT = Path(__file__).resolve().parents[1]


def test_launcher_freezes_before_gui_import(monkeypatch):
    import multiprocessing
    calls = []
    monkeypatch.setattr(multiprocessing,'freeze_support',lambda:calls.append('freeze'))
    module = types.ModuleType('qcbed.gui.__main__')
    def main():
        assert calls == ['freeze']
        calls.append('gui')
        return 0
    module.main = main
    monkeypatch.setitem(sys.modules,'qcbed.gui.__main__',module)
    monkeypatch.setattr(sys,'argv',['QCBED'])
    import pytest
    with pytest.raises(SystemExit) as result:
        runpy.run_path(str(ROOT/'packaging/macos/launcher.py'),run_name='__main__')
    assert result.value.code == 0
    assert calls == ['freeze','gui']


def test_packaging_dependencies_remain_optional():
    metadata = tomllib.loads((ROOT/'pyproject.toml').read_text())['project']
    assert not any('pyinstaller' in d.lower() or 'pyside' in d.lower() for d in metadata['dependencies'])
    assert any(d.startswith('PyInstaller==') for d in metadata['optional-dependencies']['package'])
    assert (ROOT/'src/qcbed/gui/__main__.py').exists()
    assert (ROOT/'src/qcbed/run.py').exists()
