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


def test_modern_export_reports_unsupported_without_resolving_structure(tmp_path):
    from examples.b7a_refinement_layout import real_config
    from qcbed.gui.compatibility import compatibility_report,export_document
    from qcbed.gui.document import EditorDocument
    from dataclasses import replace
    import pytest
    config=real_config()
    config=replace(config,crystallography=replace(config.crystallography,
        structure=replace(config.crystallography.structure,path=str(tmp_path/'absent.cif'))))
    report=compatibility_report(config,base_directory=tmp_path)
    assert report['errors'][0].code=='crystallography_export_unsupported'
    with pytest.raises(ValueError,match='provenance is not representable'):
        export_document(EditorDocument(config,tmp_path/'config.yml'),tmp_path/'FOR007.dat')
    assert not (tmp_path/'FOR007.dat').exists()


def test_compact_initialization_stops_before_optimizer(tmp_path,monkeypatch):
    import importlib.util
    from test_modern_grouped_run import modern_config
    spec=importlib.util.spec_from_file_location('packaging_modern_smoke',ROOT/'packaging/macos/modern_smoke.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    def forbidden(*args,**kwargs):raise AssertionError('Optimizer/evaluation must not execute')
    monkeypatch.setattr('qcbed.anneal_local_refinement.anneal',forbidden)
    monkeypatch.setattr('qcbed.anneal_local_refinement.amoeba',forbidden)
    monkeypatch.setattr('qcbed.grouped_run.evaluate_grouped',forbidden)
    result=module.initialize_only(modern_config(tmp_path),tmp_path)
    assert result['evaluations']==0 and result['optimizer_entered'] is False
