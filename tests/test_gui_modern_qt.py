from dataclasses import replace
from pathlib import Path
import time
import pytest
pytest.importorskip('PySide6')
from PySide6.QtWidgets import QApplication
from PySide6.QtCore import Qt
from qcbed import save_yaml_config,load_yaml_config
from qcbed.gui.window import ConfigWindow
from qcbed.gui.run_controller import RunController,readiness
from qcbed.gui.document import EditorDocument
from qcbed.config_editing import Id
from test_modern_grouped_run import modern_config


@pytest.fixture(scope='module')
def app():return QApplication.instance() or QApplication([])


def test_modern_edit_save_reload_and_diagnostics(app,tmp_path):
    c=modern_config(tmp_path);path=tmp_path/'modern.yml'
    save_yaml_config(c,path,base_directory=tmp_path)
    w=ConfigWindow(path)
    assert w.site_group_table.rowCount()==1
    assert w.site_group_table.item(0,0).text()=='Al_fcc'
    assert w.refinement_table.rowCount()==0
    from PySide6.QtWidgets import QLabel
    assert any('strict_b5_compatible: True' in label.text() for label in w.atomic_expert.findChildren(QLabel))
    assert not w.atomic_expert.toggle.isChecked()
    w.site_group_table.item(0,3).setText('.9')
    w.apply_site_groups.click()
    w.widgets[('crystallography','depth_angstrom')].setText('5')
    w.widgets[('experiment','gradient_option')].setCurrentText('Conventional')
    assert w.apply_edits()
    assert w.document.config.crystallography.site_groups[0].debye_waller_b==.9
    assert w.document.config.experiment.gradient_option==0
    w.widgets[('crystallography','cbed_basis','g1_hkl')].setText('0,0,0')
    assert not w.apply_edits()
    assert w.issues.topLevelItemCount()>0
    assert w.document.config.crystallography.cbed_basis.g1_hkl==(2,0,0)
    w.refresh()
    w.workflow.accept(0,w.document.config);assert w.validate_stage(1)
    assert w.validate_c_part('simulation')
    row=next(r for r in range(w.refinement_table.rowCount()) if w.refinement_table.item(r,0).data(Qt.ItemDataRole.UserRole).category.value=='depth')
    w.refinement_table.item(row,8).setText('2')
    assert w.apply_edits()
    new=tmp_path/'saved.yml';assert w.save_document(new)
    assert load_yaml_config(new)==w.document.config
    assert w.document.config.crystallography.site_groups[0].site_indices==(0,1,2,3)
    w.document.saved_config=w.document.config;w.close()


def test_modern_run_cancel_results_resume_qt(app,tmp_path):
    c=modern_config(tmp_path);path=tmp_path/'modern.yml';save_yaml_config(c,path,base_directory=tmp_path)
    output=Path(c.run.output_directory);output.mkdir()
    preserved={'notes.txt':b'user notes','final_out.dat':b'previous result'}
    for name,data in preserved.items():(output/name).write_bytes(data)
    w=ConfigWindow(path)
    for stage in range(3):w.workflow.accept(stage,w.document.config)
    w.validate_document()
    assert w.validate_output_directory();w.check_run_readiness()
    assert w.run_action.isEnabled() and not w.resume_action.isEnabled()
    w.workflow_tabs.setCurrentIndex(3)
    assert w.run_dashboard.run.isEnabled() and not w.run_dashboard.resume.isEnabled()
    w.run_dashboard.run.click()
    assert w.workflow_tabs.currentIndex()==3 and w.workflow_tabs.isTabEnabled(3)
    assert not w.workflow_tabs.isTabEnabled(2)
    assert w.run_dashboard.cancel.isEnabled() and not w.run_dashboard.run.isEnabled()
    deadline=time.monotonic()+60
    cancelled=False
    while w.controller.active:
        app.processEvents();w.poll_run()
        if not cancelled and w.results_view.results is not None and w.results_view.results.summary['latest_valid_evaluation']:
            w.run_dashboard.cancel.click();cancelled=True
            assert not w.run_dashboard.cancel.isEnabled()
        assert time.monotonic()<deadline
        time.sleep(.005)
    assert cancelled
    # A tiny run may finish before the cancellation reaches the worker; both are safe terminal outcomes.
    assert w.controller.state in ('Cancelled','Completed')
    assert w.run_dashboard.results.isEnabled()
    w.run_dashboard.results.click()
    assert w.results_actions['Show Results'].isEnabled()
    assert w.results_view.results.summary['latest_valid_evaluation'] is not None
    assert w.results_view.results.comparison is not None
    assert w.resume_action.isEnabled()
    assert w.run_dashboard.resume.isEnabled()
    w.run_dashboard.resume.click()
    while w.controller.active:
        app.processEvents();w.poll_run();time.sleep(.005)
        assert time.monotonic()<deadline
    assert w.controller.state=='Completed'
    assert w.results_view.results.summary['retained_final']
    assert w.results_view.modern_history.rowCount()>0
    assert w.results_view.beam_details.rowCount()==9
    assert all((output/name).read_bytes()==data for name,data in preserved.items())
    assert len(list(output.glob('qcbed-run-*')))==1
    assert w.results_directory.parent==output
    w.document.saved_config=w.document.config;w.close()


def test_structure_browse_cif_and_failed_read_preserves_model(app, tmp_path, monkeypatch):
    from pathlib import Path
    import os
    import numpy as np
    from qcbed.crystallography import ase_adapter
    from PySide6.QtWidgets import QFileDialog
    config = modern_config(tmp_path)
    path = tmp_path / 'modern.yml'
    save_yaml_config(config, path, base_directory=tmp_path)
    cif = Path(__file__).resolve().parents[1] / 'reference_data' / 'Al_FCC.cif'
    converted = []
    original = ase_adapter.atomic_structure_from_ase
    def convert(atoms):
        result = original(atoms)
        converted.append(result)
        return result
    monkeypatch.setattr(ase_adapter, 'atomic_structure_from_ase', convert)
    window = ConfigWindow(path)
    key = ('crystallography', 'structure', 'path')
    assert not window.widgets[key].isReadOnly()
    def choose(*args):
        assert '*.cif' in args[3] and 'POSCAR' in args[3]
        return str(cif), ''
    monkeypatch.setattr(QFileDialog, 'getOpenFileName', choose)
    window.browse(window.widgets[key], 'file')
    assert window.document.config.crystallography.structure.path == str(cif)
    assert window.document.config.crystallography.structure.format is None
    assert converted
    np.testing.assert_allclose(converted[-1].direct_lattice.lengths, [4.03893]*3)
    assert len(converted[-1].sites) == 4
    assert {site.symbol for site in converted[-1].sites} == {'Al'}
    groups=window.document.config.crystallography.site_groups
    assert len(groups)==1 and groups[0].element=='Al'
    assert groups[0].site_indices==(0,1,2,3) and groups[0].debye_waller_b is None
    assert window.site_group_table.item(0,3).text()==''
    assert groups[0].id=='Al_fcc'
    window.site_group_table.item(0,3).setText('0.86')
    assert window.document.config.crystallography.site_groups[0].debye_waller_b is None
    window.apply_site_groups.click()
    assert window.document.config.crystallography.site_groups[0].debye_waller_b==.86
    window.workflow.accept(0,window.document.config)
    window.preparation.dirty=False
    assert window.validate_stage(1),window.stage_status[1].text()
    assert window.save_document(tmp_path/'applied.yml')
    assert window.save_document(tmp_path/'applied-copy.yml')
    assert load_yaml_config(tmp_path/'applied-copy.yml').crystallography.site_groups[0].debye_waller_b==.86
    assert not window.document.config.crystallography.debye_waller
    assert not window.document.config.refinement.parameters
    assert window.document.summary()['compile'].startswith('Valid')
    window.document.config = replace(window.document.config,
        crystallography=replace(window.document.config.crystallography,
            site_groups=config.crystallography.site_groups))
    window.refresh()
    assert '4.03893' in window.modern_summary.text()
    assert 'atom_count: 4' in window.modern_summary.text()
    assert "species: ['Al']" in window.modern_summary.text()
    assert window.document.summary()['compile'] == 'Valid (committed snapshot)'
    assert 'Compile validity: Valid' in window.summary.toPlainText()
    assert len(window.atomic_axes) == 3
    assert [ax.get_title() for ax in window.atomic_axes] == [
        'Front / Beam — along +z_beam', 'Top — along +y_beam', 'Left — along +x_beam']
    assert all(len(ax.collections[0].get_offsets()) == 32 for ax in window.atomic_axes)
    assert window.site_group_table.rowCount() == 1
    assert window.site_group_table.item(0, 0).text() == 'Al_fcc'
    assert window.site_group_table.item(0, 2).text() == '0, 1, 2, 3'
    assert window.document.config.crystallography.site_groups == config.crystallography.site_groups
    from qcbed.crystallography import CrystalOrientation, ZoneAxis
    from qcbed.crystallography.projection import project_atomic_structure
    def check_views(axis):
        structure = converted[-1]
        oriented = project_atomic_structure(structure,
            CrystalOrientation(structure.direct_lattice, ZoneAxis(axis)), repeats=(2,2,2))
        xyz = np.array([(v.x,v.y,v.depth) for v in oriented])
        for ax, indices in zip(window.atomic_axes, ((0,1),(0,2),(2,1))):
            np.testing.assert_allclose(ax.collections[0].get_offsets(), xyz[:,indices])
            assert ax.get_aspect() == 1
            assert len(ax.lines) == 12
        np.testing.assert_allclose(np.array(window.atomic_3d_axes.collections[0]._offsets3d).T, xyz)
        return xyz
    initial = check_views((0,0,1))
    assert np.ptp(initial[:,0]) == pytest.approx(np.ptp(initial[:,1]))
    # Visualization remains available even if the changed axis requires new indexing.
    window.document.config = replace(window.document.config,
        crystallography=replace(window.document.config.crystallography, zone_axis=(1,1,0)))
    window.refresh()
    assert not np.allclose(check_views((1,1,0)), initial)
    before = window.document.config
    errors = []
    monkeypatch.setattr(window, 'report_error', errors.append)
    bad = tmp_path / 'bad.cif'
    bad.write_text('not a structure')
    assert not window.load_structure_path(str(bad))
    assert errors and 'previous model is unchanged' in errors[0]
    assert window.document.config == before
    # Typed relative paths also load through ASE and resolve against the document.
    window.widgets[key].setText(os.path.relpath(cif, tmp_path))
    window.widgets[key].editingFinished.emit()
    assert window.document.config == before
    window.document.saved_config = window.document.config
    window.pending = False
    window.close()
