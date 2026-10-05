"""Focused modern state transitions; no refinement campaign."""
from dataclasses import asdict,replace
import json
from pathlib import Path
import pytest
pytest.importorskip('PySide6')
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication,QLineEdit,QComboBox
from qcbed import save_yaml_config,load_yaml_config,available_refinement_targets
from qcbed.configuration import CategoryBound,RelativeBound,ParameterCategory as Category
from qcbed.gui.document import EditorDocument
from qcbed.gui.window import ConfigWindow
from test_modern_grouped_run import modern_config


@pytest.fixture
def window(tmp_path):
    app=QApplication.instance() or QApplication([])
    w=ConfigWindow();w.document=EditorDocument(modern_config(tmp_path),tmp_path/'project.yml');w.refresh()
    yield w
    w.document.saved_config=w.document.config;w.pending=False;w.site_group_pending=False;w.close();w.deleteLater();app.processEvents()


def setup(w):
    w.workflow.accept(0,w.document.config)
    assert w.validate_stage(1),w.stage_status[1].text()
    w.workflow_tabs.setCurrentIndex(2)
    assert w.validate_c_part('simulation'),w.stage_status[2].text()


def test_fresh_thickness_and_loaded_value(window,tmp_path):
    from qcbed.gui.preparation_model import new_preparation_document
    c=window.document.config
    descriptor=dict(format='qcbed-prepared-input-v1',scan=asdict(c.scan),experiment=asdict(c.experiment),
        beams=[dict(id=b.id,lattice_coordinates=b.lattice_coordinates) for b in c.crystallography.beams],scan_basis_transform=[[1,0],[0,1]])
    path=tmp_path/'prepared.json';path.write_text(json.dumps(descriptor))
    fresh=new_preparation_document(path)
    assert fresh.crystallography.depth_angstrom==1000.
    window.document=EditorDocument(fresh,tmp_path/'fresh.yml');window.refresh()
    assert window.widgets[('crystallography','depth_angstrom')].text()=='1000.0'
    # Loaded physical values do not inherit the fresh-project default.
    save_yaml_config(c,tmp_path/'loaded.yml',base_directory=tmp_path)
    window.open_document(tmp_path/'loaded.yml')
    assert window.document.config.crystallography.depth_angstrom==4.
    window.widgets[('crystallography','depth_angstrom')].setText('1000')
    setup(window)
    row=next(r for r in range(window.refinement_table.rowCount()) if window.refinement_table.item(r,0).data(Qt.ItemDataRole.UserRole).category==Category.DEPTH)
    assert window.refinement_table.item(row,3).text()=='1000.0'


def test_basis_survives_site_apply_validation_and_save(window,tmp_path):
    w=window;c=w.document.config
    basis=replace(c.crystallography.cbed_basis,g1_hkl=(1,0,0),g2_hkl=(0,1,0))
    w.document=EditorDocument(replace(c,crystallography=replace(c.crystallography,cbed_basis=basis)),tmp_path/'basis.yml');w.refresh()
    w.workflow.accept(0,w.document.config);w.workflow_tabs.setCurrentIndex(1)
    w.widgets[('crystallography','cbed_basis','g1_hkl')].setText('2,0,0')
    w.widgets[('crystallography','cbed_basis','g2_hkl')].setText('0,2,0')
    w.widgets[('crystallography','depth_angstrom')].setText('12')
    w.site_group_table.item(0,3).setText('.9');w.apply_site_groups.click()
    assert w.document.config.crystallography.cbed_basis.g1_hkl==(2,0,0)
    assert w.document.config.crystallography.cbed_basis.g2_hkl==(0,2,0)
    assert w.widgets[('crystallography','depth_angstrom')].text()=='12'
    assert w.validate_stage(1),w.stage_status[1].text()
    # B validation must preserve the independent pending C thickness edit too.
    assert w.widgets[('crystallography','depth_angstrom')].text()=='12'
    w.workflow_tabs.setCurrentIndex(2);assert w.validate_c_part('simulation')
    assert w.document.config.crystallography.depth_angstrom==12.
    assert w.validate_stage(2),w.stage_status[2].text()
    for index in (0,1,2,3,1):w.workflow_tabs.setCurrentIndex(index)
    assert w.save_document(tmp_path/'saved.yml')
    loaded=load_yaml_config(tmp_path/'saved.yml')
    assert loaded.crystallography.cbed_basis.g1_hkl==(2,0,0)
    assert loaded.crystallography.cbed_basis.g2_hkl==(0,2,0)
    assert loaded.crystallography.site_groups[0].debye_waller_b==.9


def test_configured_bounds_visible_before_selection_without_yaml_bloat(window,tmp_path):
    w=window;c=w.document.config
    rules=replace(c.refinement.bounds,specific=(),absolute=(CategoryBound(Category.DEPTH,1.,20.),CategoryBound(Category.DEBYE_WALLER,.5,1.5)),relative=(RelativeBound(Category.SCAN_COORDINATE,.2,.05),))
    w.document=EditorDocument(replace(c,refinement=replace(c.refinement,parameters=(),bounds=rules)),tmp_path/'rules.yml');w.refresh();setup(w)
    table=w.refinement_table
    for row in range(table.rowCount()):
        target=table.item(row,0).data(Qt.ItemDataRole.UserRole)
        if target.category==Category.DEPTH:
            assert (table.item(row,8).text(),table.item(row,9).text())==('1.0','20.0')
            table.item(row,1).setCheckState(Qt.CheckState.Checked)
        elif target.category==Category.DEBYE_WALLER:
            assert (table.item(row,8).text(),table.item(row,9).text())==('0.5','1.5')
        elif target.category==Category.SCAN_COORDINATE:
            start=float(table.item(row,3).text());width=max(abs(start)*.2,.05)
            assert float(table.item(row,8).text())==pytest.approx(start-width)
        elif target.category==Category.VOLTAGE:
            assert table.item(row,8).text()==''
            assert 'No canonical' in table.item(row,8).toolTip()
    assert w.validate_c_part('refinement'),w.stage_status[2].text()
    assert not w.document.config.refinement.bounds.specific
    assert w.save_document(tmp_path/'rules-saved.yml')
    assert load_yaml_config(tmp_path/'rules-saved.yml')==w.document.config
    table=w.refinement_table
    row=next(r for r in range(table.rowCount()) if table.item(r,0).data(Qt.ItemDataRole.UserRole).category==Category.DEPTH)
    table.item(row,9).setText('15')
    assert w.validate_c_part('refinement')
    assert w.document.config.refinement.bounds.specific[0].upper==15.


def test_dashboard_readiness_and_single_configuration(window,tmp_path):
    w=window;d=w.run_dashboard
    assert w.workflow_tabs.tabText(3)=='D. Run & Results'
    w.workflow_tabs.setCurrentIndex(3);w.show();QApplication.processEvents()
    assert d.run.isVisible() and not d.run.isEnabled() and not d.resume.isEnabled()
    assert 'Atomic Model must be validated' in d.reason.text()
    assert d.findChildren(QLineEdit).count(w.widgets[('run','output_directory')])==1
    assert w.widgets[('run','output_directory')] not in w.c_tabs.findChildren(QLineEdit)
    setup(w);assert w.validate_stage(2),w.stage_status[2].text()
    assert w.validate_output_directory();w.check_run_readiness()
    w.workflow_tabs.setCurrentIndex(3)
    assert d.run.isEnabled() and not d.resume.isEnabled() and not d.cancel.isEnabled()
    Path(w.document.config.run.output_directory).mkdir(exist_ok=True)
    d.update_state(check_external=True)
    assert d.status['Output directory'].text()=='Valid'
    d.navigate(1);assert w.c_tabs.currentIndex()==1 and w.workflow_tabs.currentIndex()==2
    w.widgets[('crystallography','depth_angstrom')].setText('7')
    w.workflow_tabs.setCurrentIndex(3)
    assert not d.run.isEnabled() and 'pending edits' in d.reason.text()
    assert w.workflow.receipts[:2]!=[None,None]


def test_cif_loading_preserves_thickness_and_basis_drafts(window):
    w=window
    w.widgets[('crystallography','depth_angstrom')].setText('1234')
    w.widgets[('crystallography','cbed_basis','g1_hkl')].setText('2,0,0')
    w.widgets[('crystallography','cbed_basis','g2_hkl')].setText('0,2,0')
    assert w.load_structure_path(str(Path(__file__).resolve().parents[1]/'reference_data/Al_FCC.cif'))
    assert w.widgets[('crystallography','depth_angstrom')].text()=='1234'
    w.site_group_table.item(0,3).setText('.86');w.apply_site_groups.click()
    assert w.widgets[('crystallography','depth_angstrom')].text()=='1234'
    assert w.document.config.crystallography.cbed_basis.g1_hkl==(2,0,0)
    assert w.document.config.crystallography.cbed_basis.g2_hkl==(0,2,0)


def test_absorption_bounds_use_configured_rules_only(window):
    from qcbed.gui.hard_bounds import configured_bounds
    c=window.document.config;g=c.crystallography.site_groups[0]
    c=replace(c,crystallography=replace(c.crystallography,site_groups=(replace(g,absorption=(.1,1.,1.)),)))
    infos=[i for i in available_refinement_targets(c) if i.available and i.category==Category.ABSORPTION]
    assert len(infos)==3
    assert all(configured_bounds(c,i)[:2]==(None,None) for i in infos)
    c=replace(c,refinement=replace(c.refinement,bounds=replace(c.refinement.bounds,absolute=(CategoryBound(Category.ABSORPTION,0.,1.5),))))
    assert all(configured_bounds(c,i)[:2]==(0.,1.5) for i in infos)
