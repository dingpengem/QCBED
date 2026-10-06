from dataclasses import replace
from pathlib import Path
import pytest
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication,QPushButton
from qcbed.gui.window import ConfigWindow
from qcbed.gui.document import EditorDocument
from test_projected_elastic import config


@pytest.fixture
def window(tmp_path):
    app=QApplication.instance() or QApplication([])
    c=config(tmp_path,'real',reflections=((2,0,0),(0,2,0)))
    w=ConfigWindow();w.document=EditorDocument(c,tmp_path/'project.yml');w.refresh();w.show()
    w.workflow.accept(0,c)
    yield w
    w.site_group_pending=False;w.pending=False;w.document.saved_config=w.document.config;w.close()


def test_structure_site_groups_navigation_and_basis(window,tmp_path):
    w=window
    w.report_error=lambda error:pytest.fail(str(error))
    p=tmp_path/'model.cif';p.write_bytes((Path(__file__).parents[1]/'reference_data/Al_FCC.cif').read_bytes())
    assert w.load_structure_path(str(p))
    assert w.navigation.currentItem().text(0)=='Structure'
    w.widgets[('crystallography','cbed_basis','g1_hkl')].setText('2,0,0')
    w.widgets[('crystallography','cbed_basis','g2_hkl')].setText('0,2,0')
    w.next_site_groups.click();assert w.navigation.currentItem().text(0)=='Site groups'
    w.site_group_table.item(0,3).setText('.9');assert 'Modified' in w.site_group_status.text()
    w.apply_site_groups.click();assert w.navigation.currentItem().text(0)=='Site groups'
    w.validate_atomic_button.click()
    assert w.workflow.receipts[1] is not None
    assert w.document.config.crystallography.cbed_basis.g1_hkl==(2,0,0)
    assert w.document.config.crystallography.cbed_basis.g2_hkl==(0,2,0)


def test_targets_stay_in_place_and_explicit_next(window):
    w=window;assert w.validate_stage(1);w.workflow_tabs.setCurrentIndex(2)
    assert w.validate_c_part('simulation');assert w.c_tabs.currentIndex()==0
    w.c_next['simulation'].click();assert w.c_tabs.currentIndex()==1
    p=w.elastic_panel;indices=[i for i,r in enumerate(p.rows) if r['record']]
    p.link_rows(indices,'Pair');assert w.c_tabs.currentIndex()==1
    assert 'Pair' in w.elastic_panel.live_groups.text()
    assert 'Draft optimizer layout' in w.elastic_panel.dimension_summary.text()
    w.validate_refinement_button.click();assert w.workflow.refinement_receipt is not None
    assert w.c_tabs.currentIndex()==1
    assert 'Validated optimizer layout' in w.elastic_panel.dimension_summary.text()
    assert 'Pair' in w.elastic_panel.resolved.toPlainText()
    w.c_next['refinement'].click();assert w.c_tabs.currentIndex()==2
    w.workflow_tabs.setCurrentIndex(3)
    next(b for b in w.run_dashboard.findChildren(QPushButton) if b.text()=='Go to Refinement Targets').click()
    assert w.c_tabs.currentIndex()==1
    p=w.elastic_panel;i=indices[0];old=p.table.item(i,5).text()
    p.table.item(i,6).setText('-99');w.refresh()
    assert w.c_tabs.currentIndex()==1
    assert w.elastic_panel.table.item(i,6).text()=='-99' and w.elastic_panel.table.item(i,5).text()==old
    assert w.workflow.refinement_receipt is None
    assert 'Next step:' in w.guidance.next.text()


def test_explicit_clear_discards_coefficient_drafts(window):
    w=window;assert w.validate_stage(1);assert w.validate_c_part('simulation')
    w.workflow_tabs.setCurrentIndex(2);w.c_tabs.setCurrentIndex(1)
    p=w.elastic_panel;i=next(i for i,r in enumerate(p.rows) if r['record'])
    p.table.item(i,5).setText('19')
    next(b for b in p.findChildren(QPushButton) if b.text()=='Clear saved coefficient targets').click()
    assert not w.document.config.crystallography.projected_elastic
    assert not w.elastic_panel.dirty and w.c_tabs.currentIndex()==1
    assert all(w.elastic_panel.table.item(i,0).checkState()==Qt.CheckState.Unchecked for i in range(w.elastic_panel.table.rowCount()))


def test_upstream_b_revalidation_retains_pending_target_values(window):
    w=window;assert w.validate_stage(1);assert w.validate_c_part('simulation')
    w.workflow_tabs.setCurrentIndex(2);w.c_tabs.setCurrentIndex(1)
    p=w.elastic_panel;i=next(i for i,r in enumerate(p.rows) if r['hkl']==(2,0,0) and r['record'])
    p.table.item(i,5).setText('19.5');p.table.item(i,6).setText('0');p.table.item(i,7).setText('40')
    w.workflow_tabs.setCurrentIndex(1);w.navigation.setCurrentItem(w.section_items['Site groups'])
    w.site_group_table.item(0,3).setText('.95');w.apply_site_groups.click()
    assert w.workflow.simulation_receipt is None
    assert w._target_draft_cache[2]
    assert w.validate_stage(1)
    w.workflow_tabs.setCurrentIndex(2);w.c_tabs.setCurrentIndex(0)
    assert w.validate_c_part('simulation')
    p=w.elastic_panel;i=next(i for i,r in enumerate(p.rows) if r['hkl']==(2,0,0) and r['record'])
    assert p.table.item(i,5).text()=='19.5'
    assert p.table.item(i,6).text()=='0' and p.table.item(i,7).text()=='40'
    assert p.table.item(i,0).checkState()==Qt.CheckState.Checked


def test_revert_action_discards_target_drafts(window):
    w=window;assert w.validate_stage(1);assert w.validate_c_part('simulation')
    w.workflow_tabs.setCurrentIndex(2);w.c_tabs.setCurrentIndex(1)
    p=w.elastic_panel;i=next(i for i,r in enumerate(p.rows) if r['record'])
    initial=p.table.item(i,5).text();p.table.item(i,5).setText('17')
    next(a for a in w.edit_actions if a.text()=='Revert drafts').trigger()
    assert w.elastic_panel.table.item(i,5).text()==initial
    assert not w.elastic_panel.dirty and not w.pending
    assert w.c_tabs.currentIndex()==1
