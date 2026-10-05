"""Widget-driven real-input usability acceptance; explicit pre-evaluation stop."""
import json
import sys
import time
from pathlib import Path
from types import SimpleNamespace


def main(manifest):
    from PySide6.QtWidgets import QApplication
    from PySide6.QtCore import Qt
    from qcbed.gui.window import ConfigWindow
    from qcbed import load_yaml_config,prepare_run_snapshot,run_snapshot
    out=Path(manifest['output']);out.mkdir(parents=True,exist_ok=False)
    app=QApplication.instance() or QApplication([]);w=ConfigWindow();w.show();app.processEvents()
    checks=[]
    try:
        if manifest.get('require_frozen',True):assert getattr(sys,'frozen',False)
        assert w.document is None
        assert not list(out.rglob('*.yml'))
        assert not w.run_dashboard.results.isEnabled()
        checks.append('blank New Data Preparation project exists without YAML; Open Results disabled')
        if manifest.get('example_project'):
            from qcbed.gui.preparation_model import save_raw_project
            from qcbed.yaml_config import read_yaml_document
            example=Path(manifest['example_project'])
            save_raw_project(read_yaml_document(example),out/'project.yml',base_directory=example.parent)
            w.open_document(out/'project.yml')
        v=w.preparation;v.raw_path=out/'project.yml'
        assert not v.advanced_settings.toggle.isChecked()
        v.inputs['pattern'].setText(manifest['pattern']);v.load_pattern()
        def drag(a,b):
            v.begin(SimpleNamespace(xdata=a[0],ydata=a[1],inaxes=True))
            v.end(SimpleNamespace(xdata=b[0],ydata=b[1],inaxes=True));v.accept_preview()
        drag(*manifest['diameter_endpoints'])
        v.tool.setCurrentIndex(1);r=manifest['symmetry_roi'];drag(r[:2],r[2:])
        v.signal.setCurrentText('gradient');v.run.click()
        deadline=time.monotonic()+120
        while v.active and time.monotonic()<deadline:app.processEvents();time.sleep(.02)
        assert not v.active and w.document is not None,v.status.text()
        assert w.document.config.crystallography.depth_angstrom==1000.
        w.data_panel.select_zone();assert w.validate_stage(0),w.stage_status[0].text()
        assert not v.advanced_settings.toggle.isChecked()
        checks.append('live real-pattern A without Advanced')
        w.workflow_tabs.setCurrentIndex(1)
        assert w.load_structure_path(manifest['structure'])
        if not manifest.get('example_project'):
            assert w.site_group_table.item(0,3).text()==''
        w.widgets[('crystallography','cbed_basis','g1_hkl')].setText('2,0,0')
        w.widgets[('crystallography','cbed_basis','g2_hkl')].setText('0,2,0')
        w.site_group_table.item(0,3).setText('.86');w.apply_site_groups.click()
        assert w.validate_stage(1),w.stage_status[1].text()
        assert w.document.config.crystallography.cbed_basis.g1_hkl==(2,0,0)
        assert w.document.config.crystallography.cbed_basis.g2_hkl==(0,2,0)
        checks.append('real CIF and explicit Debye–Waller B')
        w.workflow_tabs.setCurrentIndex(2)
        boundary=w.widgets[('crystallography','slice_boundaries_fractional')]
        assert boundary.text()=='0.0, 1.0' and not boundary.isReadOnly()
        assert not w.numerical_advanced.toggle.isChecked()
        app.processEvents();assert boundary.isVisible()
        assert 'one slice' in w.slice_boundary_help.text()
        assert w.validate_c_part('simulation'),w.stage_status[2].text()
        assert w.workflow.simulation_receipt is not None and not w.run_action.isEnabled()
        checks.append('visible 0, 1 accepted as one stored slice with Advanced collapsed')
        w.widgets[('crystallography','depth_angstrom')].setText(str(manifest.get('study_depth',840)))
        assert w.validate_c_part('simulation'),w.stage_status[2].text()
        assert w.recommended_numerics.isChecked()
        assert not w.numerical_advanced.toggle.isChecked()
        w.c_tabs.setCurrentIndex(0);app.processEvents();w.grab().save(str(out/'simulation.png'))
        t=w.refinement_table
        for row in range(t.rowCount()):
            category=t.item(row,0).data(Qt.ItemDataRole.UserRole).category.value
            if category not in ('depth','debye_waller'):continue
            t.item(row,1).setCheckState(Qt.CheckState.Checked)
            t.item(row,8).setText('600' if category=='depth' else '.5')
            t.item(row,9).setText('1200' if category=='depth' else '1.5')
        assert not w.advanced_optimizer_toggle.isChecked()
        assert w.validate_c_part('refinement'),w.stage_status[2].text()
        assert len(w.document.config.refinement.parameters)==2
        assert w.c_tabs.count()==3
        w.c_tabs.setCurrentIndex(2);app.processEvents()
        for name,value in manifest.get('optimizer_overrides',{}).items():
            w.widgets[('refinement','optimizer',name)].setText(str(value))
        w.widgets[('run','output_directory')].setText('')
        assert w.validate_stage(2),w.stage_status[2].text()
        assert w.workflow.ready and not w.run_action.isEnabled()
        assert 'Output directory not selected' in w.run_dashboard.reason.text()
        w.workflow_tabs.setCurrentIndex(3)
        assert not w.execution_advanced.toggle.isChecked()
        w.execution_advanced.set_expanded(True);app.processEvents()
        mode=w.widgets[('execution','execution_backend')];batch=w.widgets[('execution','batch_size')]
        assert mode.isVisible() and batch.isVisible()
        assert w.execution_advanced.isAncestorOf(mode) and not w.c_tabs.isAncestorOf(mode)
        mode.setCurrentText('batched');batch.setText('64')
        assert w.apply_edits()
        assert w.workflow.ready and w.document.config.execution.batch_size==64
        w.workflow_tabs.setCurrentIndex(3);w.execution_advanced.set_expanded(True);app.processEvents()
        for key in (('execution','execution_backend'),('execution','batch_size'),
                    ('run','max_total_evaluations'),('run','max_wall_time_seconds')):
            widget=w.widgets[key]
            deadline=time.monotonic()+3
            while time.monotonic()<deadline:
                w.run_dashboard.content_scroll.ensureWidgetVisible(widget);app.processEvents()
                if not widget.visibleRegion().isEmpty():break
                time.sleep(.01)
            assert not widget.visibleRegion().isEmpty(),key
        w.grab().save(str(out/'advanced-execution.png'))
        w.execution_advanced.set_expanded(False)
        w.widgets[('run','output_directory')].setText(str(out/'run'))
        assert w.validate_output_directory()
        assert (out/'run').is_dir()
        (out/'run'/'user-notes.txt').write_text('preserve me')
        assert w.validate_stage(2),w.stage_status[2].text()
        w.workflow_tabs.setCurrentIndex(3);app.processEvents()
        assert w.run_dashboard.run.isVisible() and w.run_dashboard.run.isEnabled()
        assert not w.run_dashboard.resume.isEnabled() and not w.run_dashboard.cancel.isEnabled()
        assert w.run_dashboard.status['Output directory'].text()=='Valid'
        assert w.run_dashboard.status['Resume'].text()=='Not available — no previous run found'
        w.run_dashboard.content_scroll.verticalScrollBar().setValue(0);app.processEvents()
        w.grab().save(str(out/'dashboard.png'))
        checks.append('fresh 1000 Å thickness, preserved g1/g2, visible D dashboard and guarded actions')
        checks.append('three C tabs valid without output; D exclusively owns output and advanced execution')
        checks.append('packaged Advanced exposes Batched and batch size 64; no duplicate C controls')
        assert w.save_document(out/'project.yml',overwrite=bool(manifest.get('example_project')))
        assert load_yaml_config(out/'project.yml')==w.document.config
        assert load_yaml_config(out/'project.yml').execution.batch_size==64
        from qcbed.yaml_config import read_yaml_document
        editable=read_yaml_document(out/'project.yml')
        assert set(editable)=={'schema_version','data_prep','atomic_model','qcbed'}
        assert 'filename' not in editable['qcbed'].get('output',{})
        assert len(list(out.rglob('*.yml')))==1 and not list(out.rglob('*.yaml'))
        for forbidden in ('internal_hk','effective_depth','checkpoint_identity','optimizer_history','validation_receipt'):
            assert forbidden not in json.dumps(editable).lower()
        saved=w.document.config
        table=w.refinement_table
        saved_rows={table.item(row,0).data(Qt.ItemDataRole.UserRole):
            (table.item(row,1).checkState(),table.item(row,8).text(),table.item(row,9).text())
            for row in range(table.rowCount())}
        w.open_document(out/'project.yml');app.processEvents()
        assert w.document.config==saved
        assert float(w.widgets[('crystallography','depth_angstrom')].text())==manifest.get('study_depth',840)
        assert w.widgets[('crystallography','slice_boundaries_fractional')].text()=='0.0, 1.0'
        assert w.widgets[('execution','execution_backend')].currentText()=='batched'
        assert w.widgets[('execution','batch_size')].text()=='64'
        assert w.widgets[('run','output_directory')].text()==str(out/'run')
        assert float(w.site_group_table.item(0,3).text())==.86
        assert w.document.config.crystallography.cbed_basis.g1_hkl==(2,0,0)
        assert w.document.config.crystallography.cbed_basis.g2_hkl==(0,2,0)
        assert w.widgets[('crystallography','cbed_basis','g1_hkl')].text()=='2, 0, 0'
        assert w.widgets[('crystallography','cbed_basis','g2_hkl')].text()=='0, 2, 0'
        for name,value in manifest.get('optimizer_overrides',{}).items():
            assert float(w.widgets[('refinement','optimizer',name)].text())==float(value)
            assert getattr(w.document.config.refinement.optimizer,name)==value
        for stage in range(2):
            assert w.validate_stage(stage),w.stage_status[stage].text()
        assert w.validate_c_part('simulation')
        table=w.refinement_table
        restored_rows={table.item(row,0).data(Qt.ItemDataRole.UserRole):
            (table.item(row,1).checkState(),table.item(row,8).text(),table.item(row,9).text())
            for row in range(table.rowCount())}
        assert restored_rows==saved_rows
        assert w.validate_c_part('refinement')
        assert w.validate_stage(2)
        assert not w.run_action.isEnabled()
        assert w.validate_output_directory()
        w.check_run_readiness()
        assert w.run_action.isEnabled() and not w.resume_action.isEnabled()
        checks.append('saved YAML reopened in GUI: physical, slicing, execution and output values restored; D revalidated')
        checks.append('single compact YAML; target selections/bounds and explicit optimizer overrides restored')
        from qcbed.gui.output_paths import modern_execution_config
        snapshot=prepare_run_snapshot(modern_execution_config(w.document.config,out),base_directory=out)
        result=run_snapshot(snapshot,cancel_requested=lambda:True)
        assert result.status=='cancelled'
        assert not list((out/'run').rglob('history.csv'))
        assert (out/'run'/'user-notes.txt').read_text()=='preserve me'
        w.check_run_readiness()
        assert w.resume_action.isEnabled(),w.resume_action.toolTip()
        assert w.open_results(snapshot.compiled.configuration.run.output_directory)
        assert w.run_dashboard.results.isEnabled()
        assert w.results_view.results.summary['recorded_evaluations']==0
        w.run_dashboard.results.click();app.processEvents()
        assert w.workspace.currentWidget() is w.results_view
        checks.append('Open Results enabled for initialization checkpoint; zero completed evaluations')
        checks.append('explicit writable-directory validation, preserved unrelated file, separate Resume readiness')
        checks.append('canonical initialization with explicit stop before evaluation')
        w.workspace.setCurrentIndex(0)
        w.workflow_tabs.setCurrentIndex(2);w.c_tabs.setCurrentIndex(2);app.processEvents();w.grab().save(str(out/'optimizer.png'))
        report=dict(status='PASS',frozen=bool(getattr(sys,'frozen',False)),checks=checks,
            optimizer_evaluations=0,project=str(out/'project.yml'),
            limitation='One-slice support retains strict transverse-periodicity checks; mesh and slicing convergence remain specimen-dependent.')
        (out/'result.json').write_text(json.dumps(report,indent=2))
    finally:
        if w.document:w.document.saved_config=w.document.config
        w.pending=False;w.preparation.dirty=False;w.close()
    return 0
