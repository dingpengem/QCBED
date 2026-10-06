"""Opt-in frozen application acceptance harness; no pytest or developer imports.

QCBED --smoke-test /absolute/manifest.json uses explicit external input paths and
writes only into a NEW output directory. Normal application startup never uses it.
"""
import json
import os
from pathlib import Path
import platform
import sys
import time
import traceback


def main(manifest_path):
    manifest = json.loads(Path(manifest_path).read_text())
    if manifest.get('startup_only'):
        from PySide6.QtWidgets import QApplication
        from qcbed.gui.window import ConfigWindow
        app=QApplication.instance() or QApplication([])
        window=ConfigWindow();window.show();app.processEvents()
        assert window.workflow_tabs.count()==4 and window.document is None
        window.close();app.processEvents()
        Path(manifest['report']).write_text(json.dumps(dict(status='PASS',startup_only=True,optimizer_evaluations=0)))
        return 0
    if manifest.get('usability'):
        from usability_smoke import main as usability_main
        return usability_main(manifest)
    if manifest.get('optimizer_editor'):
        from optimizer_smoke import main as optimizer_main
        return optimizer_main(manifest)
    if manifest.get('modern'):
        from modern_smoke import main as modern_main
        return modern_main(manifest)
    output = Path(manifest['output'])
    if not output.is_absolute(): raise ValueError('Smoke output must be absolute')
    output.mkdir(parents=True,exist_ok=False)
    from PySide6.QtCore import QTimer
    from PySide6.QtWidgets import QApplication
    from qcbed.gui.window import ConfigWindow
    from qcbed.gui.document import EditorDocument
    from qcbed import load_yaml_config, compile_config, resolve_config_paths
    app = QApplication([])
    window = ConfigWindow()
    window.show()
    results = dict(frozen=bool(getattr(sys,'frozen',False)),executable=sys.executable,
                   architecture=platform.machine(),cwd=os.getcwd(),checks=[],status='running')
    started = time.monotonic()
    stage = 'begin'
    reference = None
    failing = False
    def record(name):
        results['checks'].append(name)
        (output/'acceptance.json').write_text(json.dumps(results,indent=2)+'\n')
    def tick():
        nonlocal stage, reference, failing
        try:
            if failing:
                if not window.controller.active: app.exit(1)
                return
            if time.monotonic()-started > 180: raise TimeoutError('Frozen acceptance timeout')
            if stage == 'begin':
                import numpy, scipy.linalg, scipy.fft, h5py, matplotlib, psutil, yaml
                import PySide6
                results['versions'] = {m.__name__:getattr(m,'__version__','loaded') for m in
                    (numpy,scipy,h5py,matplotlib,psutil,yaml,PySide6)}
                assert platform.machine() == 'arm64'
                assert scipy.linalg.det(numpy.eye(2)) == 1.
                with h5py.File(output/'probe.h5','w') as f: f['data'] = [1,2]
                assert psutil.Process().memory_info().rss > 0
                record('bundled_runtime_imports')
                window.open_document(manifest['reference_yaml'])
                window.validate_document()
                assert 'Compile validity: Valid' in window.summary.toPlainText()
                record('open_reference_and_compile')
                window.widgets[('title',)].setText('Frozen acceptance edit')
                assert window.save_document(output/'saved.yml')
                reference = window.document.config
                window.open_document(output/'saved.yml')
                assert window.document.config == reference
                record('edit_save_as_reload_equivalence')
                if manifest.get('reopen_only'):
                    finish(); return
                window.open_document(manifest['run_yaml'])
                config = window.document.config
                from dataclasses import replace
                config = replace(config,run=replace(config.run,
                    output_directory=str(output/'run'),output_filename=str(output/'out.dat')))
                window.document = EditorDocument(config,output/'run.yml')
                window.refresh(); window.validate_document()
                assert window.run_action.isEnabled()
                window.start_run()
                assert window.controller.active
                stage = 'run'
            elif window.controller.active or window.timer.isActive(): return
            elif stage == 'run':
                assert window.controller.state == 'Completed', window.run_log.toPlainText()
                results['run'] = window.controller.result
                assert window.controller.result['checkpoint'] and window.controller.result['final_state']
                assert 'classification' in window.run_log.toPlainText()
                (output/'run-gui-log.txt').write_text(window.run_log.toPlainText())
                record('spawn_run_progress_checkpoints_outputs')
                window.show_results()
                assert window.results_view.results is not None
                assert window.results_view.results.history
                assert window.results_view.results.comparison is not None
                record('legacy_results_history_patterns')
                assert window.resume_action.isEnabled()
                window.start_run(True); stage = 'resume'
            elif stage == 'resume':
                assert window.controller.state == 'Completed', window.run_log.toPlainText()
                assert window.controller.result['summary']['final_objective'] == results['run']['summary']['final_objective']
                record('compatible_resume')
                from dataclasses import replace
                config = window.document.config
                config = replace(config,run=replace(config.run,
                    output_directory=str(output/'cancelled'),output_filename=str(output/'cancelled.dat')))
                window.document = EditorDocument(config,output/'cancel.yml')
                window.refresh(); window.start_run(); window.cancel_run()
                assert not window.cancel_action.isEnabled()
                stage = 'cancel'
            elif stage == 'cancel':
                assert window.controller.state == 'Cancelled', window.run_log.toPlainText()
                record('cooperative_cancel_and_restore')
                from qcbed.legacy_config import MigrationPolicy
                from qcbed.configuration import BoundsConfig, RunConfig
                from qcbed.execution import ForwardExecutionConfig
                from qcbed.gui.compatibility import export_document, compatibility_report
                policy = MigrationPolicy(123,ForwardExecutionConfig('serial'),
                    RunConfig(str(output/'imported-run'),str(output/'imported.dat')),'legacy','raise',BoundsConfig())
                window.import_legacy_document(manifest['for007'],manifest['experimental'],policy)
                before = window.document.config
                assert window.document.needs_save_as
                assert not compatibility_report(before,base_directory=output)['errors']
                export_document(window.document,output/'FOR007.dat',filename_limits=(4096,4096))
                export_document(window.document,output/'bundle',bundle=True)
                assert window.document.config == before
                assert (output/'bundle/compatibility.json').exists()
                record('legacy_import_export_bundle')
                finish()
        except BaseException:
            results.update(status='failed',traceback=traceback.format_exc())
            (output/'acceptance.json').write_text(json.dumps(results,indent=2)+'\n')
            failing = True
            window.cancel_run()
            if not window.controller.active: app.exit(1)
    def finish():
        results['status'] = 'passed'
        results['elapsed_seconds'] = time.monotonic()-started
        window.grab().save(str(output/'window.png'))
        record('app_gui_complete')
        window.pending = False
        window.document.saved_config = window.document.config
        timer.stop()
        window.close()
        app.quit()
    timer = QTimer()
    timer.setInterval(50)
    timer.timeout.connect(tick)
    timer.start()
    return app.exec()
