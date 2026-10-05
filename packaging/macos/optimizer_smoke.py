"""Opt-in frozen C-editor smoke; never invokes a run or objective."""
import json
from pathlib import Path
import sys


def main(manifest):
    from PySide6.QtWidgets import QApplication
    from PySide6.QtCore import Qt
    from qcbed.gui.window import ConfigWindow
    from qcbed import load_yaml_config
    out=Path(manifest['output']);out.mkdir(parents=True,exist_ok=False)
    app=QApplication([]);w=ConfigWindow(Path(manifest['config']));w.show();app.processEvents()
    try:
        assert getattr(sys,'frozen',False)
        assert w.validate_stage(0),w.stage_status[0].text()
        assert w.validate_stage(1),w.stage_status[1].text()
        w.workflow_tabs.setCurrentIndex(2)
        assert w.validate_c_part('simulation'),w.stage_status[2].text()
        assert [w.c_tabs.tabText(i) for i in range(w.c_tabs.count())]==['Simulation Setup','Refinement Targets','Optimizer']
        assert all(w.refinement_table.isColumnHidden(c) for c in (4,5,6))
        assert not w.advanced_optimizer.isVisible()
        assert w.validate_c_part('refinement'),w.stage_status[2].text()
        assert w.save_document(out/'saved.yml')
        assert load_yaml_config(out/'saved.yml')==w.document.config
        w.c_tabs.setCurrentIndex(2);w.advanced_optimizer_toggle.click();app.processEvents()
        assert w.advanced_optimizer.isVisible()
        (out/'result.json').write_text(json.dumps(dict(status='PASS',frozen=True,
            checks=['three C tabs','A/B/setup validation','targets without visible legacy metadata',
                    'refinement commit','compact save/reload','advanced expansion'],optimizer_evaluations=0),indent=2))
    finally:
        w.document.saved_config=w.document.config;w.pending=False;w.close()
    return 0
