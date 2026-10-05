"""Workflow checks use prepared fixtures and never launch an optimizer."""
from dataclasses import replace
from types import SimpleNamespace
import json
import numpy as np
import pytest
pytest.importorskip('PySide6')
from PySide6.QtWidgets import QApplication
from qcbed.gui.window import ConfigWindow
from qcbed.gui.document import EditorDocument
from qcbed.gui.data_model import attach_overview,provenance
from test_modern_grouped_run import modern_config
from test_gui_data_model import overview_files


@pytest.fixture
def window(tmp_path):
    app=QApplication.instance() or QApplication([])
    c=modern_config(tmp_path);pattern,manifest=overview_files(c,tmp_path)
    m=json.loads(manifest.read_text())
    m['reciprocal_frame'].update(central_disk_center_px=[3.25,4.75],basis_vector_1_px=[5.,0.],basis_vector_2_px=[0.,5.])
    m['extraction_config']=dict(selection_width=4,selection_height=4,final_binning=2)
    for i,b in enumerate(m['beam_mapping']):b['lattice_coordinates']=list(c.crystallography.beams[i].lattice_coordinates)
    manifest.write_text(json.dumps(m));c=attach_overview(c,tmp_path,pattern,manifest)
    w=ConfigWindow();w.document=EditorDocument(c,tmp_path/'config.yml');w.refresh()
    yield w
    w.document.saved_config=w.document.config;w.pending=False;w.close();w.deleteLater();app.processEvents()


def validate_all(w):
    w.data_panel.select_zone()
    assert w.validate_stage(0)
    assert w.validate_stage(1)
    assert w.validate_c_part('simulation')
    assert w.validate_stage(2)


def test_tabs_order_snapshots_and_invalidation(window):
    w=window
    assert [w.workflow_tabs.tabText(i) for i in range(3)]==['A. Data Prep','B. Atomic Model','C. Main QCBED Parameters']
    assert not w.validate_stage(1) and not w.validate_stage(0)
    validate_all(w);assert w.validate_output_directory();w.check_run_readiness();assert w.run_action.isEnabled()
    before=w.data_panel.data
    w.data_panel.beams.setCurrentIndex(1);w.data_panel.mask()
    assert w.workflow.receipts==[None,None,None] and not w.run_action.isEnabled()
    np.testing.assert_array_equal(before.intensities,w.data_panel.data.intensities)
    np.testing.assert_array_equal(before.variances,w.data_panel.data.variances)
    assert w.data_panel.data.weights[0,1]==0
    assert w.validate_stage(0);assert w.workflow.snapshot==w.document.config
    assert w.validate_stage(1);assert w.validate_c_part('simulation');assert w.validate_stage(2)
    g=w.document.config.crystallography.site_groups[0]
    c=replace(w.document.config,crystallography=replace(w.document.config.crystallography,site_groups=(replace(g,debye_waller_b=.9),)))
    w.document.commit_structure(c);w.refresh()
    assert w.workflow.receipts[0] is not None and w.workflow.receipts[1:]==[None,None]
    w.data_panel.beams.setCurrentIndex(1);w.data_panel.mask(restore=True);assert w.data_panel.data.weights[0,1]==1
    w.data_panel.mask();w.data_panel.mask(reset=True)
    np.testing.assert_array_equal(before.weights,w.data_panel.data.weights)


def test_interactive_feature_marker_annotations_and_save(window,tmp_path):
    w=window;v=w.data_panel;v.beams.setCurrentIndex(0);v.select_feature.setChecked(True)
    old=w.document.config.scan
    v.pick_overview(SimpleNamespace(xdata=4.,ydata=5.,inaxes=v.overview_figure.axes[0]))
    ax=w.data_panel.overview_figure.axes[0]
    assert len(ax.patches)==9 and len(ax.texts)==9
    assert '(0, 0)' in ax.texts[0].get_text()
    assert ax.lines[0].get_label()=='Zone axis'
    assert tuple([ax.lines[0].get_xdata()[0],ax.lines[0].get_ydata()[0]])==pytest.approx((5.125,5.375))
    assert any(l.get_label()=='Central / lattice origin' for l in ax.lines)
    assert w.document.config.scan!=old
    assert 'manual' in provenance(w.document.config)['zone_axis_selection']['method']
    saved=w.document.config;assert w.save_document(tmp_path/'save.yml')
    import yaml
    assert set(yaml.safe_load((tmp_path/'save.yml').read_text()))=={'schema_version','data_prep','atomic_model','qcbed'}
    w.open_document(tmp_path/'save.yml');assert w.document.config==saved
    assert not w.workflow.ready


def test_existing_outputs_fresh_ready_resume_missing_and_backup(window,tmp_path):
    from pathlib import Path
    from qcbed.gui.readiness import readiness
    from qcbed.gui.output_paths import backup_outputs
    w=window;out=Path(w.document.config.run.output_directory);out.mkdir()
    (out/'old.txt').write_text('retained')
    validate_all(w);assert w.validate_output_directory();w.check_run_readiness();assert w.run_action.isEnabled() and not w.resume_action.isEnabled()
    assert not [i for i in readiness(w.document.config,tmp_path) if i.severity=='error']
    (out/'checkpoint_latest.json').write_text('{}')
    assert any(i.code=='resume_incompatible' for i in readiness(w.document.config,tmp_path,resume=True))
    moves=backup_outputs(w.document.config,tmp_path)
    assert not out.exists() and (moves[0][1]/'old.txt').read_text()=='retained'


def test_run_guard_does_not_launch_before_validation(window,monkeypatch):
    calls=[];monkeypatch.setattr(window.controller,'start',lambda *a,**k:calls.append(k))
    window.start_run();assert calls==[]
