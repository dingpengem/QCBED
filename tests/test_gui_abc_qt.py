import pytest
pytest.importorskip('PySide6')
from PySide6.QtWidgets import QApplication
from qcbed import save_yaml_config
from qcbed.gui.window import ConfigWindow
from qcbed.gui.data_model import attach_overview
from test_modern_grouped_run import modern_config
from test_gui_data_model import overview_files


def test_abc_navigation_mask_views_save_and_fresh_paths(tmp_path):
    app=QApplication.instance() or QApplication([])
    c=modern_config(tmp_path);pattern,manifest=overview_files(c,tmp_path)
    c=attach_overview(c,tmp_path,pattern,manifest)
    path=tmp_path/'config.yml';save_yaml_config(c,path,base_directory=tmp_path)
    w=ConfigWindow(path)
    assert [w.navigation.topLevelItem(i).text(0) for i in range(3)]==['A. Data Prep','B. Atomic Model','C. Main QCBED Parameters']
    for i in range(3):
        item=w.navigation.topLevelItem(i);w.navigation.setCurrentItem(item)
        assert w.pages.currentIndex()==item.data(0,256)
    assert [a.get_title() for a in w.atomic_axes]==['Front / Beam — along +z_beam','Top — along +y_beam','Left — along +x_beam']
    assert w.site_group_table.columnCount()==8
    v=w.data_panel
    assert v.canvas.minimumHeight()>=300 and v.overview_canvas.minimumHeight()>=400
    assert len(v.axes)==3 and len(v.overview_figure.axes[0].patches)==9
    assert v.overview_figure.axes[0].lines[0].get_xdata()[0]==3.25
    v.beams.setCurrentIndex(1);v.mask()
    assert w.data_panel.data.weights[0,1]==0
    w.undo();assert w.data_panel.data.weights[0,1]==1
    w.redo();assert w.data_panel.data.weights[0,1]==0
    assert w.save_document(tmp_path/'other.yml')
    w.open_document(tmp_path/'other.yml');assert w.data_panel.data.weights[0,1]==0
    from pathlib import Path
    Path(w.document.config.run.output_directory).mkdir()
    w.validate_document();assert not w.run_action.isEnabled() and not w.resume_action.isEnabled()
    w.new_run_paths();assert not Path(w.document.config.run.output_directory).exists()
    assert not w.run_action.isEnabled() and not w.resume_action.isEnabled()  # Explicit A/B/C validation required.
    w.document.saved_config=w.document.config;w.close()
