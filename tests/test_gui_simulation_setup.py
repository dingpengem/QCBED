from dataclasses import replace
import pytest
pytest.importorskip('PySide6')
from PySide6.QtWidgets import QApplication
from qcbed.gui.window import ConfigWindow
from qcbed.gui.document import EditorDocument,parse_text
from qcbed import load_yaml_config
from test_modern_grouped_run import modern_config


@pytest.mark.parametrize('text',['0,1','0, 1','[0.0, 1.0]'])
def test_one_slice_parsing(text):
    assert parse_text(text,'boundaries')==(0.,1.)


def test_one_slice_visible_advanced_mesh_non001_and_save(tmp_path):
    import numpy as np
    from ase import Atoms
    from ase.io import write
    from qcbed.models import MeshConfig
    app=QApplication.instance() or QApplication([])
    c=modern_config(tmp_path)
    path=tmp_path/'orthorhombic.extxyz'
    write(path,Atoms('Al',cell=np.diag([3.,4.,5.]),scaled_positions=[[.1,.2,.3]],pbc=True))
    crystal=replace(c.crystallography,structure=replace(c.crystallography.structure,path=str(path)),
        zone_axis=(1,0,0),cbed_basis=replace(c.crystallography.cbed_basis,g1_hkl=(0,1,0),g2_hkl=(0,0,1)),
        mesh=MeshConfig(8,6),slice_boundaries_fractional=(0.,1.),
        site_groups=(replace(c.crystallography.site_groups[0],site_indices=(0,)),))
    c=replace(c,crystallography=crystal,microscope=replace(c.microscope,zone_axis=(1,0,0)))
    w=ConfigWindow();w.document=EditorDocument(c,tmp_path/'config.yml');w.refresh();w.show()
    try:
        w.workflow.accept(0,c);assert w.validate_stage(1)
        w.workflow_tabs.setCurrentIndex(2);w.c_tabs.setCurrentIndex(0);app.processEvents()
        key=('crystallography','slice_boundaries_fractional')
        assert not w.numerical_advanced.toggle.isChecked()
        assert w.widgets[key].isVisible()
        assert w.widgets[key].visibleRegion().boundingRect().width()>100
        assert not w.widgets[('crystallography','mesh','mesh_x')].isVisible()
        assert 'one slice' in w.slice_boundary_help.text()
        assert '4 × 5 Å' in w.transverse_mesh_summary.text()
        assert '0.5 / 0.833333 Å' in w.transverse_mesh_summary.text()
        w.widgets[key].setText('0,1')
        assert w.validate_c_part('simulation'),w.stage_status[2].text()
        assert w.validate_stage(2),w.stage_status[2].text()
        assert w.save_document(tmp_path/'saved.yml')
        assert load_yaml_config(tmp_path/'saved.yml').crystallography.slice_boundaries_fractional==(0.,1.)
        assert w.document.config.crystallography.depth_angstrom==4.
    finally:
        w.document.saved_config=w.document.config;w.pending=False;w.close()
