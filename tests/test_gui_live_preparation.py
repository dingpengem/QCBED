from pathlib import Path
from types import SimpleNamespace
import time
import numpy as np
import pytest
pytest.importorskip('PySide6')
from PySide6.QtWidgets import QApplication
from qcbed.gui.window import ConfigWindow
from qcbed.gui.preparation_model import prepare_files,apply_preparation,checksum
from qcbed.gui.data_model import read_data,select_zone_axis
from qcbed import save_yaml_config,load_yaml_config
from test_live_preparation import scene,selection
from test_modern_grouped_run import modern_config


@pytest.fixture
def window():
    app=QApplication.instance() or QApplication([]);w=ConfigWindow()
    yield w
    if w.document:w.document.saved_config=w.document.config
    w.pending=False;w.close();w.deleteLater();app.processEvents()


def drag(view,a,b):
    view.begin(SimpleNamespace(xdata=a[0],ydata=a[1],inaxes=True))
    view.end(SimpleNamespace(xdata=b[0],ydata=b[1],inaxes=True))


def test_live_worker_new_window_select_mask_validate(window,tmp_path):
    view=window.preparation;assert view.mode.currentText()=='New Data Preparation'
    path=tmp_path/'image.npy';np.save(path,scene())
    view.raw_path=tmp_path/'project.yml';view.inputs['pattern'].setText(str(path));view.load_pattern()
    drag(view,(70,90),(110,90));assert view.geometry is None
    view.accept_preview();assert view.geometry.radius_x_px==20
    view.tool.setCurrentIndex(1);drag(view,(86.8,87.2),(92.7,93.1))
    assert view.roi is None;view.accept_preview();assert view.roi==(86,86,94,94)
    view.binning.setValue(1)
    assert view.dirty
    assert not window.validate_stage(0)
    assert 'run the current preparation' in window.stage_status[0].text()
    view.run.click();assert view.active and view.process.pid!=__import__('os').getpid()
    app=QApplication.instance();deadline=time.monotonic()+60
    while view.active and time.monotonic()<deadline:
        app.processEvents();time.sleep(.01)
    assert not view.active,view.status.text()
    assert window.document is not None,view.status.text()
    assert window.document.config.crystallography.slice_boundaries_fractional==(0.,1.)
    assert not window.numerical_advanced.toggle.isChecked()
    assert window.data_panel.data is not None
    assert not window.validate_stage(0)  # visual origin confirmation remains required
    window.data_panel.select_zone()
    assert window.validate_stage(0),window.stage_status[0].text()
    from qcbed.gui.data_model import restored_overview
    recorded=restored_overview(window.document.config,window.document.base_directory)[1]['live_preparation']['selection']
    assert (recorded['selection_width'],recorded['selection_height'])==(40,40)
    assert recorded['symmetry_crop_half_size_px']==[4,4]
    assert recorded['symmetry_roi']==[86,86,94,94]
    before=read_data(window.document.config,window.document.base_directory)
    window.data_panel.mask()
    after=read_data(window.document.config,window.document.base_directory)
    np.testing.assert_array_equal(before.intensities,after.intensities)
    np.testing.assert_array_equal(before.variances,after.variances)
    assert window.validate_stage(0)
    assert not window.validate_stage(1)  # no borrowed atomic model


def test_aperture_threshold_preview_and_override(window,tmp_path):
    v=window.preparation;np.save(tmp_path/'p.npy',scene());v.inputs['pattern'].setText(str(tmp_path/'p.npy'));v.load_pattern()
    y,x=np.indices((180,180));image=100.*(((x-90)/20)**2+((y-90)/15)**2<1)
    np.save(tmp_path/'ap.npy',image);v.inputs['calibration_aperture'].setText(str(tmp_path/'ap.npy'));v.load_profile()
    assert v.geometry is None;v.accept_preview()
    assert (v.geometry.radius_x_px,v.geometry.radius_y_px)==(19.5,14.5)
    assert not v.inputs['disk_profile'].text()
    drag(v,(70,90),(110,90));v.accept_preview();assert v.geometry.radius_x_px==20
    assert v.diameter_endpoints is not None
    v.load_pattern();assert v.geometry.radius_x_px==20  # Continuing the same draft does not discard its geometry.
    np.save(tmp_path/'p.npy',scene()+1)  # A changed source must reset derived state.
    v.load_pattern();assert v.geometry is None and v.diameter_endpoints is None and v.geometry_input is None


def test_live_save_reload_raw_references_and_prepared_replay(tmp_path):
    image=tmp_path/'raw.npy';np.save(image,scene())
    descriptor=prepare_files(image,selection(),tmp_path/'preparation')
    c=apply_preparation(modern_config(tmp_path),descriptor)
    c=select_zone_axis(c,tmp_path,beam=0)
    saved=tmp_path/'config.yml';save_yaml_config(c,saved,base_directory=tmp_path)
    import yaml
    raw=yaml.safe_load(saved.read_text())
    assert raw['data_prep']['pattern']==str(image) and 'preparation' in raw['data_prep']
    assert load_yaml_config(saved)==c
    np.testing.assert_array_equal(read_data(c,tmp_path).intensities,read_data(load_yaml_config(saved),tmp_path).intensities)


def test_raw_yaml_opens_without_prepared_dataset(window,tmp_path):
    from qcbed.gui.preparation_model import save_raw_project
    image=tmp_path/'raw.npy';np.save(image,scene())
    project=tmp_path/'raw.yml';save_raw_project({'schema_version':2,'data_prep':{'pattern':'raw.npy'},'qcbed':{'signal_mode':'gradient'}},project,base_directory=tmp_path)
    window.open_document(project)
    assert window.document is None and window.preparation.image.shape==(180,180)
    assert window.preparation.geometry is None and window.preparation.roi is None
    assert window.preparation.signal.currentText()=='gradient'


def test_derived_controls_and_binning(window):
    from qcbed.gui.preparation_model import derive_preparation_dimensions
    from qcbed.data_preparation.models import DiskGeometry
    v=window.preparation
    assert not hasattr(v,'dimensions') and v.binning.value()==2
    assert derive_preparation_dimensions(DiskGeometry(38.,38.),(111,105,155,151),2)==(76,76,(22,23))
    assert derive_preparation_dimensions(DiskGeometry(20.,16.),(10,10,18,22),2)==(40,32,(4,6))
    assert derive_preparation_dimensions(DiskGeometry(19.5,14.5),(10,10,18,22),2)==(40,30,(4,6))
    assert derive_preparation_dimensions(DiskGeometry(19.5,14.5),(10,10,18,22),1)==(39,29,(4,6))
    assert derive_preparation_dimensions(DiskGeometry(20.,20.),(10,10,19,22),2)==(40,40,(5,6))
    assert derive_preparation_dimensions(DiskGeometry(20.1,20.),(10,10,18,22),2)==(42,40,(4,6))
    v.binning.setValue(1);assert v.binning.value()==1


def test_replay_retains_old_crop_exactly(window,tmp_path):
    image=tmp_path/'raw.npy';np.save(image,scene())
    descriptor=prepare_files(image,selection(),tmp_path/'preparation')
    window.preparation.install(descriptor)
    v=window.preparation
    assert v.roi==selection().symmetry_roi
    assert v.binning.value()==selection().final_binning
    assert '(8, 8)' in v.derived.text()  # old explicit crop, not newly derived (3,3)


def test_fractional_geometry_expands_array_only(window):
    from qcbed.gui.preparation_model import derive_preparation_dimensions
    from qcbed.data_preparation.models import DiskGeometry
    v=window.preparation
    drag(v,(0,0),(75.6,0));v.accept_preview()
    assert v.geometry==DiskGeometry(37.8,37.8)
    v.tool.setCurrentIndex(1);drag(v,(10,10),(18,22));v.accept_preview()
    assert derive_preparation_dimensions(v.geometry,v.roi,2)==(76,76,(4,6))
    geometry=DiskGeometry(38.1,37.8)
    v.geometry=geometry;v.update_derived()
    assert '78 × 76' in v.derived.text()
    assert v.geometry is geometry and (geometry.radius_x_px,geometry.radius_y_px)==(38.1,37.8)
    assert derive_preparation_dimensions(geometry,v.roi,3)==(78,78,(4,6))
    for invalid in (0,-1,True,2.5):
        with pytest.raises((ValueError,TypeError)):derive_preparation_dimensions(geometry,v.roi,invalid)


def test_arbitrary_roi_preview_boundary_and_replay(window,tmp_path):
    from qcbed.gui.preparation_model import finalize_symmetry_roi
    from qcbed.data_preparation.models import DiskGeometry
    from qcbed.data_preparation.origin_selection import prepare_roi,SelectableBeam
    roi=finalize_symmetry_roi((1.2,2.3,8.1,11.4),(30,30))
    assert roi==(4,5,12,15)
    assert finalize_symmetry_roi(roi,(30,30))==roi
    # All trial crops, not merely the overlay, stay in the image.
    beam=SelectableBeam('center',(0,0),(15.,15.),'reference')
    prepare_roi(roi,beam,DiskGeometry(30.,30.),(30,30),crop_half_size_px=(4,5))
    with pytest.raises(ValueError,match='slightly smaller region'):
        finalize_symmetry_roi((0,0,29,29),(30,30))
    v=window.preparation;v.image=scene();v.geometry=DiskGeometry(20.,20.)
    v.tool.setCurrentIndex(1);drag(v,(87.2,87.1),(92.2,92.6))
    assert v.roi is None
    v.accept_preview();assert v.roi==(87,87,93,93)
    image=tmp_path/'p.npy';np.save(image,scene())
    from dataclasses import replace
    settings=replace(selection(),symmetry_roi=v.roi,symmetry_crop_half_size_px=(3,3))
    descriptor=prepare_files(image,settings,tmp_path/'prep')
    v.install(descriptor)
    assert v.roi==(87,87,93,93) and '(3, 3)' in v.derived.text()
    assert 'A8' not in v.derived.text() and 'A8' not in v.status.text()


def test_vertical_inspector_weight_controls_and_persistence(window,tmp_path):
    from qcbed.gui.document import EditorDocument
    path=tmp_path/'raw.npy';np.save(path,scene())
    descriptor=prepare_files(path,selection(),tmp_path/'prep')
    window.document=EditorDocument(modern_config(tmp_path),tmp_path/'project.yml')
    window.refresh();window.preparation.install(descriptor);window.show()
    app=QApplication.instance();app.processEvents()
    v=window.data_panel
    assert window.navigation.currentItem()==window.section_items['Data Prep']
    assert window.preparation_scroll.maximumHeight()>600
    assert [b.text() for b in v.mask_buttons]==['Exclude','Restore Selected','Reset Manual Mask']
    assert [ax.get_title() for ax in v.axes]==['Pattern','Sigma','Final Weight']
    assert v.axes[2].images[0].get_interpolation() == 'nearest'
    assert v.axes[2].images[0].get_clim() == (0., 1.)
    assert window.document.config.experiment.path in v.summary.text()
    assert all(ax.get_subplotspec().get_gridspec().get_geometry()==(3,1) for ax in v.axes)
    target=v.overview[1]['beam_mapping'][1]['center_px']
    v.pick_overview(SimpleNamespace(xdata=target[0],ydata=target[1]))
    assert v.beams.currentIndex()==1
    before=read_data(window.document.config,tmp_path)
    for ax,array in zip(v.axes,(before.intensities,np.sqrt(before.variances),before.weights)):
        np.testing.assert_array_equal(np.asarray(ax.images[0].get_array()).ravel(),array[:,1])
    layers={}
    for label in ('Geometric Mask','Processing Validity','Manual Mask','Final Weight'):
        v.mask_layer.setCurrentText(label)
        layers[label]=np.array(v.axes[2].images[0].get_array())
    np.testing.assert_array_equal(layers['Geometric Mask']*layers['Processing Validity']*
        layers['Manual Mask'],layers['Final Weight'])
    index=int(np.flatnonzero(before.weights[:,1])[0]);y,x=divmod(index,before.scan_columns)
    def click_pixel():
        panel=window.data_panel
        panel.pick_pixel(SimpleNamespace(inaxes=panel.axes[0],xdata=x,ydata=y))
        return panel
    click_pixel().mask_buttons[0].click()
    after=read_data(window.document.config,tmp_path)
    np.testing.assert_array_equal(after.intensities,before.intensities)
    np.testing.assert_array_equal(after.variances,before.variances)
    assert after.weights[index,1]==0
    v=window.data_panel
    assert np.asarray(v.axes[0].images[0].get_array())[y,x]==before.intensities[index,1]
    assert np.asarray(v.axes[2].images[0].get_array())[y,x]==0
    assert f'Active points: {np.count_nonzero(after.weights)}' in v.summary.text()
    assert tuple(spin.value() for spin in v.box)==(x,y,x+1,y+1)
    v.mask_buttons[1].click()
    np.testing.assert_array_equal(read_data(window.document.config,tmp_path).weights,before.weights)
    click_pixel().mask_buttons[0].click();window.data_panel.mask_buttons[2].click()
    np.testing.assert_array_equal(read_data(window.document.config,tmp_path).weights,before.weights)
    click_pixel().mask_buttons[0].click();window.data_panel.select_zone()
    assert window.validate_stage(0)
    snapshot=window.workflow.preparation_descriptor
    assert window.save_document(tmp_path/'saved.yml')
    assert window.save_document(tmp_path/'save-as.yml')
    window.open_document(tmp_path/'save-as.yml');window.preparation.install(str(snapshot))
    restored=read_data(window.document.config,tmp_path)
    np.testing.assert_array_equal(restored.weights,after.weights)
    np.testing.assert_array_equal(restored.intensities,before.intensities)
    np.testing.assert_array_equal(restored.variances,before.variances)


def test_large_resizable_workspace_and_navigation(window):
    from PySide6.QtCore import Qt
    window.resize(1400,1100);window.show();QApplication.instance().processEvents()
    v=window.preparation
    assert window.preparation_split.orientation()==Qt.Orientation.Vertical
    assert window.preparation_scroll.maximumHeight()>600
    assert v.canvas.minimumWidth()>=520 and v.canvas.minimumHeight()>=520
    assert v.controls.maximumWidth()<=380
    assert v.toolbar is not None
    v.toolbar.pan();drag(v,(10,10),(40,40))
    assert v.preview is None and v.press is None
    v.toolbar.pan();drag(v,(10,10),(40,40));assert v.preview is not None
    window.preparation_split.setSizes([650,200]);QApplication.instance().processEvents()
    assert window.preparation_split.widget(1)==window.editor_split
    window.workflow_tabs.setCurrentIndex(1);window.workflow_tabs.setCurrentIndex(0)
    assert window.editor_split.parent()==window.preparation_split


def test_preparation_entry_points_without_internal_attachment(window,tmp_path):
    from PySide6.QtWidgets import QPushButton,QLabel
    from qcbed.gui.document import EditorDocument
    window.document=EditorDocument(modern_config(tmp_path),tmp_path/'project.yml');window.refresh()
    assert not hasattr(window.data_panel,'attach')
    assert window.preparation.replay.text()=='Load Previous Preparation…'
    assert window.preparation.mode.itemText(0)=='New Data Preparation'
    for widget in window.findChildren(QPushButton)+window.findChildren(QLabel):
        assert 'patch a' not in widget.text().lower()
        assert 'manifest' not in widget.text().lower()
    assert 'Load Previous Preparation' in window.data_panel.overview_figure.axes[0].texts[0].get_text()


def test_modern_has_no_separate_experimental_panel(window,tmp_path):
    from qcbed.gui.document import EditorDocument,editable_fields
    window.document=EditorDocument(modern_config(tmp_path),tmp_path/'project.yml');window.refresh()
    assert [window.navigation.topLevelItem(i).text(0) for i in range(3)]==['A. Data Prep','B. Atomic Model','C. Main QCBED Parameters']
    assert 'Experimental data' not in window.section_items
    assert window.section_items['Scan'].parent().text(0)=='A. Data Prep'
    assert ('experiment','path') not in window.widgets
    for key in (('microscope','energy_kv'),('experiment','gradient_option')):
        spec=next(f for f in editable_fields(window.document.config) if f.key==key)
        assert window.section_items[spec.section].parent().text(0)=='C. Main QCBED Parameters'
    from qcbed import load_yaml_config
    legacy=load_yaml_config(Path(__file__).resolve().parents[1]/'examples/qcbed_config.yml')
    window.document=EditorDocument(legacy,tmp_path/'legacy.yml');window.refresh()
    assert 'Experimental data' in window.section_items
    assert ('experiment','path') in window.widgets


def test_new_preparation_detaches_previous_arrays_and_geometry(window,tmp_path):
    from qcbed.gui.document import EditorDocument
    window.document=EditorDocument(modern_config(tmp_path),tmp_path/'old.yml')
    window.refresh()
    p=window.preparation
    p.mode.setCurrentIndex(1)
    p.geometry=object();p.roi=(1,2,3,4);p.preview=('old',)
    p.mode.setCurrentIndex(0)
    assert window.document is None and window.data_panel is None
    assert p.geometry is None and p.roi is None and p.preview is None
    assert p.raw_document is None and p.raw_path is None
    assert window.workflow.snapshot is None and window.workflow.preparation_descriptor is None
    assert window.workflow.receipts==[None,None,None]
    assert not window.run_action.isEnabled()
    assert 'scan' not in p.retained_bc and 'experiment' not in p.retained_bc
    assert p.retained_bc['crystallography'].beams==()


def test_stage_local_atomic_validation_without_c_boundaries(window,tmp_path):
    from dataclasses import replace
    from qcbed.gui.document import EditorDocument
    from qcbed.gui.workflow import keys
    image=tmp_path/'raw.npy';np.save(image,scene())
    descriptor=prepare_files(image,selection(),tmp_path/'prep')
    c=modern_config(tmp_path)
    c=replace(c,crystallography=replace(c.crystallography,slice_boundaries_fractional=()))
    window.document=EditorDocument(c,tmp_path/'project.yml')
    window.refresh();window.preparation.install(descriptor)
    window.data_panel.select_zone()
    assert window.validate_stage(0)
    a_receipt=window.workflow.receipts[0]
    window.load_structure_path(str(Path(__file__).resolve().parents[1]/'reference_data/Al_FCC.cif'))
    assert window.site_group_table.item(0,0).text()=='Al_fcc'
    window.site_group_table.item(0,3).setText('0.86')
    window.apply_site_groups.click()
    assert window.validate_stage(1),window.stage_status[1].text()
    assert window.stage_status[1].text()=='B. Atomic Model: Valid'
    assert window.stage_buttons[2].isEnabled()
    b_receipt=window.workflow.receipts[1]
    assert not window.validate_c_part('simulation')
    assert window.stage_status[2].text().startswith('Slice boundaries must start at 0, end at 1, and contain at least 2 strictly increasing values.')
    assert window.workflow.receipts[:2]==[a_receipt,b_receipt]
    # Numerical C changes cannot invalidate the atomic-model receipt.
    edited=replace(window.document.config,crystallography=replace(window.document.config.crystallography,
        slice_boundaries_fractional=(0.,.5,1.),depth_angstrom=7.))
    window.workflow.sync(edited)
    assert window.workflow.receipts[1]==b_receipt
