"""Focused Qt structural transactions, modal additions and explicit deletions."""
import os
os.environ.setdefault('QT_QPA_PLATFORM','offscreen')
import pytest
pytest.importorskip('PySide6')
from PySide6.QtCore import QTimer, Qt
from PySide6.QtWidgets import QApplication, QDialog, QDialogButtonBox, QMessageBox, QTableWidget
from qcbed import Id, get_object, available_refinement_targets, load_yaml_config
from qcbed.configuration import ParameterCategory as C
from qcbed.gui.structure_widgets import StructureDialog, RecordForm
from test_gui_qt import app, window
from test_gui_structure import document, region_path


@pytest.fixture
def dialog(app,document):
    dialog = StructureDialog(document)
    dialog.show()
    app.processEvents()
    yield dialog
    dialog.reject()
    dialog.deleteLater()
    app.processEvents()


def add_record(dialog,path,changes):
    dialog.select_path(path)
    errors = []
    def timeout():
        modal = QApplication.activeModalWidget()
        if modal is not None:
            from PySide6.QtWidgets import QLabel
            errors.append([label.text() for label in modal.findChildren(QLabel)])
            modal.reject()
    def fill():
        modal = QApplication.activeModalWidget()
        assert modal is not None
        form = modal.findChild(RecordForm)
        for name,text in changes.items():
            widget = form.controls[name]
            if hasattr(widget,'setText'): widget.setText(text)
            else: widget.setCurrentText(text)
        modal.findChild(QDialogButtonBox).button(QDialogButtonBox.StandardButton.Ok).click()
    QTimer.singleShot(0,fill)
    watchdog = QTimer()
    watchdog.setSingleShot(True)
    watchdog.timeout.connect(timeout)
    watchdog.start(1000)
    dialog.buttons['Add'].click()
    watchdog.stop()
    assert not errors, errors


def test_species_add_thermal_reorder_remove(dialog,document,monkeypatch):
    old = document.config
    add_record(dialog,('sample','species'),{'label':'Cr','atomic_number':'24'})
    species = dialog.draft.config.sample.species[-1]
    assert species.label == 'Cr'
    assert document.config == old
    assert 'thermal' in dialog.issues.text().lower() or 'debye' in dialog.issues.text().lower()
    add_record(dialog,region_path(old)+('debye_waller',),{'value_angstrom_squared':'.005123456789012345'})
    dialog.select_path(('sample','species',Id(species.id)))
    dialog.move(-1)
    dialog.commit()
    assert document.config.sample.species[0].id == species.id
    assert sum(t.category == C.DEBYE_WALLER for t in available_refinement_targets(document.config)) == 2
    second = StructureDialog(document)
    second.select_path(('sample','species',Id(species.id)))
    prompts = []
    def yes(*args): prompts.append(args[2]); return QMessageBox.StandardButton.Yes
    monkeypatch.setattr(QMessageBox,'question',yes)
    second.delete()
    assert 'Thermal entry uses this species' in prompts[0]
    second.commit()
    assert document.config == old
    second.deleteLater()


def test_region_add_duplicate_slice_order_atoms_and_reflections(dialog,document,tmp_path):
    original = document.config
    add_record(dialog,('sample','regions'),{'depth_angstrom':'7.25','symmetry':'INVE'})
    region = dialog.draft.config.sample.regions[-1]
    path = ('sample','regions',Id(region.id))
    assert region.depth_angstrom == 7.25
    add_record(dialog,path+('slices',),{'record_type':'S'})
    slices = get_object(dialog.draft.config,path+('slices',))
    last_path = path+('slices',Id(slices[-1].id))
    dialog.select_path(last_path+('geometry',))
    dialog.form.controls['thickness_angstrom'].setText('3.125')
    dialog.stage()
    dialog.select_path(last_path)
    dialog.move(-1)
    assert get_object(dialog.draft.config,path).slices[0].id == slices[-1].id
    add_record(dialog,last_path+('atoms',),{'x_fractional':'-0.125','y_fractional':'1.25','starting_occupancy':'.8123456789012345'})
    atom = get_object(dialog.draft.config,last_path).atoms[0]
    assert atom.y_fractional == 1.25
    dialog.select_path(last_path+('atoms',Id(atom.id)))
    dialog.duplicate()
    atoms = get_object(dialog.draft.config,last_path).atoms
    assert atoms[0].id != atoms[1].id
    add_record(dialog,path+('reflections',),{'h':'1','k':'1','elastic':'.99','absorptive':'.95'})
    reflection = get_object(dialog.draft.config,path).reflections[0]
    dialog.select_path(path+('reflections',Id(reflection.id)))
    dialog.duplicate()
    before_sort = dialog.draft.config
    dialog.select_path(path+('reflections',))
    table = dialog.scroll.widget()
    assert isinstance(table,QTableWidget)
    table.sortItems(0,Qt.SortOrder.DescendingOrder)
    assert dialog.draft.config == before_sort
    dialog.select_path(path)
    dialog.move(-1)
    dialog.commit()
    assert document.config.sample.regions[1] == original.sample.regions[0]
    assert len(document.config.sample.regions[0].reflections) == 2
    document.save(tmp_path/'complex.yml')
    assert load_yaml_config(tmp_path/'complex.yml') == document.config


def test_invalid_text_history_and_undo_redo(dialog,document,window):
    original = document.config
    path = region_path(original)
    dialog.select_path(path)
    dialog.form.controls['depth_angstrom'].setText('incomplete-')
    dialog.commit()
    assert document.config == original
    assert not document.undo_stack
    assert 'incomplete' in dialog.form.controls['depth_angstrom'].text()
    dialog.form.controls['depth_angstrom'].setText('8.5')
    dialog.commit()
    edited = document.config
    window.document = document
    window.refresh()
    window.undo()
    assert document.config == original
    window.redo()
    assert document.config == edited
    window.controller.process = object()  # The same active-state guard used by real workers.
    window.undo()
    assert document.config == edited
    window.controller.process = None


def test_absorption_linking_and_bounds(dialog,document,monkeypatch):
    from qcbed.configuration import ParameterBinding,ParameterTarget
    absorption = ('sample','absorption')
    dialog.select_path(absorption)
    dialog.form.controls['enabled'].setChecked(True)
    dialog.stage()
    add_record(dialog,absorption+('coefficients',),{'a':'.1','b':'.2','c':'.3'})
    species = document.config.sample.species[0].id
    target = ParameterTarget(C.ABSORPTION,(species,'c'))
    dialog.draft.add(('refinement','parameters'),ParameterBinding(target,(1.1,.9),.01))
    dialog.rebuild()
    dialog.select_path(absorption+('coefficients',Id(species,'species')))
    dialog.form.controls['c_mode'].setCurrentText('linked_to_b')
    monkeypatch.setattr(QMessageBox,'question',lambda *args: QMessageBox.StandardButton.Yes)
    dialog.stage()
    assert not dialog.form.controls['c'].isEnabled()
    assert not dialog.draft.config.refinement.parameters
    add_record(dialog,('refinement','bounds','absolute'),{'lower':'1','upper':'12'})
    add_record(dialog,('refinement','bounds','relative'),{'fraction':'.2','absolute_fallback':'.05'})
    region = document.config.sample.regions[0]
    add_record(dialog,('refinement','bounds','specific'),{'target':f'Region {region.id} → depth','lower':'2','upper':'8'})
    dialog.select_path(('refinement','bounds'))
    dialog.form.controls['policy'].setCurrentText('reject')
    dialog.stage()
    dialog.commit()
    coefficient = document.config.sample.absorption.coefficients[0]
    assert coefficient.c_mode == 'linked_to_b'
    assert coefficient.c == .3
    assert document.config.refinement.bounds.absolute[0].upper == 12.
    assert document.config.refinement.bounds.relative[0].fraction == .2
    assert document.config.refinement.bounds.specific[0].lower == 2.
