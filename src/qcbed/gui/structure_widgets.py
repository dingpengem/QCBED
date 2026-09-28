"""Qt structural transaction editor; drafts never replace a valid document early."""
from dataclasses import fields, is_dataclass, replace
from enum import Enum
from types import UnionType
from typing import get_args, get_origin, get_type_hints, Literal
from uuid import uuid4

from PySide6.QtCore import Qt, QSignalBlocker
from PySide6.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QFormLayout,
    QTreeWidget, QTreeWidgetItem, QWidget, QSplitter, QPushButton, QLabel,
    QLineEdit, QComboBox, QScrollArea, QMessageBox, QDialogButtonBox, QCheckBox, QTableWidget, QTableWidgetItem)
from qcbed import Id, get_object, available_refinement_targets, validate_fields, validate_config
from qcbed.configuration import (SpeciesConfig, RegionConfig, StoredSliceConfig, AtomConfig,
    ReflectionConfig, ThermalFactorConfig, AbsorptionSpeciesConfig, TargetBound,
    CategoryBound, RelativeBound, ParameterCategory)
from qcbed.models import SliceGeometry
from .structure import StructuralDraft, selector, removal_plan
from .document import display_value, parse_text, target_label


class RecordForm(QWidget):
    """Typed scalar fields only; nested objects remain navigable records."""
    def __init__(self, config, record, path, catalogue, parent=None):
        super().__init__(parent)
        self.record, self.path = record, path
        self.controls = {}
        self.parsers = {}
        self.dirty = False
        layout = QFormLayout(self)
        hints = get_type_hints(type(record))
        for field in fields(record):
            name = field.name
            value = getattr(record, name)
            annotation = hints[name]
            if name.endswith('_raw') or name == 'legacy_raw': continue
            if is_dataclass(value) and name != 'target': continue
            if isinstance(value, tuple) and (not value or is_dataclass(value[0])): continue
            if name == 'id':
                layout.addRow('Stable ID (read-only)', QLabel(value))
                continue
            choices = None
            if name == 'species': choices = [(s.id, s.id) for s in config.sample.species]
            elif name == 'target': choices = [(target_label(config, t), t.target) for t in catalogue]
            elif get_origin(annotation) is Literal: choices = [(str(v), v) for v in get_args(annotation)]
            elif isinstance(annotation, type) and issubclass(annotation, Enum): choices = [(v.value, v) for v in annotation]
            if choices is not None:
                if not any(v == value for _,v in choices): choices.append((str(value), value))
                control = QComboBox()
                for label, choice in choices: control.addItem(label, choice)
                control.setCurrentIndex(next(i for i,(_,v) in enumerate(choices) if v == value))
                parser = lambda c=control, options=choices: options[c.currentIndex()][1]
                control.currentIndexChanged.connect(self.mark_dirty)
            elif annotation is bool:
                control = QCheckBox()
                control.setChecked(value)
                parser = lambda c=control: c.isChecked()
                control.toggled.connect(self.mark_dirty)
            else:
                control = QLineEdit(display_value(value))
                def parse(c=control, a=annotation):
                    text = c.text()
                    if get_origin(a) is UnionType:
                        if not text.strip(): return None
                        a = next(t for t in get_args(a) if t is not type(None))
                    if get_origin(a) is tuple:
                        import re
                        element = get_args(a)[0]
                        return tuple(parse_text(v, 'int' if element is int else 'real')
                                     for v in re.split(r'[\s,]+',text.strip()))
                    return parse_text(text, 'int' if a is int else 'real' if a is float else 'text')
                parser = parse
                control.textChanged.connect(self.mark_dirty)
            if name == 'c' and getattr(record, 'c_mode', None) == 'linked_to_b':
                control.setToolTip('Stored C is retained; linked mode derives its physical value from B.')
                control.setEnabled(False)
            self.controls[name] = control
            self.parsers[name] = parser
            layout.addRow(name.replace('_', ' '), control)

    def mark_dirty(self, *args): self.dirty = True

    def changes(self):
        values = {name: parse() for name,parse in self.parsers.items()}
        supplied = {field.name: values.get(field.name, getattr(self.record,field.name)) for field in fields(self.record)}
        issues = validate_fields(type(self.record), supplied, section=self.path[0] if self.path else 'sample', object_path=self.path)
        if issues: raise ValueError('\n'.join(f'{i.field}: {i.message}' for i in issues))
        return values


class StructureDialog(QDialog):
    def __init__(self, document, parent=None):
        super().__init__(parent)
        self.document = document
        self.draft = StructuralDraft(document.config)
        self.catalogue = available_refinement_targets(document.config)
        self.form = None
        self.path = None
        self.setWindowTitle('Scientific structure — transactional edits')
        self.resize(1100, 800)
        layout = QVBoxLayout(self)
        note = QLabel('Stage fields to update this draft. Commit validates the entire transaction; Cancel discards it. '
                      'Up/Down explicitly changes scientific order. Tree sorting is display only. IDs cannot be renamed.')
        note.setWordWrap(True)
        note.setMaximumHeight(65)
        layout.addWidget(note)
        split = QSplitter()
        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(['Scientific record', 'Identity / order'])
        self.tree.setSortingEnabled(True)
        self.tree.currentItemChanged.connect(self.select)
        split.addWidget(self.tree)
        self.scroll = QScrollArea()
        self.scroll.setWidgetResizable(True)
        split.addWidget(self.scroll)
        layout.addWidget(split, 1)
        split.setSizes([430,650])
        self.tree.setColumnWidth(0,280)
        self.tree.setColumnWidth(1,140)
        self.tree.sortItems(0,Qt.SortOrder.AscendingOrder)
        row = QHBoxLayout()
        self.buttons = {}
        for label, callback in [('Stage fields', self.stage), ('Add', self.add), ('Duplicate', self.duplicate),
                                ('Delete', self.delete), ('Scientific up', lambda:self.move(-1)),
                                ('Scientific down', lambda:self.move(1))]:
            button = QPushButton(label)
            button.clicked.connect(lambda checked=False, cb=callback: self.invoke(cb))
            row.addWidget(button)
            self.buttons[label] = button
        layout.addLayout(row)
        self.issues = QLabel()
        self.issues.setWordWrap(True)
        self.issues.setMaximumHeight(150)
        self.issues.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(self.issues)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Save | QDialogButtonBox.StandardButton.Cancel)
        buttons.button(QDialogButtonBox.StandardButton.Save).setText('Commit transaction')
        buttons.accepted.connect(self.commit)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self.rebuild()

    def invoke(self, callback):
        try: callback()
        except (ValueError, TypeError) as exc: self.issues.setText(str(exc))

    def rebuild(self, select_path=None):
        self.form = None
        with QSignalBlocker(self.tree):
            self.tree.clear()
            def node(parent, name, value, path):
                identity = getattr(value,'id',getattr(value,'species',''))
                item = QTreeWidgetItem(parent, [name, identity if isinstance(identity,str) else ''])
                item.setData(0, Qt.ItemDataRole.UserRole, path)
                if is_dataclass(value):
                    for field in fields(value):
                        child = getattr(value, field.name)
                        if is_dataclass(child) and field.name != 'target' or isinstance(child, tuple) and (not child or is_dataclass(child[0])):
                            node(item, field.name, child, path+(field.name,))
                elif isinstance(value, tuple):
                    for i, child in enumerate(value):
                        if is_dataclass(child): node(item, f'{i+1:04d} · {getattr(child,"id",getattr(child,"species",str(getattr(child,"target",getattr(child,"category","record")))))}',child,path+(selector(child),))
                if path == select_path: self.tree.setCurrentItem(item)
            node(self.tree, 'Sample / species / regions', self.draft.config.sample, ('sample',))
            node(self.tree, 'Bounds', self.draft.config.refinement.bounds, ('refinement','bounds'))
            node(self.tree, 'Bindings (dependency resolution)', self.draft.config.refinement.parameters, ('refinement','parameters'))
            self.tree.expandToDepth(1)
        self.select(self.tree.currentItem())
        self.validate_draft()

    def validate_draft(self):
        issues = validate_config(self.draft.config)
        self.issues.setText('\n'.join(f'{i.phase.value} · {i.section} · {i.object_path} · {i.field}: {i.message}' for i in issues) or 'Draft valid. Commit creates one undoable snapshot.')
        try: self.catalogue = available_refinement_targets(self.draft.config)
        except (ValueError, TypeError): pass  # Keep last valid catalogue during incomplete transactions.

    def select(self, item, previous=None):
        if self.form and self.form.dirty:
            try: self.stage(rebuild=False)
            except (ValueError, TypeError) as exc:
                self.issues.setText(str(exc))
                with QSignalBlocker(self.tree): self.tree.setCurrentItem(previous)
                return
        self.path = item.data(0,Qt.ItemDataRole.UserRole) if item else None
        value = get_object(self.draft.config,self.path) if self.path else None
        self.form = RecordForm(self.draft.config,value,self.path,self.catalogue) if is_dataclass(value) else None
        old = self.scroll.takeWidget()
        if old: old.deleteLater()
        if isinstance(value, tuple) and value:
            records = value
            columns = [f.name for f in fields(records[0]) if not isinstance(getattr(records[0],f.name),tuple) and not is_dataclass(getattr(records[0],f.name))]
            table = QTableWidget(len(records),len(columns))
            table.setHorizontalHeaderLabels(columns)
            for row, record in enumerate(records):
                for col, field in enumerate(columns):
                    cell = QTableWidgetItem(display_value(getattr(record,field)))
                    cell.setFlags(cell.flags() & ~Qt.ItemFlag.ItemIsEditable)
                    cell.setData(Qt.ItemDataRole.UserRole,self.path+(selector(record),))
                    table.setItem(row,col,cell)
            table.horizontalHeader().setSortIndicator(-1,Qt.SortOrder.AscendingOrder)
            table.setSortingEnabled(True)
            table.resizeColumnsToContents()
            table.cellDoubleClicked.connect(lambda row,col: self.select_path(table.item(row,col).data(Qt.ItemDataRole.UserRole)))
            self.scroll.setWidget(table)
        else:
            self.scroll.setWidget(self.form or QLabel('Select a record to edit, or an ordered collection to add a record.'))
        for name in ('Duplicate', 'Scientific up', 'Scientific down'):
            self.buttons[name].setEnabled(hasattr(value, 'id'))
        self.buttons['Delete'].setEnabled(bool(self.path and isinstance(self.path[-1],Id)))
        self.buttons['Add'].setEnabled(isinstance(value,tuple))

    def select_path(self, path):
        from PySide6.QtWidgets import QTreeWidgetItemIterator
        iterator = QTreeWidgetItemIterator(self.tree)
        while iterator.value():
            item = iterator.value()
            if item.data(0,Qt.ItemDataRole.UserRole) == path:
                self.tree.setCurrentItem(item)
                return
            iterator += 1

    def stage(self, *, rebuild=True):
        if self.form and self.form.dirty:
            changes = self.form.changes()
            candidate = StructuralDraft(self.draft.config)
            if changes.get('c_mode') == 'linked_to_b':
                from qcbed import set_refinement
                from qcbed.configuration import ParameterTarget
                target = ParameterTarget(ParameterCategory.ABSORPTION, (changes['species'],'c'))
                if any(b.target == target for b in candidate.config.refinement.parameters):
                    if QMessageBox.question(self,'Linked C dependency',
                            f'Link C to B and explicitly disable the independent C refinement binding for {changes["species"]}?',
                            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                            QMessageBox.StandardButton.No) != QMessageBox.StandardButton.Yes:
                        raise ValueError('Linking cancelled: disable the independent C binding or keep C independent.')
                    candidate.config = set_refinement(candidate.config,target,enabled=False)
            candidate.update(self.path, **changes)
            self.draft.config = candidate.config
            self.form.dirty = False
        if rebuild: self.rebuild(self.path)

    def default_record(self, path):
        fresh = lambda prefix: prefix+'-'+uuid4().hex[:10]
        config = self.draft.config
        species = config.sample.species[0].id if config.sample.species else ''
        kind = path[-1]
        if kind == 'species': return SpeciesConfig(fresh('species'),'New species',1)
        if kind == 'regions': return RegionConfig(fresh('region'),1.,
            (StoredSliceConfig(fresh('slice'),SliceGeometry(1.,cell_a_angstrom=1.,cell_b_angstrom=1.,lattice_gamma_degrees=90.),()),),
            tuple(ThermalFactorConfig(s.id,0.) for s in config.sample.species))
        if kind == 'slices': return StoredSliceConfig(fresh('slice'),SliceGeometry(1.,cell_a_angstrom=1.,cell_b_angstrom=1.,lattice_gamma_degrees=90.),())
        if kind == 'atoms': return AtomConfig(fresh('atom'),species,0.,0.,1.)
        if kind == 'reflections': return ReflectionConfig(fresh('reflection'),0,0,1.,1.)
        if kind in ('debye_waller','coefficients'):
            present = {v.species for v in get_object(config,path)}
            species = next((s.id for s in config.sample.species if s.id not in present),species)
            return ThermalFactorConfig(species,0.) if kind == 'debye_waller' else AbsorptionSpeciesConfig(species,0.,0.,0.)
        if kind == 'specific':
            if not self.catalogue: raise ValueError('No available physical targets')
            return TargetBound(self.catalogue[0].target)
        if kind == 'absolute': return CategoryBound(ParameterCategory.DEPTH)
        if kind == 'relative': return RelativeBound(ParameterCategory.DEPTH,.1,.1)
        raise ValueError('Enable new bindings using the main refinement target catalogue after committing structure.')

    def add(self):
        record = self.default_record(self.path)
        dialog = QDialog(self)
        dialog.setWindowTitle('Add record — review starting values')
        layout = QVBoxLayout(dialog)
        form = RecordForm(self.draft.config,record,self.path,self.catalogue)
        layout.addWidget(form)
        preview = QLabel(repr(record))
        preview.setWordWrap(True)
        layout.addWidget(preview)
        if isinstance(record,RegionConfig): layout.addWidget(QLabel('Creates the displayed unit-cell slice and zero thermal entry for each species; edit these in the draft before committing.'))
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        layout.addWidget(buttons)
        def accept():
            try:
                value = replace(record,**form.changes())
                self.draft.add(self.path,value)
            except (ValueError,TypeError) as exc:
                preview.setText(str(exc)); return
            dialog.accept()
        buttons.accepted.connect(accept)
        buttons.rejected.connect(dialog.reject)
        if dialog.exec() == QDialog.DialogCode.Accepted: self.rebuild(self.path)

    def duplicate(self):
        self.stage(rebuild=False)
        self.draft.duplicate(self.path)
        self.rebuild(self.path[:-1])

    def delete(self):
        self.stage(rebuild=False)
        plan = removal_plan(self.draft.config,self.path)
        value = get_object(self.draft.config,self.path)
        message = 'Remove this record and its contained subtree?\n'+repr(value)+'\n\nExternal references to remove explicitly:\n'
        message += '\n'.join(f'{d.object_path} / {d.field}: {d.message}' for d in plan) or 'None'
        if QMessageBox.question(self,'Confirm scientific deletion',message,
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No) != QMessageBox.StandardButton.Yes: return
        self.draft.delete(self.path,remove_dependencies=True)
        self.rebuild(self.path[:-1])

    def move(self, offset):
        self.stage(rebuild=False)
        self.draft.move(self.path,offset)
        self.rebuild(self.path)

    def commit(self):
        try:
            self.stage(rebuild=False)
            if not self.document.commit_structure(self.draft.config):
                self.validate_draft(); return
        except (ValueError,TypeError) as exc:
            self.issues.setText(str(exc)); return
        self.accept()
