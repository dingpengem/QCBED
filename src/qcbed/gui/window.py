"""Optional Qt Widgets presentation of the frontend-independent editor document."""
from dataclasses import fields, is_dataclass
from pathlib import Path

from PySide6.QtCore import Qt, QTimer
import json
import time
from PySide6.QtGui import QAction
from PySide6.QtWidgets import (QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QFormLayout, QLineEdit, QComboBox, QPushButton, QTreeWidget, QTreeWidgetItem,
    QStackedWidget, QSplitter, QPlainTextEdit, QTableWidget, QTableWidgetItem,
    QLabel, QScrollArea, QFileDialog, QMessageBox, QTabWidget)

from qcbed import available_refinement_targets, get_object
from .document import EditorDocument, editable_fields, display_value, target_label, diagnostic_presentation


from .run_controller import RunController, readiness


class ConfigWindow(QMainWindow):
    def __init__(self, path=None):
        super().__init__()
        self.controller = RunController()
        self.document = None
        self.pending = False
        self.widgets = {}
        self.resize(1250, 850)
        toolbar = self.addToolBar('Configuration')
        self.edit_actions = []
        for title, callback in [('Open', self.choose_open), ('Save', self.save_current),
                ('Save As', self.choose_save_as), ('Apply edits', self.apply_edits),
                ('Validate', self.validate_document), ('Revert drafts', self.refresh),
                ('Structure / Bounds', self.edit_structure), ('Undo', self.undo), ('Redo', self.redo)]:
            action = QAction(title, self)
            action.triggered.connect(lambda checked=False, cb=callback: cb())
            toolbar.addAction(action)
            self.edit_actions.append(action)
        legacy_menu = self.menuBar().addMenu('Compatibility / Legacy')
        for label, callback in [('Import Legacy Project…', self.import_legacy),
                ('Export FOR007…', lambda: self.export_legacy(False)),
                ('Export Legacy Bundle…', lambda: self.export_legacy(True)),
                ('View Compatibility Report…', self.view_compatibility)]:
            action = legacy_menu.addAction(label)
            action.triggered.connect(lambda checked=False, cb=callback: cb())
            self.edit_actions.append(action)
        self.run_action = toolbar.addAction('Run', lambda: self.start_run(False))
        self.resume_action = toolbar.addAction('Resume', lambda: self.start_run(True))
        self.cancel_action = toolbar.addAction('Cancel', self.cancel_run)
        self.cancel_action.setEnabled(False)
        self.run_action.setEnabled(False)
        self.resume_action.setEnabled(False)
        self.navigation = QTreeWidget()
        self.navigation.setHeaderLabel('Configuration')
        self.pages = QStackedWidget()
        self.navigation.currentItemChanged.connect(self.navigate)
        top = QSplitter()
        top.addWidget(self.navigation)
        top.addWidget(self.pages)
        top.setStretchFactor(1, 1)
        self.issues = QTreeWidget()
        self.issues.setHeaderLabels(['Severity', 'Phase', 'Location', 'Message'])
        self.issues.itemClicked.connect(self.navigate_issue)
        self.summary = QPlainTextEdit()
        self.summary.setReadOnly(True)
        bottom = QSplitter()
        bottom.addWidget(self.issues)
        bottom.addWidget(self.summary)
        split = QSplitter(Qt.Orientation.Vertical)
        split.addWidget(top)
        split.addWidget(bottom)
        self.run_status = QLabel('Idle — Validate to check Run and Resume readiness')
        self.run_status.setWordWrap(True)
        self.run_log = QPlainTextEdit()
        self.run_log.setReadOnly(True)
        self.run_log.document().setMaximumBlockCount(3000)
        run_panel = QWidget()
        layout = QVBoxLayout(run_panel)
        layout.addWidget(self.run_status)
        layout.addWidget(self.run_log)
        split.addWidget(run_panel)
        split.setSizes([450, 220, 160])
        self.timer = QTimer(self)
        self.timer.setInterval(100)
        self.timer.timeout.connect(self.poll_run)
        from .results_view import ResultsView
        self.results_view = ResultsView()
        self.workspace = QTabWidget()
        self.workspace.addTab(split, 'Configuration')
        self.workspace.addTab(self.results_view, 'Results')
        self.setCentralWidget(self.workspace)
        self.results_directory = None
        self.follow_current_results = True
        self.last_results_refresh = 0.
        results_menu = self.menuBar().addMenu('Results')
        self.results_actions = {}
        for title, callback in (
                ('Show Results', self.show_results),
                ('Open Existing Run…', self.choose_results),
                ('Open HTML Report', lambda: self.results_view.open_artifact('report.html')),
                ('Open Output Folder', lambda: self.results_view.open_artifact(''))):
            self.results_actions[title] = results_menu.addAction(title, callback)
        self.results_actions['Open HTML Report'].setEnabled(False)
        self.results_actions['Open Output Folder'].setEnabled(False)
        follow_button = QPushButton('Show Current / Latest Run')
        follow_button.clicked.connect(self.follow_run_results)
        self.results_view.layout().insertWidget(0, follow_button)
        self.statusBar().showMessage('Open a canonical YAML configuration. Validate before Run or Resume.')
        self.setWindowTitle('QCBED configuration editor')
        if path is not None:
            self.open_document(path)

    def choose_results(self):
        path = QFileDialog.getExistingDirectory(self, 'Open canonical run output directory')
        if path:
            self.open_results(path)

    def show_results(self):
        self.refresh_results(force=True)
        self.workspace.setCurrentWidget(self.results_view)

    def follow_run_results(self):
        self.follow_current_results = True
        self.refresh_results(force=True)
        self.workspace.setCurrentWidget(self.results_view)

    def update_results_actions(self):
        for title, name in (('Open HTML Report', 'report.html'), ('Open Output Folder', '')):
            self.results_actions[title].setEnabled(
                self.results_view.results is not None and
                (self.results_view.results.directory / name).exists())

    def refresh_results(self, *, force=False):
        if not force and time.monotonic() - self.last_results_refresh < 1.:
            return
        self.last_results_refresh = time.monotonic()
        status = None
        if self.follow_current_results and self.controller.snapshot is not None:
            self.results_directory = Path(self.controller.snapshot.compiled.configuration.run.output_directory)
            status = 'Cancel pending' if self.controller.cancel_pending and self.controller.active else self.controller.state
        if self.results_directory is None:
            return
        try:
            selection = self.results_view.selection if self.results_view.results is not None and self.results_view.results.directory == self.results_directory else None
            self.results_view.load(self.results_directory, evaluation=selection, status=status)
        except ValueError as exc:
            # Checkpoint/event files may not exist yet. Never replace a valid
            # display with in-flight or incomplete state.
            if self.results_view.results is None or self.results_view.results.directory != self.results_directory:
                self.results_view.clear()
                self.results_view.summary.setPlainText(f'{status or "Unavailable"} — no valid completed evaluation yet.\n{exc}')
        self.update_results_actions()

    def open_results(self, path):
        try:
            self.results_view.load(path)
        except ValueError as exc:
            self.report_error(exc)
            return False
        self.results_directory = self.results_view.results.directory
        self.follow_current_results = False
        self.update_results_actions()
        self.workspace.setCurrentWidget(self.results_view)
        return True

    def open_document(self, path):
        document = EditorDocument.open(path)
        self.document = document
        self.refresh()

    def mark_pending(self, *args):
        self.pending = True
        self.invalidate_readiness()
        self.statusBar().showMessage('Draft edits — Apply validates and commits an immutable configuration snapshot.')

    def refresh(self):
        if self.document is None: return
        self.invalidate_readiness()
        section = self.navigation.currentItem().text(0) if self.navigation.currentItem() else 'Project'
        self.navigation.clear()
        while self.pages.count():
            page = self.pages.widget(0)
            self.pages.removeWidget(page)
            page.deleteLater()
        self.widgets = {}
        self.section_items = {}
        forms = {}
        for name in ('Project', 'Microscope', 'Regions', 'Scan', 'Experimental data',
                     'Refinement', 'Execution', 'Outputs', 'All fields (read-only)'):
            item = QTreeWidgetItem([name])
            item.setData(0, Qt.ItemDataRole.UserRole, self.pages.count())
            self.navigation.addTopLevelItem(item)
            self.section_items[name] = item
            body = QWidget()
            form = QFormLayout(body)
            scroll = QScrollArea()
            scroll.setWidgetResizable(True)
            scroll.setWidget(body)
            self.pages.addWidget(scroll)
            forms[name] = form
        for spec in editable_fields(self.document.config):
            value = getattr(get_object(self.document.config, spec.path), spec.name)
            if spec.kind == 'backend':
                widget = QComboBox()
                widget.addItems(['serial', 'batched'])
                widget.setCurrentText(value)
                widget.currentTextChanged.connect(self.mark_pending)
            else:
                widget = QLineEdit(display_value(value))
                widget.textChanged.connect(self.mark_pending)
            widget.setObjectName('.'.join(str(s) for s in spec.key))
            self.widgets[spec.key] = widget
            if spec.kind in ('file', 'directory', 'save_file'):
                container = QWidget()
                row = QHBoxLayout(container)
                row.setContentsMargins(0, 0, 0, 0)
                row.addWidget(widget)
                button = QPushButton('Browse…')
                button.clicked.connect(lambda checked=False, w=widget, k=spec.kind: self.browse(w, k))
                row.addWidget(button)
                forms[spec.section].addRow(spec.label, container)
            else:
                forms[spec.section].addRow(spec.label, widget)
        forms['Execution'].addRow(QLabel('Serial requires a blank batch size. Batched requires an explicit positive size.'))
        self.build_refinement(forms['Refinement'])
        tree = QTreeWidget()
        tree.setHeaderLabels(['Canonical field / scientific order', 'Value'])
        def populate(parent, name, value):
            item = QTreeWidgetItem(parent, [name, '' if is_dataclass(value) or isinstance(value, tuple) else display_value(value)])
            if is_dataclass(value):
                for field in fields(value): populate(item, field.name, getattr(value, field.name))
            elif isinstance(value, tuple):
                for i, entry in enumerate(value): populate(item, getattr(entry, 'id', str(i)), entry)
        populate(tree, 'QcbedConfig', self.document.config)
        tree.expandToDepth(1)
        forms['All fields (read-only)'].addRow(tree)
        self.navigation.setCurrentItem(self.section_items.get(section, self.section_items['Project']))
        self.pending = False
        self.document.issues = self.document.committed_issues
        self.show_diagnostics()
        self.setWindowTitle(f'QCBED — {self.document.path.name}' + (' *' if self.document.dirty else ''))

    def build_refinement(self, form):
        config = self.document.config
        categories = config.refinement.active_categories
        label = QLabel('Binding selection and active refinement are distinct. Current category filter (read-only): '
                       + (', '.join(str(getattr(c, 'value', c)) for c in categories) if categories is not None else 'all')
                       + '\nNew bindings require explicit D1, D2 and PLAMDA. Table sorting changes display only.')
        label.setWordWrap(True)
        form.addRow(label)
        catalogue = available_refinement_targets(config)
        table = self.refinement_table = QTableWidget(len(catalogue), 8)
        table.setHorizontalHeaderLabels(['Physical target', 'Binding', 'Active', 'Starting value', 'D1', 'D2', 'PLAMDA', 'Availability / units'])
        bindings = {b.target: b for b in config.refinement.parameters}
        for row, info in enumerate(catalogue):
            binding = bindings.get(info.target)
            values = [target_label(config, info), '', 'Yes' if info.active else 'No', repr(info.starting_value),
                *(tuple(repr(v) for v in (*binding.range_factors, binding.plamda)) if binding else ('', '', '')),
                info.reason or f'Available ({info.units})']
            for col, value in enumerate(values):
                item = QTableWidgetItem(value)
                item.setData(Qt.ItemDataRole.UserRole, info.target)
                if col not in (4, 5, 6) or not info.available:
                    item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                if col == 1:
                    item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
                    item.setCheckState(Qt.CheckState.Checked if info.selected else Qt.CheckState.Unchecked)
                    if not info.available and not info.selected:
                        item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEnabled)
                item.setToolTip(info.reason or info.label)
                table.setItem(row, col, item)
        table.setSortingEnabled(True)
        table.resizeColumnsToContents()
        table.setColumnWidth(0, 390)
        table.setMinimumHeight(350)
        table.itemChanged.connect(self.mark_pending)
        form.addRow(table)

    def drafts(self):
        text = {key: widget.currentText() if isinstance(widget, QComboBox) else widget.text()
                for key, widget in self.widgets.items()}
        table = self.refinement_table
        refinement = [(table.item(row, 0).data(Qt.ItemDataRole.UserRole),
            table.item(row, 1).checkState() == Qt.CheckState.Checked,
            *(table.item(row, col).text() for col in (4, 5, 6))) for row in range(table.rowCount())]
        return text, refinement

    def apply_edits(self):
        if self.document is None: return False
        if not self.pending: return True
        if self.document.apply(*self.drafts()):
            self.refresh()
            self.statusBar().showMessage('Valid edits committed. File readiness requires Validate.')
            return True
        self.show_diagnostics()
        self.statusBar().showMessage('Draft errors — committed configuration is unchanged.')
        return False

    def validate_document(self):
        if self.apply_edits():
            self.document.validate(check_files=True)
            self.show_diagnostics()
            self.check_run_readiness()

    def show_diagnostics(self):
        self.issues.clear()
        for issue in self.document.issues:
            severity, phase, message = diagnostic_presentation(issue)
            item = QTreeWidgetItem([severity, phase,
                f'{issue.section}: {issue.field or ""}', message])
            item.setData(0, Qt.ItemDataRole.UserRole, issue)
            self.issues.addTopLevelItem(item)
        for col in range(3): self.issues.resizeColumnToContents(col)
        summary = self.document.summary()
        labels = {'configuration': 'Configuration validity', 'compile': 'Compile validity',
                  'execution': 'Execution readiness', 'files': 'External file/output readiness'}
        lines = [f'{labels.get(key, key.replace("_", " ").capitalize())}: {value}'
                 for key, value in summary.items() if key != 'bounds']
        lines.append('\nResolved hard bounds (committed snapshot):')
        from .structure import resolved_bound_rows
        rows = resolved_bound_rows(self.document.config)
        for label, lower, upper, provenance in rows:
            lines.append(f'{label}: {lower!r} … {upper!r} ({provenance})')
        if not rows: lines.append('No active parameters')
        self.summary.setPlainText('\n'.join(lines))

    def navigate(self, item, previous=None):
        if item: self.pages.setCurrentIndex(item.data(0, Qt.ItemDataRole.UserRole))

    def navigate_issue(self, item, column=0):
        issue = item.data(0, Qt.ItemDataRole.UserRole)
        section = {'project':'Project', 'microscope':'Microscope', 'sample':'Regions',
            'scan':'Scan', 'experiment':'Experimental data', 'refinement':'Refinement',
            'execution':'Execution', 'run':'Outputs'}.get(issue.section)
        if section: self.navigation.setCurrentItem(self.section_items[section])
        widget = self.widgets.get(issue.object_path + (issue.field,))
        if widget:
            widget.setFocus()
            if isinstance(widget, QLineEdit): widget.selectAll()
        for row in range(self.refinement_table.rowCount()):
            if self.refinement_table.item(row, 0).data(Qt.ItemDataRole.UserRole) in issue.object_path:
                self.refinement_table.selectRow(row)
                self.refinement_table.scrollToItem(self.refinement_table.item(row, 0))

    def browse(self, widget, kind):
        start = str(self.document.base_directory)
        if kind == 'directory': path = QFileDialog.getExistingDirectory(self, 'Output directory', start)
        elif kind == 'file': path, _ = QFileDialog.getOpenFileName(self, 'Experimental data', start)
        else: path, _ = QFileDialog.getSaveFileName(self, 'Output filename (no file written)', start)
        if path: widget.setText(path)

    def save_document(self, path=None, *, overwrite=False):
        if not self.apply_edits(): return False
        self.document.save(path, overwrite=overwrite)
        self.refresh()
        self.statusBar().showMessage(f'Saved {self.document.path}')
        return True

    def report_error(self, exc):
        QMessageBox.warning(self, 'Configuration error', str(exc))

    def save_current(self):
        if self.document is None: return False
        if self.document.needs_save_as: return self.choose_save_as()
        try: return self.save_document(overwrite=True)
        except (OSError, ValueError, TypeError) as exc: self.report_error(exc); return False

    def choose_save_as(self):
        if self.document is None: return False
        path, _ = QFileDialog.getSaveFileName(self, 'Save canonical YAML as', str(self.document.path),
            'YAML (*.yml *.yaml)', options=QFileDialog.Option.DontConfirmOverwrite)
        if not path: return False
        overwrite = Path(path).exists()
        if overwrite and QMessageBox.question(self, 'Replace YAML?', f'Replace {path}?') != QMessageBox.StandardButton.Yes: return False
        try: return self.save_document(path, overwrite=overwrite)
        except (OSError, ValueError, TypeError) as exc: self.report_error(exc); return False

    def confirm_discard(self):
        if not self.document or not (self.pending or self.document.dirty): return True
        answer = QMessageBox.question(self, 'Unsaved edits', 'Save edits before continuing?',
            QMessageBox.StandardButton.Save | QMessageBox.StandardButton.Discard | QMessageBox.StandardButton.Cancel)
        if answer == QMessageBox.StandardButton.Save: return self.save_current()
        return answer == QMessageBox.StandardButton.Discard

    def choose_open(self):
        if not self.confirm_discard(): return
        path, _ = QFileDialog.getOpenFileName(self, 'Open canonical YAML',
            str(self.document.base_directory) if self.document else str(Path.home()), 'YAML (*.yml *.yaml)')
        if path:
            try: self.open_document(path)
            except (OSError, ValueError, TypeError) as exc: self.report_error(exc)

    def import_legacy(self):
        if self.controller.active or not self.confirm_discard(): return
        from .compatibility_widgets import ImportDialog
        dialog = ImportDialog(self)
        if not dialog.exec(): return
        try:
            self.import_legacy_document(*dialog.result_values)
        except (OSError,ValueError,TypeError) as exc: self.report_error(exc)

    def import_legacy_document(self, source, experimental, policy):
        if self.controller.active: raise RuntimeError('Cannot import during an active run')
        from .compatibility import import_document
        document = import_document(source,experimental,policy)
        self.document = document
        self.refresh()
        self.statusBar().showMessage('Imported canonical document — unsaved. Save As YAML to retain it.')

    def view_compatibility(self):
        if self.controller.active or self.document is None: return
        from .compatibility import compatibility_report
        from .compatibility_widgets import ReportDialog
        try:
            report = compatibility_report(self.document.config,base_directory=self.document.base_directory)
            ReportDialog(report,parent=self).exec()
        except (OSError,ValueError,TypeError) as exc: self.report_error(exc)

    def export_legacy(self, bundle=False):
        if self.controller.active or self.document is None: return
        if self.pending:
            self.statusBar().showMessage('Apply or revert drafts before export; export does not commit edits.')
            return
        from .compatibility import compatibility_report, export_document
        from .compatibility_widgets import ExportDialog, ReportDialog
        options = ExportDialog(self.document.base_directory,bundle=bundle,parent=self)
        if not options.exec(): return
        destination, limits = options.result_values
        try:
            report = compatibility_report(self.document.config,base_directory=self.document.base_directory,
                                          bundle=bundle,filename_limits=limits)
            if not ReportDialog(report,exporting=True,parent=self).exec(): return
            export_document(self.document,destination,bundle=bundle,filename_limits=limits)
            self.statusBar().showMessage(f'Optional legacy artifacts exported to {destination}; canonical document unchanged.')
        except (OSError,ValueError,TypeError) as exc: self.report_error(exc)

    def edit_structure(self):
        if self.controller.active or not self.apply_edits(): return
        from .structure_widgets import StructureDialog
        dialog = StructureDialog(self.document,self)
        if dialog.exec(): self.refresh()

    def undo(self):
        if self.controller.active or self.document is None: return
        if self.pending:
            self.statusBar().showMessage('Apply or revert incomplete drafts before Undo.')
            return
        if self.document.undo(): self.refresh()

    def redo(self):
        if self.controller.active or self.document is None: return
        if self.pending:
            self.statusBar().showMessage('Apply or revert incomplete drafts before Redo.')
            return
        if self.document.redo(): self.refresh()

    def invalidate_readiness(self):
        self.run_action.setEnabled(False)
        self.resume_action.setEnabled(False)
        for action in (self.run_action, self.resume_action):
            action.setToolTip('Validate the current configuration first')

    def check_run_readiness(self):
        if self.controller.active or self.pending: return
        reasons = []
        for resume, action in ((False, self.run_action), (True, self.resume_action)):
            issues = readiness(self.document.config, self.document.base_directory, resume=resume)
            errors = [diagnostic_presentation(i)[2] for i in issues if i.severity == 'error']
            action.setEnabled(not errors)
            action.setToolTip('\n'.join(errors) if errors else 'Ready')
            reasons.append(('Resume' if resume else 'Run') + ': ' + ('; '.join(errors) if errors else 'Ready'))
        self.run_log.appendPlainText('\n'.join(reasons))

    def start_run(self, resume=False):
        if self.controller.active or self.document is None or self.pending: return
        try:
            self.controller.start(self.document.config, base_directory=self.document.base_directory, resume=resume)
        except (OSError, ValueError, RuntimeError) as exc:
            self.run_log.appendPlainText(str(exc))
            self.run_status.setText('Failed to launch — ' + str(exc))
            self.invalidate_readiness()
            return
        self.follow_current_results = True
        self.results_directory = Path(self.controller.snapshot.compiled.configuration.run.output_directory)
        self.refresh_results(force=True)
        self.run_log.clear()
        self.run_log.appendPlainText('Resume' if resume else 'Fresh run')
        self.run_status.setText('Running — editor locked; immutable configuration and data snapshot')
        self.pages.setEnabled(False)
        for action in self.edit_actions: action.setEnabled(False)
        self.invalidate_readiness()
        self.cancel_action.setEnabled(True)
        self.timer.start()

    def cancel_run(self):
        if self.controller.request_cancel():
            self.cancel_action.setEnabled(False)
            self.run_log.appendPlainText('Cancellation requested; waiting for a safe boundary and checkpoint finalization.')
            self.run_status.setText('Running — cancellation pending; editor remains locked')
            self.refresh_results(force=True)

    def poll_run(self):
        for kind, payload in self.controller.poll():
            self.run_log.appendPlainText(payload if isinstance(payload, str) else json.dumps(payload, default=str))
        self.refresh_results()
        if self.controller.active:
            self.cancel_action.setEnabled(self.controller.state == 'Running' and not self.controller.cancel_pending)
            elapsed = time.monotonic() - self.controller.started
            pending = 'cancellation pending; waiting for safe completion' if self.controller.cancel_pending else 'immutable snapshot'
            self.run_status.setText(f'Running — {pending}; elapsed {elapsed:.1f} s')
            return
        self.timer.stop()
        self.cancel_action.setEnabled(False)
        self.pages.setEnabled(True)
        for action in self.edit_actions: action.setEnabled(True)
        result = self.controller.result
        if result:
            summary = result['summary']
            self.run_status.setText(f"{self.controller.state} — objective {summary.get('final_objective', 'unavailable')}; "
                f"evaluations {summary.get('evaluations', 'unavailable')}; elapsed {result['elapsed_seconds']:.2f} s; "
                f"output {summary.get('output_directory')}; checkpoint {result['checkpoint']}; final state {result['final_state']}")
        else:
            self.run_status.setText(self.controller.state + ' — see worker traceback/log below')
        self.check_run_readiness()
        self.refresh_results(force=True)
        self.statusBar().showMessage('Results updated — choose Results → Show Results to inspect.')

    def closeEvent(self, event):
        if self.controller.active:
            self.statusBar().showMessage('A run is active. Request Cancel and wait for safe completion before closing.')
            event.ignore()
            return
        if self.confirm_discard(): event.accept()
        else: event.ignore()
