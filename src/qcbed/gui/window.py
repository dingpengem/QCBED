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
    QLabel, QScrollArea, QFileDialog, QMessageBox, QTabWidget, QCheckBox)

from qcbed import available_refinement_targets, get_object
from .document import EditorDocument, editable_fields, display_value, target_label, diagnostic_presentation


from .run_controller import RunController, readiness


class ConfigWindow(QMainWindow):
    def __init__(self, path=None):
        super().__init__()
        self.controller = RunController()
        self.document = None
        from .workflow import Workflow
        self.workflow = Workflow()
        self.pending = False
        self.output_draft = False
        self.output_receipt = None
        self.widgets = {}
        self.resize(1250, 850)
        toolbar = self.addToolBar('Configuration')
        self.edit_actions = []
        for title, callback in [('Load Configuration…', self.choose_open), ('Save', self.save_current),
                ('Save As', self.choose_save_as), ('Apply edits', self.apply_edits),
                ('Check Configuration (Expert)', self.validate_document), ('Revert drafts', self.revert_drafts),
                ('Structure / Bounds', self.edit_structure), ('Undo', self.undo), ('Redo', self.redo)]:
            action = QAction(title, self)
            action.triggered.connect(lambda checked=False, cb=callback: cb())
            toolbar.addAction(action)
            self.edit_actions.append(action)
        file_menu=self.menuBar().addMenu('File')
        for action in self.edit_actions[:3]:file_menu.addAction(action)
        self.project_path_label=QLabel('Project configuration: Unsaved')
        self.statusBar().addPermanentWidget(self.project_path_label)
        legacy_menu = self.menuBar().addMenu('Compatibility / Legacy')
        for label, callback in [('Import Legacy Project…', self.import_legacy),
                ('Export FOR007…', lambda: self.export_legacy(False)),
                ('Export Legacy Bundle…', lambda: self.export_legacy(True)),
                ('View Compatibility Report…', self.view_compatibility)]:
            action = legacy_menu.addAction(label)
            action.triggered.connect(lambda checked=False, cb=callback: cb())
            self.edit_actions.append(action)
        self.run_action = toolbar.addAction('Fresh Run', lambda: self.start_run(False))
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
        self.editor_split = top
        self.workflow_tabs = QTabWidget()
        self.stage_buttons=[];self.stage_status=[]
        for i,name in enumerate(('A. Data Prep','B. Atomic Model','C. Main QCBED Parameters')):
            page=QWidget();layout=QVBoxLayout(page)
            status=QLabel('Pending');status.setWordWrap(True);layout.addWidget(status)
            button=QPushButton(('Validate Data Prep','Validate Atomic Model','Validate Optimizer')[i])
            button.clicked.connect(lambda checked=False,index=i:self.validate_stage(index))
            layout.addWidget(button);self.stage_buttons.append(button);self.stage_status.append(status)
            self.workflow_tabs.addTab(page,name)
        self.workflow_tabs.widget(0).layout().addWidget(top)
        self.workflow_tabs.currentChanged.connect(self.switch_stage)
        self.issues = QTreeWidget()
        self.issues.setHeaderLabels(['Severity', 'Phase', 'Location', 'Message'])
        self.issues.itemClicked.connect(self.navigate_issue)
        self.summary = QPlainTextEdit()
        self.summary.setReadOnly(True)
        bottom = QSplitter()
        bottom.addWidget(self.issues)
        bottom.addWidget(self.summary)
        split = QSplitter(Qt.Orientation.Vertical)
        split.addWidget(self.workflow_tabs)
        from .sections import ExpandableSection
        self.diagnostics_section=ExpandableSection('Diagnostics and compiled summary')
        self.diagnostics_section.content.addWidget(bottom)
        split.addWidget(self.diagnostics_section)
        self.run_status = QLabel('Idle — Validate to check Run and Resume readiness')
        self.run_status.setWordWrap(True)
        self.run_log = QPlainTextEdit()
        self.run_log.setReadOnly(True)
        self.run_log.document().setMaximumBlockCount(3000)
        run_panel = QWidget()
        layout = QVBoxLayout(run_panel)
        layout.addWidget(self.run_status)
        from .sections import ExpandableSection
        self.log_section=ExpandableSection('Expert / Legacy settings — execution log')
        self.log_section.content.addWidget(self.run_log);layout.addWidget(self.log_section)
        split.addWidget(run_panel)
        split.setSizes([650, 35, 80])
        self.timer = QTimer(self)
        self.timer.setInterval(100)
        self.timer.timeout.connect(self.poll_run)
        from .results_view import ResultsView
        self.results_view = ResultsView()
        self.workspace = QTabWidget()
        from .workflow_guidance import WorkflowGuidance
        self.guidance=WorkflowGuidance(self)
        configuration_page=QWidget();configuration_layout=QVBoxLayout(configuration_page)
        configuration_layout.addWidget(self.guidance);configuration_layout.addWidget(split)
        self.workspace.addTab(configuration_page, 'Configuration')
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
        from .preparation_view import NewPreparationView
        self.preparation=NewPreparationView(self)
        prep_scroll=QScrollArea();prep_scroll.setWidgetResizable(True);prep_scroll.setWidget(self.preparation)
        self.preparation_scroll=prep_scroll
        self.preparation_split=QSplitter(Qt.Orientation.Vertical)
        self.preparation_split.setChildrenCollapsible(False)
        self.preparation_split.addWidget(prep_scroll);self.preparation_split.addWidget(self.editor_split)
        self.preparation_split.setStretchFactor(0,3);self.preparation_split.setStretchFactor(1,1)
        self.editor_split.setMinimumHeight(180)
        self.preparation_split.setSizes([800,250])
        self.preparation.mode.currentIndexChanged.connect(
            lambda index:self.preparation_split.setSizes([140,900] if index else [800,250]))
        self.workflow_tabs.widget(0).layout().insertWidget(1,self.preparation_split,1)
        self.workflow_tabs.widget(0).layout().removeWidget(self.stage_buttons[0])
        self.workflow_tabs.widget(0).layout().addWidget(self.stage_buttons[0])
        from .run_dashboard import RunDashboard
        self.run_dashboard=RunDashboard(self)
        self.workflow_tabs.addTab(self.run_dashboard,'D. Run & Results')
        for action in (self.run_action,self.resume_action,self.cancel_action):
            action.changed.connect(self.run_dashboard.update_state)
        self.run_dashboard.update_state()
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
        if hasattr(self,'run_dashboard'):self.run_dashboard.update_state()

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

    def update_project_path(self):
        doc=self.document
        path=doc.path if doc is not None and not doc.needs_save_as else None
        if path is None and hasattr(self,'preparation'):
            path=getattr(self.preparation,'associated_project_path',None)
        self.project_path_label.setText('Project configuration: '+(str(path) if path else 'Unsaved'))

    def can_save_project(self):
        if self.document is None:return False
        if self.document.config.crystallography is not None:
            self.workflow.sync(self.document.config)
            if self.pending or self.output_draft or getattr(self,'site_group_pending',False) or not self.workflow.ready:
                self.statusBar().showMessage('Apply edits and validate A → B → C before saving the project.')
                return False
        return True

    def open_document(self, path):
        if hasattr(self,'preparation') and self.preparation.active:return
        from qcbed.yaml_config import read_yaml_document
        raw=read_yaml_document(path)
        prep=raw.get('data_prep',{}) if isinstance(raw,dict) else {}
        if isinstance(prep,dict) and 'pattern' in prep and not any(k in prep for k in ('prepared_dataset','preparation')):
            self.preparation.open_raw(raw,Path(path).resolve());return
        document = EditorDocument.open(path)
        self.document = document
        if hasattr(self,'preparation'):self.preparation.mode.setCurrentIndex(1)
        from .workflow import Workflow
        self.workflow=Workflow()
        self.refresh()

    def mark_pending(self, *args, stage_override=None):
        self.pending = True
        widget=self.widgets.get(('crystallography','depth_angstrom'))
        if widget is not None and self.document and self.document.config.crystallography:
            table=getattr(self,'refinement_table',None)
            if table is not None:
                try:
                    table.blockSignals(True)
                    for row in range(table.rowCount()):
                        target=table.item(row,0).data(Qt.ItemDataRole.UserRole)
                        if target.category.value=='depth':
                            value=self.document.config.crystallography.depth_angstrom
                            text=repr(value)
                            try:changed=float(widget.text())!=value
                            except ValueError:changed=True
                            if changed:text+=' (uncommitted draft: '+widget.text()+')'
                            table.item(row,3).setText(text)
                    table.blockSignals(False)
                except RuntimeError:pass
        if self.document and self.document.config.crystallography:
            sender=self.sender()
            setup_keys={s.key for s in editable_fields(self.document.config)
                        if s.section in ('Microscope','Numerical / Simulation cell')}
            if sender is not None and any(self.widgets.get(k) is sender for k in setup_keys):
                self.workflow.simulation_receipt=None
                self.workflow.refinement_receipt=None
                if hasattr(self,'c_tabs'):
                    self.validate_refinement_button.setEnabled(False)
                    self.refinement_table.hide()
                    self.refinement_notice.setText('Validate Simulation Setup to generate refinement targets.')
        stage=min(self.workflow_tabs.currentIndex(),2)
        d_edit=False
        if self.document:
            spec=next((s for s in editable_fields(self.document.config) if self.widgets.get(s.key) is self.sender()),None)
            if spec:
                stage=1 if spec.section=='Structure' else 0 if spec.section in ('Scan','Experimental data') else 2
                d_edit=bool(self.document.config.crystallography and spec.section in ('Project','Execution','Outputs'))
        if self.sender() is getattr(self,'refinement_table',None):stage=2
        if stage_override is not None:stage=stage_override
        if not d_edit:self.workflow.invalidate(stage)
        self.update_stages()
        self.invalidate_readiness()
        self.statusBar().showMessage('Draft edits — Apply validates and commits an immutable configuration snapshot.')

    def revert_drafts(self):
        # Intentional discard is different from rebuilding after upstream edits.
        self.pending=False;self._target_draft_cache=None
        if getattr(self,'elastic_panel',None) is not None:self.elastic_panel.dirty=False
        self.refresh()

    def refresh(self):
        if self.document is None: return
        from .target_drafts import capture,restore
        target_drafts=capture(self)
        view_state=None
        crystal=self.document.config.crystallography
        if crystal and getattr(self,'atomic_view_identity',None)==(crystal.structure.path,crystal.zone_axis):
            view_state=([(a.get_xlim(),a.get_ylim()) for a in self.atomic_axes],self.atomic_3d_axes.elev,self.atomic_3d_axes.azim)
        self.workflow.sync(self.document.config)
        self.invalidate_readiness()
        c_index=self.c_tabs.currentIndex() if hasattr(self,'c_tabs') else 0
        section = self.navigation.currentItem().text(0) if self.navigation.currentItem() else 'Data Prep'
        self.navigation.clear()
        while self.pages.count():
            page = self.pages.widget(0)
            self.pages.removeWidget(page)
            page.deleteLater()
        self.run_dashboard.clear_output_editor()
        self.output_draft = False
        self.widgets = {}
        self.section_items = {}
        forms = {}
        modern=self.document.config.crystallography is not None
        names=('Data Prep','Project','Microscope',* (('Structure','Site groups','Capabilities') if modern else ('Regions',)),
               'Scan',*(('Numerical / Simulation cell',) if modern else ()),*(('Experimental data',) if not modern else ()), 'Refinement','Optimizer / Run',
               'Execution','Outputs','All fields (read-only)')
        groups={}
        for label in ('A. Data Prep','B. Atomic Model','C. Main QCBED Parameters'):
            parent=QTreeWidgetItem([label]);self.navigation.addTopLevelItem(parent);groups[label]=parent
        for name in names:
            item = QTreeWidgetItem([name])
            item.setData(0, Qt.ItemDataRole.UserRole, self.pages.count())
            group=('A. Data Prep' if name in ('Data Prep','Experimental data') or (modern and name=='Scan') else
                   'B. Atomic Model' if name in ('Structure','Site groups','Regions','Capabilities') else 'C. Main QCBED Parameters')
            parent=groups[group];parent.addChild(item)
            if parent.childCount()==1:parent.setData(0,Qt.ItemDataRole.UserRole,self.pages.count())
            parent.setExpanded(True)
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
            if spec.kind in ('backend','signal','workflow'):
                widget = QComboBox()
                choices=['serial','batched'] if spec.kind=='backend' else ['Conventional','Gradient'] if spec.kind=='signal' else ['ANNEAL → AMOEBA','ANNEAL','AMOEBA']
                if spec.kind=='workflow':choices+=['PARAMIN','SIMPARA']
                widget.addItems(choices)
                shown=choices[value] if spec.kind=='signal' else {'ANSI':'ANNEAL → AMOEBA','ANNE':'ANNEAL','SIMP':'AMOEBA','PARA':'PARAMIN','SPAR':'SIMPARA'}.get(value,value) if spec.kind=='workflow' else value
                if shown not in choices:widget.addItem(shown)
                widget.setCurrentText(shown)
                widget.currentTextChanged.connect(self.mark_pending)
            else:
                widget = QLineEdit(display_value(value))
                widget.textChanged.connect(self.output_path_edited if modern and spec.key == ('run','output_directory') else self.mark_pending)
                if spec.key == ('crystallography', 'structure', 'path'):
                    widget.editingFinished.connect(lambda w=widget: self.load_structure_path(w.text()))
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
                if modern and spec.key == ('run','output_directory'):
                    self.run_dashboard.set_output_editor(container)
                else:
                    forms[spec.section].addRow(spec.label, container)
            else:
                forms[spec.section].addRow(spec.label, widget)
        forms['Execution'].addRow(QLabel('Serial requires a blank batch size. Batched requires an explicit positive size.'))
        if modern:
            note=QLabel('Modern grouped refinement supports Batched execution. Loaded settings are retained when Advanced is closed. '
                'Batch size controls tilt scheduling and memory, not physical parameters. '
                'Worker/process count is not supported by this canonical run interface.');note.setWordWrap(True)
            forms['Execution'].addRow(note)
            cadence=QLabel('Every completed evaluation (fixed for modern grouped runs)');cadence.setWordWrap(True)
            forms['Execution'].addRow('Checkpoint / result cadence',cadence)
            note=QLabel('Legacy checkpoint_every and snapshot_every values are preserved in YAML; modern grouped runs use the fixed cadence above.');note.setWordWrap(True)
            forms['Execution'].addRow(note)
        if modern:
            from .modern_widgets import populate_modern
            populate_modern(self,forms)
            from .sections import collapse_form,ExpandableSection,move_form_row
            collapse_form(forms['Scan'],'Expert / Legacy settings — recorded scan geometry')
            format_section=ExpandableSection('Advanced settings / Structure format')
            format_body=QWidget();format_form=QFormLayout(format_body)
            move_form_row(forms['Structure'],self.widgets[('crystallography','structure','format')],format_form)
            format_section.content.addWidget(format_body);forms['Structure'].addRow(format_section)
        from .data_view import DataPrepView
        self.data_panel=DataPrepView(self);forms['Data Prep'].addRow(self.data_panel)
        if not modern:
            new_run=QPushButton('Choose New Run Paths');new_run.clicked.connect(self.new_run_paths)
            forms['Outputs'].addRow(new_run)
        self.build_refinement(forms['Refinement'])
        if modern:
            self.build_c_tabs(forms)
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
        self.navigation.setCurrentItem(self.section_items.get(section, self.section_items['Data Prep']))
        self.pending = False
        self.document.issues = self.document.committed_issues
        self.show_diagnostics()
        self.update_stages()
        self.switch_stage(self.workflow_tabs.currentIndex())
        self.update_project_path()
        if modern:
            self.c_tabs.setCurrentIndex(c_index)
            restore(self,target_drafts)
            if view_state and hasattr(self,'atomic_axes'):
                for ax,(x,y) in zip(self.atomic_axes,view_state[0]):ax.set_xlim(x);ax.set_ylim(y)
                self.atomic_3d_axes.view_init(elev=view_state[1],azim=view_state[2])
        self.setWindowTitle(f'QCBED — {self.document.path.name}' + (' *' if self.document.dirty else ''))

    def field_drafts(self):
        """Capture text independently of validation receipts before a partial commit."""
        return {key:w.currentText() if isinstance(w,QComboBox) else w.text()
                for key,w in self.widgets.items()}

    def restore_field_drafts(self, drafts):
        from .document import parse_text
        specs={s.key:s for s in editable_fields(self.document.config)}
        for key,text in drafts.items():
            widget=self.widgets.get(key)
            if widget is None:continue
            spec=specs[key]
            try:
                if parse_text(text,spec.kind)==getattr(get_object(self.document.config,spec.path),spec.name):continue
            except (ValueError,TypeError,OverflowError):pass
            if isinstance(widget,QComboBox):widget.setCurrentText(text)
            else:widget.setText(text)

    def output_is_validated(self):
        return bool(self.document and not self.output_draft and
            self.output_receipt==(id(self.document),self.document.config.run.output_directory))

    def output_path_edited(self):
        self.output_draft=True
        self.output_receipt=None
        self.run_dashboard.output_validation='Not checked — Validate Output Directory'
        self.invalidate_readiness()

    def validate_output_directory(self):
        if self.controller.active or self.document is None:return False
        from dataclasses import replace
        from .output_paths import validate_output_directory
        try:
            path=validate_output_directory(self.widgets[('run','output_directory')].text(),
                self.document.base_directory,create=True)
        except (ValueError,OSError) as exc:
            self.run_dashboard.output_validation=str(exc)
            self.invalidate_readiness()
            return False
        config=replace(self.document.config,run=replace(self.document.config.run,
            output_directory=str(path),output_filename=str(path/'final_out.dat')))
        self.document._record(config)
        widget=self.widgets[('run','output_directory')]
        widget.blockSignals(True);widget.setText(str(path));widget.blockSignals(False)
        self.output_draft=False
        self.output_receipt=(id(self.document),self.document.config.run.output_directory)
        self.run_dashboard.output_validation='Valid'
        for name in ('issues','committed_issues'):
            setattr(self.document,name,tuple(i for i in getattr(self.document,name)
                if i.code not in ('output_not_writable','output_exists','output_replacement_unsafe')))
        self.update_project_path()
        self.run_dashboard.update_state()
        # Full runtime checks remain separate from the lightweight directory probe.
        QTimer.singleShot(0,self,self.check_run_readiness)
        return True

    def new_run_paths(self):
        if self.controller.active or self.pending:return
        from datetime import datetime
        from uuid import uuid4
        from qcbed import update_object
        old=Path(self.document.config.run.output_directory)
        new=old.parent/(old.name+'_'+datetime.now().strftime('%Y%m%d_%H%M%S')+'_'+uuid4().hex[:8])
        config=update_object(self.document.config,('run',),output_directory=str(new),
            output_filename=str(new/('final_out.dat' if self.document.config.crystallography else Path(self.document.config.run.output_filename).name)))
        config=update_object(config,('legacy_export',),output_filename_raw=None)
        if self.document.commit_structure(config):self.refresh();self.validate_document()

    def build_refinement(self, form):
        config = self.document.config
        categories = config.refinement.active_categories
        label = QLabel(('Choose physical quantities and hard bounds. Current category filter: ' if config.crystallography else 'Binding selection and active refinement are distinct. Current category filter (read-only): ')
                       + (', '.join(str(getattr(c, 'value', c)) for c in categories) if categories is not None else 'all')
                       + ('' if config.crystallography else '\nNew bindings require explicit D1, D2 and PLAMDA. Table sorting changes display only.'))
        label.setWordWrap(True)
        form.addRow(label)
        if config.crystallography:
            hint=QLabel('Select any nonempty subset of available physical targets; requested depth is optional. '
                'Configured hard-bound rules are shown before selection. Where no canonical default exists, enter physical hard bounds. Recommended search settings are supplied automatically. '
                'Starting values always come from the committed physical model.')
            hint.setWordWrap(True);label.setToolTip(hint.text())
        try:catalogue = () if self.document.thermal_incomplete(config) or (config.crystallography and self.workflow.simulation_receipt is None) else available_refinement_targets(config)
        except (OSError,ValueError,TypeError):catalogue=()
        modern=config.crystallography is not None
        table = self.refinement_table = QTableWidget(len(catalogue), 10 if modern else 8)
        table.setHorizontalHeaderLabels(['Physical target', 'Binding', 'Active', 'Committed starting value', 'D1 (range factor)', 'D2 (range factor)', 'PLAMDA (optimizer metadata)', 'Availability / units']+(['Hard lower bound','Hard upper bound'] if modern else []))
        resolved={}
        if modern:
            from qcbed.crystallography_config import compile_config_from_sources
            from qcbed.crystallography_refinement import physical_target_layout
            try:resolved={r['target_id']:r for r in physical_target_layout(compile_config_from_sources(config))}
            except (ValueError,TypeError,OSError):pass
        bindings = {b.target: b for b in config.refinement.parameters}
        for row, info in enumerate(catalogue):
            binding = bindings.get(info.target)
            values = [target_label(config, info)+(f' ({info.units})' if modern and info.units else ''), '', ('Active' if modern else 'Yes') if info.active else (('Available' if info.available else 'Unsupported — '+str(info.reason)) if not info.selected and modern else ('No — not selected' if not info.selected else 'No — excluded by category filter')), repr(info.starting_value),
                *(tuple(repr(v) for v in (*binding.range_factors, binding.plamda)) if binding else ('', '', '')),
                info.reason or f'Available ({info.units})']
            if modern:
                from qcbed.modern_yaml import target_name
                identity=target_name(info.target)
                limits=resolved.get(identity,{})
                if info.active and (limits.get('lower') is None or limits.get('upper') is None):
                    values[2]='Incomplete — bounds required'
                from .hard_bounds import configured_bounds
                lower,upper,bound_source=configured_bounds(config,info)
                values += [display_value(lower),display_value(upper)]
                if info.available and lower is None and upper is None and not info.selected:
                    values[2]='Available — no canonical hard-bound default'
            for col, value in enumerate(values):
                item = QTableWidgetItem(value)
                item.setData(Qt.ItemDataRole.UserRole, info.target)
                if col not in ((4,5,6,8,9) if modern else (4,5,6)) or not info.available:
                    item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                if col == 1:
                    item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
                    item.setCheckState(Qt.CheckState.Checked if info.selected else Qt.CheckState.Unchecked)
                    if not info.available and not info.selected:
                        item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEnabled)
                tips={4:'Dimensionless scaling factor D1; PARAMIN TOP factor. Not a hard bound.',
                      5:'Dimensionless scaling factor D2; PARAMIN BOT factor. Not a hard bound.',
                      6:'AMOEBA: scaled simplex step; physical step = PLAMDA × (D2−D1) × initial value. SIMPARA: absolute physical uncertainty width. ANNEAL and PARAMIN do not use PLAMDA.'}
                item.setToolTip(tips.get(col,info.reason or info.label))
                if modern and col in (8,9):item.setToolTip(bound_source+'; physical hard bounds, not D1/D2/PLAMDA')
                table.setItem(row, col, item)
        table.setSortingEnabled(True)
        table.resizeColumnsToContents()
        table.setColumnWidth(0, 230 if modern else 390)
        table.setMinimumHeight(230 if modern else 350)
        table.itemChanged.connect(self.refinement_draft_changed)
        form.addRow(table)
        if modern:
            for col in (4,5,6,7):table.setColumnHidden(col,True)
            table.horizontalHeaderItem(1).setText('Refine')
            table.horizontalHeaderItem(2).setText('Status')
            from .projected_elastic import ProjectedElasticPanel
            self.elastic_panel=ProjectedElasticPanel(self)
            form.addRow(self.elastic_panel)
            return
        report=QPushButton('Execution readiness details')
        report.clicked.connect(self.execution_details)
        form.addRow(report)

    def build_c_tabs(self, forms):
        """C owns physical/refinement/optimizer settings; D owns execution."""
        from .sections import ExpandableSection,move_form_row
        self.c_tabs=QTabWidget()
        self.c_guidance={};self.c_next={}
        move_form_row(forms['Numerical / Simulation cell'],self.widgets[('crystallography','depth_angstrom')],forms['Microscope'])
        boundary=self.widgets[('crystallography','slice_boundaries_fractional')]
        move_form_row(forms['Numerical / Simulation cell'],boundary,forms['Microscope'])
        boundary_help=('Normalized slice boundaries along the simulation-cell propagation direction. '
                       '0, 1 → one slice; 0, 0.5, 1 → two slices. '
                       'Start at 0, end at 1, and use strictly increasing values.')
        boundary.setToolTip(boundary_help)
        self.slice_boundary_help=QLabel(boundary_help);self.slice_boundary_help.setWordWrap(True)
        self.slice_boundary_help.hide()  # Full convention is available in the field tooltip.
        mesh_help=('Number of numerical sampling points across the transverse periodic multislice '
                   'simulation cell perpendicular to the selected zone axis. Controls projected-potential/'
                   'wavefield sampling; not the CBED detector pixel count. Columns (X) sample P; rows (Y) '
                   'sample Q. These lattice directions can be oblique.')
        for name in ('mesh_x','mesh_y'):
            self.widgets[('crystallography','mesh',name)].setToolTip(mesh_help)
        sample=getattr(self,'compiled_simulation_for_display',None)
        text='Simulation-cell dimensions unavailable — validate Simulation Setup with a complete model.'
        if sample is not None:
            g=sample.model.slices[0];mesh=sample.model.mesh
            text=(f'Committed transverse simulation cell: |P| × |Q| = {g.cell_a_angstrom:g} × {g.cell_b_angstrom:g} Å; '
                  f'γ = {g.lattice_gamma_degrees:g}°. Sampling steps along P/Q: '
                  f'{g.cell_a_angstrom/mesh.mesh_x:g} / {g.cell_b_angstrom/mesh.mesh_y:g} Å. '
                  'For an orthogonal cell these are dx/dy; revalidate after mesh edits.')
        self.transverse_mesh_summary=QLabel(text);self.transverse_mesh_summary.setWordWrap(True)
        forms['Numerical / Simulation cell'].addRow(self.transverse_mesh_summary)
        optimizer_advanced=ExpandableSection('Advanced optimizer settings')
        body=QWidget();advanced_form=QFormLayout(body)
        for spec in editable_fields(self.document.config):
            if spec.path==('refinement','optimizer') and spec.name not in ('method','ftol','maxfsimp','maxfpara'):
                move_form_row(forms['Optimizer / Run'],self.widgets[spec.key],advanced_form)
        optimizer_advanced.content.addWidget(body);forms['Optimizer / Run'].addRow(optimizer_advanced)
        self.optimizer_advanced_settings=optimizer_advanced
        groups=(('Simulation Setup',('Microscope','Numerical / Simulation cell')),
                ('Refinement Targets',('Refinement',)),
                ('Optimizer',('Optimizer / Run',)))
        self.d_sections={}
        for name in ('Project','Execution','All fields (read-only)'):
            section=ExpandableSection('Expert / Legacy settings — canonical state' if name=='All fields (read-only)' else 'Advanced execution settings' if name=='Execution' else name)
            body=forms[name].parentWidget()
            for i in range(self.pages.count()):
                old=self.pages.widget(i)
                if isinstance(old,QScrollArea) and old.widget() is body:
                    old.takeWidget();break
            section.content.addWidget(body)
            self.run_dashboard.execution_settings.addWidget(section)
            self.d_sections[name]=section
            if name=='Execution':self.execution_advanced=section
        for title,names in groups:
            body=QWidget();layout=QVBoxLayout(body)
            key={'Simulation Setup':'simulation','Refinement Targets':'refinement','Optimizer':'optimizer'}[title]
            self.c_guidance[key]=QLabel();self.c_guidance[key].setWordWrap(True)
            for name in names:
                form=forms[name]
                # Release the original section scroll area's ownership before
                # moving its content; otherwise it still resizes the C content.
                for i in range(self.pages.count()):
                    old=self.pages.widget(i)
                    if isinstance(old,QScrollArea) and old.widget() is form.parentWidget():
                        old.takeWidget();break
                form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
                form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
                if name=='Numerical / Simulation cell':
                    self.recommended_numerics=QCheckBox('Use recommended numerical settings')
                    self.recommended_numerics.setChecked(True)
                    layout.addWidget(self.recommended_numerics)
                    note=self.numerical_profile_note=QLabel('Retains the current configured mesh/cutoff; check transverse mesh and slicing convergence for your specimen.')
                    note.setWordWrap(True);self.recommended_numerics.setToolTip(note.text());note.hide()
                    self.numerical_advanced=ExpandableSection('Advanced numerical settings')
                    self.numerical_advanced.content.addWidget(form.parentWidget());layout.addWidget(self.numerical_advanced)
                    self.recommended_numerics.toggled.connect(lambda on:self.numerical_advanced.set_expanded(not on))
                else:
                    layout.addWidget(QLabel(name))
                    layout.addWidget(form.parentWidget())
            if title=='Simulation Setup':
                button=QPushButton('Validate Simulation Setup');button.setObjectName('validateSimulationSetup')
                button.clicked.connect(lambda:self.validate_c_part('simulation'))
            elif title=='Refinement Targets':
                self.refinement_notice=QLabel(('Refinement validated.' if self.workflow.refinement_receipt is not None else 'Refinement requires re-validation.') if self.workflow.simulation_receipt is not None else
                    'Validate Simulation Setup to generate refinement targets.')
                self.refinement_notice.setWordWrap(True);layout.addWidget(self.refinement_notice)
                button=self.validate_refinement_button=QPushButton('Validate Refinement Targets');button.setObjectName('validateRefinement')
                button.setEnabled(self.workflow.simulation_receipt is not None)
                button.clicked.connect(lambda:self.validate_c_part('refinement'))
                self.refinement_table.setVisible(self.workflow.simulation_receipt is not None)
                header=self.refinement_table.horizontalHeader()
                for visual,logical in enumerate((0,1,3,8,9,4,5,6,2,7)):
                    header.moveSection(header.visualIndex(logical),visual)
            elif title=='Optimizer':
                from .optimizer_panel import build
                build(self,layout,forms['Optimizer / Run'])
                button=QPushButton('Validate Optimizer')
                button.clicked.connect(lambda:self.validate_stage(2))
            from .workflow_guidance import button as nav_button
            from .responsive_actions import ResponsiveActions
            actions=[]
            if key!='simulation':actions.append(nav_button(self,'Back: '+('Simulation Setup' if key=='refinement' else 'Refinement Targets'),2,'Simulation Setup' if key=='refinement' else 'Refinement Targets'))
            next_title={'simulation':'Refinement Targets','refinement':'Optimizer','optimizer':'Run & Results'}[key]
            self.c_next[key]=nav_button(self,'Next: '+next_title,3 if key=='optimizer' else 2,None if key=='optimizer' else next_title)
            actions.extend((button,self.c_next[key]))
            layout.addStretch()
            # The scrollable content owns its minimum height; actions never share
            # that geometry or get covered by an expanding table/details panel.
            from PySide6.QtWidgets import QLayout
            layout.setSizeConstraint(QLayout.SizeConstraint.SetMinAndMaxSize)
            scroll=QScrollArea();scroll.setWidgetResizable(True);scroll.setWidget(body)
            page=QWidget();page_layout=QVBoxLayout(page);page_layout.setContentsMargins(0,0,0,0)
            page_layout.addWidget(scroll,1)
            page_layout.addWidget(self.c_guidance[key])
            footer=ResponsiveActions(actions);footer.setObjectName('cBottomActions')
            page_layout.addWidget(footer)
            self.c_tabs.addTab(page,title)
        page_index=self.pages.addWidget(self.c_tabs)
        root=self.navigation.topLevelItem(2)
        root.takeChildren()
        root.setData(0,Qt.ItemDataRole.UserRole,page_index)
        root.setData(0,Qt.ItemDataRole.UserRole+1,0)
        for tab,(title,names) in enumerate(groups):
            item=QTreeWidgetItem([title]);root.addChild(item)
            item.setData(0,Qt.ItemDataRole.UserRole,page_index)
            item.setData(0,Qt.ItemDataRole.UserRole+1,tab)
            self.section_items[title]=item
            for name in names:self.section_items[name]=item
        self.c_tabs.setTabEnabled(1,True)
        self.stage_buttons[2].setText('Validate Optimizer')
        self.stage_buttons[2].setVisible(False)

    def validate_c_part(self, part):
        from .workflow import simulation_key
        if self.controller.active or self.document is None:return False
        try:
            if any(r is None for r in self.workflow.receipts[:2]):
                raise ValueError('Validate A and B. Atomic Model first.')
            if part=='refinement' and self.workflow.simulation_receipt is None:
                raise ValueError('Validate Simulation Setup to generate refinement targets.')
            if part=='refinement' and getattr(self,'elastic_panel',None) is not None:
                if not self.elastic_panel.apply(refresh=False):
                    self.stage_status[2].setText(self.elastic_panel.message.text())
                    self.invalidate_readiness()
                    return False
            texts, bindings, bounds=self.drafts()
            sections=('Microscope','Numerical / Simulation cell') if part=='simulation' else ()
            from .document import parse_text
            other_drafts={}
            for spec in editable_fields(self.document.config):
                if spec.section in sections:continue
                old=getattr(get_object(self.document.config,spec.path),spec.name)
                try:unchanged=parse_text(texts[spec.key],spec.kind)==old
                except (ValueError,TypeError,OverflowError):unchanged=False
                if not unchanged:other_drafts[spec.key]=texts[spec.key]
            if not self.document.apply(texts,bindings,bounds,stage=part):
                self.show_diagnostics()
                if part=='simulation' and hasattr(self,'numerical_advanced'):self.numerical_advanced.set_expanded(True)
                if part=='refinement':self.refinement_notice.setText('Needs attention — '+ '; '.join(i.message for i in self.document.issues if i.severity=='error'))
                self.stage_status[2].setText("; ".join(i.message for i in self.document.issues if i.severity=="error"))
                return False
            self.workflow.sync(self.document.config)
            if part=='simulation':
                self.workflow.simulation_receipt=simulation_key(self.document.config)
                self.workflow.refinement_receipt=None
            else:
                from .workflow import refinement_key
                self.workflow.refinement_receipt=refinement_key(self.document.config)
            self.refresh()
            for key,text in other_drafts.items():
                widget=self.widgets[key]
                if isinstance(widget,QComboBox):widget.setCurrentText(text)
                else:widget.setText(text)
            self.c_tabs.setCurrentIndex(0 if part=='simulation' else 1)
            self.guidance.update_state()
            if part=='refinement':self.refinement_notice.setText('✓ Refinement Targets Valid — Next: Optimizer')
            self.stage_status[2].setText('Simulation Setup valid — validate refinement.' if part=='simulation' else
                '✓ Refinement Targets Valid — Next: Optimizer.')
            self.statusBar().showMessage(self.stage_status[2].text())
            return True
        except (ValueError,TypeError,OSError) as exc:
            self.stage_status[2].setText(str(exc));self.invalidate_readiness();return False

    def execution_details(self):
        from PySide6.QtWidgets import QDialog,QVBoxLayout,QPlainTextEdit
        from .execution_details import describe_execution
        text=describe_execution(self.document,self.refinement_table)
        dialog=QDialog(self);dialog.setWindowTitle('Execution readiness details')
        layout=QVBoxLayout(dialog);view=QPlainTextEdit();view.setReadOnly(True)
        view.setPlainText(text);layout.addWidget(view);dialog.resize(950,700)
        dialog.exec()

    def refinement_draft_changed(self,item):
        table=self.refinement_table
        row=item.row()
        modern=self.document.config.crystallography is not None
        if modern and item.column()==1 and item.checkState()==Qt.CheckState.Checked:
            from .optimizer_panel import fill_defaults
            fill_defaults(self,row)
        table.blockSignals(True)
        selected=table.item(row,1).checkState()==Qt.CheckState.Checked
        missing=[name for col,name in ((4,'D1'),(5,'D2'),(6,'PLAMDA'))
                 if not table.item(row,col).text().strip()]
        required=' and '.join((', '.join(missing[:-1]),missing[-1])) if len(missing)>1 else (missing[0] if missing else '')
        status=('Incomplete refinement binding — '+required+(' is required.' if len(missing)==1 else ' are required.')) if missing else 'Draft — Validate Refinement required'
        if selected and table.columnCount()>8:
            try:
                lo,hi=float(table.item(row,8).text()),float(table.item(row,9).text())
                initial=float(table.item(row,3).text())
                valid=lo < hi and lo <= initial <= hi
            except ValueError:valid=False
            if not valid:status+=' Hard bounds must be ordered and contain the starting value.'
        if modern:
            if not selected:status='Available'
            elif missing:status='Incomplete — advanced search settings required'
            elif not valid:
                status='Incomplete — bounds required' if not table.item(row,8).text().strip() or not table.item(row,9).text().strip() else 'Invalid — bounds must be ordered and contain the starting value'
            else:status='Ready — validate targets to activate'
        table.item(row,2).setText(status if selected or modern else 'Not selected')
        self.workflow.refinement_receipt=None
        table.blockSignals(False)
        self.mark_pending()
        if modern:
            from .optimizer_panel import refresh
            refresh(self)

    def drafts(self):
        if self.document.config.crystallography:
            from .optimizer_panel import resolve_unused
            resolve_unused(self)
        text = {key: widget.currentText() if isinstance(widget, QComboBox) else widget.text()
                for key, widget in self.widgets.items()}
        for spec in editable_fields(self.document.config):
            if spec.path==('refinement','optimizer') and not self.widgets[spec.key].isEnabled():
                text[spec.key]=display_value(getattr(get_object(self.document.config,spec.path),spec.name))
        table = self.refinement_table
        refinement = [(table.item(row, 0).data(Qt.ItemDataRole.UserRole),
            table.item(row, 1).checkState() == Qt.CheckState.Checked,
            *(table.item(row, col).text() for col in (4, 5, 6))) for row in range(table.rowCount())]
        bounds=[]
        if self.document.config.crystallography is not None:
            from .hard_bounds import configured_bounds
            catalogue={i.target:i for i in available_refinement_targets(self.document.config)} if table.rowCount() else {}
            bounds=[(table.item(r,0).data(Qt.ItemDataRole.UserRole),table.item(r,8).text(),table.item(r,9).text())
                for r in range(table.rowCount()) if table.item(r,8).text() or table.item(r,9).text() or table.item(r,1).checkState()==Qt.CheckState.Checked]
            # Retain category rules without serializing a specific copy per row.
            def inherited_unchanged(row):
                target,lo,hi=row
                if any(b.target==target for b in self.document.config.refinement.bounds.specific):return False
                lower,upper,_=configured_bounds(self.document.config,catalogue[target])
                return (lo,hi)==(display_value(lower),display_value(upper)) and lower is not None and upper is not None
            bounds=[row for row in bounds if not inherited_unchanged(row)]
        return text, refinement, bounds

    def apply_edits(self, *, stage=None):
        if stage not in (0,1,'simulation','execution') and getattr(self,'elastic_panel',None) is not None:
            if not self.elastic_panel.apply(refresh=False):return False
        if getattr(self,'site_group_pending',False):
            self.stage_status[1].setText('Apply Site-Group Parameters before validating or saving.')
            return False
        if self.document is None: return False
        if not self.pending and (not self.output_draft or stage=='execution'): return True
        fields_before=self.field_drafts() if stage is not None else {}
        if self.document.apply(*self.drafts(),stage=stage):
            self.refresh()
            self.restore_field_drafts(fields_before)
            self.statusBar().showMessage('Valid edits committed. File readiness requires Validate.')
            return True
        self.show_diagnostics()
        self.statusBar().showMessage('Draft errors — committed configuration is unchanged.')
        return False

    def validate_document(self):
        if self.apply_edits():
            self.document.validate(check_files=self.document.config.crystallography is None,execution=self.document.config.crystallography is None)
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
                 for key, value in summary.items() if key not in ('bounds','resolved_layout')]
        if summary.get('resolved_layout'):
            lines.append('\n'+summary['resolved_layout'])
        lines.append('\nResolved hard bounds (committed snapshot):')
        from .structure import resolved_bound_rows
        try:rows = () if self.document.thermal_incomplete(self.document.config) else resolved_bound_rows(self.document.config)
        except (OSError,ValueError,TypeError):rows=()
        for label, lower, upper, provenance in rows:
            lines.append(f'{label}: {lower!r} … {upper!r} ({provenance})')
        if not rows: lines.append('No active parameters')
        if self.document.files_checked:
            from .execution_details import external_checks
            lines.append('\nExternal checks:')
            lines.extend(external_checks(self.document))
        for issue in self.document.issues:
            if issue.severity=='error':
                lines.append(f'{issue.section}: {issue.message}')
        self.summary.setPlainText('\n'.join(lines))
        if hasattr(self,'run_dashboard'):self.run_dashboard.update_state()

    def navigate(self, item, previous=None):
        if item:
            self.pages.setCurrentIndex(item.data(0, Qt.ItemDataRole.UserRole))
            tab=item.data(0,Qt.ItemDataRole.UserRole+1)
            if tab is not None and hasattr(self,'c_tabs'):self.c_tabs.setCurrentIndex(tab)
            root=item.parent() or item
            index=self.navigation.indexOfTopLevelItem(root)
            if index>=0 and index!=self.workflow_tabs.currentIndex():self.workflow_tabs.setCurrentIndex(index)

    def navigate_issue(self, item, column=0):
        issue = item.data(0, Qt.ItemDataRole.UserRole)
        section = {'project':'Project', 'microscope':'Microscope', 'sample':'Regions',
            'scan':'Scan', 'experiment':'Experimental data', 'refinement':'Refinement',
            'execution':'Execution', 'run':'Outputs', 'crystallography':'Structure'}.get(issue.section)
        if section=='Experimental data' and self.document.config.crystallography is not None:section='Data Prep'
        if self.document.config.crystallography and section in ('Project','Execution','Outputs'):
            self.workflow_tabs.setCurrentIndex(3)
            if issue.field!='output_directory':
                self.d_sections['Project' if section=='Project' else 'Execution'].set_expanded(True)
        elif section: self.navigation.setCurrentItem(self.section_items[section])
        widget = self.widgets.get(issue.object_path + (issue.field,))
        if widget:
            widget.setFocus()
            if isinstance(widget, QLineEdit): widget.selectAll()
        for row in range(self.refinement_table.rowCount()):
            if self.refinement_table.item(row, 0).data(Qt.ItemDataRole.UserRole) in issue.object_path:
                self.refinement_table.selectRow(row)
                self.refinement_table.scrollToItem(self.refinement_table.item(row, 0))

    def load_structure_path(self, path):
        if not path.strip(): return False
        current = self.document.config.crystallography.structure.path
        if path == str(current): return True
        drafts=self.field_drafts()
        drafts.pop(('crystallography','structure','path'),None)
        drafts.pop(('crystallography','structure','format'),None)
        try:
            self.document.load_structure_file(path)
        except (ValueError, OSError, TypeError) as exc:
            self.report_error(f'Could not load the structure file. The previous model is unchanged.\n{exc}')
            return False
        self.refresh()
        self.restore_field_drafts(drafts)
        self.workflow_tabs.setCurrentIndex(1)
        self.navigation.setCurrentItem(self.section_items['Structure'])
        self.statusBar().showMessage('Structure loaded. Confirm zone axis and g1/g2, then Next: Configure Site Groups.')
        return True

    def browse(self, widget, kind):
        start = str(self.document.base_directory)
        if widget is self.widgets.get(('crystallography', 'structure', 'path')):
            path, _ = QFileDialog.getOpenFileName(self, 'Select structure file', start,
                'Structure files (*.cif *.extxyz *.xyz *.vasp POSCAR CONTCAR);;All files (*)')
            if path: self.load_structure_path(path)
            return
        if kind == 'directory': path = QFileDialog.getExistingDirectory(self, 'Output directory', start)
        elif kind == 'file': path, _ = QFileDialog.getOpenFileName(self, 'Select input file', start)
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
        if not self.can_save_project():return False
        if self.document.needs_save_as: return self.choose_save_as()
        try: return self.save_document(overwrite=True)
        except (OSError, ValueError, TypeError) as exc: self.report_error(exc); return False

    def choose_save_as(self):
        if not self.can_save_project():return False
        path, _ = QFileDialog.getSaveFileName(self, 'Save canonical YAML as', str(self.document.path),
            'YAML (*.yml *.yaml)', options=QFileDialog.Option.DontConfirmOverwrite)
        if not path: return False
        overwrite = Path(path).exists()
        if overwrite and QMessageBox.question(self, 'Replace YAML?', f'Replace {path}?') != QMessageBox.StandardButton.Yes: return False
        try: return self.save_document(path, overwrite=overwrite)
        except (OSError, ValueError, TypeError) as exc: self.report_error(exc); return False

    def confirm_discard(self):
        if not self.document or not (self.pending or self.output_draft or getattr(self,"site_group_pending",False) or self.document.dirty): return True
        answer = QMessageBox.question(self, 'Unsaved edits', 'Save edits before continuing?',
            QMessageBox.StandardButton.Save | QMessageBox.StandardButton.Discard | QMessageBox.StandardButton.Cancel)
        if answer == QMessageBox.StandardButton.Save: return self.save_current()
        return answer == QMessageBox.StandardButton.Discard

    def choose_open(self):
        if not self.preparation.allow_replacement():return
        if not self.confirm_discard(): return
        path, _ = QFileDialog.getOpenFileName(self, 'Load Configuration',
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
        if hasattr(self,'preparation'):self.preparation.mode.setCurrentIndex(1)
        from .workflow import Workflow
        self.workflow=Workflow()
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
        if self.document.config.crystallography is not None:
            self.navigation.setCurrentItem(self.section_items['Structure'])
            return
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
        if hasattr(self,'run_dashboard'):self.run_dashboard.update_state()

    def check_run_readiness(self):
        if self.controller.active or self.pending or self.output_draft: return
        if self.document:self.workflow.sync(self.document.config)
        if self.document is None or not self.workflow.ready or self.document.thermal_incomplete(self.document.config):
            self.invalidate_readiness()
            return
        if self.document.config.crystallography and not self.output_is_validated():
            self.invalidate_readiness()
            return
        reasons = []
        for resume, action in ((False, self.run_action), (True, self.resume_action)):
            issues = readiness(self.document.config, self.document.base_directory, resume=resume)
            errors = [diagnostic_presentation(i)[2] for i in issues if i.severity == 'error']
            if not self.workflow.ready:errors.insert(0,'Validate A → B → C before execution')
            action.setEnabled(not errors)
            action.setToolTip('\n'.join(errors) if errors else 'Ready')
            reasons.append(('Resume' if resume else 'Run') + ': ' + ('; '.join(errors) if errors else 'Ready'))
        self.run_log.appendPlainText('\n'.join(reasons))
        if hasattr(self,'run_dashboard'):self.run_dashboard.update_state(check_external=True)

    def start_run(self, resume=False):
        if self.controller.active or self.document is None or self.pending or self.output_draft: return
        self.workflow.sync(self.document.config)
        if not self.workflow.ready or (self.document.config.crystallography and not self.output_is_validated()):
            self.update_stages();self.invalidate_readiness();return
        replace_outputs=False
        if not resume and self.document.config.crystallography is None:
            from .output_paths import existing_outputs
            paths=existing_outputs(self.document.config,self.document.base_directory)
            if paths:
                answer=QMessageBox.question(self,'Replace previous outputs?',
                    'Fresh Run will use the selected paths. Existing outputs will be moved to sibling backups first:\n'+'\n'.join(map(str,paths)))
                if answer!=QMessageBox.StandardButton.Yes:return
                replace_outputs=True
        try:
            self.controller.start(self.document.config, base_directory=self.document.base_directory, resume=resume, replace_outputs=replace_outputs)
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
        self.workflow_tabs.setCurrentIndex(3)
        for index in range(3):self.workflow_tabs.setTabEnabled(index,False)
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
        for index in range(3):self.workflow_tabs.setTabEnabled(index,True)
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
        self.run_dashboard.update_state(check_external=True)
        self.statusBar().showMessage('Results updated — choose Results → Show Results to inspect.')

    def closeEvent(self, event):
        if self.preparation.active:
            self.statusBar().showMessage('Data preparation is active; wait for the worker to finish.')
            event.ignore();return
        if self.controller.active:
            self.statusBar().showMessage('A run is active. Request Cancel and wait for safe completion before closing.')
            event.ignore()
            return
        if self.confirm_discard(): event.accept()
        else: event.ignore()

    def switch_stage(self,index):
        if index==3:
            self.editor_split.hide()
            self.check_run_readiness()
            self.run_dashboard.update_state(check_external=True)
            return
        if index==0:self.preparation_split.addWidget(self.editor_split)
        else:self.workflow_tabs.widget(index).layout().addWidget(self.editor_split)
        for i in range(self.navigation.topLevelItemCount()):
            self.navigation.topLevelItem(i).setHidden(i!=index)
        item=self.navigation.currentItem()
        root=(item.parent() or item) if item else None
        if root is None or self.navigation.indexOfTopLevelItem(root)!=index:
            root=self.navigation.topLevelItem(index)
            if root:
                preferred=self.section_items.get(('Simulation Setup','Refinement Targets','Optimizer')[self.c_tabs.currentIndex()]) if index==2 and hasattr(self,'c_tabs') else None
                self.navigation.setCurrentItem(preferred or root.child(0))
        # Reparented editor widgets must be explicitly shown on the selected tab.
        self.editor_split.show()

    def update_stages(self):
        if hasattr(self,'guidance'):self.guidance.update_state()
        for i in range(3):
            valid=self.workflow.receipts[i] is not None
            if not valid or not self.stage_status[i].text().startswith(('B. Atomic Model: Valid' if i==1 else chr(65+i)+' validated')):
                self.stage_status[i].setText(chr(65+i)+(' validated' if valid else ' pending'))
            self.stage_buttons[i].setEnabled(not self.controller.active and all(v is not None for v in self.workflow.receipts[:i]))

    def validate_stage(self,index):
        if index==2 and self.document and self.document.config.crystallography:
            if self.workflow.simulation_receipt is None:
                self.stage_status[2].setText('Validate Simulation Setup before Validate C.')
                return False
            if self.workflow.refinement_receipt is None and not self.validate_c_part('refinement'):return False
        if index==2 and self.document and self.document.config.crystallography:
            from .document import parse_text
            key=('crystallography','slice_boundaries_fractional')
            widget=self.widgets.get(key)
            text=widget.text() if widget else ', '.join(map(str,self.document.config.crystallography.slice_boundaries_fractional))
            try:parse_text(text,'boundaries')
            except ValueError as exc:
                self.stage_status[index].setText(str(exc))
                self.invalidate_readiness()
                if widget:widget.setFocus()
                return False
        if self.preparation.active or (index==0 and self.preparation.dirty):
            self.stage_status[0].setText('A processing — wait for preparation to complete' if self.preparation.active
                else 'A pending — run the current preparation, inspect results, then Validate A')
            return False
        if self.document is None or self.controller.active or not self.apply_edits(stage='execution' if index==2 and self.document.config.crystallography else index):return False
        self.workflow.sync(self.document.config)
        try:
            if any(v is None for v in self.workflow.receipts[:index]):raise ValueError('Validate previous stages first')
            detail=''
            if index==0:
                from .workflow import validate_data
                detail=validate_data(self.document.config,self.document.base_directory)
            elif index==1:
                if self.document.thermal_incomplete(self.document.config):
                    raise ValueError(self.document.thermal_message(self.document.config))
                if self.document.config.crystallography:
                    from .workflow import validate_atomic_model
                    validate_atomic_model(self.document.config)
                    issues=()
                else:
                    from qcbed import validate_config,ValidationPhase
                    issues=validate_config(self.document.config,phases=(ValidationPhase.FIELD,ValidationPhase.COMPILE))
                errors=[i.message for i in issues if i.severity=='error']
                if errors:raise ValueError('; '.join(errors))
            else:
                issues=self.document.validate(check_files=self.document.config.crystallography is None,execution=self.document.config.crystallography is None)
                self.show_diagnostics()
                errors=[i.message for i in issues if i.severity=='error']
                if errors:raise ValueError('; '.join(errors))
            self.workflow.accept(index,self.document.config,base_directory=self.document.base_directory)
            self.update_stages()
            self.stage_status[index].setText('B. Atomic Model: Valid' if index==1 else chr(65+index)+' validated'+(' — '+detail if detail else ''))
            self.check_run_readiness()
            return True
        except (ValueError,OSError,TypeError,KeyError) as exc:
            self.workflow.invalidate(index);self.update_stages()
            self.stage_status[index].setText(chr(65+index)+' pending — '+str(exc))
            self.invalidate_readiness();return False
