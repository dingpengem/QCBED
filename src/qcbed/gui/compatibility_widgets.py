"""Optional Qt forms for one-time legacy migration and explicit export review."""
from pathlib import Path
import json
from PySide6.QtWidgets import (QDialog, QFormLayout, QVBoxLayout, QHBoxLayout,
    QLineEdit, QComboBox, QLabel, QDialogButtonBox, QWidget, QPushButton,
    QFileDialog, QPlainTextEdit, QScrollArea)
from qcbed.legacy_config import MigrationPolicy
from qcbed.configuration import RunConfig, BoundsConfig
from qcbed.execution import ForwardExecutionConfig
from .document import parse_text


class ImportDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Import Legacy Project — explicit migration policy')
        self.resize(760,740)
        self.result_values = None
        outer = QVBoxLayout(self)
        note = QLabel('Legacy files do not define Python execution policies. Supply every choice below. '
                      'Import creates an unsaved canonical document; it does not write files or start a run. '
                      'Relative paths use the FOR007 directory.')
        note.setWordWrap(True)
        outer.addWidget(note)
        body = QWidget()
        form = QFormLayout(body)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(body)
        outer.addWidget(scroll)
        self.controls = {}
        fields = [('source','FOR007 input','file'),('experimental','Header-bearing experimental input','file'),
            ('seed','Random seed','int'),('execution','Execution backend',('serial','batched')),
            ('batch_size','Batch size (required for batched; blank for serial)','optional_int'),
            ('objective','Objective policy',('legacy','masked')),
            ('invalid','Invalid-trial policy',('raise','reject')),
            ('bounds','Hard bounds policy',('none',)),
            ('output_directory','Canonical run output directory','directory'),
            ('output_filename','Canonical run output filename','output'),
            ('snapshot_every','Snapshot interval (0 disables)','int'),
            ('checkpoint_every','Checkpoint interval','int'),
            ('max_evaluations','Evaluation budget (0 unlimited)','int'),
            ('max_seconds','Wall-time seconds (0 unlimited)','real')]
        for key,label,kind in fields:
            if isinstance(kind,tuple):
                widget = QComboBox()
                widget.addItem('Choose explicitly…',None)
                for value in kind: widget.addItem(value,value)
            else: widget = QLineEdit()
            self.controls[key] = widget
            if kind in ('file','directory','output'):
                row = QWidget(); layout = QHBoxLayout(row); layout.setContentsMargins(0,0,0,0)
                layout.addWidget(widget)
                button = QPushButton('Browse…')
                button.clicked.connect(lambda checked=False,w=widget,k=kind:self.browse(w,k))
                layout.addWidget(button)
                form.addRow(label,row)
            else: form.addRow(label,widget)
        self.error = QLabel(); self.error.setWordWrap(True); outer.addWidget(self.error)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.button(QDialogButtonBox.StandardButton.Ok).setText('Import canonical document')
        buttons.accepted.connect(self.submit)
        buttons.rejected.connect(self.reject)
        outer.addWidget(buttons)

    def browse(self,widget,kind):
        start = str(Path.home())
        if kind == 'directory': path = QFileDialog.getExistingDirectory(self,'Choose directory',start)
        elif kind == 'file': path,_ = QFileDialog.getOpenFileName(self,'Choose legacy input',start)
        else: path,_ = QFileDialog.getSaveFileName(self,'Choose output path (not written)',start)
        if path: widget.setText(path)

    def values(self):
        value = {}
        for key,control in self.controls.items():
            value[key] = control.currentData() if isinstance(control,QComboBox) else control.text()
            if key != 'batch_size' and (value[key] is None or value[key] == ''):
                raise ValueError(f'{key}: an explicit choice/value is required')
        source = Path(value['source'])
        if not source.is_absolute(): raise ValueError('Choose an absolute FOR007 input path')
        base = source.parent
        path = lambda name: str((base/value[name]).resolve())
        batch = parse_text(value['batch_size'],'optional_int')
        if value['execution'] == 'batched' and batch is None:
            raise ValueError('Batched migration requires an explicit batch size')
        snapshot = parse_text(value['snapshot_every'],'int')
        evaluations = parse_text(value['max_evaluations'],'int')
        seconds = parse_text(value['max_seconds'],'real')
        if snapshot < 0 or evaluations < 0 or seconds < 0: raise ValueError('Intervals/budgets cannot be negative')
        policy = MigrationPolicy(parse_text(value['seed'],'int'),ForwardExecutionConfig(value['execution'],batch),
            RunConfig(path('output_directory'),path('output_filename'),snapshot or None,
                      parse_text(value['checkpoint_every'],'int'),evaluations or None,seconds or None),
            value['objective'],value['invalid'],BoundsConfig(value['bounds']))
        return source,Path(path('experimental')),policy

    def submit(self):
        try: self.result_values = self.values()
        except (ValueError,TypeError) as exc:
            self.error.setText(str(exc)); return
        self.accept()


class ExportDialog(QDialog):
    def __init__(self, base, *, bundle, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Export Legacy Bundle' if bundle else 'Export FOR007')
        self.base, self.bundle = base,bundle
        self.result_values = None
        layout = QFormLayout(self)
        self.destination = QLineEdit()
        row = QWidget(); r = QHBoxLayout(row); r.setContentsMargins(0,0,0,0); r.addWidget(self.destination)
        browse = QPushButton('Browse…'); r.addWidget(browse)
        def choose():
            path,_ = QFileDialog.getSaveFileName(self,'Choose NEW bundle directory' if bundle else 'Choose NEW FOR007 file',str(base),
                                                 options=QFileDialog.Option.DontConfirmOverwrite)
            if path: self.destination.setText(path)
        browse.clicked.connect(choose)
        layout.addRow('New directory' if bundle else 'New file',row)
        note = QLabel('Exports are exclusive: choose a new destination; existing files/directories are never replaced. '
                      'A bundle writes FOR007.dat, a separate header-bearing in.dat and compatibility.json. '
                      'Standalone export requires header-validation data and verified executable filename widths.')
        note.setWordWrap(True); layout.addRow(note)
        self.exp_width,self.out_width = QLineEdit(),QLineEdit()
        if not bundle:
            layout.addRow('Verified EXPINT filename width',self.exp_width)
            layout.addRow('Verified BMNAME filename width',self.out_width)
        self.error = QLabel(); self.error.setWordWrap(True); layout.addRow(self.error)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.button(QDialogButtonBox.StandardButton.Ok).setText('Review compatibility…')
        buttons.accepted.connect(self.submit); buttons.rejected.connect(self.reject)
        layout.addRow(buttons)

    def submit(self):
        try:
            if not self.destination.text(): raise ValueError('Choose a new destination')
            destination = (self.base/self.destination.text()).resolve()
            if destination.exists(): raise FileExistsError('Destination already exists; choose a new path')
            limits = None if self.bundle else (parse_text(self.exp_width.text(),'int'),parse_text(self.out_width.text(),'int'))
            if limits and min(limits) < 1: raise ValueError('Filename widths must be positive')
            self.result_values = destination,limits
        except (OSError,ValueError,TypeError) as exc: self.error.setText(str(exc)); return
        self.accept()


class ReportDialog(QDialog):
    def __init__(self, report, *, exporting=False, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Legacy compatibility report')
        self.resize(850,650)
        layout = QVBoxLayout(self)
        self.text = QPlainTextEdit(); self.text.setReadOnly(True)
        sections = ['FAITHFULLY REPRESENTABLE (only if export validation succeeds)', *report['represented'],
            '\nOMITTED — FOR007 has no representation / compatibility caveats',
            json.dumps(report['omissions'],indent=2,default=lambda x:getattr(x,'value',str(x))),
            '\nUNSUPPORTED / UNEXPORTABLE']
        sections.extend(f'{i.code} · {i.field}: {i.message}' for i in report['errors'])
        if not report['errors']: sections.append('No representability errors for the selected mode.')
        self.text.setPlainText('\n'.join(sections)); layout.addWidget(self.text)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Cancel | QDialogButtonBox.StandardButton.Ok)
        buttons.button(QDialogButtonBox.StandardButton.Ok).setText('Export acknowledging omissions' if exporting else 'Close')
        buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(not exporting or not report['errors'])
        buttons.accepted.connect(self.accept); buttons.rejected.connect(self.reject); layout.addWidget(buttons)
