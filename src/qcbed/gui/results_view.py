"""Persistent Qt workspace displaying only durable recorded results."""
import numpy as np
from PySide6.QtCore import QUrl
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
    QPlainTextEdit, QTableWidget, QTableWidgetItem, QAbstractItemView,
    QTabWidget, QComboBox, QLabel, QCheckBox)
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from .results import load_results, completed_events


def table(headers):
    widget = QTableWidget(0, len(headers))
    widget.setHorizontalHeaderLabels(headers)
    widget.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
    return widget


def fill(widget, rows):
    widget.setRowCount(len(rows))
    for i, row in enumerate(rows):
        for j, value in enumerate(row):
            widget.setItem(i, j, QTableWidgetItem('Unavailable' if value is None else str(value)))
    widget.resizeColumnsToContents()


class ResultsView(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.results = None
        self.selection = None
        self.live_status = None
        layout = QVBoxLayout(self)
        self.provenance = QLabel('No valid completed evaluation yet.')
        self.provenance.setWordWrap(True)
        layout.addWidget(self.provenance)
        controls = QHBoxLayout()
        self.evaluations = QComboBox()
        self.evaluations.currentIndexChanged.connect(self.select_evaluation)
        controls.addWidget(QLabel('Inspect evaluation:'))
        controls.addWidget(self.evaluations)
        layout.addLayout(controls)
        self.tabs = QTabWidget()
        layout.addWidget(self.tabs)
        self.summary = QPlainTextEdit()
        self.summary.setReadOnly(True)
        self.summary.setPlainText('No completed evaluation yet. Open an existing run or start a run.')
        self.tabs.addTab(self.summary, 'Overview')
        self.parameters = table(['Physical target', 'Category', 'Initial', 'Displayed evaluation value', 'Delta', 'Units', 'Resolved bounds', 'Latest', 'Best', 'Retained'])
        self.tabs.addTab(self.parameters, 'Refined Parameters')
        history_page = QWidget()
        history_layout = QVBoxLayout(history_page)
        self.history_figure = Figure(figsize=(7, 3), layout='constrained')
        self.history_canvas = FigureCanvasQTAgg(self.history_figure)
        self.history_canvas.mpl_connect('pick_event', self.pick_evaluation)
        history_layout.addWidget(self.history_canvas)
        self.history = table(['Evaluation', 'Objective', 'Stage', 'Loop', 'Event', 'Execution segment', 'Elapsed seconds'])
        history_layout.addWidget(self.history)
        self.tabs.addTab(history_page, 'χ² History')
        self.modern_history=table(['Evaluation','YNORM','Requested depth Å','Effective depth Å','Debye–Waller B Å²'])
        physical_page=QWidget();physical_layout=QVBoxLayout(physical_page)
        self.physical_target=QComboBox();physical_layout.addWidget(self.physical_target)
        self.physical_target.currentIndexChanged.connect(self.plot_physical_history)
        self.physical_figure=Figure(figsize=(7,4),layout='constrained')
        self.physical_canvas=FigureCanvasQTAgg(self.physical_figure)
        physical_layout.addWidget(self.physical_canvas);physical_layout.addWidget(self.modern_history)
        self.tabs.addTab(physical_page,'Physical / YNORM history')
        self.beam_details=table(['Beam identity','Active samples','Squared-objective contribution','Fraction of squared objective'])
        self.tabs.addTab(self.beam_details,'Beam details')
        comparison_page = QWidget()
        comparison_layout = QVBoxLayout(comparison_page)
        self.pattern_status = QLabel('Pattern unavailable.')
        self.pattern_status.setWordWrap(True)
        comparison_layout.addWidget(self.pattern_status)
        self.beams = QComboBox()
        self.beams.currentIndexChanged.connect(self.show_beam)
        comparison_layout.addWidget(self.beams)
        self.common_scale = QCheckBox('Common Experimental/Calculated scale (display only)')
        self.common_scale.setChecked(True)
        self.common_scale.toggled.connect(lambda: self.show_beam(self.beams.currentIndex()))
        comparison_layout.addWidget(self.common_scale)
        self.pattern_figure = Figure(figsize=(9, 3), layout='constrained')
        self.pattern_canvas = FigureCanvasQTAgg(self.pattern_figure)
        comparison_layout.addWidget(self.pattern_canvas)
        self.comparison = table(['Tilt / scan point', 'Experimental', 'Calculated (scaled)', 'Difference'])
        comparison_layout.addWidget(self.comparison)
        self.tabs.addTab(comparison_page, 'Pattern Matching')
        files_page = QWidget()
        files_layout = QVBoxLayout(files_page)
        files_layout.addWidget(QLabel('Open the production HTML report in your browser. Files remain unchanged.'))
        self.file_buttons = {}
        for label, name in [('Open Output Folder', ''), ('Open HTML Report', 'report.html'),
                            ('Open final output file', 'final_out.dat'),
                            ('Open checkpoint location', 'checkpoint_latest.json'),
                            ('Open current snapshot location', 'latest_evaluation.npz'),
                            ('Open retained snapshot location', 'final_snapshot.json')]:
            button = QPushButton(label)
            button.setEnabled(False)
            button.clicked.connect(lambda checked=False, n=name: self.open_artifact(n))
            self.file_buttons[name] = button
            files_layout.addWidget(button)
        files_layout.addStretch()
        self.tabs.addTab(files_page, 'Report / Files')

    def clear(self):
        self.results = None
        self.selection = None
        self.provenance.setText('No valid completed evaluation yet.')
        self.evaluations.blockSignals(True)
        self.evaluations.clear()
        self.evaluations.blockSignals(False)
        self.parameters.setRowCount(0)
        self.history.setRowCount(0)
        self.modern_history.setRowCount(0)
        self.beam_details.setRowCount(0)
        self.physical_figure.clear();self.physical_canvas.draw_idle()
        self.beams.clear()
        self.show_beam(-1)
        self.history_figure.clear()
        self.history_canvas.draw_idle()
        for button in self.file_buttons.values():
            button.setEnabled(False)

    def load(self, directory, *, evaluation=None, status=None):
        result = load_results(directory, evaluation=evaluation)
        old_beam = self.beams.currentIndex()
        self.results, self.selection, self.live_status = result, evaluation, status
        summary = result.summary
        lines = ['PARTIAL / latest-valid' if summary['partial'] or status in ('Running', 'Cancel pending', 'Cancelled', 'Interrupted') else 'Completed — recorded results']
        lines.append(f"Latest valid completed evaluation: {summary['latest_valid_evaluation']}")
        self.provenance.setText(
            f"Latest valid completed evaluation: {summary['latest_valid_evaluation']} — "
            f"displaying {summary['displayed_evaluation']} — {status or summary['status']}"
            + (' — partial/latest-valid' if summary['partial'] or status in ('Running', 'Cancel pending', 'Cancelled', 'Interrupted') else ''))
        lines.append(f"Displayed evaluation: {summary['displayed_evaluation']}" +
                     (' (historical)' if summary['historical'] else ''))
        for label, value in (
                ('Run status', status or summary['status']),
                ('Objective / χ² (recorded objective definition)', summary['retained_objective']),
                ('YNORM', summary['retained_ynorm']), ('Stage', summary['recorded_stage']),
                ('Attempted evaluation count', summary['recorded_evaluations']),
                ('Elapsed seconds (recorded)', summary['recorded_elapsed_seconds']),
                ('Output directory', str(result.directory)),
                ('Checkpoint', summary['checkpoint_status']), ('Final-state artifact', summary['final_state_status'])):
            lines.append(f'{label}: {value}')
        if 'signal_mode' in summary:
            for key in ('signal_mode','initial_parameters','requested_depth','effective_depth','best','last_evaluated','retained_final','optimizer_result','convergence'):
                value=summary.get(key)
                if isinstance(value,dict):
                    value={k:v for k,v in value.items() if k in ('requested_depth','effective_depth','B','objective','ynorm','stage','global_evaluation_index','stop_reason','converged','convergence','evaluations','iterations')}
                lines.append(f'{key}: {value}')
        lines.extend('Note: ' + warning for warning in result.warnings)
        self.summary.setPlainText('\n'.join(lines))
        events = completed_events(result.history)
        self.evaluations.blockSignals(True)
        self.evaluations.clear()
        self.evaluations.addItem('Latest valid (follow live)', None)
        for row in events:
            index = row['global_evaluation_index']
            self.evaluations.addItem(f"Evaluation {index} — {row.get('stage', '')}", index)
        self.evaluations.setCurrentIndex(max(0, self.evaluations.findData(evaluation)))
        self.evaluations.blockSignals(False)
        fill(self.parameters, [(p['label'], p['category'], p['initial'], p['final'], p['delta'], p['units'],
                               p['bounds'] if p['bounds'] else 'Unbounded / not recorded',p.get('latest'),p.get('best'),p.get('retained')) for p in result.parameters])
        fill(self.history, [(r['global_evaluation_index'], r.get('objective'), r.get('stage'),
                            r.get('loop_index'), 'BOUND REJECTION' if r['classification'] == 'bound_rejection' else r['classification'],
                            r.get('restart_count'), r.get('elapsed_seconds')) for r in result.history])
        self.modern_history.setColumnCount(2+len(result.parameters))
        self.modern_history.setHorizontalHeaderLabels(['Evaluation','YNORM']+[p['label'] for p in result.parameters])
        fill(self.modern_history,[(r['global_evaluation_index'],r.get('ynorm'),*r.get('physical_active_vector',r.get('active_vector',[]))) for r in events])
        previous=self.physical_target.currentText()
        self.physical_target.blockSignals(True);self.physical_target.clear()
        self.physical_target.addItems([p['label'] for p in result.parameters])
        self.physical_target.setCurrentIndex(max(0,self.physical_target.findText(previous)))
        self.physical_target.blockSignals(False)
        self.plot_physical_history()
        comparison=result.comparison or {}
        fill(self.beam_details,[(str(m),b['active_count'],b['squared_objective_contribution'],b['fraction_of_total_error'])
            for m,b in zip(comparison.get('beam_metadata',()),comparison.get('per_beam',()))])
        self.plot_history(events)
        self.beams.blockSignals(True)
        self.beams.clear()
        if result.comparison is not None:
            for index, beam in enumerate(result.comparison['beams']):
                metadata=result.comparison.get('beam_metadata')
                self.beams.addItem((f"{metadata[index]['stable_id']} — original hkl {metadata[index]['original_hkl']}" if metadata else f'Beam ({beam[0]}, {beam[1]}) — ordered entry {index + 1}'))
        self.beams.setCurrentIndex(max(0, min(old_beam, self.beams.count()-1)))
        self.beams.blockSignals(False)
        self.show_beam(self.beams.currentIndex())
        for name, button in self.file_buttons.items():
            button.setEnabled((result.directory / name).exists())

    def plot_physical_history(self, *_):
        if not self.results: return
        result=self.results;index=self.physical_target.currentIndex()
        rows=[r for r in completed_events(result.history) if r.get('physical_active_vector') is not None]
        self.physical_figure.clear()
        first,second=self.physical_figure.subplots(2,1)
        if 0<=index<len(result.parameters):
            first.plot([r['global_evaluation_index'] for r in rows],[r['physical_active_vector'][index] for r in rows],'.-')
            first.set_ylabel(result.parameters[index]['label'])
        second.plot([r['global_evaluation_index'] for r in rows],[r.get('ynorm') for r in rows],'.-')
        second.set_ylabel('YNORM');second.set_xlabel('Recorded evaluation')
        self.physical_canvas.draw_idle()

    def select_evaluation(self, index):
        if self.results is not None:
            self.load(self.results.directory, evaluation=self.evaluations.itemData(index), status=self.live_status)

    def pick_evaluation(self, event):
        if event.artist is getattr(self, '_points', None) and len(event.ind):
            number = self._point_evaluations[event.ind[0]]
            self.evaluations.setCurrentIndex(self.evaluations.findData(number))

    def plot_history(self, events):
        self.history_figure.clear()
        ax = self.history_figure.add_subplot()
        self._point_evaluations = [r['global_evaluation_index'] for r in events]
        values = [r['objective'] for r in events]
        ax.plot(self._point_evaluations, values, linewidth=.8)
        self._points = ax.scatter(self._point_evaluations, values, s=15, picker=5, label='Completed evaluation')
        previous = None
        for row in events:
            boundary = (row.get('restart_count'), row.get('loop_index'), row.get('stage'))
            if previous is not None and boundary != previous:
                ax.axvline(row['global_evaluation_index'], color='grey', linestyle=':', alpha=.6)
            previous = boundary
        rejected = [r['global_evaluation_index'] for r in self.results.history if r['classification'] == 'bound_rejection']
        if rejected:
            ax.scatter(rejected, [.96]*len(rejected), transform=ax.get_xaxis_transform(),
                       marker='x', color='red', label='Bound rejection (no objective; marker at top)')
        if rejected:
            left, right = ax.get_xlim()
            ax.set_xlim(min(left, min(rejected)-.1), max(right, max(rejected)+.1))
        ax.axvline(self.results.summary['displayed_evaluation'] or 0, color='tab:blue', alpha=.25)
        ax.set_xlabel('Evaluation number — dotted lines: stage/loop/resume boundaries')
        ax.set_ylabel('Recorded objective / χ²')
        ax.legend(fontsize='small')
        self.history_canvas.draw_idle()

    def show_beam(self, index):
        self.pattern_figure.clear()
        self.comparison.setRowCount(0)
        if self.results is None or self.results.comparison is None or index < 0:
            self.pattern_status.setText('No stored pattern for this evaluation. Propagation is never recomputed.')
            self.pattern_canvas.draw_idle()
            return
        data = self.results.comparison
        experimental, calculated = data['experimental'][:, index], data['scaled'][:, index]
        residual = calculated-experimental
        fill(self.comparison, [(i, x, y, y-x) for i, (x, y) in enumerate(zip(experimental, calculated))])
        shape = data.get('scan_shape')
        if shape is None or len(shape) != 2 or int(np.prod(shape)) != len(experimental):
            self.pattern_status.setText('Stored scan dimensions unavailable; numeric values shown without inventing geometry.')
            self.pattern_canvas.draw_idle()
            return
        shape = tuple(int(v) for v in shape)
        coordinates = None
        if 'scan' in data:
            from qcbed.scan import LegacyScanDefinition, generate_scan_tilts
            try:
                coordinates = generate_scan_tilts(LegacyScanDefinition(**data['scan'])).reshape(*shape, 2)
            except (ValueError, TypeError):
                coordinates = None
        self.pattern_status.setText('Stored scaled objective-space intensity (includes GRADIENT if used). '
            'Residual = calculated − experimental; excluded support is NaN/hidden; colors affect display only. '
            + ('Axes: undilated reciprocal H/K.' if coordinates is not None else 'Axes: stored scan column/row; physical geometry unavailable.'))
        finite=np.concatenate((experimental[np.isfinite(experimental)],calculated[np.isfinite(calculated)]))
        lo,hi=(float(finite.min()),float(finite.max())) if finite.size else (0.,1.)
        finite_residual=residual[np.isfinite(residual)]
        limit=max(float(np.max(np.abs(finite_residual))) if finite_residual.size else 0.,np.finfo(float).eps)
        for ax, values, title in zip(self.pattern_figure.subplots(1, 3),
                                    (experimental, calculated, residual),
                                    ('Experimental', 'Calculated', 'Residual')):
            options = dict(cmap='coolwarm' if title == 'Residual' else 'viridis')
            if title == 'Residual':
                options.update(vmin=-limit, vmax=limit)
            elif self.common_scale.isChecked():
                options.update(vmin=lo, vmax=hi)
            if coordinates is not None and min(shape) > 1:
                artist = ax.pcolormesh(coordinates[..., 0], coordinates[..., 1], values.reshape(shape),
                                       shading='nearest', **options)
                ax.set_xlabel('H'); ax.set_ylabel('K')
            else:
                artist = ax.imshow(values.reshape(shape), origin='lower', **options)
                ax.set_xlabel('Scan column'); ax.set_ylabel('Scan row')
            ax.set_title(title)
            self.pattern_figure.colorbar(artist, ax=ax)
        self.pattern_canvas.draw_idle()

    def open_artifact(self, name):
        if self.results is None or name not in self.file_buttons:
            return False
        path = self.results.directory / name
        if not path.exists():
            return False
        if name.endswith(('.json', '.npz')):
            path = path.parent
        return QDesktopServices.openUrl(QUrl.fromLocalFile(str(path)))
