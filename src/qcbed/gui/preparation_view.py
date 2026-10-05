"""Interactive raw-input selection. Scientific computation lives in Patch A."""
import multiprocessing
from pathlib import Path
import traceback
import uuid
from PySide6.QtCore import QTimer
from PySide6.QtWidgets import (QWidget,QVBoxLayout,QHBoxLayout,QLabel,QPushButton,QLineEdit,
    QFileDialog,QComboBox,QSpinBox)
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg,NavigationToolbar2QT
from matplotlib.patches import Ellipse,Rectangle
import numpy as np
from .preparation_model import load_image,prepare_files,apply_preparation,checksum,derive_preparation_dimensions,finalize_symmetry_roi
from qcbed.data_preparation.models import DiskGeometry,SignalMode
from qcbed.data_preparation.live import LiveSelection


def preparation_worker(arguments,connection):
    try:
        result=prepare_files(**arguments,progress=lambda message:connection.send(('progress',message)))
        connection.send(('complete',result))
    except BaseException as exc:
        message=str(exc)
        if message.startswith(('Insufficient detector room','ROI must lie inside','Fractional crop bounds')):
            message='The selected region is too close to the image boundary. Please select a slightly smaller region around the centrosymmetric feature.'
        elif message.startswith('ROI corners must lie inside'):
            message='Please select a smaller region around the centrosymmetric feature within one diffraction disk.'
        connection.send(('error',message,traceback.format_exc()))
    finally:connection.close()


class NewPreparationView(QWidget):
    def __init__(self,window):
        super().__init__();self.window=window;self.image=None;self.geometry=None;self.roi=None;self.preview=None;self.process=None;self.press=None;self.dirty=False;self.raw_document=None;self.raw_path=None;self.raw_roi=None;self.diameter_endpoints=None;self.geometry_source=None;self.geometry_input=None
        root=QHBoxLayout(self);self.controls=QWidget();self.controls.setMaximumWidth(380)
        layout=QVBoxLayout(self.controls);root.addWidget(self.controls)
        self.image_panel=QWidget();image_layout=QVBoxLayout(self.image_panel);root.addWidget(self.image_panel,1)
        self.mode=QComboBox();self.mode.addItems(['New Data Preparation','Load Previous Preparation']);layout.addWidget(self.mode)
        self.status=QLabel('Import a full CBED. Optional aperture: calibration and threshold geometry. Without one, draw a disk diameter.');self.status.setWordWrap(True);layout.addWidget(self.status)
        from .sections import ExpandableSection
        self.advanced_settings=ExpandableSection('Advanced settings / Data Processing')
        self.inputs={}
        for key,label in [('pattern','Full CBED'),('calibration_aperture','Aperture (optional)'),('disk_profile','Disk-profile override (optional)')]:
            row=QHBoxLayout();row.addWidget(QLabel(label));edit=QLineEdit();self.inputs[key]=edit;row.addWidget(edit)
            button=QPushButton('Browse…');button.clicked.connect(lambda checked=False,k=key:self.browse(k));row.addWidget(button)
            (self.advanced_settings.content if key=='disk_profile' else layout).addLayout(row)
            edit.textChanged.connect(self.invalidate)
        self.tool=QComboBox();self.tool.addItems(['Draw Disk Diameter','Select Centrosymmetric Feature']);layout.addWidget(self.tool)
        self.figure=Figure(figsize=(6,4),layout='constrained');self.canvas=FigureCanvasQTAgg(self.figure);self.canvas.setMinimumSize(520,520)
        self.toolbar=NavigationToolbar2QT(self.canvas,self.image_panel);image_layout.addWidget(self.toolbar);image_layout.addWidget(self.canvas,1)
        self.canvas.mpl_connect('button_press_event',self.begin);self.canvas.mpl_connect('button_release_event',self.end)
        self.accept=QPushButton('Accept geometry / Snap ROI preview');self.accept.clicked.connect(self.accept_preview);layout.addWidget(self.accept)
        self.derived=QLabel('Beam extraction size / symmetry crop: derived after geometry and ROI confirmation')
        self.derived.setWordWrap(True);layout.addWidget(self.derived)
        row=QHBoxLayout()
        row.addWidget(QLabel('Binning'));self.binning=QSpinBox();self.binning.setRange(1,100000);self.binning.setValue(2)
        row.addWidget(self.binning);self.binning.valueChanged.connect(self.invalidate);self.advanced_settings.content.addLayout(row);layout.addWidget(self.advanced_settings)
        self.signal=QComboBox();self.signal.addItems(['conventional','gradient']);layout.addWidget(self.signal);self.signal.currentIndexChanged.connect(self.invalidate)
        self.run=QPushButton('Run Data Preparation');self.run.clicked.connect(self.start);layout.addWidget(self.run)
        self.replay=QPushButton('Load Previous Preparation…');self.replay.clicked.connect(self.load_previous);layout.addWidget(self.replay)
        self.mode.currentIndexChanged.connect(self.switch_mode);self.replay.hide()
        self.timer=QTimer(self);self.timer.setInterval(100);self.timer.timeout.connect(self.poll)
        self.inputs['pattern'].editingFinished.connect(self.load_pattern)
        self.inputs['disk_profile'].editingFinished.connect(self.load_profile)
        self.inputs['calibration_aperture'].editingFinished.connect(self.load_profile)

    def reset_preparation(self):
        """Detach all previous A inputs; keep only explicitly separated B/C settings."""
        from dataclasses import replace
        doc=self.window.document
        if doc is not None:
            self.associated_project_path=None if doc.needs_save_as else doc.path
        if doc is not None and doc.config.crystallography is not None:
            c=doc.config
            self.retained_bc=dict(title=c.title,microscope=c.microscope,
                crystallography=replace(c.crystallography,beams=(),scan_basis_transform=((1,0),(0,1))),
                refinement=c.refinement,execution=c.execution,run=c.run)
        self.window.document=None
        self.image=None;self.geometry=None;self.roi=None;self.preview=None
        self.raw_roi=None;self.diameter_endpoints=None;self.geometry_source=None;self.geometry_input=None
        self.raw_document=None;self.raw_path=None
        while self.window.pages.count():
            page=self.window.pages.widget(0);self.window.pages.removeWidget(page);page.deleteLater()
        self.window.navigation.clear()
        self.window.data_panel=None
        self.window.summary.clear();self.window.issues.clear()
        self.window.pending=False
        self.invalidate();self.draw();self.window.update_project_path()

    def allow_replacement(self):
        if self.active:
            self.status.setText('Wait for the current preparation worker to finish.')
            return False
        if not self.dirty:return True
        from PySide6.QtWidgets import QMessageBox
        answer=QMessageBox.question(self,'Replace unfinished preparation?',
            'Discard the current preparation draft and replace it?',
            QMessageBox.StandardButton.Discard | QMessageBox.StandardButton.Cancel,
            QMessageBox.StandardButton.Cancel)
        return answer==QMessageBox.StandardButton.Discard

    def switch_mode(self,index):
        if index==0:
            if not self.allow_replacement():
                self.mode.blockSignals(True);self.mode.setCurrentIndex(1);self.mode.blockSignals(False)
                return
            self.reset_preparation()
            for widget in self.inputs.values():
                widget.blockSignals(True);widget.clear();widget.blockSignals(False)
        self.replay.setVisible(bool(index))
        self.image_panel.setVisible(not index)
        for i in range(2,self.controls.layout().count()-1):
            item=self.controls.layout().itemAt(i)
            if item.widget():item.widget().setVisible(not index)
            elif item.layout():
                for j in range(item.layout().count()):
                    if item.layout().itemAt(j).widget():item.layout().itemAt(j).widget().setVisible(not index)
        if index==1:
            self.dirty=False
            self.status.setText('Replay uses existing prepared arrays and provenance. New preparation never reads those values.')

    def open_raw(self,raw,path):
        from qcbed.modern_yaml import _mapping
        _mapping(raw,('schema_version','data_prep','atomic_model','qcbed'),('data_prep',),path='raw project')
        if type(raw.get('schema_version',2)) is not int or raw.get('schema_version',2)!=2:raise ValueError('Expected raw project schema version 2')
        prep=_mapping(raw['data_prep'],('pattern','aperture','disk_profile'),('pattern',),path='data_prep')
        self.reset_preparation()
        self.retained_bc=None
        self.raw_document=raw;self.raw_path=path
        self.associated_project_path=path
        self.window.update_project_path()
        self.signal.setCurrentText(raw.get('qcbed',{}).get('signal_mode','conventional'))
        self.window.document=None
        while self.window.pages.count():
            page=self.window.pages.widget(0);self.window.pages.removeWidget(page);page.deleteLater()
        self.window.navigation.clear();self.window.workflow.invalidate(0)
        self.mode.blockSignals(True);self.mode.setCurrentIndex(0);self.mode.blockSignals(False)
        for key,field in (('pattern','pattern'),('aperture','calibration_aperture'),('disk_profile','disk_profile')):
            self.inputs[field].setText(str((path.parent/prep[key]).resolve()) if key in prep else '')
        self.load_pattern()
        if prep.get('aperture') or prep.get('disk_profile'):self.load_profile()
        self.window.workflow_tabs.setCurrentIndex(0)

    def save_raw(self):
        try:
            from .preparation_model import save_raw_project
            from qcbed.modern_yaml import encode_modern
            path,_=QFileDialog.getSaveFileName(self,'Save raw-input project','','YAML (*.yml *.yaml)')
            if not path:return
            destination=Path(path).resolve()
            raw=self.raw_document or (encode_modern(self.window.document.config,destination) if self.window.document else {'schema_version':2})
            raw=dict(raw,data_prep={key:self.inputs[field].text() for key,field in
                (('pattern','pattern'),('aperture','calibration_aperture'),('disk_profile','disk_profile')) if self.inputs[field].text()})
            save_raw_project(raw,destination,base_directory=self.raw_path.parent if self.raw_path else (self.window.document.base_directory if self.window.document else destination.parent),overwrite=destination.exists())
            self.raw_document=raw;self.raw_path=destination
            self.status.setText('Raw-input draft saved; no prepared data or geometry invented.')
        except (ValueError,OSError) as exc:self.status.setText(str(exc))

    @property
    def active(self):return self.process is not None

    def update_derived(self):
        if self.geometry is None or self.roi is None:
            self.derived.setText('Beam extraction size / symmetry crop: waiting for confirmed disk geometry and interactive ROI')
            return
        try:
            w,h,crop=derive_preparation_dimensions(self.geometry,self.roi,self.binning.value())
            self.derived.setText(f'Derived beam extraction (expanded to binning): {w} × {h} px; symmetry crop half size: {crop} px')
        except ValueError as exc:self.derived.setText(str(exc))

    def invalidate(self,*args):
        self.update_derived()
        self.dirty=True
        self.window.workflow.invalidate(0);self.window.invalidate_readiness()
        self.window.stage_status[0].setText('A pending — new preparation inputs changed')

    def browse(self,key):
        path,_=QFileDialog.getOpenFileName(self,key,'','Images (*.dat *.txt *.npy *.tif *.tiff *.png);;All files (*)')
        if path:
            self.inputs[key].setText(path)
            if key=='pattern':self.load_pattern()
            if key in ('disk_profile','calibration_aperture'):self.load_profile()

    def load_pattern(self):
        if self.image is not None:
            try:
                if checksum(self.inputs['pattern'].text())==self.pattern_checksum:return
            except OSError:pass
        raw,path=self.raw_document,self.raw_path
        self.reset_preparation()
        self.raw_document,self.raw_path=raw,path
        try:
            self.image=load_image(self.inputs['pattern'].text());self.pattern_checksum=checksum(self.inputs['pattern'].text());self.geometry=None;self.roi=None;self.preview=None;self.raw_roi=None;self.diameter_endpoints=None;self.geometry_source=None;self.geometry_input=None;self.draw()
            self.status.setText('Draw a diameter or load a matched disk profile; then select a rough symmetry ROI. No previous geometry is reused.')
        except (ValueError,OSError) as exc:self.image=None;self.status.setText(str(exc))

    def load_profile(self):
        try:
            from qcbed.data_preparation.disk_detection import estimate_disk_geometry
            g=estimate_disk_geometry(load_image(self.inputs['disk_profile'].text() or self.inputs['calibration_aperture'].text())).geometry
            self.profile_image=load_image(self.inputs['disk_profile'].text() or self.inputs['calibration_aperture'].text())
            self.preview=('profile',g);self.draw();self.status.setText(f'Threshold geometry estimate: rx={g.radius_x_px:g}, ry={g.radius_y_px:g} px. Confirm geometry with Accept.')
        except (ValueError,OSError) as exc:self.status.setText(str(exc))

    def begin(self,event):
        if self.toolbar.mode:return
        if not self.active and event.inaxes and event.xdata is not None:self.press=(event.xdata,event.ydata)

    def end(self,event):
        if self.toolbar.mode:
            self.press=None;return
        if self.press is None or event.xdata is None or self.active:return
        a=self.press;self.press=None;b=(event.xdata,event.ydata)
        try:
            if self.tool.currentIndex()==0:
                radius=float(np.linalg.norm(np.array(b)-a)/2)
                self.preview=('manual',DiskGeometry(radius,radius),a,b)
                self.status.setText(f'Drawn diameter={2*radius:.12g} px; radius={radius:.12g} px. Accept or redraw; values are not rounded.')
            else:
                raw=(min(a[0],b[0]),min(a[1],b[1]),max(a[0],b[0]),max(a[1],b[1]))
                self.preview=('roi',finalize_symmetry_roi(raw,None if self.image is None else self.image.shape),raw)
                self.status.setText(f'ROI {raw}; integer preview {self.preview[1]}. Explicit Accept/Snap required.')
            self.draw()
        except ValueError as exc:self.status.setText(str(exc))

    def accept_preview(self):
        if not self.preview:return
        kind,value,*rest=self.preview
        if kind=='roi':
            self.roi=value;self.raw_roi=rest[0]
        else:
            self.geometry=value
            self.geometry_source=kind
            if kind=='manual':
                self.diameter_endpoints=tuple(rest);self.inputs['disk_profile'].clear()
            else:
                self.diameter_endpoints=None
                source=self.inputs['disk_profile'].text() or self.inputs['calibration_aperture'].text()
                self.geometry_input=(source,checksum(source))
        self.preview=None;self.invalidate();self.draw()
        self.status.setText(f'Geometry: {self.geometry}; symmetry ROI: {self.roi}. Extraction follows disk diameter; symmetry crop follows the accepted ROI. See derived dimensions below.')

    def draw(self):
        self.figure.clear();ax=self.figure.add_subplot()
        if self.image is not None:ax.imshow(self.image,origin='lower')
        if self.roi:ax.add_patch(Rectangle(self.roi[:2],self.roi[2]-self.roi[0],self.roi[3]-self.roi[1],fill=False,color='cyan'))
        if self.preview:
            kind,value,*rest=self.preview
            if kind=='manual':
                a,b=rest;ax.plot([a[0],b[0]],[a[1],b[1]],'y-');ax.add_patch(Ellipse(tuple((np.array(a)+b)/2),2*value.radius_x_px,2*value.radius_y_px,fill=False,color='yellow'))
            elif kind=='profile':
                ax.clear();ax.imshow(self.profile_image,origin='lower')
                from qcbed.data_preparation.disk_detection import estimate_disk_geometry
                estimate=estimate_disk_geometry(self.profile_image);y,x=np.nonzero(self.profile_image>estimate.final_threshold)
                ax.add_patch(Ellipse((float(x.mean()),float(y.mean())),2*value.radius_x_px,2*value.radius_y_px,fill=False,color='yellow'))
                ax.set_title('Threshold geometry preview (profile frame)')
            elif kind=='roi':ax.add_patch(Rectangle(value[:2],value[2]-value[0],value[3]-value[1],fill=False,color='yellow'))
        ax.set(xlabel='Detector x / column',ylabel='Detector y / row');self.canvas.draw_idle()

    def start(self):
        if self.active or self.window.controller.active:return
        try:
            if self.image is None:raise ValueError('Import a full CBED pattern first.')
            if checksum(self.inputs['pattern'].text())!=self.pattern_checksum:raise ValueError('Full pattern changed; reload and confirm selections again.')
            if self.geometry_source=='profile' and (not self.geometry_input or (self.inputs['disk_profile'].text() or self.inputs['calibration_aperture'].text())!=self.geometry_input[0] or checksum(self.geometry_input[0])!=self.geometry_input[1]):raise ValueError('Geometry source changed; estimate and confirm it again.')
            if self.geometry is None:raise ValueError('Confirm disk-profile geometry or draw and accept a disk diameter.')
            if self.roi is None:raise ValueError('Select and Accept/Snap a centrosymmetric feature ROI.')
            q=self.binning.value()
            w,h,(cx,cy)=derive_preparation_dimensions(self.geometry,self.roi,q)
            if self.window.document is None and self.raw_path is None:
                folder=QFileDialog.getExistingDirectory(self,'Choose project folder for preparation cache')
                if not folder:return
                self.raw_path=Path(folder).resolve()/'project.yml'
            base=self.window.document.base_directory if self.window.document else self.raw_path.parent
            base.mkdir(parents=True,exist_ok=True)
            profile=self.inputs['disk_profile'].text() or None
            selection=LiveSelection(self.geometry if getattr(self,'geometry_source',None)=='manual' else None,self.roi,w,h,q,SignalMode(self.signal.currentText()),(cx,cy),raw_symmetry_roi=self.raw_roi,diameter_endpoints_px=self.diameter_endpoints)
            args=dict(pattern=self.inputs['pattern'].text(),selection=selection,output=str(base/'.qcbed-data'/('preparation-'+uuid.uuid4().hex)),
                calibration_aperture=self.inputs['calibration_aperture'].text() or None,disk_profile=profile)
            context=multiprocessing.get_context('spawn');self.connection,child=context.Pipe(duplex=False)
            self.process=context.Process(target=preparation_worker,args=(args,child));self.process.start();child.close()
            self.status.setText('Preparing the current draft…')
            self.window.stage_status[0].setText('A processing — current preparation draft')
            self.setEnabled(False)
            for action in self.window.edit_actions:action.setEnabled(False)
            self.window.pages.setEnabled(False)
            self.run.setEnabled(False);self.window.workflow.invalidate(0);self.window.invalidate_readiness();self.timer.start()
        except (ValueError,OSError) as exc:self.status.setText(str(exc))

    def poll(self):
        finished=False
        while self.connection.poll():
            try:message=self.connection.recv()
            except EOFError:break
            if message[0]=='progress':self.status.setText(message[1])
            elif message[0]=='error':
                self.status.setText('Preparation failed: '+message[1]);self.window.run_log.appendPlainText(message[2]);finished=True
            elif message[0]=='complete':
                self.install(message[1]);finished=True
        if finished or not self.process.is_alive():
            if not finished:self.status.setText('Preparation worker interrupted before publishing a result.')
            self.process.join(timeout=.1);self.process=None;self.connection.close();self.timer.stop();self.run.setEnabled(True);self.setEnabled(True);self.window.pages.setEnabled(True)
            for action in self.window.edit_actions:action.setEnabled(True)

    def install(self,path):
        try:
            if self.window.document is None:
                from .preparation_model import new_preparation_document
                from .document import EditorDocument
                if self.raw_document and 'atomic_model' in self.raw_document and 'qcbed' in self.raw_document:
                    from qcbed.modern_yaml import decode_modern
                    raw=dict(self.raw_document,data_prep=dict(preparation=path),qcbed=dict(self.raw_document['qcbed'],signal_mode=self.signal.currentText()))
                    config=decode_modern(raw,self.raw_path.parent)
                else:
                    config=new_preparation_document(path)
                    if getattr(self,'retained_bc',None):
                        from dataclasses import replace
                        bc=dict(self.retained_bc)
                        bc['crystallography']=replace(bc['crystallography'],
                            beams=config.crystallography.beams)
                        config=replace(config,**bc)
                associated=getattr(self,'associated_project_path',None)
                self.window.document=EditorDocument(config,associated or Path(path).parent/'project.yml')
                self.window.document.needs_save_as=associated is None
                self.window.document.saved_config=None
            config=apply_preparation(self.window.document.config,path)
            from .data_model import restored_overview
            overview=restored_overview(config,Path(path).resolve().parent)
            if overview:
                selection=overview[1].get('live_preparation',{}).get('selection')
                if selection:
                    # Replay displays recorded values; never rederive older preparations.
                    self.roi=tuple(selection['symmetry_roi'])
                    self.geometry=DiskGeometry(**overview[1]['physical_disk_geometry'])
                    self.binning.setValue(selection['final_binning'])
                    self.derived.setText(f"Recorded beam extraction: {selection['selection_width']} × {selection['selection_height']} px; symmetry crop half size: {tuple(selection['symmetry_crop_half_size_px'])} px")
            self.dirty=False
            self.window.document.commit_preparation(config);self.window.refresh();self.mode.setCurrentIndex(1)
            self.window.navigation.setCurrentItem(self.window.section_items['Data Prep'])
            self.status.setText('Live preparation complete. Extraction uses nearest detector pixels (up to half-pixel registration offsets); physical geometry is unchanged. Inspect beam boxes, confirm the refined origin, inspect/mask Pattern / Sigma / Weight, then Validate A.')
        except (ValueError,OSError) as exc:self.status.setText(str(exc))

    def load_previous(self):
        if not self.allow_replacement():return
        path,_=QFileDialog.getOpenFileName(self,'Load Previous Preparation','','Saved preparation (*.json)')
        if path:self.install(path)
