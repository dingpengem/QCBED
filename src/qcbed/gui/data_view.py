"""Linked experimental signal/uncertainty/support inspection; no preprocessing."""
import numpy as np
from PySide6.QtWidgets import QWidget,QVBoxLayout,QHBoxLayout,QLabel,QPushButton,QComboBox,QSpinBox,QLayout,QCheckBox
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg,NavigationToolbar2QT
from matplotlib.patches import Rectangle
from matplotlib.widgets import RectangleSelector
from .data_model import read_data,edit_mask,restored_overview,provenance


class DataPrepView(QWidget):
    def __init__(self,window):
        super().__init__();self.window=window;self.data=None;self.overview=None;self.physical_labels={}
        layout=QVBoxLayout(self);layout.setSizeConstraint(QLayout.SizeConstraint.SetMinimumSize);self.summary=QLabel();self.summary.setWordWrap(True);layout.addWidget(self.summary)
        layout.addWidget(QLabel('Selected beam'));self.beams=QComboBox();layout.addWidget(self.beams)
        self.mask_layer=QComboBox()
        self.mask_layer.addItems(['Final Weight','Geometric Mask','Processing Validity','Manual Mask'])
        layout.addWidget(self.mask_layer)
        hint=QLabel('Final Weight enters chi-square. Geometry uses the existing beam-radius and binned coverage rules. '
            'Processing Validity shows additional numerical/gradient exclusions within that support '
            '(outside it, processing is neutral). Manual Mask records user exclusions. '
            'Even conventional mode retains positive-intensity validity checks.')
        hint.setWordWrap(True);layout.addWidget(hint)
        self.mask_layer.currentIndexChanged.connect(self.draw)
        self.figure=Figure(figsize=(5,10),layout='constrained');self.canvas=FigureCanvasQTAgg(self.figure);self.canvas.setMinimumHeight(850);layout.addWidget(self.canvas)
        controls=QHBoxLayout();self.box=[]
        for label in ('x start','y start','x stop (exclusive)','y stop (exclusive)'):
            controls.addWidget(QLabel(label));spin=QSpinBox();controls.addWidget(spin);self.box.append(spin)
        from .sections import ExpandableSection
        self.advanced_settings=ExpandableSection('Advanced settings / Exact pixel selection')
        self.advanced_settings.content.addLayout(controls);layout.addWidget(self.advanced_settings)
        actions=QHBoxLayout();layout.addLayout(actions)
        self.mask_buttons=[]
        for title,restore,reset in (('Exclude',False,False),('Restore Selected',True,False),('Reset Manual Mask',False,True)):
            b=QPushButton(title);b.clicked.connect(lambda checked=False,r=restore,s=reset:self.mask(restore=r,reset=s));actions.addWidget(b);self.mask_buttons.append(b)
        self.select_feature=QCheckBox('Select symmetry center on full pattern (selected reference beam)')
        layout.addWidget(self.select_feature)
        confirm=QPushButton('Confirm recorded zone-axis position')
        confirm.clicked.connect(lambda:self.select_zone())
        layout.addWidget(confirm)
        self.overview_figure=Figure(figsize=(6,4),layout='constrained');self.overview_canvas=FigureCanvasQTAgg(self.overview_figure);self.overview_canvas.setMinimumHeight(550);layout.addWidget(self.overview_canvas)
        self.overview_toolbar=NavigationToolbar2QT(self.overview_canvas,self)
        layout.insertWidget(1,self.overview_toolbar)
        layout.removeWidget(self.overview_canvas);layout.insertWidget(2,self.overview_canvas)
        layout.removeWidget(self.select_feature);layout.insertWidget(3,self.select_feature)
        layout.removeWidget(confirm);layout.insertWidget(4,confirm)
        self.overview_canvas.mpl_connect('button_press_event',self.pick_overview)
        self.canvas.mpl_connect('button_press_event',self.pick_pixel)
        try:
            doc=window.document;self.data=read_data(doc.config,doc.base_directory)
            try:self.overview=restored_overview(doc.config,doc.base_directory)
            except (ValueError,OSError,KeyError,TypeError) as exc:self.summary.setText('Overview unavailable: '+str(exc)+'\n')
            if doc.config.crystallography:
                from qcbed import compile_config_from_sources
                if window.workflow.receipts[1] is not None:
                    self.physical_labels={b.stable_id:b.original_hkl for b in compile_config_from_sources(doc.config).crystallography.indexed_beams}
            d=self.data;mode='Gradient' if d.gradient_option else 'Conventional'
            binary=np.all((d.weights==0)|(d.weights==1))
            self.summary.setText(self.summary.text()+f'{d.intensities.shape[1]} ordered beams; (tilt, beam) input {d.intensities.shape}; scan {d.scan_rows}×{d.scan_columns}.\n'
                f'Source: {doc.config.experiment.path}\n'
                f'Weight values: {np.unique(d.weights).tolist()}; display: row-major (y, x).\n'
                f'{mode} signal loaded; other mode requires its matching prepared file. Active points: {np.count_nonzero(d.weights)}.\n'
                'Pattern = configured experimental signal; sigma = sqrt(configured variance), not an inferred noise model.\n'
                'Weight = binary χ² inclusion. Exclude changes only weight; Restore uses original support. Pixel x=column, y=row.\n'
                +('' if binary else 'Nonbinary legacy weights: inspection only; masking disabled.'))
            for b in self.mask_buttons:b.setEnabled(bool(binary))
            for i,hk in enumerate(d.beam_indices_hk):
                label=self.overview[1]['beam_mapping'][i]['stable_beam_id'] if self.overview else str(hk)
                self.beams.addItem(f'{i+1}: {label}')
            for spin,limit in zip(self.box,(d.scan_columns,d.scan_rows,d.scan_columns,d.scan_rows)):spin.setRange(0,limit)
            self.box[2].setValue(1);self.box[3].setValue(1)
            self.beams.currentIndexChanged.connect(self.draw);self.draw()
        except (ValueError,OSError,TypeError,KeyError) as exc:
            self.summary.setText('Experimental inspection unavailable: '+str(exc))
            for b in self.mask_buttons:b.setEnabled(False)

    def draw(self,*args):
        if self.data is None:return
        d=self.data;i=self.beams.currentIndex();self.figure.clear();self.axes=self.figure.subplots(3,1)
        self.selectors=[]
        sigma=np.sqrt(np.where(d.variances>=0,d.variances,np.nan))
        layer=self.mask_layer.currentText()
        weight=d.weights[:,i].reshape(d.scan_rows,d.scan_columns)
        layer_image=weight
        unavailable=False
        if layer!='Final Weight':
            meta=provenance(self.window.document.config)
            if layer=='Manual Mask':
                layer_image=np.ones_like(weight)
                for point,beam in meta.get('excluded_indices',[]):
                    if beam==i:layer_image.ravel()[point]=0
                unavailable=not bool(meta)
            else:
                mapping=self.overview[1]['beam_mapping'][i] if self.overview else {}
                layers=(self.overview[1].get('live_preparation',{}).get('mask_layers',{}) if self.overview else {})
                record=layers.get(mapping.get('stable_beam_id'))
                unavailable=record is None
                if record is not None:layer_image=np.asarray(record['geometric' if layer=='Geometric Mask' else 'processing'])
        for ax,a,name in zip(self.axes,(d.intensities,sigma,d.weights),('Pattern','Sigma','Weight')):
            options = dict(interpolation='nearest')
            if name == 'Weight' and np.all((a == 0) | (a == 1)):
                options.update(cmap='gray', vmin=0., vmax=1.)
            if name=='Weight' and unavailable:
                ax.text(.5,.5,'Layer unavailable in this older preparation.\nReprepare to record it.',ha='center')
                ax.set_title(layer);continue
            ax.imshow(layer_image if name=='Weight' else a[:,i].reshape(d.scan_rows,d.scan_columns,order='C'),
                      origin='lower', **options)
            ax.set_title(layer if name=='Weight' else name);ax.set_xlabel('Column x');ax.set_ylabel('Row y')
            self.selectors.append(RectangleSelector(ax,self.select_box,useblit=False,button=[1],minspanx=.1,minspany=.1))
        self.canvas.draw_idle();self.overview_figure.clear();ax=self.overview_figure.add_subplot()
        if self.overview:
            image,m=self.overview;ax.imshow(image,origin='lower');mapping=m['beam_mapping']
            for index,beam in enumerate(mapping):
                for roi in m['rois']:
                    if roi['beam_id']!=beam['stable_beam_id']:continue
                    x,y=roi['center_px'];g=roi['geometry'];rx,ry=g['radius_x_px'],g['radius_y_px']
                    ax.add_patch(Rectangle((x-rx,y-ry),2*rx,2*ry,fill=False,edgecolor='yellow' if index==i else 'cyan'))
                    ax.text(x,y,self.beam_label(beam,index),color='yellow' if index==i else 'white',fontsize=7,ha='center')
            ax.plot(*m['reciprocal_frame']['zone_axis_origin_px'],'mo',mfc='none',markersize=12,label='Zone axis')
            ax.plot(*m['reciprocal_frame']['zone_axis_origin_px'],'m+',markersize=12)
            central=m['reciprocal_frame'].get('central_disk_center_px')
            if central is not None:ax.plot(*central,'cx',label='Central / lattice origin')
            ax.legend();ax.set_title('Full CBED — recorded ROI bounding boxes');ax.set_xlabel('Detector x / column');ax.set_ylabel('Detector y / row')
        else:ax.text(.5,.5,'Full CBED / zone-axis metadata unavailable.\nUse New Data Preparation or Load Previous Preparation.',ha='center');ax.set_axis_off()
        self.overview_canvas.draw_idle()

    def pick_pixel(self,event):
        if self.data is None or event.inaxes not in self.axes or event.xdata is None:return
        x,y=int(np.floor(event.xdata+.5)),int(np.floor(event.ydata+.5))
        if 0<=x<self.data.scan_columns and 0<=y<self.data.scan_rows:
            for spin,value in zip(self.box,(x,y,x+1,y+1)):spin.setValue(value)
            for ax in self.axes:
                for line in list(ax.lines):line.remove()
                ax.plot(x,y,'r+',markersize=10)
            self.canvas.draw_idle()

    def select_box(self,a,b):
        if None in (a.xdata,a.ydata,b.xdata,b.ydata):return
        x0,x1=sorted((a.xdata,b.xdata));y0,y1=sorted((a.ydata,b.ydata))
        for spin,value in zip(self.box,(int(np.ceil(x0)),int(np.ceil(y0)),int(np.floor(x1))+1,int(np.floor(y1))+1)):spin.setValue(value)

    def pick_overview(self,event):
        if self.overview_toolbar.mode:return
        if not self.overview or event.xdata is None:return
        if self.select_feature.isChecked():
            self.select_zone((event.xdata,event.ydata));return
        matches=[]
        for i,b in enumerate(self.overview[1]['beam_mapping']):
            for r in self.overview[1]['rois']:
                if r['beam_id']==b['stable_beam_id']:
                    x,y=r['center_px'];g=r['geometry']
                    if abs(event.xdata-x)<=g['radius_x_px'] and abs(event.ydata-y)<=g['radius_y_px']:matches.append(((event.xdata-x)**2+(event.ydata-y)**2,i))
        if matches:self.beams.setCurrentIndex(min(matches)[1])

    def commit(self,config):
        i=self.beams.currentIndex();box=tuple(spin.value() for spin in self.box)
        if self.window.document.commit_preparation(config):
            self.window.refresh();self.window.data_panel.beams.setCurrentIndex(i)
            for spin,value in zip(self.window.data_panel.box,box):spin.setValue(value)

    def mask(self,*,restore=False,reset=False):
        if self.window.controller.active or self.window.pending:return
        try:
            d=self.window.document;self.commit(edit_mask(d.config,d.base_directory,beam=self.beams.currentIndex(),box=None if reset else tuple(s.value() for s in self.box),restore=restore,reset=reset))
        except (ValueError,OSError,TypeError,KeyError) as exc:self.summary.setText(str(exc))

    def beam_label(self,beam,index):
        label=str(tuple(beam['lattice_coordinates'])) if beam.get('lattice_coordinates') is not None else beam['stable_beam_id']
        hkl=self.physical_labels.get(beam['stable_beam_id'])
        if hkl is not None:
            label+='  '+str(hkl)
        elif beam.get('physical_indices') is not None:label+='  '+str(beam['physical_indices'])
        return label

    def select_zone(self,feature=None):
        if self.window.controller.active or self.window.pending:return
        try:
            from .data_model import select_zone_axis
            d=self.window.document
            self.commit(select_zone_axis(d.config,d.base_directory,beam=self.beams.currentIndex(),feature=feature))
        except (ValueError,OSError,KeyError,TypeError) as exc:self.summary.setText(str(exc))
