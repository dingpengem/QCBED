"""Dedicated physical E_g editor; no independent scattering calculation."""
from dataclasses import replace
from cmath import phase
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QWidget,QVBoxLayout,QLabel,QTableWidget,QTableWidgetItem,QComboBox,QPushButton,QSizePolicy,QGroupBox,QGridLayout,QLineEdit,QInputDialog,QMessageBox,QAbstractItemView,QPlainTextEdit
from qcbed.configuration import ProjectedElasticConfig
from qcbed.projected_elastic import catalogue,components,coordinates,complex_value
from qcbed.crystallography_config import compile_config_from_sources


def display(z): return f'{z.real:.17g} {z.imag:+.17g}i'


class ProjectedElasticPanel(QWidget):
    def __init__(self, window):
        super().__init__();self.window=window;self.document_identity=window.document;self.dirty=False;self.resets=set()
        self.setSizePolicy(QSizePolicy.Policy.Expanding,QSizePolicy.Policy.Preferred)
        layout=QVBoxLayout(self)
        from .sections import ExpandableSection
        from .responsive_actions import ResponsiveActions
        from PySide6.QtWidgets import QLayout
        layout.setSizeConstraint(QLayout.SizeConstraint.SetMinAndMaxSize)
        layout.addWidget(QLabel('Projected Elastic Coefficients'))
        self.help_section=ExpandableSection('Details — coefficient convention and sharing')
        label=QLabel('Projected Elastic Coefficients\nProjected elastic coefficient E_g⁽ˢ⁾: the slice-resolved transverse Fourier coefficient of the elastic potential in the zone-axis-oriented multislice simulation cell. '
            'Reflections are displayed using their original crystallographic hkl identities. '
            'IAM E_g is calculated from the atomic model for each physical slice. '
            'E_g is separate from absorption A_g; these are not conventional 3D unit-cell structure factors. '
            'Friedel partners are linked. Complex bounds use real, imaginary order; no universal bounds are supplied.')
        label.setWordWrap(True);self.help_section.content.addWidget(label)
        self.table=QTableWidget(0,11);layout.addWidget(self.table)
        self.table.setHorizontalHeaderLabels(['Refine','Reflection','Slice','IAM E_g','Mode','Starting value','Lower','Upper','Status','Reset to IAM','Parameter sharing'])
        self.details=QGroupBox('Selected coefficient details');details_layout=QGridLayout(self.details)
        self.details_text=QLabel('Select a reflection to inspect its components.');self.details_text.setWordWrap(True)
        details_layout.addWidget(self.details_text,0,0,1,3)
        details_layout.addWidget(QLabel('Lower'),1,1);details_layout.addWidget(QLabel('Upper'),1,2)
        self.bound_labels=[];self.bound_fields={}
        for component in range(2):
            name=QLabel();self.bound_labels.append(name);details_layout.addWidget(name,component+2,0)
            for col in (6,7):
                field=QLineEdit();self.bound_fields[component,col]=field
                field.setAccessibleName(f'{"Real" if component==0 else "Imaginary"} {"lower" if col==6 else "upper"} bound')
                details_layout.addWidget(field,component+2,col-5)
                field.textEdited.connect(lambda text,component=component,col=col:self.edit_bound(component,col,text))
        self.row_details=ExpandableSection('Selected coefficient details')
        self.row_details.content.addWidget(self.details)
        reset_selected=QPushButton('Reset selected to IAM')
        reset_selected.clicked.connect(lambda:self.reset(self.table.currentRow()) if self.table.currentRow()>=0 and not self.rows[self.table.currentRow()]['reason'] else None)
        self.row_details.content.addWidget(reset_selected)
        layout.addWidget(self.row_details)
        buttons=[]
        for col,(title,action) in enumerate((('Share One Refinement Parameter',self.link_selected),('Make Independent',self.unlink_selected))):
            button=QPushButton(title);button.clicked.connect(action);buttons.append(button)
        self.sharing_instructions=QLabel('1. Check Refine for the targets.  2. Select two or more compatible rows (Ctrl/Cmd or Shift).\n3. Share One Refinement Parameter.  4. Check the group below.  5. Validate Refinement Targets.')
        self.sharing_instructions.setWordWrap(True);layout.insertWidget(1,self.sharing_instructions)
        self.sharing_semantics=QLabel('Shared targets use common optimizer coordinates and stay tied during refinement. Use Make Independent to separate them again.')
        self.sharing_semantics.setWordWrap(True);layout.insertWidget(2,self.sharing_semantics)
        layout.insertWidget(3,ResponsiveActions(buttons))
        self.sharing_example=ExpandableSection('Example — sharing (200) and (020)')
        example=QLabel('Real-only (200) and (020) share one Re(E_g⁽ˢ⁾) coordinate.\nIndependent: 2 parameters → Shared: 1 parameter.\nWith sample depth: Independent total: 3 → Shared total: 2.\nReal-only: shared real component; each imaginary part stays fixed.\nComplex: shared real and imaginary components (two coordinates).\nMagnitude-only: shared magnitude; each phase stays fixed.\nShared refinement parameter: user-selected constraint.\nFriedel partner: automatic physical constraint, always preserved.')
        example.setWordWrap(True);self.sharing_example.content.addWidget(example);layout.addWidget(self.sharing_example)
        self.group_summary=QLabel();self.group_summary.setWordWrap(True);layout.insertWidget(4,self.group_summary)
        help_text=QLabel('Parameter sharing reduces the number of independent optimizer variables by making selected targets use the same refinement coordinate throughout refinement. Shared real component: Re stays equal, each Im remains fixed. Shared complex coefficient: one common real and one common imaginary parameter. Shared magnitude: equal magnitude, each IAM phase remains fixed. Same slice required. Mandatory Friedel partners are separate and do not add optimizer parameters. Select rows with Ctrl/Cmd or Shift.')
        help_text.setWordWrap(True);self.help_section.content.addWidget(help_text);layout.addWidget(self.help_section)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.dimension_summary=QLabel();self.dimension_summary.setWordWrap(True);layout.insertWidget(5,self.dimension_summary)
        self.live_groups=QLabel();self.live_groups.setWordWrap(True)
        self.group_details=ExpandableSection('Shared parameter groups — details')
        self.group_details.content.addWidget(self.live_groups);layout.addWidget(self.group_details)
        self.message=QLabel();self.message.setWordWrap(True);layout.insertWidget(6,self.message)
        self.apply_button=QPushButton('Apply Projected Elastic Coefficients');self.apply_button.clicked.connect(self.apply)
        self.help_section.content.addWidget(self.apply_button)
        clear=QPushButton('Clear saved coefficient targets');clear.clicked.connect(self.clear_targets);self.help_section.content.addWidget(clear)
        try:self.rows=catalogue(compile_config_from_sources(window.document.config)) if window.workflow.simulation_receipt is not None else ()
        except (ValueError,TypeError,OSError) as exc:self.rows=();self.message.setText(str(exc))
        self.table.setRowCount(len(self.rows))
        for i,row in enumerate(self.rows):
            record=row['record'];mode=row['mode'];row['_editor_mode']=mode
            values=['',str(row['hkl']),str(row['interval']),display(row['iam']),'',self.start_text(row['initial'],mode),
                ', '.join('' if v is None else str(v) for v in record.lower[:len(components(mode))]) if record else '',
                ', '.join('' if v is None else str(v) for v in record.upper[:len(components(mode))]) if record else '',
                row['reason'] or ('Stale IAM' if row['stale'] else 'Available')]
            for j,value in enumerate(values):
                item=QTableWidgetItem(value)
                item.setToolTip(value)
                if j not in (5,6,7):item.setFlags(item.flags()&~Qt.ItemFlag.ItemIsEditable)
                if j==0:
                    item.setFlags(item.flags()|Qt.ItemFlag.ItemIsUserCheckable)
                    item.setCheckState(Qt.CheckState.Checked if record and record.refine else Qt.CheckState.Unchecked)
                if row['reason']:item.setFlags(item.flags()&~Qt.ItemFlag.ItemIsEnabled)
                self.table.setItem(i,j,item)
            self.table.setItem(i,10,QTableWidgetItem(('Shared: '+record.equality_group) if record and record.equality_group else 'Independent'))
            self.table.item(i,10).setFlags(self.table.item(i,10).flags() & ~Qt.ItemFlag.ItemIsEditable)
            combo=QComboBox();combo.addItem('Real only','real');combo.addItem('Complex (Re + Im)','complex');combo.addItem('Magnitude only','magnitude')
            combo.setCurrentIndex(combo.findData(mode));combo.setEnabled(not row['reason']);self.table.setCellWidget(i,4,combo)
            combo.currentIndexChanged.connect(lambda _,i=i:self.mode_changed(i))
            reset=QPushButton('Reset to IAM');reset.setEnabled(not row['reason']);reset.clicked.connect(lambda _,i=i:self.reset(i));self.table.setCellWidget(i,9,reset)
        self.table.currentCellChanged.connect(lambda *_:self.refresh_details())
        self.table.itemChanged.connect(self.changed)
        self.table.resizeColumnsToContents();self.table.setMinimumHeight(230)
        self.table.setColumnHidden(3,True);self.table.setColumnHidden(9,True)
        for col,width in ((1,100),(2,90),(4,160),(5,135),(6,90),(7,90),(8,110),(10,130)):self.table.setColumnWidth(col,width)
        self.table.setToolTip('Hard bounds use physical components. Complex bounds: real, imaginary. Select a row and expand Selected coefficient details for IAM and Friedel information.')
        self.window.refinement_table.itemChanged.connect(self.update_dimensions)
        self.update_dimensions()
        self.resolved=QPlainTextEdit();self.resolved.setReadOnly(True);self.resolved.setMinimumHeight(220)
        self.layout_details=ExpandableSection('Resolved Refinement Parameter Layout')
        self.layout_details.content.addWidget(self.resolved);layout.addWidget(self.layout_details)
        if window.workflow.refinement_receipt is not None:
            from .resolved_layout import resolved_layout
            try:self.resolved.setPlainText(resolved_layout(compile_config_from_sources(window.document.config))[1])
            except (ValueError,TypeError,OSError) as exc:self.resolved.setPlainText('Layout unavailable: '+str(exc))
        else:self.resolved.setPlainText('Resolved Refinement Parameter Layout — Validate Refinement Targets to resolve the current selection.')

    def refresh_preserving_drafts(self):
        w=self.window;fields=w.field_drafts();table=w.refinement_table
        drafts={table.item(i,0).data(Qt.ItemDataRole.UserRole):[(col,table.item(i,col).text(),table.item(i,col).checkState()) for col in (1,4,5,6,8,9)] for i in range(table.rowCount())}
        pending=w.pending;selected={(self.rows[i]['hkl'],self.rows[i]['interval']) for i in self.selected_rows()}
        w.refresh();w.restore_field_drafts(fields);table=w.refinement_table;table.blockSignals(True)
        for i in range(table.rowCount()):
            for col,text,check in drafts.get(table.item(i,0).data(Qt.ItemDataRole.UserRole),[]):
                table.item(i,col).setText(text)
                if col==1:table.item(i,col).setCheckState(check)
        table.blockSignals(False)
        from PySide6.QtCore import QItemSelectionModel
        for i,row in enumerate(w.elastic_panel.rows):
            if (row['hkl'],row['interval']) in selected:w.elastic_panel.table.selectionModel().select(w.elastic_panel.table.model().index(i,0),QItemSelectionModel.SelectionFlag.Select|QItemSelectionModel.SelectionFlag.Rows)
        w.elastic_panel.update_dimensions()
        if pending:w.mark_pending(stage_override=2)

    def selected_rows(self):
        return sorted({index.row() for index in self.table.selectionModel().selectedRows()})

    def link_rows(self,indices,name,*,use_first=False,preview=False):
        from qcbed.projected_elastic import coordinates,complex_value
        import numpy as np
        if len(indices)<2:raise ValueError('Select at least two coefficient rows')
        for i in indices:
            if self.rows[i]['reason']:raise ValueError(self.rows[i]['reason'])
            self.table.item(i,0).setCheckState(Qt.CheckState.Checked)
        records=list(self.collect());keys={(self.rows[i]['hkl'],self.rows[i]['interval']) for i in indices}
        c=self.window.document.config
        candidate=replace(c,crystallography=replace(c.crystallography,projected_elastic=tuple(records)))
        rows=catalogue(compile_config_from_sources(candidate));chosen=[r for r in rows if (r['hkl'],r['interval']) in keys]
        first=next(r for r in chosen if (r['hkl'],r['interval'])==(self.rows[indices[0]]['hkl'],self.rows[indices[0]]['interval']))
        if any(r['mode']!=first['mode'] or r['interval']!=first['interval'] for r in chosen):
            raise ValueError('Shared targets must have the same refinement mode and physical slice interval')
        different=any(not np.allclose(coordinates(r['initial'],r['mode']),coordinates(first['initial'],first['mode']),rtol=1e-12,atol=1e-12) for r in chosen)
        if different and not use_first:raise ValueError('Different starting values: explicitly choose Use first selected; values will not be averaged')
        lookup={(r['hkl'],r['interval']):r for r in chosen}
        for j,record in enumerate(records):
            if (record.hkl,record.interval) not in keys:continue
            row=lookup[(record.hkl,record.interval)]
            if different:
                z=complex_value(coordinates(first['initial'],first['mode']),record.mode,row['iam'] if record.mode=='magnitude' else row['initial'])
                record=replace(record,initial=(z.real,z.imag),reference=row['reference'])
            records[j]=replace(record,equality_group=name)
        candidate=replace(c,crystallography=replace(c.crystallography,projected_elastic=tuple(records)))
        compile_config_from_sources(candidate)  # validates intersected bounds and common coordinates
        if preview:return candidate
        self.commit_sharing(candidate,'Sharing changed optimizer dimension')

    def commit_sharing(self,candidate,title):
        w=self.window
        before=len(compile_config_from_sources(w.document.config).packed_starting_vector)
        after=len(compile_config_from_sources(candidate).packed_starting_vector)
        w.document._record(candidate)
        w.workflow.refinement_receipt=None;w.workflow.invalidate(2)
        self.dirty=False;w.invalidate_readiness()
        self.refresh_preserving_drafts()
        message=f'{title}: {before} → {after}. Validate Refinement Targets to accept this layout.'
        w.elastic_panel.message.setText(message)
        w.stage_status[2].setText(message);w.statusBar().showMessage(message)

    def dimension_counts(self,records=None):
        categories=self.window.document.config.refinement.active_categories
        ordinary=0;table=self.window.refinement_table
        for i in range(table.rowCount()):
            target=table.item(i,0).data(Qt.ItemDataRole.UserRole)
            if table.item(i,1).checkState()==Qt.CheckState.Checked and (categories is None or target.category in categories):ordinary+=1
        physical=before=after=ordinary;groups={}
        if categories is not None and 'projected_elastic' not in categories:return physical,before,after
        if records is None:
            items=[(r['record'].equality_group if r['record'] else None,self.table.cellWidget(i,4).currentData()) for i,r in enumerate(self.rows) if not r['reason'] and self.table.item(i,0).checkState()==Qt.CheckState.Checked]
        else:items=[(r.equality_group,r.mode) for r in records if r.refine]
        for i,(group,mode) in enumerate(items):
            n=len(components(mode));physical+=1;before+=n;groups[('group',group) if group else ('target',i)]=n
        return physical,before,after+sum(groups.values())

    def update_dimensions(self,*_):
        from .resolved_layout import resolved_layout
        try:
            compiled=compile_config_from_sources(self.window.document.config)
            counts,_=resolved_layout(compiled)
            text=f"Active physical targets: {counts['active_physical_targets']}\nOptimizer parameters: {len(compiled.packed_starting_vector)}"
            text=text.replace('\n','  |  ')
        except (ValueError,TypeError,OSError) as exc:text='Optimizer parameters unavailable — '+str(exc)
        draft=self.dirty or self.window.pending or self.window.workflow.refinement_receipt is None
        if draft:
            physical,_,dimension=self.dimension_counts()
            text=f'Draft optimizer layout (estimate)\nActive physical targets: {physical}  |  Optimizer parameters: {dimension}'
        else:text='Validated optimizer layout\n'+text
        groups={}
        for i,row in enumerate(self.rows):
            record=row['record']
            if record and record.equality_group:
                groups.setdefault(record.equality_group,[]).append((i,row))
        text+=f'  |  Shared parameter groups: {len(groups)}'
        lines=['Shared Refinement Parameters']
        for name,members in groups.items():
            lines.append(name)
            for i,row in members:
                mode=self.table.cellWidget(i,4).currentText()
                lines.append(f"  {row['hkl']}, slice {row['interval']} — {mode}; start {self.table.item(i,5).text()}; bounds [{self.table.item(i,6).text()}, {self.table.item(i,7).text()}]")
            modes={self.table.cellWidget(i,4).currentData() for i,_ in members}
            lines.append('  '+({'real':'Shared Re(E_g); own fixed Im','complex':'Shared Re(E_g) and Im(E_g)','magnitude':'Shared |E_g|; own IAM phase'}[next(iter(modes))] if len(modes)==1 else 'Incomplete: incompatible modes — choose a common mode or Make Independent'))
        if not groups:lines.append('Independent — no user-defined shared groups')
        self.live_groups.setText('\n'.join(lines))
        self.group_summary.setText('\n'.join(f"Shared {name}: "+', '.join(str(row['hkl']) for _,row in members)+' — '+', '.join(sorted({self.table.cellWidget(i,4).currentText() for i,_ in members})) for name,members in groups.items()) or 'Independent — no user-selected shared groups')
        self.dimension_summary.setText(text)

    def sharing_confirmation(self,indices,candidate):
        records=candidate.crystallography.projected_elastic
        keys={(self.rows[i]['hkl'],self.rows[i]['interval']) for i in indices}
        chosen=[r for r in records if (r.hkl,r.interval) in keys];mode=chosen[0].mode
        before=self.dimension_counts()[2];after=self.dimension_counts(records)[2]
        semantics={'real':'one real optimizer parameter (Shared real component). Each imaginary component remains fixed',
            'complex':'one real and one imaginary optimizer parameter (Shared complex coefficient)',
            'magnitude':'one magnitude optimizer parameter (Shared magnitude). Each IAM phase remains fixed'}[mode]
        return f'These {len(chosen)} targets will share {semantics}.\nOptimizer dimensionality: {before} → {after}.\nThe shared coordinate(s) remain equal throughout refinement; mandatory Friedel relations are unchanged.'

    def link_selected(self):
        indices=self.selected_rows()
        if len(indices)<2:self.message.setText('Select at least two coefficient rows');return
        name,ok=QInputDialog.getText(self,'Shared parameter group','Shared group name:',text='Group 1')
        if not ok:return
        states=[self.table.item(i,0).checkState() for i in indices]
        prior_dirty=self.dirty;prior_pending=self.window.pending;prior_message=self.message.text();prior_resolved=self.resolved.toPlainText()
        use_first=False
        try:
            try:candidate=self.link_rows(indices,name,preview=True)
            except ValueError as exc:
                if not str(exc).startswith('Different starting values'):raise
                answer=QMessageBox.question(self,'Choose shared starting value',str(exc)+'\nUse the first selected row’s starting component(s)?',QMessageBox.StandardButton.Yes|QMessageBox.StandardButton.Cancel)
                if answer!=QMessageBox.StandardButton.Yes:return
                use_first=True;candidate=self.link_rows(indices,name,use_first=True,preview=True)
            answer=QMessageBox.question(self,'Share refinement parameters',self.sharing_confirmation(indices,candidate),QMessageBox.StandardButton.Yes|QMessageBox.StandardButton.Cancel)
            if answer==QMessageBox.StandardButton.Yes:
                self.link_rows(indices,name,use_first=use_first)
                return
        except ValueError as exc:self.message.setText(str(exc))
        finally:
            if self.window.elastic_panel is self:
                self.table.blockSignals(True)
                for i,state in zip(indices,states):self.table.item(i,0).setCheckState(state)
                self.table.blockSignals(False)
                self.dirty=prior_dirty;self.window.pending=prior_pending;self.update_dimensions();self.resolved.setPlainText(prior_resolved)
                if self.message.text()=='Draft coefficient edits — Apply or validate to commit.':self.message.setText(prior_message)

    def unlink_selected(self):
        names={self.rows[i]['record'].equality_group for i in self.selected_rows() if self.rows[i]['record'] is not None}
        records=tuple(replace(r,equality_group=None) if r.equality_group in names else r for r in self.collect())
        c=self.window.document.config;candidate=replace(c,crystallography=replace(c.crystallography,projected_elastic=records))
        self.commit_sharing(candidate,'Optimizer dimension')

    def refresh_details(self):
        i=self.table.currentRow()
        if i<0:
            self.details.setEnabled(False)
            return
        row=self.rows[i];mode=self.table.cellWidget(i,4).currentData()
        self.details.setEnabled(not row['reason'])
        def describe(name,z):
            angle=f'{phase(z):.10g} rad' if z else 'undefined (zero magnitude)'
            return f'{name}: Re={z.real:.12g}; Im={z.imag:.12g}; |E_g|={abs(z):.12g}; phase={angle}'
        lines=[f"{row['hkl']}, slice {row['interval']}; mandatory Friedel partner {row['linked_hkl']}",describe('IAM',row['iam'])]
        try:
            text=self.table.item(i,5).text().strip()
            z=complex(text.replace(' ','').replace('i','j')) if mode=='complex' else complex_value((float(text),),mode,row['initial'] if mode=='real' else row['iam'])
            lines.append(describe('Starting value',z))
        except ValueError:lines.append('Starting value: incomplete or invalid')
        self.details_text.setText('\n'.join(lines))
        for component in range(2):
            visible=component==0 or mode=='complex'
            self.bound_labels[component].setText('Magnitude' if mode=='magnitude' else ('Real' if component==0 else 'Imaginary'))
            self.bound_labels[component].setVisible(visible)
            for col in (6,7):
                field=self.bound_fields[component,col];field.setVisible(visible)
                parts=self.table.item(i,col).text().split(',')
                field.setText(parts[component].strip() if component<len(parts) else '')

    def edit_bound(self,component,col,text):
        i=self.table.currentRow()
        if i<0:return
        count=len(components(self.table.cellWidget(i,4).currentData()))
        parts=self.table.item(i,col).text().split(',')
        parts=(parts+['']*count)[:count];parts[component]=text
        self.table.item(i,col).setText(', '.join(parts))

    def clear_targets(self):
        self.dirty=False  # Explicit Clear discards coefficient drafts as well as saved targets.
        cache=getattr(self.window,'_target_draft_cache',None)
        if cache:self.window._target_draft_cache=(cache[0],cache[1],{})
        c=self.window.document.config
        self.window.document._record(replace(c,crystallography=replace(c.crystallography,projected_elastic=())))
        self.window.workflow.refinement_receipt=None
        self.window.refresh()

    @staticmethod
    def start_text(value,mode):
        return display(value) if mode=='complex' else repr(float(coordinates(value,mode)[0]))

    def changed(self,*_):
        if hasattr(self,'resolved'):self.resolved.setPlainText('Draft changes — Validate Refinement Targets to resolve the layout.')
        self.window.workflow.refinement_receipt=None
        self.dirty=True;self.window.mark_pending(stage_override=2);self.message.setText('Draft coefficient edits — Apply or validate to commit.');self.refresh_details();self.update_dimensions()

    def mode_changed(self,i):
        row=self.rows[i];mode=self.table.cellWidget(i,4).currentData();old=row.get('_editor_mode',row['mode'])
        text=self.table.item(i,5).text()
        try:
            value=complex(text.replace(' ','').replace('i','j')) if old=='complex' else complex_value((float(text),),old,row['iam'] if old=='magnitude' else row['initial'])
        except ValueError:value=row['initial']
        self.table.blockSignals(True)
        self.table.item(i,5).setText(self.start_text(value,mode))
        for col in (6,7):
            parts=self.table.item(i,col).text().split(',')
            if mode=='complex':new=', '.join((parts+[''])[:2])
            elif mode=='real' and old=='complex':new=parts[0]
            elif mode==old:new=parts[0]
            else:new=''  # Magnitude and signed real bounds are not interchangeable.
            self.table.item(i,col).setText(new)
        self.table.blockSignals(False);row['_editor_mode']=mode;self.changed()

    def reset(self,i):
        row=self.rows[i];self.resets.add(i)
        self.table.item(i,5).setText(self.start_text(row['iam'],self.table.cellWidget(i,4).currentData()))
        self.table.item(i,8).setText('Reset to current IAM (apply to commit)');self.changed()

    def collect(self):
        records=[]
        for i,row in enumerate(self.rows):
            old=row['record'];mode=self.table.cellWidget(i,4).currentData()
            selected=self.table.item(i,0).checkState()==Qt.CheckState.Checked
            if row['reason']:
                if old:records.append(old)
                continue
            text=self.table.item(i,5).text().strip()
            reset=i in self.resets
            edited=text!=self.start_text(row['initial'],row['mode']) or mode!=row['mode']
            def limits(col):
                parts=self.table.item(i,col).text().strip().split(',')
                if len(parts)==1 and not parts[0]:return (None,None)
                if len(parts)!=len(components(mode)):raise ValueError('Bounds require one value, or real, imaginary for Complex mode')
                v=tuple(float(p) for p in parts)
                return v if len(v)==2 else (v[0],None)
            if reset: initial=None;reference=None
            elif edited:
                value=complex(text.replace(' ','').replace('i','j')) if mode=='complex' else complex_value((float(text),),mode,row['initial'] if mode=='real' else row['iam'])
                initial=(value.real,value.imag);reference=row['reference']
            else:initial=old.initial if old else None;reference=old.reference if old else None
            lo,hi=limits(6),limits(7)
            if selected or old or edited:
                records.append(ProjectedElasticConfig(row['hkl'],row['interval'],mode,selected,initial,reference,lo,hi,old.step_size if old else .01,old.equality_group if old else None))
        return tuple(records)

    def apply(self,checked=False,*,refresh=True):
        if not self.dirty:return True
        try:
            c=self.window.document.config
            candidate=replace(c,crystallography=replace(c.crystallography,projected_elastic=self.collect()))
            compile_config_from_sources(candidate)
            self.window.document._record(candidate);self.window.workflow.refinement_receipt=None;self.dirty=False
            if refresh:
                drafts=self.window.field_drafts();self.window.refresh();self.window.restore_field_drafts(drafts)
            return True
        except (ValueError,TypeError,OSError) as exc:self.message.setText(str(exc));return False
