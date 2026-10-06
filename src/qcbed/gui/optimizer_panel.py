"""Presentation of existing search metadata, independent of physical target selection."""
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QCheckBox,QTableWidget,QTableWidgetItem,QLabel
from qcbed.optimizer_defaults import recommended_binding,binding_settings

TIPS=(
 'Dimensionless range factor D1; PARAMIN TOP. AMOEBA uses it in scaled coordinates. Not a physical hard bound. Signed factors are preserved; PARAMIN requires nonzero factors.',
 'Dimensionless range factor D2; PARAMIN BOT. AMOEBA requires D1 ≠ D2 and a nonzero starting value. Signed factors are preserved. Not a physical hard bound.',
 'AMOEBA: dimensionless simplex displacement; physical step = PLAMDA × (D2−D1) × initial value; signed/zero values retain backend semantics. SIMPARA: nonnegative absolute uncertainty width in target units. Not a hard bound. Unused by ANNEAL/PARAMIN.')


def fill_defaults(window,row,*,replace=False):
    table=window.refinement_table
    target=table.item(row,0).data(Qt.ItemDataRole.UserRole)
    try:values=binding_settings(recommended_binding(target))
    except ValueError:return
    table.blockSignals(True)
    for col,value in zip((4,5,6),values):
        if replace or not table.item(row,col).text().strip():table.item(row,col).setText(repr(value))
    table.blockSignals(False)


def build(window,layout,form):
    window.recommended_optimizer=QCheckBox('Use recommended optimizer settings')
    defaults=True
    for b in window.document.config.refinement.parameters:
        try:defaults &= b==recommended_binding(b.target)
        except ValueError:defaults=False
    window.recommended_optimizer.setChecked(defaults)
    layout.addWidget(window.recommended_optimizer)
    note=QLabel('Recommended per-parameter search metadata comes from the documented reference profiles. '
        'High-level workflow settings remain explicit below. Loaded overrides are preserved; selecting recommended settings replaces them.')
    note.setWordWrap(True);window.recommended_optimizer.setToolTip(note.text());note.hide()
    from .sections import ExpandableSection
    expert=ExpandableSection('Expert / Legacy settings — per-target search metadata')
    window.advanced_optimizer_toggle=expert.toggle
    window.advanced_optimizer=expert.body
    window.optimizer_metadata=QTableWidget(0,4)
    window.optimizer_metadata.setHorizontalHeaderLabels(['Physical parameter','Search-range factor (D1)','Search-range factor (D2)','Perturbation parameter (PLAMDA)'])
    window.optimizer_metadata.setMinimumHeight(250);expert.content.addWidget(window.optimizer_metadata)
    layout.addWidget(expert)
    window.optimizer_defaults_status=QLabel('Advanced optimizer settings: using recommended defaults' if defaults else
        'Advanced optimizer settings: explicit overrides loaded')
    layout.addWidget(window.optimizer_defaults_status)
    window.optimizer_metadata.itemChanged.connect(lambda item:edit(window,item))
    window.recommended_optimizer.toggled.connect(lambda on:recommend(window,on))
    window.widgets[('refinement','optimizer','method')].currentTextChanged.connect(lambda _:refresh(window))
    refresh(window)


def refresh(window):
    if not hasattr(window,'optimizer_metadata'):return
    table=window.refinement_table;advanced=window.optimizer_metadata
    method=window.widgets[('refinement','optimizer','method')].currentText()
    amoeba='AMOEBA' in method;anneal='ANNEAL' in method;para=method=='PARAMIN';spar=method=='SIMPARA'
    elastic=bool(window.document.config.crystallography and any(r.refine for r in window.document.config.crystallography.projected_elastic))
    used=(False,False,amoeba or para or spar) if elastic else (amoeba or para,amoeba or para,amoeba or spar)
    advanced.blockSignals(True)
    rows=[r for r in range(table.rowCount()) if table.item(r,1).checkState()==Qt.CheckState.Checked]
    advanced.setRowCount(len(rows))
    for dest,source in enumerate(rows):
        target=table.item(source,0).data(Qt.ItemDataRole.UserRole)
        for col in range(4):
            text=table.item(source,0 if col==0 else col+3).text()
            if col and not used[col-1]:text='Not used by the current optimizer'
            item=QTableWidgetItem(text);item.setData(Qt.ItemDataRole.UserRole,target)
            if not col or not used[col-1] or window.recommended_optimizer.isChecked():item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            if col:item.setToolTip(('With E_g targets: physical search step for AMOEBA/PARAMIN; absolute uncertainty width for SIMPARA. D1/D2 are retained but unused.' if elastic else TIPS[col-1]))
            advanced.setItem(dest,col,item)
    advanced.resizeColumnsToContents();advanced.blockSignals(False)
    enabled={'maxtemp':anneal or spar,'mintemp':anneal or spar,'coolrate':anneal or spar,
        'ftol':amoeba or para,'chimin':amoeba,'maxfsimp':anneal or amoeba or spar,
        'maxfpara':para,'seed':anneal or spar or para}
    for name,on in enabled.items():
        widget=window.widgets.get(('refinement','optimizer',name))
        if widget:
            widget.setEnabled(on);widget.setToolTip(('AMOEBA: soft search-evaluation budget, excluding N+1 initialization vertices; checked with > at cycle boundaries, so total evaluations can exceed this number. ANNEAL: legacy iteration counter (up to budget+1). SIMPARA: soft evaluation budget. Use D Advanced maximum completed iterations for an exact additional iteration cap.' if name=='maxfsimp' else '') if on else 'Not used by the current optimizer')


def recommend(window,on):
    if on:
        for row in range(window.refinement_table.rowCount()):
            if window.refinement_table.item(row,1).checkState()==Qt.CheckState.Checked:
                fill_defaults(window,row,replace=True)
    window.optimizer_defaults_status.setText('Advanced optimizer settings: using recommended defaults' if on else 'Advanced optimizer settings: explicit overrides')
    window.workflow.refinement_receipt=None
    window.workflow.invalidate(2);window.pending=True;window.invalidate_readiness()
    refresh(window)


def edit(window,item):
    if item.column()==0:return
    target=item.data(Qt.ItemDataRole.UserRole)
    table=window.refinement_table
    row=next(r for r in range(table.rowCount()) if table.item(r,0).data(Qt.ItemDataRole.UserRole)==target)
    table.item(row,item.column()+3).setText(item.text())


def resolve_unused(window):
    """Unused fields retain explicit values; empty inactive inputs need no user entry."""
    method=window.widgets[('refinement','optimizer','method')].currentText()
    used=('AMOEBA' in method or method=='PARAMIN',)*2+('AMOEBA' in method or method=='SIMPARA',)
    table=window.refinement_table
    table.blockSignals(True)
    try:
        for row in range(table.rowCount()):
            if table.item(row,1).checkState()!=Qt.CheckState.Checked:continue
            target=table.item(row,0).data(Qt.ItemDataRole.UserRole)
            try:values=binding_settings(recommended_binding(target))
            except ValueError:continue
            for col,on,value in zip((4,5,6),used,values):
                if not on and not table.item(row,col).text().strip():table.item(row,col).setText(repr(value))
    finally:table.blockSignals(False)
