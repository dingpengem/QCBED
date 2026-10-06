"""Preserve editor drafts by physical identity when widgets are rebuilt."""
from PySide6.QtCore import Qt


def capture(w):
    panel=getattr(w,'elastic_panel',None)
    if panel is None or panel.document_identity is not w.document:return None
    cached=getattr(w,'_target_draft_cache',None)
    ordinary=dict(cached[1]) if cached and cached[0] is w.document else {}
    if w.pending:
        table=w.refinement_table
        for i in range(table.rowCount()):
            ordinary[table.item(i,0).data(Qt.ItemDataRole.UserRole)]=[(j,table.item(i,j).text(),table.item(i,j).checkState()) for j in (1,4,5,6,8,9)]
    elastic=dict(cached[2]) if cached and cached[0] is w.document else {}
    if panel.dirty:
        for i,r in enumerate(panel.rows):
            elastic[(r['hkl'],r['interval'])]=(panel.table.item(i,0).checkState(),panel.table.cellWidget(i,4).currentData(),[panel.table.item(i,j).text() for j in (5,6,7)])
    return ordinary,elastic


def restore(w,snapshot):
    if not snapshot:return
    ordinary,elastic=snapshot;table=w.refinement_table;table.blockSignals(True)
    restored_ordinary=set();ordinary_changed=False
    def same(a,b):
        if a==b:return True
        try:return float(a)==float(b)
        except ValueError:return False
    for i in range(table.rowCount()):
        identity=table.item(i,0).data(Qt.ItemDataRole.UserRole)
        if identity in ordinary:restored_ordinary.add(identity)
        for j,text,check in ordinary.get(identity,[]):
            ordinary_changed |= table.item(i,j).checkState()!=check if j==1 else not same(table.item(i,j).text(),text)
            table.item(i,j).setText(text)
            if j==1:table.item(i,j).setCheckState(check)
    table.blockSignals(False)
    p=w.elastic_panel;p.table.blockSignals(True);restored=set()
    for i,r in enumerate(p.rows):
        identity=(r['hkl'],r['interval'])
        if identity not in elastic:continue
        check,mode,texts=elastic[identity];combo=p.table.cellWidget(i,4);combo.blockSignals(True);combo.setCurrentIndex(combo.findData(mode));combo.blockSignals(False)
        r['_editor_mode']=mode
        p.table.item(i,0).setCheckState(check)
        for j,text in zip((5,6,7),texts):p.table.item(i,j).setText(text)
        restored.add(identity)
    p.table.blockSignals(False)
    if elastic:
        p.dirty=True;w.pending=True
        p.message.setText('Preserved target drafts — validate against the updated model.' if len(restored)==len(elastic) else 'Target drafts retained while Simulation Setup is unvalidated; review any unresolved hkl/slice identities after regeneration.')
    w._target_draft_cache=(w.document,{k:v for k,v in ordinary.items() if k not in restored_ordinary},{k:v for k,v in elastic.items() if k not in restored})
    if ordinary_changed or w._target_draft_cache[1]:w.pending=True
    p.update_dimensions()
