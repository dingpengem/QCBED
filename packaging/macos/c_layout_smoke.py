"""Geometry assertions shared by source-widget and frozen acceptance checks."""
import json
from PySide6.QtCore import QPoint
from PySide6.QtWidgets import QPushButton,QScrollArea


def settle(app):
    for _ in range(8):app.processEvents()


def assert_controls(w,app):
    settle(app)
    page=w.c_tabs.currentWidget()
    from qcbed.gui.responsive_actions import ResponsiveActions
    footer=page.findChild(ResponsiveActions,'cBottomActions')
    assert footer is not None and footer.isVisible()
    controls=footer.findChildren(QPushButton)
    rects=[]
    for button in controls:
        assert button.isVisible(),button.text()
        assert button.width()>=button.minimumSizeHint().width(),button.text()
        assert button.height()>=button.minimumSizeHint().height(),button.text()
        assert footer.rect().contains(button.geometry()),(button.text(),button.geometry(),footer.rect())
        rect=button.rect();rect.moveTopLeft(button.mapTo(page,QPoint(0,0)))
        assert page.rect().contains(rect),(button.text(),rect,page.rect())
        assert all(not rect.intersects(other) for other in rects),button.text()
        rects.append(rect)
    scroll=page.findChild(QScrollArea)
    assert scroll.geometry().bottom()<footer.geometry().top()
    assert scroll.horizontalScrollBar().maximum()==0,'C content must fit the viewport; tables scroll internally'
    return scroll


def audit(w,app,out):
    from qcbed.gui.responsive_actions import ResponsiveActions
    dimensions=[]
    w.workflow_tabs.setCurrentIndex(2)
    for width,height in ((960,720),(1250,850),(1700,950)):
        w.resize(width,height);settle(app)
        for tab,name in enumerate(('simulation','targets','optimizer')):
            w.c_tabs.setCurrentIndex(tab);scroll=assert_controls(w,app)
            section=(w.numerical_advanced,w.elastic_panel.layout_details,w.optimizer_advanced_settings)[tab]
            assert not section.toggle.isChecked()
            if tab==1:
                p=w.elastic_panel
                assert p.table.isColumnHidden(3) and p.table.isColumnHidden(9)
                assert p.sharing_instructions.isVisible() and p.group_summary.isVisible()
                scroll.ensureWidgetVisible(p.sharing_instructions);settle(app)
                for actions in p.findChildren(ResponsiveActions):
                    scroll.ensureWidgetVisible(actions);settle(app)
                    for b in actions.items:
                        assert actions.rect().contains(b.geometry())
                        visible=b.rect();visible.moveTopLeft(b.mapTo(scroll.viewport(),QPoint(0,0)))
                        assert scroll.viewport().rect().contains(visible),b.text()
                    assert not actions.items[0].geometry().intersects(actions.items[1].geometry())
                dimensions.append(p.width())
            w.grab().save(str(out/f'{name}-{width}.png'))
            section.set_expanded(True);settle(app);assert_controls(w,app)
            assert scroll.verticalScrollBar().maximum()>=0
            section.set_expanded(False);settle(app);assert_controls(w,app)
            assert len(w.c_tabs.widget(tab).findChildren(ResponsiveActions,'cBottomActions'))==1
    assert dimensions[0]<dimensions[1]<dimensions[2],dimensions
    return dimensions


def run(w,app,out):
    widths=audit(w,app,out)
    report=json.loads((out/'result.json').read_text())
    report['checks']+=['C responsive controls: narrow/normal/wide; expanded/collapsed details; no action collision or clipping','Coefficient content expands with viewport; bottom validation remains accessible','Visible sharing instructions, immediate group summary and expandable mode-specific example']
    report['coefficient_panel_widths']=widths
    (out/'result.json').write_text(json.dumps(report,indent=2))
    return 0
