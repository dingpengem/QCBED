"""Consistent display-only expandable sections; never modify canonical values."""
from PySide6.QtWidgets import QWidget,QVBoxLayout,QPushButton,QFormLayout,QScrollArea
from PySide6.QtCore import QTimer


class ExpandableSection(QWidget):
    def __init__(self,title,parent=None,*,expanded=False):
        super().__init__(parent)
        self.title=title
        layout=QVBoxLayout(self);layout.setContentsMargins(0,0,0,0)
        self.toggle=QPushButton();self.toggle.setCheckable(True)
        layout.addWidget(self.toggle)
        self.body=QWidget();self.content=QVBoxLayout(self.body);self.content.setContentsMargins(8,0,0,0)
        layout.addWidget(self.body)
        self.toggle.toggled.connect(self._show)
        self.set_expanded(expanded)

    def _show(self,on):
        self.body.setVisible(on)
        self.toggle.setText(('▼ ' if on else '▶ ')+self.title)
        if on:QTimer.singleShot(0,self._reveal)

    def _reveal(self):
        try:
            parent=self.parentWidget()
            while parent is not None:
                if isinstance(parent,QScrollArea):
                    parent.ensureWidgetVisible(self.body,20,20);return
                parent=parent.parentWidget()
        except RuntimeError:
            pass # A refresh may have replaced the old immutable-document view.


    def set_expanded(self,on):
        self.toggle.setChecked(on);self._show(on)


def move_form_row(source,widget,destination):
    row=source.takeRow(widget)
    label=row.labelItem.widget() if row.labelItem else None
    field=row.fieldItem.widget()
    if label:destination.addRow(label,field)
    else:destination.addRow(field)


def collapse_form(form,title,*,expanded=False):
    section=ExpandableSection(title,expanded=expanded)
    body=QWidget();destination=QFormLayout(body)
    while form.rowCount():
        row=form.takeRow(0)
        label=row.labelItem.widget() if row.labelItem else None
        field=row.fieldItem.widget() if row.fieldItem else None
        if field:
            if label:destination.addRow(label,field)
            else:destination.addRow(field)
    section.content.addWidget(body);form.addRow(section)
    return section
