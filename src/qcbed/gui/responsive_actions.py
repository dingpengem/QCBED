"""Layout-only action rows which stack when their labels no longer fit."""
from PySide6.QtCore import QSize
from PySide6.QtWidgets import QWidget,QBoxLayout,QLayout,QSizePolicy


class ResponsiveActions(QWidget):
    def __init__(self,widgets):
        super().__init__()
        self.items=list(widgets)
        self.row=QBoxLayout(QBoxLayout.Direction.LeftToRight,self)
        self.row.setContentsMargins(0,0,0,0)
        self.row.setSpacing(6)
        self.row.setSizeConstraint(QLayout.SizeConstraint.SetNoConstraint)
        for widget in self.items:self.row.addWidget(widget)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,QSizePolicy.Policy.Minimum)

    def minimumSizeHint(self):
        return QSize(max((w.minimumSizeHint().width() for w in self.items),default=0),self.row.minimumSize().height())

    def resizeEvent(self,event):
        super().resizeEvent(event)
        needed=sum(w.minimumSizeHint().width() for w in self.items)+6*max(0,len(self.items)-1)
        direction=QBoxLayout.Direction.TopToBottom if self.width()<needed else QBoxLayout.Direction.LeftToRight
        if self.row.direction()!=direction:
            self.row.setDirection(direction)
            self.updateGeometry()
