"""Navigation-only guidance; validation receipts remain the authority."""
from PySide6.QtWidgets import QWidget,QVBoxLayout,QHBoxLayout,QPushButton,QLabel


def go(w,stage,section=None):
    w.workflow_tabs.setCurrentIndex(stage)
    if section and section in w.section_items:w.navigation.setCurrentItem(w.section_items[section])


def button(w,text,stage,section=None):
    b=QPushButton(text);b.clicked.connect(lambda:go(w,stage,section));return b


class WorkflowGuidance(QWidget):
    def __init__(self,w):
        super().__init__();self.w=w;self.seen=set();layout=QVBoxLayout(self);row=QHBoxLayout();layout.addLayout(row)
        self.cards=[]
        for i,title in enumerate(('A. Data Prep','B. Atomic Model','C. Main Parameters','D. Run & Results')):
            b=button(w,title,i);row.addWidget(b);self.cards.append(b)
        self.next=QLabel('Next step: Define disk geometry and run Data Preparation.');self.next.setWordWrap(True);layout.addWidget(self.next)
    def update_state(self):
        w=self.w;f=w.workflow;valid=[r is not None for r in f.receipts]+[w.run_action.isEnabled()]
        self.seen.update(i for i,v in enumerate(valid) if v)
        for i,(b,title,v) in enumerate(zip(self.cards,('A. Data Prep','B. Atomic Model','C. Main Parameters','D. Run & Results'),valid)):
            state='Valid' if v else 'Not started' if w.document is None else 'Needs revalidation' if i in self.seen else 'Blocked' if i and not all(valid[:i]) else 'Ready to validate'
            if (i==0 and getattr(getattr(w,'preparation',None),'active',False)) or (i==3 and w.controller.active):state='In progress'
            color='#286247' if v else '#8a5d19' if state in ('Needs revalidation','Ready to validate') else '#666666'
            b.setText(('✓ ' if v else '○ ')+title+' — '+state)
            b.setStyleSheet(f'QPushButton {{border: {2 if w.workflow_tabs.currentIndex()==i else 1}px solid {color}; color: {color}; padding: 7px;}}')
        if w.document is None:next_step='Define disk geometry and run Data Preparation.'
        elif not valid[0]:next_step='Inspect the prepared data and Validate Data Prep.'
        elif not valid[1]:next_step='Confirm structure, zone axis and g1/g2; configure and Apply Site Groups, then Validate Atomic Model.'
        elif f.simulation_receipt is None:next_step='Configure and Validate Simulation Setup.'
        elif f.refinement_receipt is None:next_step='Select refinement targets and Validate Refinement Targets.'
        elif not valid[2]:next_step='Review optimizer settings and Validate Optimizer.'
        elif not w.output_is_validated():next_step='Select and Validate Output Directory in D.'
        else:next_step='Run is ready.' if valid[3] else 'Check D readiness for remaining runtime inputs.'
        self.next.setText('Next step: '+next_step)
        for key,label in getattr(w,'c_guidance',{}).items():
            ok={'simulation':f.simulation_receipt is not None,'refinement':f.refinement_receipt is not None,'optimizer':valid[2]}[key]
            label.setText(('✓ Valid' if ok else '○ Needs validation')+' — '+{'simulation':'C1 Simulation Setup','refinement':'C2 Refinement Targets','optimizer':'C3 Optimizer'}[key])
        for key,b in getattr(w,'c_next',{}).items():b.setEnabled({'simulation':f.simulation_receipt is not None,'refinement':f.refinement_receipt is not None,'optimizer':valid[2]}[key])
