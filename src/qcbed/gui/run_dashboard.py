"""Workflow dashboard and the single output-directory editor."""
from pathlib import Path
from PySide6.QtWidgets import QWidget,QVBoxLayout,QHBoxLayout,QLabel,QPushButton,QFormLayout,QScrollArea


class RunDashboard(QWidget):
    def __init__(self, window):
        super().__init__()
        self.window=window
        outer=QVBoxLayout(self);outer.setContentsMargins(0,0,0,0)
        self.content_scroll=QScrollArea();self.content_scroll.setWidgetResizable(True)
        outer.addWidget(self.content_scroll)
        content=QWidget();self.content_scroll.setWidget(content)
        layout=QVBoxLayout(content)
        title=QLabel('Run & Results');layout.addWidget(title)
        layout.addWidget(QLabel('Execution & Output'))
        self.output_form=QFormLayout();layout.addLayout(self.output_form)
        self.output_editor=None
        self.output_validation='Not checked — Validate Output Directory'
        status_body=QWidget();form=QFormLayout(status_body)
        scroll=QScrollArea();scroll.setWidgetResizable(True);scroll.setWidget(status_body)
        scroll.setMinimumHeight(180);layout.addWidget(scroll,1)
        self.status={}
        for key in ('Output directory','Run','Resume','A. Data Prep','B. Atomic Model',
                    'C Simulation Setup','Refinement Targets','Optimizer','External inputs',
                    'Active physical targets','Optimizer parameters','Compile validity','Execution readiness'):
            self.status[key]=QLabel('Invalid — validate configuration');self.status[key].setWordWrap(True)
            form.addRow(key,self.status[key])
        from .workflow_guidance import button as nav_button
        for title,stage,section in [('Go to Data Prep',0,None),('Go to Atomic Model',1,'Structure'),('Go to Simulation Setup',2,'Simulation Setup'),('Go to Refinement Targets',2,'Refinement Targets'),('Go to Optimizer',2,'Optimizer')]:
            layout.addWidget(nav_button(window,title,stage,section))
        controls=QHBoxLayout();layout.addLayout(controls)
        self.run=QPushButton('Run QCBED');self.resume=QPushButton('Resume QCBED')
        self.cancel=QPushButton('Cancel');self.results=QPushButton('Open Results')
        for button in (self.run,self.resume,self.cancel,self.results):
            button.setMinimumHeight(40);controls.addWidget(button)
        self.run.clicked.connect(lambda:window.start_run(False))
        self.resume.clicked.connect(lambda:window.start_run(True))
        self.cancel.clicked.connect(window.cancel_run)
        self.results.clicked.connect(window.show_results)
        self.reason=QLabel();self.reason.setWordWrap(True);layout.addWidget(self.reason)
        self.run_state=QLabel();self.run_state.setWordWrap(True);layout.addWidget(self.run_state)
        links=QHBoxLayout();layout.addLayout(links)
        for title,index in (('Go to Simulation Setup',0),('Go to Refinement Targets',1),('Go to Optimizer',2)):
            button=QPushButton(title);button.clicked.connect(lambda checked=False,i=index:self.navigate(i));links.addWidget(button)
        check=self.readiness_button=QPushButton('Validate Execution / Run Readiness')
        check.clicked.connect(window.validate_document);layout.addWidget(check)
        self.execution_settings=QVBoxLayout();layout.addLayout(self.execution_settings)
        layout.addStretch()
        self.external=[]

    def clear_output_editor(self):
        self.window.d_sections={}
        while self.output_form.rowCount():self.output_form.removeRow(0)
        while self.execution_settings.count():
            item=self.execution_settings.takeAt(0);item.widget().hide();item.widget().deleteLater()
        self.output_editor=None
        self.output_validation='Not checked — Validate Output Directory'

    def set_output_editor(self,container):
        self.output_editor=container
        self.output_form.addRow('Output directory',container)
        button=QPushButton('Validate Output Directory')
        button.clicked.connect(self.window.validate_output_directory)
        container.layout().addWidget(button)

    def navigate(self,index):
        self.window.workflow_tabs.setCurrentIndex(2)
        if hasattr(self.window,'c_tabs'):self.window.c_tabs.setCurrentIndex(index)

    def update_state(self, *, check_external=False):
        w=self.window
        if hasattr(w,'guidance'):w.guidance.update_state();doc=w.document;flow=w.workflow
        if check_external and doc:
            from .execution_details import external_checks
            self.external=external_checks(doc)
        def status(key,valid,reason):
            self.status[key].setText('✓ Valid' if valid else '○ Blocked — '+reason)
            self.status[key].setStyleSheet('padding: 6px; border-left: 3px solid '+('#39805b' if valid else '#a77838')+';')
        status('A. Data Prep',flow.receipts[0] is not None,'Data Prep must be validated')
        status('B. Atomic Model',flow.receipts[1] is not None,'Atomic Model must be validated')
        modern=doc is not None and doc.config.crystallography is not None
        status('C Simulation Setup',flow.simulation_receipt is not None if modern else flow.ready,'Validate Simulation Setup')
        status('Refinement Targets',flow.refinement_receipt is not None if modern else flow.ready,'Validate Refinement Targets')
        status('Optimizer',flow.receipts[2] is not None,'Validate C')
        inputs=[s for s in self.external if s.startswith(('Prepared experimental','Preparation metadata','Structure/CIF'))]
        failures=[s for s in inputs if ': BLOCKED' in s]
        self.status['External inputs'].setText('\n'.join(failures) if failures else 'Valid — prepared data, metadata and structure checked' if inputs else 'Not checked — load a project')
        self.status['External inputs'].setToolTip('\n'.join(inputs))
        editor=w.widgets.get(('run','output_directory'))
        selected=editor.text().strip() if editor is not None else (doc.config.run.output_directory.strip() if doc else '')
        if not selected:self.output_validation='Not selected'
        elif w.output_is_validated():self.output_validation='Valid'
        elif self.output_validation in ('Valid','Not selected'):
            self.output_validation='Not checked — Validate Output Directory'
        self.status['Output directory'].setText(self.output_validation)
        if self.output_editor:self.output_editor.setEnabled(not w.controller.active)
        self.readiness_button.setEnabled(bool(doc) and not w.controller.active)
        for section in getattr(w,'d_sections',{}).values():section.setEnabled(not w.controller.active)
        summary=doc.summary() if doc else {}
        for key in ('Active physical targets','Optimizer parameters'):
            value=summary.get(key.lower().replace(' ','_'))
            self.status[key].setText(str(value)+' committed' if value is not None else 'Unavailable — validate configuration')
        count=summary.get('optimizer_parameters',0)
        errors=[i.message for i in doc.issues if i.severity=='error'] if doc else []
        status('Compile validity',bool(doc and flow.ready and not errors),'Validate current committed configuration')
        self.run.setEnabled(w.run_action.isEnabled() and not w.controller.active)
        self.resume.setEnabled(w.resume_action.isEnabled() and not w.controller.active)
        self.cancel.setEnabled(w.cancel_action.isEnabled())
        self.results.setEnabled(w.results_directory is not None or w.results_view.results is not None)
        reasons=[]
        if not doc:reasons.append('Load a project or complete Data Prep')
        elif not w.controller.active:
            if modern and not w.output_is_validated():
                reasons.append('Validate Output Directory' if selected else 'Output directory not selected')
            if w.pending or getattr(w,'site_group_pending',False):reasons.append('Apply pending edits and revalidate')
            if flow.receipts[0] is None:reasons.append('Data Prep must be validated')
            if flow.receipts[1] is None:reasons.append('Atomic Model must be validated')
            if modern and flow.simulation_receipt is None:reasons.append('Validate Simulation Setup')
            if not count:reasons.append('No refinement targets selected')
            if modern and flow.refinement_receipt is None:reasons.append('Validate Refinement Targets')
            if flow.receipts[2] is None:reasons.append('Validate C')
            reasons.extend(errors)
            table=getattr(w,'refinement_table',None)
            if modern and table is not None:
                from PySide6.QtCore import Qt
                for row in range(table.rowCount()):
                    if table.item(row,1).checkState()==Qt.CheckState.Checked:
                        state=table.item(row,2).text()
                        if state.startswith(('Invalid','Incomplete')):
                            reasons.append(table.item(row,0).text()+': '+state)
            for action,label in ((w.run_action,'Run'),):
                if not action.isEnabled() and action.toolTip() and action.toolTip()!='Validate the current configuration first':
                    reasons.append(label+': '+action.toolTip())
        self.reason.setText('\n'.join(dict.fromkeys(reasons)) or ('Run active — settings locked' if w.controller.active else 'Ready.'))
        self.status['Execution readiness'].setText('Running' if w.controller.active else 'Valid — Run available' if self.run.isEnabled() else 'Blocked — see reasons below')
        self.status['Run'].setText('Ready' if self.run.isEnabled() else 'Running' if w.controller.active else 'Not ready — see reasons below')
        resume_text=w.resume_action.toolTip()
        if doc and modern and selected and not w.output_draft:
            from .output_paths import modern_execution_config
            try:
                previous=modern_execution_config(doc.config,doc.base_directory,resume=True)
                if not (Path(previous.run.output_directory)/'checkpoint_latest.json').is_file():
                    resume_text='Not available — no previous run found'
            except (OSError,ValueError):pass
        self.status['Resume'].setText('Ready' if self.resume.isEnabled() else resume_text.removeprefix('Resume: '))
        text=w.run_status.text()
        result=w.results_view.results
        if result and 'stop_reason' in result.summary:
            summary=result.summary
            text+=f"\nIterations: {summary['iterations'] if summary['iterations'] is not None else 'Not recorded'}\nObjective evaluations: {summary['objective_evaluations']}\nStop reason: {summary['stop_reason']}"
        self.run_state.setText(text)
