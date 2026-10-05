"""Read-only explanation of current widget drafts and committed readiness."""
from qcbed.config_validation import ValidationPhase


def external_checks(document):
    from .data_model import read_data,provenance,restored_overview
    c=document.config
    rows=[]
    def check(label,operation):
        try:
            operation()
            rows.append(label+': PASS')
        except (OSError,ValueError,TypeError,KeyError) as exc:
            rows.append(label+': BLOCKED — '+str(exc))
    check('Prepared experimental data '+c.experiment.path,lambda:read_data(c,document.base_directory))
    def assets():
        meta=provenance(c)
        if meta:restored_overview(c,document.base_directory)
    check('Preparation metadata/cache integrity (when recorded)',assets)
    if c.crystallography:
        from qcbed.crystallography_config import load_config_structure
        check('Structure/CIF '+c.crystallography.structure.path,lambda:load_config_structure(c))
    from .output_paths import validate_output_directory,check_replacement
    check('Output directory',lambda:validate_output_directory(c.run.output_directory,document.base_directory))
    if c.crystallography is None:
        check('Legacy output replacement',lambda:check_replacement(c,document.base_directory))
    for issue in document.issues:
        if issue.phase==ValidationPhase.FILES and issue.severity=='error':
            rows.append(f'{issue.section}.{issue.field}: BLOCKED — {issue.message}')
    rows.append('Raw CBED: not a refinement runtime input; required when recomputing preparation.')
    rows.append('Output filename: internally managed as final_out.dat; no user input required.')
    return rows


def describe_execution(document,table):
    from PySide6.QtCore import Qt
    from qcbed import available_refinement_targets
    c=document.config
    catalogue={v.target:v for v in available_refinement_targets(c)}
    bindings={b.target:b for b in c.refinement.parameters}
    lines=['Grouped Run accepts any nonempty subset of supported complete physical bindings.',
           'Each active target needs ordered physical hard bounds containing its starting value; target-specific domains still apply.','']
    for row in range(table.rowCount()):
        target=table.item(row,0).data(Qt.ItemDataRole.UserRole)
        info=catalogue[target];binding=bindings.get(target)
        selected=table.item(row,1).checkState()==Qt.CheckState.Checked
        if target.category.value not in ('depth','debye_waller') and not selected and binding is None:continue
        missing=[name for col,name in ((4,'D1'),(5,'D2'),(6,'PLAMDA'),(8,'lower bound'),(9,'upper bound'))
                 if col<table.columnCount() and not table.item(row,col).text().strip()]
        reason=('not selected' if not selected else 'incomplete — '+', '.join(missing)+' required'
                if missing else 'not committed — Apply edits required' if binding is None
                else table.item(row,2).text())
        lines.extend([table.item(row,0).text(),
            f'  Binding checkbox: {selected}; committed binding: {binding is not None}',
            f'  Visible starting value: {table.item(row,3).text()}; committed: {info.starting_value}',
            '  Visible D1 / D2 / PLAMDA: '+ ' / '.join(table.item(row,col).text() or '(missing)' for col in (4,5,6)),
            '  Committed D1 / D2 / PLAMDA: '+(str((*binding.range_factors,binding.plamda)) if binding else '(none)'),
            f'  Active/status: {reason}'])
    return '\n'.join(lines+['','External checks:']+external_checks(document))
