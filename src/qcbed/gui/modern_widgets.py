"""Read-only structure/provenance views using the existing B1/B2 plotting API."""
from PySide6.QtWidgets import QLabel,QTableWidget,QTableWidgetItem,QPushButton
from PySide6.QtCore import Qt
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from qcbed.crystallography_config import compile_config_from_sources
from qcbed.modern_frontend import CAPABILITIES


def populate_modern(window,forms):
    config=window.document.config
    window.compiled_simulation_for_display=None
    help=QLabel('Use original/conventional-cell Miller indices for the reciprocal basis.\n'
        'Select supported physical refinement targets in C. Membership and occupancy are fixed inputs; explicit absorption coefficients may be refined.')
    help.setWordWrap(True);forms['Structure'].addRow(help)
    try:
        from qcbed.crystallography_config import load_config_structure
        from qcbed.crystallography import CrystalOrientation, ZoneAxis
        structure = load_config_structure(config)
        lattice = structure.direct_lattice
        orientation = CrystalOrientation(lattice, ZoneAxis(config.crystallography.zone_axis))
        diagnostic = dict(lattice_lengths_angstrom=lattice.lengths.tolist(),
            lattice_angles_degrees=lattice.angles.tolist(), atom_count=len(structure.sites),
            species=list(dict.fromkeys(site.symbol for site in structure.sites)))
        try:
            from .document import parse_text
            parse_text(', '.join(map(str,config.crystallography.slice_boundaries_fractional)),'boundaries')
            if window.document.thermal_incomplete(config):
                raise ValueError(window.document.thermal_message(config))
            bundle = compile_config_from_sources(config).crystallography
            window.compiled_simulation_for_display=bundle.simulation_sample
            diagnostic.update(bundle.diagnostics)
            diagnostic.update(requested_depth=config.crystallography.depth_angstrom,
                effective_depth=bundle.simulation_sample.effective_depth_angstrom)
        except (ValueError, OSError, TypeError) as exc:
            diagnostic['Model validation'] = str(exc)
        from .sections import ExpandableSection
        window.atomic_expert=ExpandableSection('Expert / Legacy settings — model diagnostics')
        details=QLabel('\n'.join(f'{k}: {v}' for k,v in diagnostic.items()));details.setWordWrap(True)
        window.atomic_expert.content.addWidget(details);forms['Structure'].addRow(window.atomic_expert)
        basic={k:v for k,v in diagnostic.items() if k in ('lattice_lengths_angstrom','lattice_angles_degrees','atom_count','species')}
        window.modern_summary=QLabel('\n'.join(f'{k}: {v}' for k,v in basic.items()))
        window.modern_summary.setWordWrap(True);forms['Structure'].addRow(window.modern_summary)
        from qcbed.crystallography.projection import project_atomic_structure
        import numpy as np
        from itertools import product
        from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
        # All views share the validated beam-frame transform; replication is display-only.
        sites=project_atomic_structure(structure,orientation,repeats=(2,2,2))
        coordinates=np.array([(v.x,v.y,v.depth) for v in sites])
        corners=np.array(list(product((0.,2.),repeat=3)))
        boundary=orientation.cartesian_to_beam(corners @ lattice.vectors)
        edges=[(i,j) for i in range(8) for j in range(i+1,8)
               if np.count_nonzero(corners[i]!=corners[j])==1]
        labels=('x_beam (Å)','y_beam (Å)','z_beam (Å)')
        symbols=list(dict.fromkeys(v.symbol for v in sites))
        figure=Figure(figsize=(12,4.5),layout='constrained')
        window.atomic_axes=figure.subplots(1,3)
        for ax,(first,second,title) in zip(window.atomic_axes,(
                (0,1,'Front / Beam — along +z_beam'),
                (0,2,'Top — along +y_beam'),
                (2,1,'Left — along +x_beam'))):
            for symbol in symbols:
                rows=np.array([v.symbol==symbol for v in sites])
                ax.scatter(coordinates[rows,first],coordinates[rows,second],label=symbol)
            for i,j in edges:
                ax.plot(boundary[[i,j],first],boundary[[i,j],second],color='gray',linewidth=.8)
            ax.set(xlabel=labels[first],ylabel=labels[second],title=title)
            ax.set_aspect('equal');ax.legend()
        canvas=FigureCanvasQTAgg(figure)
        canvas.setMinimumSize(800,400)
        forms['Structure'].addRow(canvas)
        figure3=Figure(figsize=(10,6),layout='constrained')
        ax=window.atomic_3d_axes=figure3.add_subplot(111,projection='3d')
        for symbol in symbols:
            rows=np.array([v.symbol==symbol for v in sites])
            ax.scatter(*coordinates[rows].T,label=symbol)
        for i,j in edges:
            ax.plot(*boundary[[i,j]].T,color='gray',linewidth=.8)
        extent=np.ptp(boundary,axis=0)
        ax.set_box_aspect(extent)
        length=float(max(extent))*.6
        ax.quiver(0,0,0,0,0,length,color='black',arrow_length_ratio=.1)
        ax.text(0,0,length,'Zone axis '+str(config.crystallography.zone_axis))
        ax.set(xlabel=labels[0],ylabel=labels[1],zlabel=labels[2],
               title='3D Atomic Model — 2 × 2 × 2 cells (visualization only)')
        ax.legend()
        canvas3=FigureCanvasQTAgg(figure3);canvas3.setMinimumSize(650,500)
        forms['Structure'].addRow(NavigationToolbar2QT(canvas3,window))
        forms['Structure'].addRow(canvas3)
    except (ValueError,OSError,TypeError) as exc:
        forms['Structure'].addRow(QLabel(str(exc)))
    groups=config.crystallography.site_groups
    table=window.site_group_table=QTableWidget(len(groups),8)
    table.setHorizontalHeaderLabels(['Stable group','Element','Original sites','Debye–Waller B (Å²)','Occupancy','Absorption (a,b,c; blank=disabled)','Refine Debye–Waller B?','Resolved Debye–Waller bounds'])
    from qcbed.crystallography.site_groups import group_target
    from qcbed.crystallography_refinement import physical_target_layout
    try:layout={} if window.document.thermal_incomplete(config) else {p['target_id']:p for p in physical_target_layout(compile_config_from_sources(config))}
    except (ValueError,OSError):layout={}
    for i,g in enumerate(groups):
        bound=layout.get(f'site_group.{g.id}.debye_waller_b',{})
        specific=next((b for b in config.refinement.bounds.specific if b.target==group_target(g.id)),None)
        limits=(bound.get('lower'),bound.get('upper')) if bound else ('Specific (inactive)',specific.lower,specific.upper) if specific else 'Inactive / resolve on selection'
        for j,value in enumerate((g.id,g.element,g.site_indices,g.debye_waller_b,g.occupancy,
                                 'Disabled' if g.absorption is None else g.absorption,
                                 'Yes' if bound else 'No', limits)):
            if j==2:value=', '.join(map(str,g.site_indices))
            if j==3 and value is None:value=''
            if j==5:value='' if g.absorption is None else ', '.join(map(str,g.absorption))
            item=QTableWidgetItem(str(value))
            if j>=6:item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            table.setItem(i,j,item)
    window.site_group_pending=False
    def draft_changed(*args):
        window.site_group_pending=True
        window.workflow.invalidate(1);window.invalidate_readiness()
        window.statusBar().showMessage('Site-group drafts — select Apply Site-Group Parameters.')
    def apply_groups():
        from qcbed.configuration import AtomicSiteGroup
        from .document import EditorDocument
        candidate=EditorDocument(window.document.config,window.document.path)
        drafts=window.field_drafts()
        try:
            for row in range(table.rowCount()):
                text=[table.item(row,j).text().strip() for j in range(6)]
                if not text[3]:raise ValueError(f'Debye–Waller B is required for site group {text[0]}.')
                group=AtomicSiteGroup(text[0],text[1],tuple(int(v.strip()) for v in text[2].split(',')),
                    float(text[3]),float(text[4]),
                    tuple(float(v.strip()) for v in text[5].split(',')) if text[5] else None)
                candidate.edit_site_group(row,group)
            # Commit B's independent orientation/basis drafts with the group edit.
            if not candidate.apply(drafts,stage=1):
                raise ValueError('; '.join(i.message for i in candidate.issues if i.severity=='error'))
        except (ValueError,TypeError,OSError) as exc:
            window.statusBar().showMessage(str(exc))
            window.stage_status[1].setText(str(exc))
            return
        window.document._record(candidate.config)
        window.document.validate()
        window.site_group_pending=False
        window.refresh()
        window.restore_field_drafts(drafts)
    table.cellChanged.connect(draft_changed)
    window.apply_site_groups=QPushButton('Apply Site-Group Parameters')
    window.apply_site_groups.clicked.connect(apply_groups)
    forms['Site groups'].addRow(window.apply_site_groups)
    table.resizeColumnsToContents();forms['Site groups'].addRow(table)
    from .sections import ExpandableSection
    from PySide6.QtWidgets import QCheckBox
    window.group_advanced=ExpandableSection('Advanced settings / Site-group membership')
    membership=QCheckBox('Show and edit member-site indices')
    for column in (2,6,7):table.setColumnHidden(column,True)
    membership.toggled.connect(lambda on:table.setColumnHidden(2,not on))
    window.group_advanced.content.addWidget(membership);forms['Site groups'].addRow(window.group_advanced)
    note=QLabel('Initial groups share element, source occupancy and source Debye–Waller B; they are not a symmetry-equivalence determination. Edit and review before validation.\n'
        'Use C → Refinement Targets to choose physical quantities and hard bounds.\n'
        'ANNEAL, AMOEBA, PARAMIN and SIMPARA support the dynamic active layout; optimizer-specific restrictions are checked before Run.')
    note.setWordWrap(True);forms['Site groups'].addRow(note)
    audit=QTableWidget(len(CAPABILITIES),7)
    audit.setHorizontalHeaderLabels(['Family','Capability','Scope','Units','Initial source','Bounds source','Binding / limitation'])
    for i,row in enumerate(CAPABILITIES):
        for j,value in enumerate(row):
            item=QTableWidgetItem(value);item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable);audit.setItem(i,j,item)
    audit.resizeColumnsToContents();audit.setMinimumHeight(330)
    capabilities=ExpandableSection('Expert / Legacy settings — capability audit')
    capabilities.content.addWidget(audit);forms['Capabilities'].addRow(capabilities)
    label=QLabel('Signal mode selects a different prepared experimental signal/objective. Use Data Prep to prepare or reopen matching data; switching mode does not transform existing data.\n'
        'Conventional and gradient objective values are not directly comparable. Gradient support is binary inclusion; variance is the statistical denominator.')
    label.setWordWrap(True);forms['Microscope'].addRow(label)
