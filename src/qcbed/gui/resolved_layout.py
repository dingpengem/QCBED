"""Read-only presentation of the actual compiled optimizer layout."""
from qcbed.projected_elastic import layout_rows


def resolved_layout(compiled):
    ordinary=[p for p in compiled.layout.parameters if p.category!='projected_elastic']
    rows=layout_rows(compiled) if compiled.crystallography is not None else ()
    groups={};owners={};independent=[]
    for p in rows:
        meta=p['projected_elastic']
        for m in meta.get('members',[meta]):owners[(tuple(m['hkl']),tuple(m['interval']))]=m
        if meta.get('equality_group'):groups.setdefault(meta['equality_group'],[]).append(p)
        else:independent.append(p)
    independent_owners={(tuple(p['projected_elastic']['hkl']),tuple(p['projected_elastic']['interval'])) for p in independent}
    counts=dict(active_physical_targets=len(ordinary)+len(owners),shared_parameter_groups=len(groups),
        independent_targets=len(ordinary)+len(independent_owners),mandatory_friedel_relationships=len(owners),optimizer_parameters=compiled.parameter_count)
    def identity(m):return f"{tuple(m['hkl'])}, slice {tuple(m['interval'])}"
    lines=['Resolved Refinement Parameter Layout (committed compiled snapshot)']
    modes={'real':'Real only — shared Re(E_g); each imaginary component remains fixed',
           'complex':'Complex — shared Re(E_g) and Im(E_g)',
           'magnitude':'Magnitude only — shared |E_g|; each IAM phase remains fixed'}
    def coordinate(p):return f"  {p['projected_elastic']['component']}: start {p['initial_value']}; hard bounds [{p['lower']}, {p['upper']}]"
    for name,parameters in groups.items():
        meta=parameters[0]['projected_elastic']
        lines.extend(['',f'Shared Parameter Group: {name}','Members:'])
        lines.extend('  '+identity(m) for m in meta['members'])
        lines.extend(['Mode: '+modes[meta['mode']],f'Optimizer dimensions contributed: {len(parameters)}'])
        lines.extend(coordinate(p) for p in parameters)
    lines.extend(['','Independent targets:'])
    if compiled.crystallography is not None:
        from qcbed.crystallography_refinement import physical_target_layout
        lines.extend('  '+p['target_id'] for p in physical_target_layout(compiled) if p['category']!='projected_elastic')
    else:
        for p in ordinary:lines.append(f'  {p.category} {p.target}')
    for p in independent:lines.append('  '+identity(p['projected_elastic'])+' / '+p['projected_elastic']['component'])
    if not ordinary and not independent:lines.append('  None')
    lines.extend(['','Mandatory Friedel relationships (not extra coordinates):'])
    lines.extend('  '+identity(m)+f" ↔ {tuple(m['linked_hkl'])}" for m in owners.values())
    lines.append('')
    for key,value in counts.items():lines.append(key.replace('_',' ').capitalize()+': '+str(value))
    lines.append(f'Final optimizer-vector dimension: {compiled.parameter_count}')
    assert len(ordinary)+len(rows)==compiled.parameter_count
    return counts,'\n'.join(lines)
