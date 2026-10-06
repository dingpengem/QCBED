"""Slice-resolved projected elastic coefficients; unchanged scattering/propagation.

E and A remain separate. Absolute E replacement is defined even at zero IAM E;
no division by an atomic reference or modification of absorptive coefficients.
"""
from dataclasses import asdict, replace
import hashlib
import json
import numpy as np

from .configuration import ProjectedElasticConfig
from .structure_factors import elastic_structure_factors
from .absorption import absorptive_structure_factors
from .phase_grating import phase_grating_from_structure_factors, reciprocal_cutoff_mask, legacy_nyquist_mask
from .crystallography.index_transform import original_hkl_to_internal_HK

CATEGORY = 'projected_elastic'
UNITS = 'source-scaled projected coefficient (V·Å interpretation)'


def owner(hkl):
    hkl = tuple(hkl)
    return tuple(-v for v in hkl) if any(hkl) and next(v for v in hkl if v) < 0 else hkl


def key(record):
    return tuple(record.hkl), tuple(record.interval)


def label(hkl, interval):
    return f'Projected elastic coefficient E_g {tuple(hkl)}, slice [{interval[0]:g}, {interval[1]:g}]'


def components(mode):
    return ('real', 'imaginary') if mode == 'complex' else (mode,)


def coordinates(value, mode):
    return (value.real, value.imag) if mode == 'complex' else (abs(value) if mode == 'magnitude' else value.real,)


def complex_value(values, mode, reference):
    if mode == 'complex': return complex(*values)
    if mode == 'real': return complex(values[0], reference.imag)
    if abs(reference) == 0:
        raise ValueError('Magnitude-only mode requires a nonzero IAM phase')
    if values[0] < 0: raise ValueError('Magnitude must be nonnegative')
    if values[0] == abs(reference):return reference
    return (reference / abs(reference)) * values[0]


def model_reference(compiled):
    c = compiled.source_configuration
    s = compiled.crystallography.original_structure
    value = dict(lattice=s.direct_lattice.vectors.tolist(), sites=[asdict(a) for a in s.sites],
        groups=[asdict(g) for g in c.crystallography.site_groups],
        thermal=[asdict(t) for t in c.crystallography.debye_waller],
        absorption=asdict(c.crystallography.absorption), voltage=c.microscope.energy_kv,
        zone=c.crystallography.zone_axis, boundaries=c.crystallography.slice_boundaries_fractional)
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def iam_grids(compiled):
    sample = compiled.crystallography.simulation_sample
    return tuple(elastic_structure_factors(sample.model.species, atoms, stored.reciprocal_grid,
        applied_debye_waller_angstrom_squared=applied.debye_waller_angstrom_squared,
        vacancy_fraction=applied.vacancy_fraction)
        for atoms, stored, applied in zip(sample.model.atoms_by_slice, sample.prepared.slice_types, sample.parameters.slices))


def catalogue(compiled):
    """All measured reflections plus explicitly configured hkl, independent of beam order."""
    crystal = compiled.crystallography
    c = compiled.source_configuration.crystallography
    records = {key(r): r for r in c.projected_elastic}
    if len(records) != len(c.projected_elastic): raise ValueError('Duplicate projected elastic hkl/slice identity')
    intervals = tuple(zip(c.slice_boundaries_fractional[:-1], c.slice_boundaries_fractional[1:]))
    for r in records.values():
        if r.interval not in intervals: raise ValueError(f'{label(r.hkl,r.interval)}: slice changed; remove/rebind the stale target explicitly')
    reflections = {owner(b.original_hkl) for b in crystal.indexed_beams} | {r.hkl for r in records.values()}
    grids = iam_grids(compiled)
    digest = model_reference(compiled)
    rows = []
    for hkl in sorted(reflections):
        hk = original_hkl_to_internal_HK(hkl, crystal.zone_cell)
        for slice_index, interval in enumerate(intervals):
            grid = crystal.simulation_sample.prepared.slice_types[slice_index].reciprocal_grid
            plus = np.argwhere((grid.indices.h == hk.H) & (grid.indices.k == hk.K))
            minus = np.argwhere((grid.indices.h == -hk.H) & (grid.indices.k == -hk.K))
            reason = None
            if not any(hkl): reason = '(000) is read-only: uniform elastic phase'
            mask = reciprocal_cutoff_mask(grid, c.dsap_inverse_angstrom) & legacy_nyquist_mask(grid)
            positions = None
            if plus.shape != (1,2) or minus.shape != (1,2):
                reason = 'Friedel pair outside the transverse mesh'
                iam = 0j
            else:
                positions = (tuple(plus[0]), tuple(minus[0]))
                iam = complex(grids[slice_index][positions[0]])
                if not all(mask[p] for p in positions): reason = 'Friedel pair excluded by cutoff/Nyquist mask'
            default = 'real' if abs(iam.imag) <= 1e-8 else 'complex'
            record = records.get((hkl,interval))
            mode = record.mode if record else default
            if mode == 'magnitude' and abs(iam) == 0: reason = 'Magnitude-only mode unavailable: IAM phase is undefined at zero'
            initial = complex(*record.initial) if record and record.initial is not None else iam
            stale = bool(record and record.initial is not None and record.reference != digest)
            rows.append(dict(hkl=hkl,interval=interval,slice_index=slice_index,positions=positions,
                iam=iam, initial=initial,mode=mode,record=record,reference=digest,stale=stale,
                reason=reason,label=label(hkl,interval),linked_hkl=tuple(-v for v in hkl)))
    return tuple(rows)


def active(row, config):
    return bool(row['record'] and row['record'].refine and
        (config.refinement.active_categories is None or CATEGORY in config.refinement.active_categories))


def parameter_groups(rows, config):
    """Optional equality groups share optimizer coordinates, never Friedel owners."""
    grouped={}
    for i,row in enumerate(rows):
        record=row['record']
        name=('group',record.equality_group) if record.equality_group else ('target',i)
        grouped.setdefault(name,[]).append((i,row))
    result=[]
    for identity,members in grouped.items():
        i,first=members[0];mode=first['mode'];n=len(components(mode))
        if identity[0]=='group':
            prefix=f'Linked refinement group {identity[1]}: '
            if len(members)<2:raise ValueError(prefix+'requires at least two physical members; unlink the remaining member')
            if any(r['mode']!=mode or r['interval']!=first['interval'] for _,r in members):
                raise ValueError(prefix+'members must have the same mode and physical slice interval')
            if any(active(r,config)!=active(first,config) for _,r in members):
                raise ValueError(prefix+'all members must share Refine activation')
            if any(not np.allclose(coordinates(r['initial'],mode),coordinates(first['initial'],mode),rtol=1e-12,atol=1e-12) for _,r in members):
                raise ValueError(prefix+'different starting values require an explicit common starting value; no averaging')
        else:prefix=first['label']+': '
        lower=tuple(max((r['record'].lower[j] for _,r in members if r['record'].lower[j] is not None),default=None) for j in range(n))
        upper=tuple(min((r['record'].upper[j] for _,r in members if r['record'].upper[j] is not None),default=None) for j in range(n))
        if identity[0]=='group':
            for v,lo,hi in zip(coordinates(first['initial'],mode),lower,upper):
                if (lo is not None and hi is not None and lo>=hi) or (lo is not None and v<lo) or (hi is not None and v>hi):
                    raise ValueError(prefix+'bound intersection is empty or excludes the chosen start; define compatible member bounds')
        result.append((i,members,lower,upper))
    return tuple(result)


def attach_coefficients(compiled):
    """Append physical components to the modern layout; never change legacy state layout."""
    from .refinement import RefinementParameterSpec, RefinementParameterLayout
    from .refinement_bounds import RefinementBound, RefinementBounds
    rows = catalogue(compiled)
    stored_rows = tuple(r for r in rows if r['record'] is not None)
    specs = list(compiled.layout.parameters)
    bounds = list(compiled.bounds.entries)
    steps = list(compiled.plamda)
    for row in stored_rows:
        if row['reason'] and (row['record'].refine or row['record'].initial is not None): raise ValueError(row['label']+': '+row['reason'])
    for record_index,members,lower,upper in parameter_groups(stored_rows,compiled.source_configuration):
        row=members[0][1];record=row['record']
        if not active(row, compiled.source_configuration): continue
        for _,member in members:
            member['initial']=complex_value(coordinates(row['initial'],row['mode']),row['mode'],member['iam'] if row['mode']=='magnitude' else member['initial'])
        for j, v in enumerate(coordinates(row['initial'], row['mode'])):
            specs.append(RefinementParameterSpec(len(specs),CATEGORY,(record_index,j),float(v),(1.05,.95)))
            if lower[j] is not None or upper[j] is not None:
                bounds.append(RefinementBound(CATEGORY,(record_index,j),lower[j],upper[j]))
            steps.append(record.step_size)
    optimizer = replace(compiled.optimizer, plamda=tuple(steps))
    result = replace(compiled, layout=RefinementParameterLayout(tuple(specs)),
        bounds=RefinementBounds(tuple(bounds)), optimizer=optimizer,
        crystallography=replace(compiled.crystallography, projected_elastic=stored_rows))
    return replace(result, crystallography=replace(result.crystallography,
        simulation_sample=apply_grids(result, stored_rows)))


def apply_grids(compiled, rows):
    sample = compiled.crystallography.simulation_sample
    changed = [r for r in rows if r['initial'] != r['iam'] and r['positions'] is not None]
    if not changed: return sample  # Exact no-change baseline, including printed-phase residuals.
    grids = list(iam_grids(compiled))
    affected = set()
    for row in changed:
        i = row['slice_index']; plus, minus = row['positions']
        grids[i][plus] = row['initial']; grids[i][minus] = row['initial'].conjugate()
        affected.add(i)
    stored = list(sample.prepared.slice_types)
    for i in sorted(affected):
        grid = stored[i].reciprocal_grid
        p = sample.parameters; a = p.slices[i]
        absorption = (absorptive_structure_factors(sample.model.species,sample.model.atoms_by_slice[i],grid,
            voltage_kv=p.applied_voltage_kv,applied_debye_waller_angstrom_squared=a.debye_waller_angstrom_squared,
            coefficients=p.absorption_coefficients,vacancy_fraction=a.vacancy_fraction,
            absorption_raw=sample.model.absorption_raw) if p.absorption_coefficients is not None else np.zeros_like(grids[i]))
        grating = phase_grating_from_structure_factors(grids[i],absorption,grid,
            voltage_kv=p.applied_voltage_kv,dsap_inverse_angstrom=p.dsap_inverse_angstrom,method=p.phase_grating_method)
        stored[i] = replace(stored[i],phase_grating=grating)
    return replace(sample, prepared=replace(sample.prepared,slice_types=tuple(stored)))


def layout_rows(compiled):
    result = []
    for i,members,lower,upper in parameter_groups(compiled.crystallography.projected_elastic,compiled.source_configuration):
        row=members[0][1]
        if not active(row,compiled.source_configuration): continue
        r = row['record']
        metadata = dict(hkl=list(row['hkl']),interval=list(row['interval']),mode=row['mode'],
            iam=[row['iam'].real,row['iam'].imag],initial_complex=[row['initial'].real,row['initial'].imag],
            reference=row['reference'],linked_hkl=list(row['linked_hkl']),stale=row['stale'])
        if r.equality_group:
            metadata['equality_group']=r.equality_group
            metadata['members']=[dict(hkl=list(m['hkl']),interval=list(m['interval']),mode=m['mode'],iam=[m['iam'].real,m['iam'].imag],initial_complex=[m['initial'].real,m['initial'].imag],linked_hkl=list(m['linked_hkl'])) for _,m in members]
        group_label=('Linked coefficient group '+r.equality_group+' {'+', '.join(str(m['hkl']) for _,m in members)+'}, slice '+str(row['interval'])) if r.equality_group else row['label']
        for j, (component,value) in enumerate(zip(components(row['mode']),coordinates(row['initial'],row['mode']))):
            result.append(dict(target_id=f'{CATEGORY}:{(r.equality_group+":") if r.equality_group else ""}{json.dumps([row["hkl"],row["interval"],component],separators=(",",":"))}',
                label=group_label+' / '+component,category=CATEGORY,scope='slice-reflection',site_group=None,
                initial_value=float(value),lower=lower[j],upper=upper[j],units=UNITS,
                binding=dict(category=CATEGORY,identifiers=[*map(str,row['hkl']),*map(str,row['interval']),component]),
                internal_target=[i,j],selected=True,projected_elastic=dict(metadata,component=component)))
    return tuple(result)


def check_ready(compiled):
    for row in compiled.crystallography.projected_elastic:
        if row['stale'] and (row['record'].initial is not None or active(row,compiled.source_configuration)):
            raise ValueError(row['label']+': stale IAM reference; Reset to IAM or recommit the starting value')
        if row['mode']=='magnitude' and row['record'].initial is not None and not row['stale']:
            expected=complex_value((abs(row['initial']),),'magnitude',row['iam'])
            if not np.isclose(row['initial'],expected,rtol=1e-14,atol=1e-14):
                raise ValueError(row['label']+': magnitude starting value must retain IAM phase')
        if active(row, compiled.source_configuration) and row['mode'] == 'magnitude' and row['record'].lower[0] is not None and row['record'].lower[0] < 0:
            raise ValueError(row['label']+': magnitude lower bound must be nonnegative')


def trial_records(compiled, vector):
    records = list(compiled.source_configuration.crystallography.projected_elastic)
    lookup = {key(r):i for i,r in enumerate(records)}
    for i,members,lower,upper in parameter_groups(compiled.crystallography.projected_elastic,compiled.source_configuration):
        if not active(members[0][1],compiled.source_configuration):continue
        values=[float(v) for p,v in zip(compiled.layout.parameters,vector) if p.category==CATEGORY and p.target[0]==i]
        for _,row in members:
            value=complex_value(values,row['mode'],row['iam'] if row['mode']=='magnitude' else row['initial'])
            records[lookup[(row['hkl'],row['interval'])]]=replace(row['record'],initial=(value.real,value.imag),reference=None)
    return tuple(records)


def results_coefficients(layout, initial, *, displayed=None, latest=None, best=None, retained=None):
    """Reconstruct physical complex values only from recorded layout/vector metadata."""
    groups = {}
    for index,p in enumerate(layout):
        meta = p.get('projected_elastic')
        if meta is not None:
            for member in meta.get('members',[meta]):
                identity = (tuple(member['hkl']),tuple(member['interval']))
                item=dict(p,projected_elastic=dict(meta,**member))
                groups.setdefault(identity,[]).append((index,item))
    result=[]
    for (hkl,interval), entries in groups.items():
        meta=entries[0][1]['projected_elastic'];reference=complex(*(meta['iam'] if meta['mode']=='magnitude' else meta['initial_complex']))
        def restore(vector):
            if vector is None:return None
            values=[vector[i] for i,_ in entries]
            return complex_value(values,meta['mode'],reference)
        states={name:restore(v) for name,v in dict(initial=initial,displayed=displayed,latest=latest,best=best,retained=retained).items()}
        value=states['displayed'] if states['displayed'] is not None else states['initial']
        result.append(dict(equality_group=meta.get('equality_group'),hkl=hkl,interval=interval,mode=meta['mode'],iam=complex(*meta['iam']),**states,
            real=value.real,imaginary=value.imag,magnitude=abs(value),phase_radians=float(np.angle(value)) if value != 0 else None,
            imaginary_fixed=meta['mode']=='real',bounds=[(p['lower'],p['upper']) for _,p in entries],
            parameter_indices=[i for i,_ in entries],linked_hkl=meta['linked_hkl']))
    return tuple(result)
