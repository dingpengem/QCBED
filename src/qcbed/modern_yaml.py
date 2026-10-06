"""A/B/C document frontend over QcbedConfig; no alternate numerical model.

Prepared-input descriptors reference existing arrays and editor provenance. They
are input documents, distinct from compiled run snapshots and result artifacts.
"""
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path

from .configuration import QcbedConfig, LegacyExportConfig

OPTIMIZER_DEFAULTS = dict(loops=1,maxfsimp=20,maxfpara=3,ftol=1e-6,chimin=1e-6,
                         maxtemp=100.,mintemp=1.,coolrate=.9,seed=2026,paraopt=0,loud=0)
OPTIMIZER_NAMES = dict(loops='cycles',maxfsimp='stage_budget',maxfpara='coordinate_search_budget',
    ftol='objective_tolerance',chimin='stop_objective',maxtemp='maximum_temperature',
    mintemp='minimum_temperature',coolrate='cooling_rate',seed='seed',
    paraopt='coordinate_search_option',loud='verbosity')
WORKFLOWS = {'ANSI':['ANNEAL','AMOEBA'],'ANNE':['ANNEAL'],'SIMP':['AMOEBA'],
             'PARA':['PARAMIN'],'SPAR':['SIMPARA']}


def _mapping(value, allowed, required=(), *, path):
    if not isinstance(value,dict):raise ValueError(f'{path}: expected mapping')
    extra=set(value)-set(allowed);missing=set(required)-set(value)
    if extra:raise ValueError(f'{path}: unknown fields {sorted(extra)}')
    if missing:raise ValueError(f'{path}: missing required fields {sorted(missing)}')
    return value


def _pairs(items):
    result={}
    for k,v in items:
        if k in result:raise ValueError(f'Duplicate descriptor key: {k}')
        result[k]=v
    return result


def _resolve(base,value):
    if not isinstance(value,str):raise ValueError('File reference must be a string')
    return str((base/value).resolve())


def target_name(target):
    """Format an existing catalogue identity; never infer or enumerate targets."""
    if target.category.value=='depth' and target.identifiers==('sample','requested_depth'):
        return 'sample.requested_depth'
    if target.category.value=='debye_waller' and target.identifiers[0]=='site_group':
        return f'site_group.{target.identifiers[1]}.debye_waller_b'
    return target.category.value+':'+json.dumps(target.identifiers,separators=(',',':'))


def _descriptor(config):
    from .config_serialization import _plain
    c=config.crystallography
    body=Path(config.experiment.path)
    side=body.with_suffix(body.suffix+'.prep.json')
    result=dict(format='qcbed-prepared-input-v1',
        experiment=_plain(asdict(config.experiment)),scan=_plain(asdict(config.scan)),
        beams=_plain([asdict(b) for b in c.beams]),scan_basis_transform=_plain(c.scan_basis_transform))
    # Exact body identity remains owned by the existing canonical runner.
    # This link preserves Patch A source, mask and manual-origin provenance.
    if side.exists():
        result['preparation_metadata']=dict(path=str(side),sha256=hashlib.sha256(side.read_bytes()).hexdigest())
    return result


def publish_prepared_descriptor(config, directory):
    """Publish the current A state as an immutable, content-addressed replay descriptor."""
    descriptor=_descriptor(config)
    content=(json.dumps(descriptor,indent=2,allow_nan=False)+'\n').encode()
    folder=Path(directory)/'.qcbed-data';folder.mkdir(exist_ok=True)
    path=folder/(hashlib.sha256(content).hexdigest()+'.prepared.json')
    try:
        with path.open('xb') as f:f.write(content)
    except FileExistsError:
        if path.read_bytes()!=content:raise ValueError('Prepared descriptor collision')
    return path


def encode_modern(config,destination):
    """Write an immutable preparation descriptor, return the compact YAML mapping."""
    from .config_serialization import _plain
    raw=_plain(asdict(config));c=raw['crystallography'];r=raw['refinement']
    path=publish_prepared_descriptor(config,destination.parent)
    descriptor=_descriptor(config)
    atomic={k:c[k] for k in ('structure','zone_axis','cbed_basis')}
    if c['site_groups']:atomic['site_groups']=[{('sites' if k=='site_indices' else k):v for k,v in g.items()} for g in c['site_groups']]
    if c['debye_waller']:atomic['thermal_factors']=c['debye_waller']
    if c['absorption']!={'enabled':False,'coefficients':[],'legacy_raw':None}:atomic['absorption']=c['absorption']
    optimizer={OPTIMIZER_NAMES[k]:v for k,v in r['optimizer'].items() if k!='method' and v!=OPTIMIZER_DEFAULTS.get(k)}
    method=r['optimizer']['method']
    if method not in WORKFLOWS:raise ValueError(f'Unsupported workflow for compact YAML: {method}')
    qcbed=dict(title=config.title,accelerating_voltage_kv=config.microscope.energy_kv,
        signal_mode={0:'conventional',1:'gradient'}[config.experiment.gradient_option],
        requested_depth_angstrom=c['depth_angstrom'],
        numerics=dict(mesh=c['mesh'],slice_boundaries_fractional=c['slice_boundaries_fractional'],
                      reciprocal_cutoff_inverse_angstrom=c['dsap_inverse_angstrom']),
        workflow=WORKFLOWS[method])
    if optimizer:qcbed['optimizer']=optimizer
    ref={}
    if c["projected_elastic"]:ref["projected_elastic"]=[{k:v for k,v in r.items() if k!="equality_group" or v is not None} for r in c["projected_elastic"]]
    if config.refinement.parameters:
        from .optimizer_defaults import recommended_binding
        ref['targets']={}
        for b in config.refinement.parameters:
            try:default=recommended_binding(b.target)
            except ValueError:default=None
            settings={}
            if default is None or b.range_factors!=default.range_factors:settings['range_factors']=list(b.range_factors)
            if default is None or b.plamda!=default.plamda:settings['step_size']=b.plamda
            ref['targets'][target_name(b.target)]=settings
        if len(ref['targets'])!=len(config.refinement.parameters):raise ValueError('Duplicate refinement target')
    bounds=config.refinement.bounds
    if bounds.specific:
        ref['bounds']={target_name(b.target):[b.lower,b.upper] for b in bounds.specific}
        if len(ref['bounds'])!=len(bounds.specific):raise ValueError('Duplicate bound target')
    for field,default in (('iref',-1.),('objective_mode','masked'),('invalid_trial_policy','reject'),('active_categories',None)):
        if r[field]!=default:ref[field]=r[field]
    if bounds.policy!='reject':ref['bounds_policy']=bounds.policy
    if bounds.absolute:ref['category_bounds']=r['bounds']['absolute']
    if bounds.relative:ref['relative_bounds']=r['bounds']['relative']
    if ref:qcbed['refinement']=ref
    if raw['execution']!={'execution_backend':'serial','batch_size':None}:qcbed['execution']=raw['execution']
    output=dict(directory=config.run.output_directory)
    defaults=dict(snapshot_every=100,checkpoint_every=25,max_total_evaluations=1000,max_wall_time_seconds=3600.,max_stage_iterations=None)
    for field,default in defaults.items():
        if raw['run'][field]!=default:output[field]=raw['run'][field]
    qcbed['output']=output
    prep=dict(prepared_dataset=str(path))
    metadata=descriptor.get('preparation_metadata')
    if metadata:
        side=json.loads(Path(metadata['path']).read_text())
        overview=side.get('overview') or {}
        if overview.get('manifest'):
            manifest=json.loads(Path(overview['manifest']).read_text())
            if 'live_preparation' in manifest:
                inputs=manifest['inputs']
                prep=dict(pattern=inputs['pattern']['path'],preparation=str(path))
                if 'calibration_aperture' in inputs:prep['aperture']=inputs['calibration_aperture']['path']
                if 'disk_profile' in inputs:prep['disk_profile']=inputs['disk_profile']['path']
    return dict(schema_version=2,data_prep=prep,atomic_model=atomic,qcbed=qcbed)


def decode_modern(value,base):
    """Decode strict input frontend, then use existing model/compiler/catalogue."""
    from .yaml_config import _decode
    from .config_serialization import _plain
    from .config_paths import resolve_config_paths
    from .config_targets import available_refinement_targets
    from .configuration import ParameterBinding,TargetBound
    _mapping(value,('schema_version','data_prep','atomic_model','qcbed'),('data_prep','atomic_model','qcbed'),path='config')
    if type(value.get('schema_version',2)) is not int or value.get('schema_version',2)!=2:raise ValueError('Expected modern YAML schema_version 2')
    prep=_mapping(value['data_prep'],('prepared_dataset','preparation','pattern','aperture','disk_profile'),path='data_prep')
    if 'prepared_dataset' in prep and 'preparation' in prep:raise ValueError('Select one preparation reference')
    reference=prep.get('preparation',prep.get('prepared_dataset'))
    if reference is None:raise ValueError('Raw project requires interactive Data Preparation before canonical compilation')
    for name in ('pattern','aperture','disk_profile'):
        if name in prep:_resolve(base,prep[name])
    file=Path(_resolve(base,reference))
    descriptor=json.loads(file.read_text(),object_pairs_hook=_pairs)
    _mapping(descriptor,('format','experiment','scan','beams','scan_basis_transform','preparation_metadata'),
        ('format','experiment','scan','beams','scan_basis_transform'),path='prepared_dataset')
    if descriptor['format']!='qcbed-prepared-input-v1':raise ValueError('Unknown prepared-input descriptor format')
    experiment=dict(descriptor['experiment'])
    experiment['path']=_resolve(file.parent,experiment['path'])
    metadata=descriptor.get('preparation_metadata')
    if metadata is not None:
        _mapping(metadata,('path','sha256'),('path','sha256'),path='preparation_metadata')
        side=Path(_resolve(file.parent,metadata['path']))
        expected=Path(experiment['path']).with_suffix(Path(experiment['path']).suffix+'.prep.json')
        if side!=expected or hashlib.sha256(side.read_bytes()).hexdigest()!=metadata['sha256']:
            raise ValueError('Prepared mask/origin metadata reference or checksum mismatch')
    if 'pattern' in prep:
        if metadata is None:raise ValueError('Raw input references require a preparation provenance record')
        mask_meta=json.loads(side.read_text())
        overview=mask_meta.get('overview') or {}
        manifest=json.loads(Path(overview['manifest']).read_text())
        inputs=manifest.get('inputs',{})
        for name,record_name in (('pattern','pattern'),('aperture','calibration_aperture'),('disk_profile','disk_profile')):
            if name in prep:
                source=inputs.get(record_name,{})
                if _resolve(base,prep[name])!=source.get('path') or hashlib.sha256(Path(source['path']).read_bytes()).hexdigest()!=source.get('sha256'):
                    raise ValueError(f'{name}: raw source differs from cached preparation; prepare again')
    atomic=_mapping(value['atomic_model'],('structure','zone_axis','cbed_basis','site_groups','thermal_factors','absorption'),
        ('structure','zone_axis','cbed_basis'),path='atomic_model')
    q=_mapping(value['qcbed'],('title','accelerating_voltage_kv','signal_mode','requested_depth_angstrom','numerics',
        'workflow','optimizer','refinement','execution','output','legacy_export'),
        ('accelerating_voltage_kv','signal_mode','requested_depth_angstrom','numerics','workflow','output'),path='qcbed')
    mode=q['signal_mode']
    if mode not in ('gradient','conventional'):raise ValueError('qcbed.signal_mode: expected gradient or conventional')
    expected_mode={'conventional':0,'gradient':1}[mode]
    if type(experiment.get('gradient_option')) is not int or experiment['gradient_option']!=expected_mode:
        raise ValueError('Signal mode differs from prepared descriptor; select matching prepared data')
    numerics=_mapping(q['numerics'],('mesh','slice_boundaries_fractional','reciprocal_cutoff_inverse_angstrom'),
        ('mesh','slice_boundaries_fractional','reciprocal_cutoff_inverse_angstrom'),path='qcbed.numerics')
    workflow=q['workflow']
    matches=[k for k,v in WORKFLOWS.items() if v==workflow]
    if not matches:raise ValueError('Unsupported qcbed.workflow sequence')
    optimizer_input=_mapping(q.get('optimizer',{}),OPTIMIZER_NAMES.values(),path='qcbed.optimizer')
    optimizer={key:optimizer_input[name] for key,name in OPTIMIZER_NAMES.items() if name in optimizer_input}
    ref=_mapping(q.get('refinement',{}),('targets','bounds','iref','objective_mode','invalid_trial_policy','active_categories',
        'bounds_policy','category_bounds','relative_bounds','projected_elastic'),path='qcbed.refinement')
    output=_mapping(q['output'],('directory','filename','output_file','snapshot_every','checkpoint_every','max_total_evaluations','max_wall_time_seconds','max_stage_iterations'),
        ('directory',),path='qcbed.output')
    directory=_resolve(base,output['directory']) if output['directory'].strip() else ''
    run={k:v for k,v in output.items() if k not in ('directory','filename','output_file')}
    run.update(output_directory=directory,output_filename=str(Path(directory)/'final_out.dat'))
    structure=atomic['structure']
    if isinstance(structure,str):structure=dict(path=structure)
    crystal=dict(structure=structure,zone_axis=atomic['zone_axis'],cbed_basis=atomic['cbed_basis'],
        site_groups=[{('site_indices' if k=='sites' else k):v for k,v in _mapping(g,('id','element','sites','debye_waller_b','occupancy','absorption'),('id','element','sites','debye_waller_b'),path='atomic_model.site_groups').items()} for g in atomic.get('site_groups',[])],debye_waller=atomic.get('thermal_factors',[]),
        absorption=atomic.get('absorption',{}),projected_elastic=ref.get('projected_elastic',[]),mesh=numerics['mesh'],
        depth_angstrom=q['requested_depth_angstrom'],slice_boundaries_fractional=numerics['slice_boundaries_fractional'],
        dsap_inverse_angstrom=numerics['reciprocal_cutoff_inverse_angstrom'],beams=descriptor['beams'],scan_basis_transform=descriptor['scan_basis_transform'])
    raw=dict(title=q.get('title','QCBED'),sample=None,microscope=dict(energy_kv=q['accelerating_voltage_kv'],zone_axis=atomic['zone_axis']),
        crystallography=crystal,scan=descriptor['scan'],experiment=experiment,execution=q.get('execution',{}),run=run,
        legacy_export=q.get('legacy_export',{}),refinement=dict(optimizer=dict(OPTIMIZER_DEFAULTS,**optimizer,method=matches[0]),parameters=[],
            iref=ref.get('iref',-1.),objective_mode=ref.get('objective_mode','masked'),invalid_trial_policy=ref.get('invalid_trial_policy','reject'),
            active_categories=ref.get('active_categories'),bounds=dict(policy=ref.get('bounds_policy','reject'),specific=[],
                absolute=ref.get('category_bounds',[]),relative=ref.get('relative_bounds',[]))))
    config=resolve_config_paths(_decode(raw,QcbedConfig,'config'),base_directory=base)
    catalogue={target_name(i.target):i.target for i in available_refinement_targets(config)}
    def target(name):
        if name not in catalogue:raise ValueError(f'Unknown physical refinement target: {name}')
        return catalogue[name]
    bindings=[];bounds=[]
    for name,v in _mapping(ref.get('targets',{}),catalogue,path='qcbed.refinement.targets').items():
        _mapping(v,('range_factors','step_size'),path=name)
        from .optimizer_defaults import recommended_binding
        default=recommended_binding(target(name)) if 'range_factors' not in v or 'step_size' not in v else None
        bindings.append(_decode(dict(target=_plain(asdict(target(name))),
            range_factors=v.get('range_factors',list(default.range_factors) if default else None),
            plamda=v.get('step_size',default.plamda if default else None)),ParameterBinding,name))
    for name,v in _mapping(ref.get('bounds',{}),catalogue,path='qcbed.refinement.bounds').items():
        if not isinstance(v,list) or len(v)!=2:raise ValueError(f'{name}: bounds require [lower,upper]')
        bounds.append(_decode(dict(target=_plain(asdict(target(name))),lower=v[0],upper=v[1]),TargetBound,name))
    return replace(config,refinement=replace(config.refinement,parameters=tuple(bindings),bounds=replace(config.refinement.bounds,specific=tuple(bounds))))
