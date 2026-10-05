"""B6a read-only input consistency gate. No propagation/objective imports or calls."""
from dataclasses import asdict
import json
import numpy as np
from .config_compiler import CompiledConfig
from .crystallography import simulation_cell_compatibility
from .crystallography.orientation import zone_law_value
from .crystallography.prepared_indexing import adapt_indexed_prepared_dataset
from .data_preparation import prepared_dataset_manifest

# Same float64 metric tolerance as B4b noncubic reciprocal-vector regressions.
ATOL=1e-14
RTOL=1e-14


class PreflightError(ValueError):
    def __init__(self,code,message):
        self.code=code
        super().__init__(message)


def _require(condition,code,message):
    if not condition: raise PreflightError(code,message)


def _close(actual,expected,code,message):
    _require(np.allclose(actual,expected,atol=ATOL,rtol=RTOL),code,message)


def require_physical_config(config):
    """No inference from legacy zone labels or the demonstration structure."""
    _require(config.crystallography is not None,'missing_physical_inputs',
        'Missing original crystal structure, original zone axis and original-cell g1_hkl/g2_hkl: '
        'an explicit crystallography configuration is required; legacy labels are not sufficient.')


def _statistics(values):
    return dict(min=float(values.min()),max=float(values.max()),median=float(np.median(values)),
        nan_count=int(np.isnan(values).sum()),inf_count=int(np.isinf(values).sum()),
        nonpositive_count=int((values<=0).sum()))


def preflight_inputs(compiled, datasets, *, baseline_manifests, input_kind,
                     physical_input_confirmation=None):
    """Check compiled B4b sample against unchanged prepared dataset(s).

    Baseline manifests must come from the preceding validated preparation stage,
    not be generated from a potentially edited input inside this function.
    Real input requires an explicit confirmation record; this checks consistency,
    not the experimental truth of an unconfirmed indexing assignment.
    Returns JSON-compatible diagnostics only, without modifying any input.
    """
    _require(input_kind in ('synthetic','confirmed_real'),'input_kind','Declare synthetic or confirmed_real input')
    if input_kind=='confirmed_real':
        _require(isinstance(physical_input_confirmation,str) and bool(physical_input_confirmation.strip()),
            'unconfirmed_physical_inputs','Real structure, zone and g1_hkl/g2_hkl require explicit confirmation provenance')
    _require(isinstance(compiled,CompiledConfig) and compiled.crystallography is not None,
        'missing_physical_inputs','Compiled B4b crystallography is required; no indexing is inferred')
    c=compiled.crystallography
    sample=c.simulation_sample
    depth=sample.depth_sample
    cell=c.zone_cell
    _require(depth.simulation_cell.zone_cell is cell,'sample_basis_provenance',
        'HK indexing and simulation sample must reference the exact same ZoneAxisCell object')
    _require(c.orientation.direct_lattice is cell.original_lattice and
        c.original_structure.direct_lattice is cell.original_lattice,'original_basis_provenance',
        'Original structure, orientation and zone cell provenance differ')
    _require(c.zone_axis==cell.zone_axis==c.orientation.zone_axis==c.reflection_basis.zone_axis,
             'zone_axis_mismatch','Zone-axis provenance differs')
    compatibility=simulation_cell_compatibility(cell)
    _require(compatibility.compatible,'strict_b5',compatibility.reason)
    _require(compiled.state.model.slices==sample.model.slices,'sample_state_mismatch','Compiled slice geometries differ from B5 sample')
    _require(tuple((s.label,s.atomic_number) for s in compiled.state.model.species)==
             tuple((s.label,s.atomic_number) for s in sample.model.species),
             'sample_state_mismatch','Compiled species differ from B5 sample')
    def atoms(model):
        return tuple(tuple((a.species_index_one_based,a.x_fractional,a.y_fractional,a.starting_occupancy,a.record_type)
                           for a in row) for row in model.atoms_by_slice)
    _require(atoms(compiled.state.model)==atoms(sample.model),'sample_state_mismatch','Compiled atoms differ from B5 sample')
    _require(len(compiled.state.regions)==1 and compiled.state.regions[0].depth_angstrom==sample.requested_depth_angstrom,
        'sample_state_mismatch','Compiled region depth differs from B5 sample')
    thermal=tuple(r.debye_waller for r in compiled.state.regions for _ in range(r.slice_count))
    _require(thermal==tuple(s.debye_waller_angstrom_squared for s in sample.parameters.slices),
        'thermal_mismatch','Compiled and prepared B factors differ')
    _require(all(s.debye_waller_b is not None and np.isfinite(s.debye_waller_b) and s.debye_waller_b>=0
                 and 0<s.occupancy<=1 for s in depth.original_structure.sites),
        'thermal_occupancy','Missing or invalid thermal/occupancy metadata')
    _require(depth.recast_structure.direct_lattice is cell.direct_lattice,
        'sample_basis_provenance','Recast atoms and zone cell provenance differ')
    def original_atoms(structure):
        return tuple((s.symbol,s.fractional_position,s.occupancy) for s in structure.sites)
    _require(original_atoms(c.original_structure)==original_atoms(depth.original_structure),
        'original_atom_mismatch','Original structure positions/species/occupancy differ from B5 source')
    _require(all(a.debye_waller_b is None or a.debye_waller_b==b.debye_waller_b
                 for a,b in zip(c.original_structure.sites,depth.original_structure.sites)),
        'thermal_mismatch','Explicit original site B differs from B5 source')
    _require(sorted(i for slab in depth.slices for i in slab.source_indices)==list(range(len(c.original_structure.sites))),
        'atom_ownership','Slice source indices must retain each original site exactly once')
    _require(len(c.original_structure.sites)==len(depth.recast_structure.sites)==
             sum(len(s.sites) for s in depth.slices),'atom_count','Atom ownership/count mismatch')
    z=c.orientation.beam_direction_cartesian
    cz=cell.direct_lattice.vectors[2]/np.linalg.norm(cell.direct_lattice.vectors[2])
    sz=depth.simulation_cell.frame[2]
    _close(z,cz,'beam_direction','B2 and B3 propagation directions differ')
    _close(z,sz,'beam_direction','B2 and B5 propagation directions differ')
    _close(depth.simulation_cell.direct_lattice.reciprocal().vectors[:2],
           cell.reciprocal_lattice.vectors[:2],'reciprocal_mismatch','Simulation and zone transverse reciprocal bases differ')
    original_reciprocal=c.original_structure.direct_lattice.reciprocal().vectors
    g1=np.array(c.reflection_basis.hkl_1)@original_reciprocal
    g2=np.array(c.reflection_basis.hkl_2)@original_reciprocal
    table=[]
    for beam in c.indexed_beams:
        n1,n2=beam.experimental_index
        hkl=beam.original_hkl
        _require(hkl==tuple(n1*a+n2*b for a,b in zip(c.reflection_basis.hkl_1,c.reflection_basis.hkl_2)),
            'basis_algebra','Beam original hkl does not follow the configured ordered experimental basis')
        law=zone_law_value(hkl,c.zone_axis)
        _require(law==0,'zone_law',f'{beam.stable_id}: original hkl {hkl} violates zone law: {law}')
        g=np.array(hkl)@original_reciprocal
        gi=np.array(beam.internal_hk.as_HKL())@cell.reciprocal_lattice.vectors
        gb=n1*g1+n2*g2
        _close(g,gi,'reciprocal_mismatch',f'{beam.stable_id}: original/internal reciprocal vectors differ')
        _close(g,gb,'basis_algebra',f'{beam.stable_id}: experimental basis algebra differs')
        parallel=float(abs(np.dot(gi,z)))
        _close(parallel,0.,'reciprocal_plane',f'{beam.stable_id}: reciprocal vector leaves beam plane')
        table.append(dict(stable_id=beam.stable_id,n1=n1,n2=n2,
            original_h=hkl[0],original_k=hkl[1],original_l=hkl[2],
            internal_H=beam.internal_hk.H,internal_K=beam.internal_hk.K,zone_law_value=law,
            reciprocal_vector_error=float(np.linalg.norm(g-gi)),
            reciprocal_component_error=float(np.max(np.abs(g-gi))),
            basis_algebra_error=float(np.linalg.norm(g-gb)),beam_parallel_component=parallel))
    _require(bool(table),'empty_beams','No configured beams')
    expected=tuple((b.stable_id,b.experimental_index) for b in c.indexed_beams)
    hk=tuple((b.internal_hk.H,b.internal_hk.K) for b in c.indexed_beams)
    _require(hk==compiled.experiment.beam_indices_hk==sample.parameters.beam_indices_hk,
             'beam_hk_mismatch','Sample/compiled experimental beam channels differ')
    _require(bool(datasets) and set(datasets)==set(baseline_manifests),
             'baseline_modes','Require matching nonempty datasets and baseline manifests')
    modes={}
    for mode,data in datasets.items():
        _require(mode==data.mode.value,'signal_mode','Dataset mode label mismatch')
        _require(tuple((b.id,n) for b,n in zip(data.beams,data.lattice_coordinates))==expected,
            'beam_order','Dataset beam IDs/order/coordinates differ from compiled configuration')
        _require(data.variances is not None,'variance_missing','Scientific variance is required')
        _require(data.intensities.shape==data.variances.shape==data.weights.shape==
            (data.scan_geometry.scan_rows*data.scan_geometry.scan_columns,len(table)),
            'array_shape','Signal/variance/weight shapes differ from scan and beams')
        active=data.weights>0
        _require(np.any(active) and np.isfinite(data.variances).all() and np.all(data.variances[active]>0),
            'invalid_variance','All active variances must be finite and positive; no floor applied')
        manifest=json.loads(json.dumps(prepared_dataset_manifest(data)))
        baseline=json.loads(json.dumps(baseline_manifests[mode]))
        for name in ('intensities','variances','weights'):
            _require(manifest['arrays'][name]==baseline['arrays'][name],'experimental_regression',
                f'{mode}: {name} differs from preparation baseline')
        _require(manifest['scan_geometry']==baseline['scan_geometry'],'a10_regression','A10 geometry differs from baseline')
        _require(manifest['beam_mapping']==baseline['beam_mapping'],'beam_order','Beam metadata differs from baseline')
        for name in ('corner_a','corner_b','corner_c','scan_rows','scan_columns'):
            _require(getattr(data.scan_geometry,name)==getattr(compiled.state.scan,name),
                     'a10_regression',f'{name} differs from compiled scan')
        indexed=c.index_dataset(data)
        _require(indexed.zone_cell is cell,'sample_basis_provenance','Adapter changed zone-cell provenance')
        _require(indexed.indexed_beams==c.indexed_beams,'beam_hk_mismatch','Adapter and compiled indexing differ')
        try:
            adapted=c.adapt_dataset(data,objective_mode=compiled.runner.objective_mode)
            direct=adapt_indexed_prepared_dataset(indexed,scan_basis_transform=c.scan_basis_transform,
                                                objective_mode=compiled.runner.objective_mode)
        except ValueError as exc:
            raise PreflightError('adapter_rejected',str(exc)) from exc
        _require(adapted.scan==direct.scan==compiled.state.scan,'scan_basis','Compiled scan and direct B4a adapter differ')
        _require(adapted.experimental.beam_indices_hk==hk,'beam_hk_mismatch','Adapter HK differs')
        for name in ('intensities','variances','weights'):
            _require(np.array_equal(getattr(adapted.experimental,name),getattr(data,name)),
                     'adapter_array_mismatch',f'Adapter changed {name}')
        with np.errstate(divide='raise',invalid='raise',over='raise'):
            np.divide(data.weights,data.variances,out=np.zeros_like(data.weights),where=active)
        per_beam=[]
        for i,beam in enumerate(data.beams):
            a=active[:,i]
            _require(np.any(a),'empty_support',f'{beam.id}: no active samples')
            per_beam.append(dict(stable_id=beam.id,active_samples=int(a.sum()),
                signal_min=float(data.intensities[a,i].min()),signal_max=float(data.intensities[a,i].max()),
                variance=_statistics(data.variances[a,i])))
        modes[mode]=dict(signal_shape=list(data.intensities.shape),variance_shape=list(data.variances.shape),
            weight_shape=list(data.weights.shape),beam_count=len(data.beams),samples_per_beam=len(data.intensities),
            active_count=int(active.sum()),support_per_beam=active.sum(axis=0).tolist(),
            active_variance=_statistics(data.variances[active]),beams=per_beam,
            scan_geometry=asdict(data.scan_geometry),signal_regression='exact',a10_regression='exact',
            beam_order_regression='exact',adapter='accepted without evaluation',
            gradient_option=adapted.experimental.gradient_option)
    return dict(status='PASS',input_kind=input_kind,physical_input_confirmation=physical_input_confirmation,
        tolerances=dict(atol=ATOL,rtol=RTOL,integer_relations='exact'),
        beam_table=table,sample=dict(c.diagnostics,
            original_lattice=c.original_structure.direct_lattice.vectors.tolist(),
            zone_lattice=cell.direct_lattice.vectors.tolist(),
            simulation_lattice=depth.simulation_cell.direct_lattice.vectors.tolist(),
            simulation_frame=depth.simulation_cell.frame.tolist(),recast_atom_count=len(depth.recast_structure.sites),
            slice_thicknesses_angstrom=[s.thickness_angstrom for s in depth.slices],
            slice_depth_intervals_fractional=[[s.lower_fraction,s.upper_fraction] for s in depth.slices],
            sequence=list(sample.parameters.slice_sequence_zero_based),
            requested_depth_angstrom=sample.requested_depth_angstrom,effective_depth_angstrom=sample.effective_depth_angstrom),
        max_errors=dict(reciprocal_vector=max(r['reciprocal_vector_error'] for r in table),
            reciprocal_component=max(r['reciprocal_component_error'] for r in table),
            basis_algebra=max(r['basis_algebra_error'] for r in table),
            beam_direction=float(max(np.linalg.norm(z-cz),np.linalg.norm(z-sz))),
            beam_parallel_component=max(r['beam_parallel_component'] for r in table)),
        hk_sample_provenance='identical ZoneAxisCell object',modes=modes,
        configured_gradient_option=compiled.experiment.gradient_option,
        propagation_executed=False,objective_executed=False,
        interpretation='Algebraic/metadata consistency only; does not independently establish experimental indexing.')
