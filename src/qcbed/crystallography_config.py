"""File-input boundary and pure crystallography-to-production configuration bridge."""
from dataclasses import dataclass, replace
from pathlib import Path
import numpy as np

from .configuration import (SampleConfig, SpeciesConfig, RegionConfig, StoredSliceConfig,
    AtomConfig, ScanConfig, ExperimentalDataConfig, CrystalScanConfig,
    CrystalExperimentalDataConfig)
from .crystallography import (AtomicStructureModel, ZoneAxis, CrystalOrientation,
    ZoneAxisCell, OriginalReflectionBasis, index_experimental_beams,
    build_zone_axis_simulation_sample, ZoneAxisSimulationSample, index_prepared_dataset,
    adapt_indexed_prepared_dataset, internal_scan_basis)
from .absorption import AbsorptionCoefficients


class CrystallographyConfigError(ValueError):
    def __init__(self, field, message):
        self.field = field
        super().__init__(f'crystallography.{field}: {message}')


@dataclass(frozen=True)
class CompiledCrystallography:
    original_structure: AtomicStructureModel
    zone_axis: ZoneAxis
    orientation: CrystalOrientation
    zone_cell: ZoneAxisCell
    reflection_basis: OriginalReflectionBasis
    simulation_sample: ZoneAxisSimulationSample
    source_path: str
    scan_basis_transform: tuple[tuple[int, int], tuple[int, int]]
    scan_config: CrystalScanConfig
    site_groups: tuple = ()
    target_aliases: tuple = ()

    @property
    def diagnostics(self):
        return dict(source_path=self.source_path,atom_count=len(self.original_structure.sites),
            cell_volume_angstrom_cubed=float(np.linalg.det(self.original_structure.direct_lattice.vectors)),
            zone_axis=self.zone_axis.uvw,g1_hkl=self.reflection_basis.hkl_1,g2_hkl=self.reflection_basis.hkl_2,
            transform=self.zone_cell.transform.tolist(),strict_b5_compatible=True,
            simulation_atom_count=sum(len(s.sites) for s in self.simulation_sample.depth_sample.slices),
            simulation_slice_count=len(self.simulation_sample.depth_sample.slices))

    @property
    def simulation_cell(self):
        return self.simulation_sample.depth_sample.simulation_cell

    @property
    def indexed_beams(self):
        return self.simulation_sample.indexed_beams

    def index_dataset(self, dataset):
        return index_prepared_dataset(dataset,self.reflection_basis,self.zone_cell)

    def adapt_dataset(self, dataset, *, objective_mode):
        for name in ("scan_rows","scan_columns","corner_a","corner_b","corner_c"):
            if getattr(dataset.scan_geometry,name) != getattr(self.scan_config,name):
                raise ValueError(f"prepared dataset {name} differs from compiled scan configuration")
        indexed = self.index_dataset(dataset)
        expected = tuple((b.stable_id,b.experimental_index) for b in self.indexed_beams)
        actual = tuple((b.stable_id,b.experimental_index) for b in indexed.indexed_beams)
        if actual != expected:
            raise ValueError('prepared dataset beam IDs/order/coordinates differ from compiled configuration')
        return adapt_indexed_prepared_dataset(indexed,
            scan_basis_transform=self.scan_basis_transform,objective_mode=objective_mode)


def load_config_structure(config):
    """Read a single local ASE frame, returning canonical QCBED objects only.

    This is deliberately outside compile_config. No implicit CWD base, wrapping,
    frame selection, thermal-factor inference or disordered-CIF approximation.
    """
    source = config.crystallography.structure
    path = Path(source.path)
    if not path.is_absolute():
        raise CrystallographyConfigError('structure.path','resolve against an explicit document base first')
    if not path.is_file():
        raise CrystallographyConfigError('structure.path',f'structure file not found: {path}')
    try:
        from ase.io import read
    except ImportError as exc:
        raise CrystallographyConfigError('structure.path',
            'ASE is required; install qcbed[crystallography] in this environment') from exc
    try:
        frames = read(path,index=':',format=source.format)
        if len(frames) != 1:
            raise ValueError('exactly one structure frame is required')
        atoms = frames[0]
        if atoms.info.get('occupancy'):
            # ASE CIF retains explicit occupancies here; refuse unsupported disorder.
            if any(len(v)!=1 or next(iter(v.values())) != 1 for v in atoms.info['occupancy'].values()):
                raise ValueError('disordered/partial CIF occupancy requires an explicit canonical structure')
        from .crystallography.ase_adapter import atomic_structure_from_ase
        return atomic_structure_from_ase(atoms)
    except Exception as exc:
        raise CrystallographyConfigError('structure.path',f'cannot read structure {path}: {exc}') from exc


def compile_config_from_sources(config):
    """I/O compilation frontend for already document-resolved paths."""
    from .config_compiler import compile_config
    return compile_config(config,structure=load_config_structure(config)
                          if config.crystallography is not None else None)


def compile_crystallography_config(config, *, structure):
    """Pure adapter; reuse B1–B5 and the existing pre-sliced compiler."""
    from .config_compiler import compile_config
    c = config.crystallography
    if config.sample is not None:
        raise CrystallographyConfigError('structure','sample and crystallography modes are mutually exclusive; set sample: null')
    if not isinstance(config.scan,CrystalScanConfig) or not isinstance(config.experiment,CrystalExperimentalDataConfig):
        raise CrystallographyConfigError('structure','crystallography mode requires scan without basis_1/basis_2 and experiment without beam_indices_hk')
    if not isinstance(structure,AtomicStructureModel):
        raise CrystallographyConfigError('structure.path',
            'unresolved structure: use compile_config_from_sources(), or pass canonical structure= to pure compile_config()')
    try:
        zone = ZoneAxis(c.zone_axis)
        if ZoneAxis(config.microscope.zone_axis) != zone:
            raise ValueError('microscope.zone_axis must match crystallography.zone_axis')
    except (ValueError,TypeError) as exc:
        raise CrystallographyConfigError('zone_axis',str(exc)) from exc
    try:
        basis = OriginalReflectionBasis(c.cbed_basis.g1_hkl,c.cbed_basis.g2_hkl,zone)
    except (ValueError,TypeError) as exc:
        raise CrystallographyConfigError('cbed_basis',str(exc)) from exc
    aliases = ()
    if c.site_groups:
        from .crystallography.site_groups import group_compiler_inputs, remap_refinement
        structure, thermal_records, absorption_config, aliases = group_compiler_inputs(config, structure)
        c = replace(c,debye_waller=thermal_records,absorption=absorption_config)
        effective_refinement = remap_refinement(config.refinement,aliases)
    else:
        effective_refinement = config.refinement
    cell = ZoneAxisCell(structure.direct_lattice,zone)
    symbols = tuple(dict.fromkeys(s.symbol for s in structure.sites))
    thermal = {t.species:t.value_angstrom_squared for t in c.debye_waller}
    if len(thermal)!=len(c.debye_waller) or set(thermal)!=set(symbols):
        raise CrystallographyConfigError('debye_waller','No complete thermal model: provide explicit site groups covering every atom with a B value, or a legacy per-element debye_waller mapping. Site-group configs do not require the legacy mapping.')
    if any(s.debye_waller_b is not None and s.debye_waller_b != thermal[s.symbol] for s in structure.sites):
        raise CrystallographyConfigError('debye_waller','explicit structure B conflicts with configured B; no overwrite')
    applied = replace(structure,sites=tuple(replace(s,debye_waller_b=thermal[s.symbol]) for s in structure.sites))
    try:
        beams = index_experimental_beams(((b.id,b.lattice_coordinates) for b in c.beams),basis,cell)
    except (ValueError,TypeError) as exc:
        raise CrystallographyConfigError('beams',str(exc)) from exc
    if c.absorption.enabled:
        coeff = {v.species:v for v in c.absorption.coefficients}
        if len(coeff)!=len(c.absorption.coefficients) or set(coeff)!=set(symbols):
            raise CrystallographyConfigError('absorption','require one coefficient record per chemical symbol')
        if any(v.c_mode != 'independent' for v in coeff.values()):
            raise CrystallographyConfigError('absorption','linked C is not supported in the crystallography bridge')
        absorption = tuple(AbsorptionCoefficients(coeff[s].a,coeff[s].b,coeff[s].c) for s in symbols)
    else:
        if c.absorption.coefficients:
            raise CrystallographyConfigError('absorption','disabled absorption requires empty coefficients')
        absorption = None
    try:
        sample = build_zone_axis_simulation_sample(applied,cell,
            slice_boundaries_fractional=c.slice_boundaries_fractional,depth_angstrom=c.depth_angstrom,
            mesh=c.mesh,indexed_beams=beams,energy_kv=config.microscope.energy_kv,
            dsap_inverse_angstrom=c.dsap_inverse_angstrom,absorption_coefficients=absorption,
            reciprocal_ordering='unshifted',phase_grating_method='direct')
    except (ValueError,TypeError) as exc:
        raise CrystallographyConfigError('simulation_sample',str(exc)) from exc
    from .crystallography.prepared_indexing import internal_scan_basis_from_crystallography
    try:
        scan_basis = internal_scan_basis_from_crystallography(basis,cell,scan_basis_transform=c.scan_basis_transform)
    except (ValueError,TypeError) as exc:
        raise CrystallographyConfigError('scan_basis_transform',str(exc)) from exc
    stored = tuple(StoredSliceConfig(f'crystal-slice-{i}',geometry,
        tuple(AtomConfig(f'site-{index}',site.symbol,*site.fractional_position[:2],site.occupancy)
              for index,site in zip(slab.source_indices,slab.sites)))
        for i,(geometry,slab) in enumerate(zip(sample.model.slices,sample.depth_sample.slices)))
    derived = replace(config,crystallography=None,refinement=effective_refinement,
        microscope=replace(config.microscope,zone_axis=zone.uvw),
        sample=SampleConfig(c.mesh,tuple(SpeciesConfig(s,s,sp.atomic_number)
            for s,sp in zip(symbols,sample.model.species)),
            (RegionConfig('crystal',c.depth_angstrom,stored,c.debye_waller),),
            c.dsap_inverse_angstrom,absorption=c.absorption),
        scan=ScanConfig(config.scan.scan_columns,config.scan.scan_rows,*scan_basis,
            config.scan.corner_a,config.scan.corner_b,config.scan.corner_c),
        experiment=ExperimentalDataConfig(config.experiment.path,
            tuple((b.internal_hk.H,b.internal_hk.K) for b in beams),config.experiment.gradient_option,
            config.experiment.header,config.experiment.format))
    # Output selection is an execution concern. The legacy adapter requires a
    # nonempty directory spelling even for a pure model preview; never persist
    # this placeholder or use it as a run destination.
    preview=derived if derived.run.output_directory.strip() else replace(derived,
        run=replace(derived.run,output_directory='.'))
    compiled = compile_config(preview)
    if preview is not derived:compiled=replace(compiled,configuration=derived)
    bundle = CompiledCrystallography(structure,zone,CrystalOrientation(structure.direct_lattice,zone),
        cell,basis,sample,c.structure.path,c.scan_basis_transform,config.scan,
        c.site_groups,aliases)
    return replace(compiled,crystallography=bundle,source_configuration=config)
