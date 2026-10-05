"""Narrow strict B5 bridge into the existing immutable propagation snapshot."""
from dataclasses import dataclass
import math
import numpy as np

from .simulation_cell import DepthResolvedSample, build_depth_resolved_sample
from .structure import _ELEMENTS
from .index_transform import original_hkl_to_internal_HK
from ..legacy_io import LegacyInput, Species, AtomSite
from ..models import SliceGeometry, MeshConfig
from ..simulation import (SamplePreparationParameters, SliceAppliedParameters,
                          PreparedMultisliceSample, prepare_multislice_sample)


@dataclass(frozen=True)
class ZoneAxisSimulationSample:
    depth_sample: DepthResolvedSample
    indexed_beams: tuple
    model: LegacyInput
    parameters: SamplePreparationParameters
    prepared: PreparedMultisliceSample
    requested_depth_angstrom: float
    effective_depth_angstrom: float


def build_zone_axis_simulation_sample(structure, zone_cell, *,
        slice_boundaries_fractional, depth_angstrom, mesh: MeshConfig,
        indexed_beams, energy_kv, dsap_inverse_angstrom, absorption_coefficients,
        reciprocal_ordering, phase_grating_method) -> ZoneAxisSimulationSample:
    """Prepare compatible depth slabs with unchanged legacy scattering kernels.

    Explicit no-vacancy, unit-dilation, unity-multiplier subset. B must be supplied
    and identical for every site of each species (including empty slabs). No
    default thermal values, species splitting or binding reinterpretation.
    Absorption coefficients follow first-appearance species order; None disables.
    Depth uses the legacy whole-slice cyclic sequence, including final overshoot.
    """
    depth = build_depth_resolved_sample(structure, zone_cell,
                       slice_boundaries_fractional=slice_boundaries_fractional)
    if isinstance(depth_angstrom, bool) or not math.isfinite(depth_angstrom) or depth_angstrom <= 0:
        raise ValueError('depth_angstrom must be finite and positive')
    beams = tuple(indexed_beams)
    for beam in beams:
        if original_hkl_to_internal_HK(beam.original_hkl, zone_cell) != beam.internal_hk:
            raise ValueError('indexed beam HK does not match sample zone cell')
    if len({b.stable_id for b in beams}) != len(beams):
        raise ValueError('indexed beam stable IDs must be unique')
    symbols = tuple(dict.fromkeys(s.symbol for s in structure.sites))
    thermal = []
    for symbol in symbols:
        values = {s.debye_waller_b for s in structure.sites if s.symbol == symbol}
        if None in values or len(values) != 1:
            raise ValueError(f'{symbol}: require explicit identical per-species B; per-site B mapping unsupported')
        thermal.append(values.pop())
    p,q,_ = depth.simulation_cell.direct_lattice.vectors
    a,b = float(np.linalg.norm(p)),float(np.linalg.norm(q))
    gamma = math.degrees(math.atan2(float(np.linalg.norm(np.cross(p,q))),float(np.dot(p,q))))
    geometries = tuple(SliceGeometry(s.thickness_angstrom,cell_a_angstrom=a,
        cell_b_angstrom=b,lattice_gamma_degrees=gamma) for s in depth.slices)
    atoms = tuple(tuple(AtomSite(symbols.index(s.symbol)+1,*s.fractional_position[:2],
        s.occupancy,'S') for s in slab.sites) for slab in depth.slices)
    model = LegacyInput('B5 strict depth sample',energy_kv,False,zone_cell.zone_axis.uvw,
        mesh,geometries,tuple(Species(s,_ELEMENTS.index(s)+1) for s in symbols),atoms,
        '', '', absorption_raw='YE' if absorption_coefficients is not None else 'NO')
    from ..depth_geometry import requested_depth_to_sample_geometry
    extent = requested_depth_to_sample_geometry(depth_angstrom,
        (s.thickness_angstrom for s in geometries))
    sequence, total = extent.sequence, extent.effective_depth
    parameters = SamplePreparationParameters(energy_kv,
        tuple(SliceAppliedParameters(tuple(thermal),0.,1.) for _ in geometries),
        dsap_inverse_angstrom,absorption_coefficients,tuple(sequence),
        tuple((b.internal_hk.H,b.internal_hk.K) for b in beams),
        reciprocal_ordering,phase_grating_method,'legacy','unity')
    prepared = prepare_multislice_sample(model,parameters)
    return ZoneAxisSimulationSample(depth,beams,model,parameters,prepared,
                                    float(depth_angstrom),total)
