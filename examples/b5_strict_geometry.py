"""Report strict B5 compatibility; no forward simulation or objective evaluation."""
import numpy as np
from qcbed.crystallography import (
    DirectLattice, ZoneAxis, ZoneAxisCell, AtomicSite, AtomicStructureModel,
    simulation_cell_compatibility, build_depth_resolved_sample,
    original_hkl_to_internal_HK,
)


def main():
    for axis in ((0, 0, 1), (1, 1, 0), (1, 1, 1)):
        lattice = DirectLattice(np.diag([4., 4., 4.]))
        cell = ZoneAxisCell(lattice, ZoneAxis(axis))
        print("Original lattice (Å):", lattice.vectors.tolist())
        print("T:", cell.transform.tolist())
        print("Recast direct:", cell.direct_lattice.vectors.tolist())
        print("Recast reciprocal (1/Å):", cell.reciprocal_lattice.vectors.tolist())
        report = simulation_cell_compatibility(cell)
        print(f'{axis}: compatible={report.compatible}; C offsets={report.longitudinal_offsets}')
        print(report.reason)
        if report.compatible:
            structure = AtomicStructureModel(lattice, (
                AtomicSite('Si', (0, 0, .1), debye_waller_b=.3),
                AtomicSite('Si', (.25, .25, .6), debye_waller_b=.3)))
            sample = build_depth_resolved_sample(
                structure, cell, slice_boundaries_fractional=(0, .5, 1))
            print('Simulation cell:', sample.simulation_cell.direct_lattice.vectors.tolist())
            print('Simulation frame:', sample.simulation_cell.frame.tolist())
            print('Atoms before/after:', len(structure.sites), sum(len(s.sites) for s in sample.slices))
            print('Supercell replication: none (unimodular integer shear)')
            errors = []
            for hkl in ((0,0,0),(1,0,0),(0,1,0),(2,-1,0)):
                hk = original_hkl_to_internal_HK(hkl,cell)
                g = np.array(hkl) @ lattice.reciprocal().vectors
                gi = np.array(hk.as_HKL()) @ sample.simulation_cell.direct_lattice.reciprocal().vectors
                error = float(np.linalg.norm(g-gi))
                errors.append(error)
                print('hkl/HK/g_original/g_internal/error:', hkl,hk,g,gi,error)
            print('Maximum reciprocal-vector error (1/Å):', max(errors))
            errors = []
            for slab in sample.slices:
                for index, site in zip(slab.source_indices,slab.sites):
                    xyz = np.array(site.fractional_position) @ sample.simulation_cell.direct_lattice.vectors
                    xyz += sample.c_period_translations[index]*sample.simulation_cell.direct_lattice.vectors[2]
                    original = np.array(structure.sites[index].fractional_position) @ lattice.vectors
                    errors.append(float(np.linalg.norm(xyz-original)))
            print('Maximum Cartesian-position error (Å):',max(errors))
            beam = cell.direct_lattice.vectors[2]
            beam = beam/np.linalg.norm(beam)
            print('Beam-direction mismatch:',float(np.linalg.norm(sample.simulation_cell.frame[2]-beam)))
            for slab in sample.slices:
                print(f'  depth fraction [{slab.lower_fraction}, {slab.upper_fraction}): '
                      f'{slab.thickness_angstrom} Å; source sites {slab.source_indices}')


if __name__ == '__main__':
    main()
