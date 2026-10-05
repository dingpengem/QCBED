"""B3 original-hkl to internal-HK algebra; no experimental beam association."""
import argparse
import json
from itertools import product
from pathlib import Path
import numpy as np
from qcbed.crystallography import (DirectLattice,AtomicSite,AtomicStructureModel,
    ZoneAxis,CrystalOrientation,build_zone_axis_cell,original_hkl_to_zone_HKL,
    zone_HKL_to_original_hkl,original_hkl_to_internal_HK,OriginalReflectionBasis,
    reflection_from_basis_coordinates,original_reflection_basis_vectors,
    recast_structure_to_zone_cell,is_zone_reflection_original,is_zolz_zone_index)


def numerical_checks():
    metrics=dict(reciprocal_vector_invariance_error=0.,cartesian_atom_position_error=0.,
                 beam_direction_mismatch=0.,beam_parallel_reciprocal_error=0.,
                 reciprocal_duality_error=0.,reciprocal_transform_error=0.)
    tested=0
    axes=((0,0,1),(1,0,0),(1,1,0),(1,1,1),(1,-1,0),(2,1,1),(3,2,1))
    positions=np.array([[-.2,1.7,.35],[.4,-.7,2.1]])
    for cell in (np.eye(3)*4,np.diag([3.,4.,5.]),np.array([[3.,0,0],[1.,4.,0],[.5,-.7,5.]])):
        d=DirectLattice(cell);b=d.reciprocal().vectors
        s=AtomicStructureModel(d,tuple(AtomicSite('Al',tuple(p)) for p in positions))
        for axis in axes:
            z=ZoneAxis(axis);q=build_zone_axis_cell(d,z);o=CrystalOrientation(d,z)
            c=q.direct_lattice.vectors[2]; bz=q.reciprocal_lattice.vectors
            recast=recast_structure_to_zone_cell(s,q)
            errors=dict(cartesian_atom_position_error=float(abs(np.array([a.fractional_position for a in recast.sites])@q.direct_lattice.vectors-positions@cell).max()),
                beam_direction_mismatch=float(abs(c/np.linalg.norm(c)-o.beam_direction_cartesian).max()),
                beam_parallel_reciprocal_error=float(abs(bz[:2]@o.beam_direction_cartesian).max()),
                reciprocal_duality_error=float(abs(bz@q.direct_lattice.vectors.T-np.eye(3)).max()),
                reciprocal_transform_error=float(abs(bz-q.inverse_transform.T@b).max()))
            for k,v in errors.items():metrics[k]=max(metrics[k],v)
            for hkl in product(range(-3,4),repeat=3):
                H=original_hkl_to_zone_HKL(hkl,q)
                if zone_HKL_to_original_hkl(H,q)!=hkl:raise ValueError('integer round trip failed')
                if is_zone_reflection_original(hkl,z)!=is_zolz_zone_index(H):raise ValueError('zone law equivalence failed')
                tested+=1
                g=np.array(H)@bz
                metrics['reciprocal_vector_invariance_error']=max(metrics['reciprocal_vector_invariance_error'],float(abs(np.array(hkl)@b-g).max()))
                if H[2]==0:
                    metrics['beam_parallel_reciprocal_error']=max(metrics['beam_parallel_reciprocal_error'],abs(float(g@o.beam_direction_cartesian)))
    return dict(**metrics,exact_round_trips_and_zone_law_checks=tested)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('tmp/b3_zone_axis_cell'))
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    d=DirectLattice(np.eye(3)*4.05)
    examples=[]
    for axis,reflections in (((0,0,1),[(1,0,0),(1,1,0),(0,0,1)]),
            ((1,1,0),[(1,-1,0),(0,0,2),(1,-1,2),(1,0,0)]),
            ((1,1,1),[(1,-1,0),(1,0,-1),(1,1,-2),(1,1,1)])):
        q=build_zone_axis_cell(d,ZoneAxis(axis));rows=[]
        for hkl in reflections:
            H=original_hkl_to_zone_HKL(hkl,q)
            rows.append(dict(original_hkl=hkl,zone_HKL=H,reconstructed_original_hkl=zone_HKL_to_original_hkl(H,q),
                internal_HK=None if H[2] else list(original_hkl_to_internal_HK(hkl,q).as_HKL()[:2])))
        entry=dict(zone_axis=q.zone_axis.uvw,original_direct=d.vectors.tolist(),
            transform=q.transform.tolist(),determinant=1,zone_direct=q.direct_lattice.vectors.tolist(),
            zone_reciprocal=q.reciprocal_lattice.vectors.tolist(),reflection_mappings=rows)
        examples.append(entry);print(json.dumps(entry,indent=2))
    basis=OriginalReflectionBasis((1,-1,0),(0,0,2),ZoneAxis((1,1,0)))
    q=build_zone_axis_cell(d,basis.zone_axis)
    mapping=[]
    for n in ((1,0),(0,1),(1,1),(-1,1)):
        hkl=reflection_from_basis_coordinates(basis,*n)
        mapping.append(dict(basis_coordinates=n,original_hkl=hkl,internal_HK0=original_hkl_to_internal_HK(hkl,q).as_HKL()))
    report=dict(examples=examples,user_original_basis=[basis.hkl_1,basis.hkl_2],
        basis_g_vectors=original_reflection_basis_vectors(d.reciprocal(),basis).tolist(),
        algebra_only_mapping=mapping,checks=numerical_checks())
    (args.output/'zone_cell_diagnostics.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(mapping=mapping,checks=report['checks']),indent=2))


if __name__=='__main__':main()
