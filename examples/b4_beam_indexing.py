"""Explicit user-basis demonstration, not an indexing assignment of real pattern.dat."""
import argparse
import csv
from dataclasses import asdict
from itertools import product
import json
from pathlib import Path
import numpy as np
from qcbed.crystallography import (DirectLattice,ZoneAxis,OriginalReflectionBasis,
    build_zone_axis_cell,index_experimental_beams,index_experimental_beam,
    original_reflection_basis_vectors)

INDICES=((0,0),(-1,0),(1,0),(0,-1),(0,1),(-1,-1),(1,1),(-1,1),(1,-1))


def metric_checks():
    cases=[('cubic',np.eye(3)*4.05,(1,1,0),((1,-1,0),(0,0,2))),
        ('orthorhombic',np.diag([3.,4.,5.]),(1,1,0),((1,-1,0),(0,0,2))),
        ('hexagonal',np.array([[3.,0,0],[-1.5,3*np.sqrt(3)/2,0],[0,0,5.]]),(0,0,1),((1,0,0),(0,1,0))),
        ('triclinic',np.array([[3.,0,0],[1.,4.,0],[.5,-.7,5.]]),(1,1,1),((1,-1,0),(0,1,-1)))]
    results=[]
    for name,cell,zone,refs in cases:
        d=DirectLattice(cell);z=ZoneAxis(zone);q=build_zone_axis_cell(d,z)
        basis=OriginalReflectionBasis(*refs,z);b=d.reciprocal();g=original_reflection_basis_vectors(b,basis)
        errors=[]
        for n in product(range(-3,4),repeat=2):
            beam=index_experimental_beam(stable_id='test',experimental_index=n,reflection_basis=basis,zone_cell=q)
            errors.append(float(np.linalg.norm(b.g_vector(beam.original_hkl)-np.array(n)@g)))
        results.append(dict(metric=name,zone=zone,basis=refs,pairs_tested=len(errors),maximum_difference_norm=max(errors)))
    return results


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('tmp/b4_beam_indexing'))
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    d=DirectLattice(np.eye(3)*4.05);z=ZoneAxis((1,1,0))
    basis=OriginalReflectionBasis((1,-1,0),(0,0,2),z);q=build_zone_axis_cell(d,z)
    beams=index_experimental_beams([(f'lattice_{n1}_{n2}',(n1,n2)) for n1,n2 in INDICES],basis,q)
    g=original_reflection_basis_vectors(d.reciprocal(),basis)
    rows=[];diagnostics=[]
    print('stable_id n1 n2 original_h original_k original_l internal_H internal_K')
    for beam in beams:
        n1,n2=beam.experimental_index;h,k,l=beam.original_hkl
        row=dict(stable_id=beam.stable_id,n1=n1,n2=n2,original_h=h,original_k=k,original_l=l,
                 internal_H=beam.internal_hk.H,internal_K=beam.internal_hk.K)
        rows.append(row);print(' '.join(str(v) for v in row.values()))
        original=d.reciprocal().g_vector(beam.original_hkl);combined=n1*g[0]+n2*g[1]
        diagnostics.append(dict(**asdict(beam),g_original=original.tolist(),g_from_basis=combined.tolist(),
            difference_norm=float(np.linalg.norm(original-combined))))
    with (args.output/'indexed_beams.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    checks=metric_checks()
    report=dict(source='Demonstration basis, not experimentally established HKLs for pattern.dat',
        correspondence='experimental g1 -> original hkl_1; experimental g2 -> original hkl_2; no swap/sign search',
        zone=z.uvw,user_original_basis=[basis.hkl_1,basis.hkl_2],original_direct_lattice=d.vectors.tolist(),
        zone_transform=q.transform.tolist(),g_vector_units='angstrom^-1 (no 2pi)',beams=diagnostics,
        metric_checks=checks,maximum_reciprocal_consistency_error=max(c['maximum_difference_norm'] for c in checks))
    (args.output/'indexed_beams.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(checks,indent=2))
    print('Maximum reciprocal consistency error:',report['maximum_reciprocal_consistency_error'])


if __name__=='__main__':main()
