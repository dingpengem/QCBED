"""Offline FCC Al zone projections and supplied-HKL geometry diagnostics."""
import argparse
import json
from pathlib import Path
import numpy as np
from ase.build import bulk
from qcbed.crystallography import (atomic_structure_from_ase,ZoneAxis,CrystalOrientation,
    plot_projected_atomic_structure,project_hkl,zone_law_value,DirectLattice)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('tmp/b2_zone_axis_projection'))
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    structure=atomic_structure_from_ase(bulk('Al','fcc',a=4.05,cubic=True))
    zones=[((0,0,1),[(1,0,0),(0,1,0),(1,1,0)]),
           ((1,1,0),[(1,-1,0),(0,0,1)]),((1,1,1),[(1,-1,0),(1,0,-1)])]
    report=[]
    for uvw,hkls in zones:
        o=CrystalOrientation(structure.direct_lattice,ZoneAxis(uvw))
        fig,_=plot_projected_atomic_structure(structure,o,(2,2,2))
        path=args.output/('zone_'+''.join(map(str,uvw))+'.png');fig.savefig(path,dpi=160)
        import matplotlib.pyplot as plt
        plt.close(fig)
        projections=[project_hkl(structure.direct_lattice.reciprocal(),o,h) for h in hkls]
        entry=dict(zone_axis=o.zone_axis.uvw,beam_direction=o.beam_direction_cartesian.tolist(),
            beam_basis=o.beam_basis.tolist(),figure=str(path),
            supplied_hkls=[dict(hkl=p.hkl,g_cartesian=p.g_cartesian.tolist(),g_beam=p.g_beam.tolist(),
                zone_law=zone_law_value(p.hkl,o.zone_axis)) for p in projections])
        report.append(entry);print(json.dumps(entry,indent=2))
    # Report measured errors for cubic, orthorhombic and triclinic checks.
    errors=[]
    vectors=np.array([[1.2,-3.4,.5],[.2,2.,-1.],[0,0,0]])
    for cell in (np.eye(3)*4,np.diag([3.,4.,5.]),np.array([[3.,0,0],[1.,4.,0],[.5,-.7,5.]])):
        d=DirectLattice(cell)
        for uvw,hkls in zones:
            o=CrystalOrientation(d,ZoneAxis(uvw));r=o.beam_basis
            errors.append(dict(orthonormality_error=float(abs(r@r.T-np.eye(3)).max()),
                determinant=float(np.linalg.det(r)),round_trip_error=float(abs(o.beam_to_cartesian(o.cartesian_to_beam(vectors))-vectors).max()),
                zone_parallel_error=max(abs(project_hkl(d.reciprocal(),o,h).beam_parallel_component) for h in hkls)))
    diagnostics=dict(examples=report,numerical_checks=errors)
    (args.output/'geometry_diagnostics.json').write_text(json.dumps(diagnostics,indent=2)+'\n')
    print('Max errors:',{key:max(e[key] for e in errors) for key in ('orthonormality_error','round_trip_error','zone_parallel_error')})
    print('Determinant range:',min(e['determinant'] for e in errors),max(e['determinant'] for e in errors))


if __name__=='__main__':main()
