"""Deterministic FCC aluminium B1 example; no external structure database."""
import argparse
from pathlib import Path
from ase.build import bulk
from qcbed.crystallography import atomic_structure_from_ase,plot_atomic_structure


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=Path('tmp/b1_crystal_structure.png'))
    args=parser.parse_args()
    atoms=bulk('Al','fcc',a=4.05,cubic=True)
    model=atomic_structure_from_ase(atoms)
    direct=model.direct_lattice; reciprocal=direct.reciprocal()
    print('Direct lattice (rows, Å):\n',direct.vectors)
    print('Reciprocal lattice (rows, Å^-1, no 2π):\n',reciprocal.vectors)
    print('Cell volume (Å^3):',direct.volume)
    for hkl in ((1,0,0),(1,1,0),(1,1,1),(2,0,0)):
        print(f'd{hkl} = {reciprocal.d_spacing(hkl):.12g} Å (geometric spacing, not an allowed-reflection claim)')
    fig,_=plot_atomic_structure(model,repeats=(2,2,2))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(args.output,dpi=160)
    print('Figure:',args.output)


if __name__=='__main__': main()
