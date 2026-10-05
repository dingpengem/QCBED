"""Explicitly revalidate reviewed A6 centers after canonical geometric A5 localization.

This diagnostic does not change automatic A6 selection or invent observations.
Each --defining-center must pass the shared geometry using the NEW central mask.
"""
import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import numpy as np
from qcbed.data_preparation import (DiskGeometry, locate_central_disk,
    LatticeDetectionConfig, LatticeCandidate, infer_lattice_seed, bootstrap_lattice)
from qcbed.data_preparation.arc_catalogue import generate_arc_catalogue, refinement_comparison
from qcbed.data_preparation.arc_geometry import arc_edge_support, _validated_center
from qcbed.data_preparation.central_arc import write_central_arc_localization
from qcbed.data_preparation.diagnostics import (load_pattern, write_seed_diagnostics,
                                               write_bootstrap_diagnostics)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pattern',type=Path,default=Path(__file__).resolve().parents[1]/'reference_data/pattern.dat')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--radii',type=float,nargs=2,required=True)
    parser.add_argument('--geometry-provenance',required=True)
    parser.add_argument('--defining-center',type=float,nargs=2,action='append',required=True)
    args=parser.parse_args()
    if len(args.defining_center)!=2: parser.error('Exactly two previously reviewed defining centers required')
    args.output.mkdir(parents=True,exist_ok=False)
    image=load_pattern(args.pattern);geometry=DiskGeometry(*args.radii)
    central=locate_central_disk(image,geometry)
    write_central_arc_localization(args.output,image,central)
    config=LatticeDetectionConfig()
    arrays,_=arc_edge_support(image,geometry,central.central_disk_center_px,config=config)
    candidates=[]
    for i,center in enumerate(args.defining_center):
        metrics=_validated_center(np.array(center),arrays,geometry,config)
        if metrics is None: raise ValueError(f'Reviewed defining center {center} fails current boundary validation')
        candidates.append(LatticeCandidate(f'revalidated_{i}',tuple(center),metrics['confidence'],**metrics))
    seed=infer_lattice_seed(tuple(candidates),central.central_disk_center_px,geometry,config=config)
    write_seed_diagnostics(args.output,image,seed)
    bootstrap=bootstrap_lattice(image,seed,geometry)
    outputs=write_bootstrap_diagnostics(args.output,image,bootstrap,geometry)
    catalogue=generate_arc_catalogue(bootstrap,image.shape,geometry)
    outputs.update(write_bootstrap_diagnostics(args.output,image,catalogue,geometry,catalogue=True))
    result=dict(refinement_comparison=refinement_comparison(bootstrap),catalogue=asdict(catalogue),counts=catalogue.status_counts,input=str(args.pattern.resolve()),sha256=hashlib.sha256(args.pattern.read_bytes()).hexdigest(),
        shape=image.shape,geometry_provenance=args.geometry_provenance,
        selection='Explicit previously reviewed centers, revalidated against current edges; not automatic candidate selection',
        central=asdict(central),seed=asdict(seed),bootstrap=asdict(bootstrap),outputs=outputs)
    (args.output/'validation_summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    lines=['# A7 arc validation', '', 'Status: '+bootstrap.status,
        'Catalogue: '+catalogue.status, 'Physical radii/provenance: '+args.geometry_provenance,
        '', '| Site | Prediction | Status | Usable circumference | Reason |', '|---|---|---|---|---|']
    lines += [f'| {s.integer_indices} | {s.predicted_center_px} | {s.status} | {s.usable_boundary_fraction:.4f} | {s.reason} |' for s in bootstrap.sites]
    lines += ['', f'Additional independent confirmations: {bootstrap.newly_confirmed_count}',
        f'Required observations: {bootstrap.required_observations}',bootstrap.message,
        'No refined origin/basis/RMS or seed-to-fit changes are available unless a global fit succeeds.',
        'A5/A6 science unchanged. No inputparas. A8 not started.']
    (args.output/'validation_summary.md').write_text('\n'.join(lines)+'\n')
    print('Central:',central.central_disk_center_px)
    print('Basis:',seed.basis_vector_1_px,seed.basis_vector_2_px,'angle:',seed.basis_angle_deg)
    for site in bootstrap.sites:
        print(site.integer_indices,site.predicted_center_px,site.status,
              None if site.candidate is None else site.candidate.center_px)
    print('New confirmations:',bootstrap.newly_confirmed_count)


if __name__=='__main__': main()
