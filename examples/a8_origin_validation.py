"""A8 manual beam/feature validation. No real reference or feature is guessed."""
import argparse
import csv
import json
import hashlib
from pathlib import Path
from types import SimpleNamespace
from qcbed.data_preparation.diagnostics import load_pattern
from qcbed.data_preparation import DiskGeometry
from qcbed.data_preparation.origin import (select_symmetry_reference,SymmetrySearchConfig,
    measure_symmetry_origin,write_symmetry_diagnostics)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pattern',type=Path,required=True)
    p.add_argument('--catalogue',type=Path,required=True,help='A7 arc beam_catalogue.csv')
    p.add_argument('--beam-id');p.add_argument('--indices',type=int,nargs=2)
    p.add_argument('--feature-offset',type=float,nargs=2,required=True)
    p.add_argument('--central-center',type=float,nargs=2,required=True)
    p.add_argument('--radii',type=float,nargs=2,required=True)
    p.add_argument('--search-half-size',type=float,nargs=2)
    p.add_argument('--crop-half-size',type=float,nargs=2)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.beam_id is None and a.indices is None: p.error('The symmetry reference beam must be specified')
    with a.catalogue.open() as stream: rows=list(csv.DictReader(stream))
    sites=tuple(SimpleNamespace(beam_id=r['beam_id'],integer_indices=(int(r['n1']),int(r['n2'])),
        predicted_center_px=(float(r['predicted_x']),float(r['predicted_y'])),status=r['status']) for r in rows)
    reference=select_symmetry_reference(SimpleNamespace(sites=sites),beam_id=a.beam_id,
        lattice_indices=a.indices,feature_offset_px=a.feature_offset)
    image=load_pattern(a.pattern)
    result=measure_symmetry_origin(image,DiskGeometry(*a.radii),a.central_center,reference,
        config=SymmetrySearchConfig(a.search_half_size,a.crop_half_size))
    print(write_symmetry_diagnostics(a.output,image,result))
    (a.output/'input_provenance.json').write_text(json.dumps(dict(
        pattern=str(a.pattern.resolve()),pattern_sha256=hashlib.sha256(a.pattern.read_bytes()).hexdigest(),
        catalogue=str(a.catalogue.resolve()),catalogue_sha256=hashlib.sha256(a.catalogue.read_bytes()).hexdigest(),
        shape=image.shape,selection='User-selected diffraction-condition beam and local feature; no automatic selection',
        preprocessing='No smoothing, intensity normalization, or detector correction performed by A8'),indent=2)+'\n')
    print('Zone-axis origin:',result.zone_axis_origin_px)


if __name__=='__main__': main()
