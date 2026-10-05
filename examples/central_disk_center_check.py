"""Compatibility comparison of historical intensity and local arc centers; not the canonical workflow."""
import argparse
import hashlib
import json
from pathlib import Path
from qcbed.data_preparation import DiskGeometry, CentralArcConfig, LatticeDetectionConfig, compare_central_disk_centers
from qcbed.data_preparation.diagnostics import load_pattern
from qcbed.data_preparation.central_arc import write_central_center_comparison


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pattern',type=Path,default=Path(__file__).resolve().parents[1]/'reference_data/pattern.dat')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--radii',type=float,nargs=2,required=True)
    parser.add_argument('--geometry-provenance',required=True)
    parser.add_argument('--search-radius',type=float,default=8.)
    parser.add_argument('--search-step',type=float,default=.5)
    parser.add_argument('--radius-tolerance',type=float,default=2.)
    args=parser.parse_args()
    image=load_pattern(args.pattern)
    result=compare_central_disk_centers(image,DiskGeometry(*args.radii),
        search_config=CentralArcConfig(args.search_radius,args.search_step),
        arc_config=LatticeDetectionConfig(radius_tolerance_px=args.radius_tolerance),geometry_provenance=args.geometry_provenance)
    print(write_central_center_comparison(args.output,image,result))
    (args.output/'input_provenance.json').write_text(json.dumps(dict(
        pattern=str(args.pattern.resolve()),sha256=hashlib.sha256(args.pattern.read_bytes()).hexdigest(),
        shape=list(image.shape),geometry_provenance=args.geometry_provenance,
        preprocessing='Provided prepared pattern; no additional correction/noise processing'),indent=2)+'\n')
    print('Status:',result.status,'intensity:',result.intensity_center_px,'arc:',result.arc_center_px,
          'delta:',result.difference_vector_px,'distance:',result.difference_magnitude_px)


if __name__=='__main__': main()
