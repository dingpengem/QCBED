"""Compare CURRENT Patch A geometry with recorded DM coordinates; no input fitting.

Reads existing diagnostics only, without rerunning or modifying A5-A9.
"""
import argparse
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import numpy as np
from qcbed.data_preparation import ExtractionConfig
from qcbed.data_preparation.extraction_diagnostics import load_a7_extraction_catalogue
from qcbed.data_preparation.scan_conversion import scan_from_catalogue

# Recorded output values only. Never used to construct any geometric input.
DM_COMPARISON = ((-.508156,.407641),(.407641,.508156),(-.407641,-.508156))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--a8-measurement',type=Path,required=True)
    parser.add_argument('--a9-summary',type=Path,required=True)
    parser.add_argument('--catalogue',type=Path,required=True)
    parser.add_argument('--central-beam-id',required=True)
    parser.add_argument('--canonical-basis-ids',nargs=2,required=True)
    parser.add_argument('--basis-transform',nargs=4,type=int,required=True,help='Row-major 2x2 integer matrix; scan columns = canonical columns @ transform')
    parser.add_argument('--binning',type=int,required=True,help='Explicit confirmed coordinate-conversion binning; A9 arrays remain unchanged')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    measurement=json.loads(args.a8_measurement.read_text())
    summary=json.loads(args.a9_summary.read_text())
    catalogue=load_a7_extraction_catalogue(args.catalogue)
    digest=hashlib.sha256(args.catalogue.read_bytes()).hexdigest()
    if summary['provenance']['catalogue_sha256']!=digest:
        raise ValueError('A9 catalogue identity does not match supplied catalogue')
    saved_config=ExtractionConfig(**summary['config'])
    if args.binning!=saved_config.final_binning:
        raise ValueError('A10 binning must match the physical A9 extraction; rerun A9 with the confirmed setting')
    config=saved_config
    result=scan_from_catalogue(catalogue,central_beam_id=args.central_beam_id,
        canonical_basis_beam_ids=args.canonical_basis_ids,
        zone_axis_origin_px=measurement['zone_axis_origin_px'],extraction_config=config,
        basis_transform=np.array(args.basis_transform).reshape(2,2))
    if tuple(measurement['central_disk_center_px'])!=result.reference_point_px:
        raise ValueError('A8 central geometry differs from selected A7 central beam')
    corners=np.array([result.scan_geometry.corner_a,result.scan_geometry.corner_b,result.scan_geometry.corner_c])
    differences=corners-np.array(DM_COMPARISON)
    report=dict(**asdict(result),handedness=result.handedness,
        validation='Current reconstructed geometry versus historical output values; not matched-input regression',
        historical_dm_corners=DM_COMPARISON,signed_differences=differences.tolist(),
        absolute_differences=np.abs(differences).tolist(),maximum_absolute_difference=float(np.abs(differences).max()),
        saved_a9_binning=saved_config.final_binning,confirmed_conversion_binning=args.binning,
        a9_array_geometry_matches=saved_config.final_binning==args.binning,
        inputs={k:dict(path=str(getattr(args,k).resolve()),sha256=hashlib.sha256(getattr(args,k).read_bytes()).hexdigest())
                for k in ('a8_measurement','a9_summary','catalogue')},
        retained_beam_ids=[b['beam_id'] for b in summary['beams'] if b['retained_by_a9']],
        uncertainty=summary['uncertainty'])
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'scan_geometry.json').write_text(json.dumps(report,indent=2)+'\n')
    text='# A10 current-geometry comparison\n\n```json\n'+json.dumps(report,indent=2)+'\n```\n'
    (args.output/'scan_geometry.md').write_text(text)
    print(json.dumps(report,indent=2))


if __name__=='__main__': main()
