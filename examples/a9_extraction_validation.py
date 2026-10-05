"""A9 diagnostic extraction with explicit selection, dimensions and uncertainty semantics."""
import argparse
import hashlib
from pathlib import Path
import numpy as np
from qcbed.data_preparation import (DiskGeometry, ExtractionConfig, SignalMode,
    WeightMode, prepare_experimental_qcbed)
from qcbed.data_preparation.diagnostics import load_pattern
from qcbed.data_preparation.extraction_diagnostics import (
    load_a7_extraction_catalogue, write_extraction_diagnostics)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pattern',type=Path,required=True)
    p.add_argument('--catalogue',type=Path,required=True)
    selection=p.add_mutually_exclusive_group(required=True)
    selection.add_argument('--beam-ids',nargs='+')
    selection.add_argument('--all-lattice-sites',action='store_true',help='Explicitly consider all non-rejected A7 predictions, including unconfirmed sites')
    p.add_argument('--minimum-beam-frame-coverage',type=float,default=.90)
    p.add_argument('--radii',nargs=2,type=float,required=True)
    size=p.add_mutually_exclusive_group(required=True)
    size.add_argument('--selection',nargs=2,type=int,metavar=('WIDTH','HEIGHT'),help='Diagnostic-only numerical crop; not physical scan geometry')
    size.add_argument('--physical-selection-from-disk',action='store_true',help='Explicit physical-diameter policy; requires integral diameters')
    p.add_argument('--diagnostic-window',nargs=2,type=int,metavar=('WIDTH','HEIGHT'),help='Optional presentation metadata only; never controls extraction or A10')
    p.add_argument('--binning',type=int,required=True)
    p.add_argument('--expansion',type=int,required=True)
    p.add_argument('--signal-mode',choices=[m.value for m in SignalMode],required=True)
    p.add_argument('--weight-mode',choices=[m.value for m in WeightMode],default='beam_mask')
    p.add_argument('--uncertainty',type=Path)
    p.add_argument('--uncertainty-kind',choices=('sigma','variance'))
    p.add_argument('--noise-valid-mask',type=Path,help='Validated A4 final Boolean .npy mask')
    p.add_argument('--uncertainty-provenance')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    noise={}
    if a.uncertainty is not None:
        if not a.uncertainty_kind or not a.uncertainty_provenance:
            p.error('Uncertainty requires explicit --uncertainty-kind and --uncertainty-provenance')
        noise[a.uncertainty_kind+'_map']=load_pattern(a.uncertainty)
        noise['uncertainty_provenance']=a.uncertainty_provenance
    elif a.uncertainty_kind or a.uncertainty_provenance:
        p.error('Uncertainty declaration requires --uncertainty')
    if a.noise_valid_mask:
        noise['noise_valid_mask']=np.load(a.noise_valid_mask,allow_pickle=False)
    image=load_pattern(a.pattern)
    catalogue=load_a7_extraction_catalogue(a.catalogue)
    geometry=DiskGeometry(*a.radii)
    ids=[s.beam_id for s in catalogue.sites if s.status != 'rejected'] if a.all_lattice_sites else a.beam_ids
    settings=dict(final_binning=a.binning,expansion_factor=a.expansion,
        minimum_beam_frame_coverage=a.minimum_beam_frame_coverage,
        diagnostic_window_shape_px=None if a.diagnostic_window is None else tuple(reversed(a.diagnostic_window)))
    config=(ExtractionConfig.from_disk_geometry(geometry,**settings) if a.physical_selection_from_disk
            else ExtractionConfig(*a.selection,**settings))
    result=prepare_experimental_qcbed(image,catalogue,geometry,selected_beam_ids=ids,
        config=config,signal_mode=SignalMode(a.signal_mode),weight_mode=WeightMode(a.weight_mode),**noise)
    provenance=dict(pattern=str(a.pattern.resolve()),pattern_sha256=hashlib.sha256(a.pattern.read_bytes()).hexdigest(),
        catalogue=str(a.catalogue.resolve()),catalogue_sha256=hashlib.sha256(a.catalogue.read_bytes()).hexdigest(),
        selected_ids=ids,shape=image.shape,
        preprocessing='Pattern used as supplied; no A3 correction, smoothing or A4 provenance inferred',
        geometry='Manually supplied physical radii; canonical A7 predicted centers unchanged')
    summary=write_extraction_diagnostics(a.output,image,result,geometry,provenance=provenance)
    print(f"Selected {summary['selected_beam_count']}; extracted {summary['extracted_beam_count']}; rejected {summary['rejected_for_coverage']}")
    for b in summary['beams']:
        print(b['beam_id'],b['signal_shape'],b['coverage_fraction'],b['signal_mean'],b['final_nonzero_weight_count'])


if __name__=='__main__': main()
