"""Live Patch A orchestration; no archived geometry, catalogues or origin inputs.

Calibration and geometry remain separate backend operations; a user-supplied
aperture is also thresholded for disk geometry unless an override is supplied. Arrays use (row,column); selections use detector (x,y) pixels.
"""
from dataclasses import dataclass, asdict
import numpy as np
from .models import DiskGeometry, SignalMode
from .noise import _pattern, characterize_noise, sigma_map_from_noise_model
from .variance import variance_map_from_sigma_map
from .disk_detection import estimate_disk_geometry
from .central_arc import locate_central_disk
from .lattice_detection import detect_lattice_candidates, infer_lattice_seed, _assign, _config
from itertools import combinations
from .lattice_bootstrap import bootstrap_lattice
from .arc_catalogue import generate_arc_catalogue
from .lattice_refinement import refine_lattice_and_generate_catalogue
from .origin_selection import SelectableBeam, nearest_beam, prepare_roi
from .origin import measure_symmetry_origin
from .extraction import ExtractionConfig, WeightMode, prepare_experimental_qcbed
from .scan_conversion import scan_from_catalogue, assemble_prepared_dataset
from ..detector.psf import calibrate_psf_from_aperture
from ..detector.correction import correct_cbed_with_otf


@dataclass(frozen=True)
class LiveSelection:
    geometry: DiskGeometry | None
    symmetry_roi: tuple[float,float,float,float]
    selection_width: int
    selection_height: int
    final_binning: int = 1
    signal_mode: SignalMode = SignalMode.CONVENTIONAL
    symmetry_crop_half_size_px: tuple[int,int] = (8,8)
    psf_seed: int = 0
    raw_symmetry_roi: tuple[float,float,float,float] | None = None
    diameter_endpoints_px: tuple[tuple[float,float],tuple[float,float]] | None = None


@dataclass(frozen=True)
class LivePreparation:
    dataset: object
    corrected_pattern: np.ndarray
    sigma: np.ndarray
    variance: np.ndarray
    validity: np.ndarray
    geometry: DiskGeometry
    symmetry: object
    provenance: dict


def _live_arc_seed(candidates, center, geometry):
    """Reuse the validated two-arc seed for sparse, mutually consistent evidence.

    Geometry is supplied upstream and never fitted here. Every sparse candidate
    must agree with the selected basis; no intensity fallback or threshold change.
    Dense observations retain the existing lattice-fitting path.
    """
    config = _config(None)
    candidates = tuple(sorted(candidates, key=lambda c: (c.center_px, c.id)))
    if not 2 < len(candidates) < config.min_inliers:
        return infer_lattice_seed(candidates, center, geometry, config=config)
    relative = np.array([c.center_px for c in candidates]) - center
    choices = []
    for pair in combinations(candidates, 2):
        try:
            seed = infer_lattice_seed(pair, center, geometry, config=config)
        except ValueError:
            continue
        matrix = np.column_stack((seed.basis_vector_1_px, seed.basis_vector_2_px))
        _, residual, inliers = _assign(matrix, relative, config.lattice_residual_tolerance_px)
        if not np.all(inliers):
            continue
        # Geometric confidence only; deterministic ties follow physical centers.
        choices.append((sum(c.confidence for c in pair), -float(sum(residual**2)), seed))
    if not choices:
        raise ValueError('Sparse arc centers do not support a consistent non-collinear lattice seed')
    return max(choices, key=lambda item: item[:2])[2]


def resolve_geometry(*, geometry=None, disk_profile=None, aperture=None):
    """Explicit geometry > optional profile override > aperture threshold.

    Manual drawing is an explicit geometry override, never inferred from the
    array dimensions. No geometry source means interactive input is required.
    """
    if geometry is not None:
        if not isinstance(geometry,DiskGeometry):raise TypeError('Expected DiskGeometry')
        return geometry,'explicit/manual geometry override'
    if disk_profile is not None:
        return estimate_disk_geometry(disk_profile).geometry,'optional diffraction disk-profile override'
    if aperture is not None:
        return estimate_disk_geometry(aperture).geometry,'aperture threshold geometry'
    raise ValueError('No aperture image supplied. Draw and confirm a disk diameter before continuing.')


def prepare_signal(pattern, *, calibration_aperture=None, disk_profile=None,
                   geometry=None, psf_seed=0, psf_config=None, noise_config=None):
    """Explicit geometry precedence; no resampling or fallback after calibration failure."""
    image=_pattern(pattern)
    geometry,geometry_source=resolve_geometry(geometry=geometry,disk_profile=disk_profile,aperture=calibration_aperture)
    if calibration_aperture is None:
        variance=np.maximum(image,0.)
        return image.copy(),np.sqrt(variance),variance,image>0,geometry,dict(
            noise_source='Poisson / intensity fallback',geometry_source=geometry_source,
            correction='none',uncertainty='variance=max(I,0); sigma=sqrt(max(I,0)); positive intensity support')
    # Existing OTF correction requires the same array grid, not the same optical
    # magnification or diffraction-disk size. Geometry above is a separate threshold
    # measurement, never inferred from the frame dimensions.
    if np.shape(calibration_aperture)!=image.shape:
        raise ValueError('Detector calibration and CBED array grids must match for existing OTF correction; automatic resampling is unsupported.')
    response=calibrate_psf_from_aperture(calibration_aperture,config=psf_config,seed=psf_seed)
    corrected=correct_cbed_with_otf(image,response)
    fit=characterize_noise(corrected,config=noise_config)
    sigma=sigma_map_from_noise_model(corrected,fit.model)
    return corrected,sigma,variance_map_from_sigma_map(sigma),fit.noise_mask,geometry,dict(
        noise_source='calibrated',geometry_source=geometry_source,psf_seed=psf_seed,
        correction='calibrate_psf_from_aperture -> correct_cbed_with_otf',
        detector_provenance=list(response.provenance),noise_model=asdict(fit.model),noise_config=asdict(fit.config),
        uncertainty='A4 sigma squared; final A4 noise mask')


def prepare_from_raw_cbed(pattern, selection, *, calibration_aperture=None, disk_profile=None,
                          progress=lambda message:None, psf_config=None, noise_config=None):
    """Call unchanged scientific functions with current live inputs only."""
    if not isinstance(selection,LiveSelection):raise TypeError('Expected LiveSelection')
    if selection.symmetry_roi is None:raise ValueError('Centrosymmetric feature has not been selected.')
    progress('Signal correction and uncertainty')
    image,sigma,variance,valid,geometry,meta=prepare_signal(pattern,calibration_aperture=calibration_aperture,
        disk_profile=disk_profile,geometry=selection.geometry,psf_seed=selection.psf_seed,
        psf_config=psf_config,noise_config=noise_config)
    cfg=ExtractionConfig(selection.selection_width,selection.selection_height,final_binning=selection.final_binning,
        physical_selection_shape_px=(selection.selection_height,selection.selection_width),
        sampling_policy='nearest_detector_grid')
    progress('Central disk and reciprocal-lattice detection')
    central=locate_central_disk(image,geometry)
    candidates=detect_lattice_candidates(image,geometry,central.central_disk_center_px)
    seed=_live_arc_seed(candidates,central.central_disk_center_px,geometry)
    meta['lattice_detection'] = dict(candidate_count=len(candidates),
        validated_arc_centers=[c.center_px for c in candidates],
        seed=asdict(seed), geometry=asdict(geometry))
    if seed.seed_method=='two_validated_arc_centers':
        bootstrap=bootstrap_lattice(image,seed,geometry)
        catalogue=generate_arc_catalogue(bootstrap,image.shape,geometry)
        meta['lattice_status']=bootstrap.message
    else:
        full=refine_lattice_and_generate_catalogue(image,seed,geometry)
        from .arc_catalogue import ArcBeamCatalogue
        from .lattice_bootstrap import BootstrapSite
        catalogue=ArcBeamCatalogue(tuple(BootstrapSite(b.lattice_indices,b.predicted_center_px,b.observed_center_px,None,
            'reference' if b.is_central else b.status.value,'live refined lattice',beam_id=b.observation.id) for b in full.beams),
            'live refined lattice',True)
        meta['lattice_status']='refined reciprocal lattice'
    sites=catalogue.sites
    meta['lattice_sites'] = [asdict(s) for s in sites]
    def site(index):
        matches=[s for s in sites if s.integer_indices==index]
        if len(matches)!=1:raise ValueError(f'Missing unique reciprocal basis beam {index}')
        return matches[0]
    center=site((0,0)).predicted_center_px
    progress('Refine Symmetry Center')
    x0,y0,x1,y1=selection.symmetry_roi
    beam,_=nearest_beam(tuple(SelectableBeam(s.beam_id,s.integer_indices,s.predicted_center_px,s.status) for s in sites),
        ((x0+x1)/2,(y0+y1)/2),selection_radius_px=max(geometry.radius_x_px,geometry.radius_y_px))
    reference,search=prepare_roi(selection.symmetry_roi,beam,geometry,image.shape,
        crop_half_size_px=selection.symmetry_crop_half_size_px)
    symmetry=measure_symmetry_origin(image,geometry,center,reference,config=search)
    if not symmetry.crop_has_contrast:raise ValueError('Selected symmetry crop has no contrast; reselect a feature.')
    progress('A9 extraction and A10 scan geometry')
    selected=tuple(s.beam_id for s in sites if s.status in ('reference','observed','defining','confirmed','weak','weak_evidence','missing','edge_clipped'))
    calibrated=calibration_aperture is not None
    validity_input={'noise_valid_mask':valid} if calibrated else {'detector_valid_mask':valid}
    extraction=prepare_experimental_qcbed(image,catalogue,geometry,selected_beam_ids=selected,config=cfg,
        signal_mode=selection.signal_mode,weight_mode=WeightMode.NOISE_CHARACTERIZATION if calibrated else WeightMode.BEAM_MASK,
        variance_map=variance,**validity_input,uncertainty_provenance=meta['noise_source']+': '+meta['uncertainty'])
    meta['extraction'] = [dict(indices=b.integer_indices, status=b.status, coverage=b.coverage_fraction,
        physical_center=b.predicted_center_px, sampled_bounds=b.bounds_xyxy,
        sampling_offset_px=(b.bounds_xyxy[0]-(b.predicted_center_px[0]-cfg.selection_width/2),
                            b.bounds_xyxy[1]-(b.predicted_center_px[1]-cfg.selection_height/2))) for b in extraction.beams]
    # Diagnostic layers only; extraction/objective arrays are not modified.
    from .extraction import _integrate
    layers={}
    for beam in extraction.beams:
        h,w=beam.roi_geometry_mask.shape
        yy,xx=np.indices((h,w))
        radial=(xx-w/2)**2+(yy-h/2)**2 <= (w/2)**2
        geometric=_integrate((beam.roi_geometry_mask & radial &
            beam.detector_in_frame_mask).astype(float),cfg,threshold=True)>0
        final=beam.final_weight>0
        if np.any(final & ~geometric):
            raise ValueError('Diagnostic geometric layer does not contain final support')
        # Processing eligibility is evaluated within geometric support; outside
        # that domain it is neutral, so geometric & processing equals base weight.
        processing=final | ~geometric
        layers[beam.beam.id]=dict(geometric=geometric.astype(int).tolist(),
            processing=processing.astype(int).tolist())
    meta['mask_layers']=layers
    meta['sampling_policy']='nearest detector grid (ties to even); pixel registration approximation; physical geometry and scan basis unchanged'
    conversion=scan_from_catalogue(catalogue,central_beam_id=site((0,0)).beam_id,
        canonical_basis_beam_ids=(site((1,0)).beam_id,site((0,1)).beam_id),
        zone_axis_origin_px=symmetry.zone_axis_origin_px,extraction_config=cfg,basis_transform=((1,0),(0,1)))
    dataset=assemble_prepared_dataset(extraction,conversion,intensity_units='input intensity units',physical_geometry=geometry,
        provenance=(meta['noise_source'],'Live catalogue; identity scan basis; physical HK belongs to B'))
    meta.update(selection=asdict(selection),symmetry=asdict(symmetry),central_localization=asdict(central),
                beam_selection='All non-rejected sites; A9 frame/support rejection retained')
    return LivePreparation(dataset,image,sigma,variance,valid,geometry,symmetry,meta)
