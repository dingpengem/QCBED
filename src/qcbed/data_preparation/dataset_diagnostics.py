"""Read-only JSON-compatible metadata for the A1 prepared dataset, not a new model."""
from dataclasses import asdict
import hashlib


def prepared_dataset_manifest(dataset):
    """Describe stable ordering, immutable numerical arrays and preparation settings.

    Large arrays remain ordinary NPY artifacts, as in the existing A9 diagnostics.
    Hashes cover contiguous float64 bytes in public (tilt,beam) C order.
    Physical HK remains null unless explicitly assigned in BeamObservation.
    """
    arrays = {}
    for name in ('intensities', 'variances', 'weights'):
        value = getattr(dataset, name)
        arrays[name] = None if value is None else dict(shape=list(value.shape),
            dtype=str(value.dtype), sha256=hashlib.sha256(value.tobytes(order='C')).hexdigest())
    mapping = []
    for i, beam in enumerate(dataset.beams):
        mapping.append(dict(stable_beam_id=beam.id,
            lattice_coordinates=list(dataset.lattice_coordinates[i]) if dataset.lattice_coordinates else None,
            physical_indices=None if beam.indices is None else list(beam.indices),
            original_extraction_block=dataset.extraction_block_indices[i] if dataset.extraction_block_indices else None,
            retained_extraction_block=i, refinement_column=i,
            center_px=list(beam.center_px),
            center_source=dataset.center_sources[i] if dataset.center_sources else None,
            nonzero_weight_pixels=int((dataset.weights[:, i] != 0).sum())))
    cfg = dataset.extraction_config
    return dict(arrays=arrays, beam_mapping=mapping,
        array_order='tilt=row*scan_columns+column; beam=ordered retained extraction block',
        beam_major_flattening='array.T.ravel(order="C")',
        signal_mode=dataset.mode.value, intensity_units=dataset.intensity_units,
        weight_mode=None if dataset.weight_mode is None else dataset.weight_mode.value,
        physical_disk_geometry=None if dataset.physical_disk_geometry is None else asdict(dataset.physical_disk_geometry),
        extraction_config=None if cfg is None else asdict(cfg),
        raw_pixels_per_beam=None if cfg is None else cfg.selection_width*cfg.selection_height,
        final_pixels_per_beam=dataset.intensities.shape[0],
        flattened_signal_size=dataset.intensities.size,
        flattened_weight_size=dataset.weights.size,
        flattened_variance_size=None if dataset.variances is None else dataset.variances.size,
        rois=[asdict(r) for r in dataset.rois], reciprocal_frame=asdict(dataset.reciprocal_frame),
        scan_geometry=asdict(dataset.scan_geometry), provenance=list(dataset.provenance))
