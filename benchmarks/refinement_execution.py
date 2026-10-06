"""Complete real-data trial timings and untimed existing CBED plot diagnostics.

Run with qcbed-python from the repository root. No optimizer or physics cache.
"""

import os
for _name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

from collections import defaultdict
from contextlib import ExitStack
from datetime import datetime
from html import escape
import json
from pathlib import Path
import sys
from time import perf_counter
from unittest.mock import patch

import numpy as np

from profile_thickness_only import thickness_only
from qcbed import ForwardExecutionConfig
from qcbed.experimental import parse_experimental_intensities
from qcbed.legacy_io import _parse_for007_records
from qcbed.refinement_io import _CaptureRecords
from qcbed.refinement_driver import parse_legacy_refinement
from qcbed.refinement import pack_refinement_parameters
import qcbed.refinement_trial as trial
from qcbed.validation import CbedComparison, _simulated
from qcbed.visualization import plot_cbed_comparison


def timed_trial(state, vector, experimental, execution):
    times = defaultdict(float)
    def wrap(name, fn):
        def timed(*args, **kwargs):
            start = perf_counter()
            try:
                return fn(*args, **kwargs)
            finally:
                times[name] += perf_counter() - start
        return timed
    with ExitStack() as stack:
        for name, category in (("build_tilt_propagators", "forward"),
                ("simulate_slice_sequence", "forward"),
                ("simulate_tilt_series_batched", "forward"),
                ("_gradient", "gradient"), ("compute_selected_objective", "ynorm_objective")):
            stack.enter_context(patch.object(trial, name, wrap(category, getattr(trial, name))))
        start = perf_counter()
        result = trial.evaluate_refinement_trial(state, vector, experimental,
            objective_mode="masked", execution=execution)
        total = perf_counter() - start
    return result, dict(total_seconds=total, forward_seconds=times["forward"],
        gradient_seconds=times["gradient"], ynorm_objective_seconds=times["ynorm_objective"],
        other_seconds=total-sum(times.values()))


def preview(destination, label, experimental, result, baseline_images=None):
    """Use the existing four-panel comparison; duplicate trial panels are labeled.

    Verify actual image arrays/masks/scales, not just inputs to the plotter.
    """
    simulated = _simulated(result, experimental)
    comparison = CbedComparison(experimental, simulated, simulated)
    rows = ["<!doctype html><meta charset='utf-8'><title>CBED trial preview</title>",
        "<style>body{font-family:sans-serif;max-width:1400px;margin:auto}img{width:100%}</style>",
        f"<h1>{escape(label)}: evaluated reference trial</h1>",
        "<p>Diagnostic only; no optimization or improvement claimed. Masked fit view: "
        "weight-zero pixels are grey in every panel. The two simulation panels show "
        "the same evaluated trial using the existing comparison plotter. Experimental "
        "and simulated values share one scale; residual = YNORM-scaled simulation minus "
        "experiment, with a symmetric zero-centered scale. No display normalization. "
        "Plotting is excluded from all trial timings.</p>"]
    images = []
    for b, (h, k) in enumerate(experimental.beam_indices_hk):
        figure = plot_cbed_comparison(comparison, beam_index=b, mask_zero_weight=True)
        figure.axes[1].set_title("Evaluated trial (YNORM-scaled)")
        figure.axes[2].set_title("Same evaluated trial")
        figure.axes[3].set_title("Trial minus experimental")
        plotted = [ax.images[0] for ax in figure.axes[:4]]
        expected_mask = (experimental.weights[:, b] == 0).reshape(experimental.scan_rows, experimental.scan_columns)
        for image in plotted:
            np.testing.assert_array_equal(np.ma.getmaskarray(image.get_array()), expected_mask)
        assert plotted[0].get_clim() == plotted[1].get_clim() == plotted[2].get_clim()
        low, high = plotted[3].get_clim()
        assert low == -high
        values = tuple((im.get_array().copy(), im.get_clim()) for im in plotted)
        if baseline_images is not None:
            for (actual, limits), (expected, expected_limits) in zip(values, baseline_images[b]):
                np.testing.assert_array_equal(actual.data, expected.data)
                np.testing.assert_array_equal(actual.mask, expected.mask)
                assert limits == expected_limits
        images.append(values)
        filename = f"{label}_beam_{b+1:02d}.png"
        figure.savefig(destination / filename, dpi=90)
        figure.clear()
        rows.append(f"<h2>Beam H={h}, K={k}</h2><img loading='lazy' src='{filename}' alt='Beam H={h}, K={k}'>")
    (destination / f"{label}.html").write_text("\n".join(rows), encoding="utf-8")
    return images


def main():
    root = Path(__file__).resolve().parents[1]
    destination = root / "tmp" / "refinement_execution" / datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    destination.mkdir(parents=True, exist_ok=False)
    source = root / "tests/fixtures/legacy_reference" / "FOR007.DAT"
    original, _ = parse_legacy_refinement(source, seed=2026)
    state, _ = thickness_only(original)
    baseline_state = repr(state)
    records = _CaptureRecords(source)
    _parse_for007_records(records)
    row = records.values["output beam count and indices"]
    beams = tuple(zip(row[1::2], row[2::2]))
    experimental = parse_experimental_intensities(source.parent / "in.dat", beam_indices_hk=beams, scan=state.scan)
    vector = pack_refinement_parameters(state)
    results = {}
    report = dict(interpreter=sys.executable, numpy=np.__version__, objective_mode="masked",
        depth_angstrom=float(vector[0]), gradient_option=experimental.gradient_option,
        output_directory=str(destination), cases={})
    for label, execution in (("serial", ForwardExecutionConfig()),
                            ("batched32", ForwardExecutionConfig("batched", 32)),
                            ("batched64", ForwardExecutionConfig("batched", 64))):
        print(f"Starting {label}", flush=True)
        result, metrics = timed_trial(state, vector, experimental, execution)
        results[label] = result
        reference = results["serial"]
        for name in ("amplitudes", "calculated_intensities", "objective_intensities"):
            a, b = getattr(result, name), getattr(reference, name)
            np.testing.assert_array_equal(a, b)
            metrics[f"max_{name}_residual"] = float(np.max(np.abs(a-b)))
        assert result.objective == reference.objective and result.scaling == reference.scaling
        assert result.physical_state == reference.physical_state and repr(state) == baseline_state
        metrics.update(objective=float(result.objective), ynorm=float(result.scaling),
            objective_difference=float(result.objective-reference.objective),
            ynorm_difference=float(result.scaling-reference.scaling))
        report["cases"][label] = metrics
        metrics["speedup"] = report["cases"]["serial"]["total_seconds"] / metrics["total_seconds"]
        (destination / "results.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
        np.savez_compressed(destination / f"{label}.npz", amplitudes=result.amplitudes,
            raw=result.calculated_intensities, objective_intensities=result.objective_intensities,
            objective=result.objective, ynorm=result.scaling)
        print(json.dumps(metrics), flush=True)
    # All timing scopes have ended before validation/visualization assembly.
    images = preview(destination, "serial", experimental, results["serial"])
    preview(destination, "batched64", experimental, results["batched64"], images)
    report["plotted_arrays_masks_and_scales_identical"] = True
    (destination / "results.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    (destination / "index.html").write_text(
        "<!doctype html><meta charset='utf-8'><h1>Reference trial backend comparison</h1>"
        "<p>All nine beams; diagnostic trials, not an optimization. Serial/batched plotted arrays, "
        "masks and color scales are exactly identical. Plotting excluded from timings.</p>"
        "<ul><li><a href='serial.html'>Serial reference preview</a></li>"
        "<li><a href='batched64.html'>Batched 64 preview</a></li>"
        "<li><a href='results.json'>Timings and residuals</a></li></ul>", encoding="utf-8")
    print(f"Preview: {destination / 'index.html'}", flush=True)


if __name__ == "__main__":
    main()
