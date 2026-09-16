"""Informational CPU timings; run with the qcbed-python interpreter.

No timing assertions. Each row is the median of three calls after one warmup;
preparation is measured separately. No process/thread scheduling is added.
"""

from statistics import median
from time import perf_counter

import numpy as np

from qcbed import (
    AtomSite, Species, MeshConfig, SliceGeometry, AbsorptionCoefficients,
    SliceAppliedParameters, SamplePreparationParameters, prepare_multislice_sample,
    simulate_one_tilt, simulate_tilt_series_serial, simulate_tilt_series_batched,
)
from qcbed.legacy_io import LegacyInput


def main():
    model = LegacyInput("benchmark", 200, False, (0, 0, 1), MeshConfig(16, 16),
        (SliceGeometry(2, .125, -.2, 20, 10, 60), SliceGeometry(3, -.1, .3, 20, 10, 60)),
        (Species("Al", 13),), ((AtomSite(1, .25, .125, .5, "S"),),
                             (AtomSite(1, .375, .25, .75, "S"),)),
        "out.dat", "in.dat", title_raw="benchmark", absorption_raw="YES")
    parameters = SamplePreparationParameters(200,
        (SliceAppliedParameters((.86,), .1, 1.1), SliceAppliedParameters((1.2,), .2, .95)),
        .5, (AbsorptionCoefficients(0, 1, 1),), (0, 1)*10,
        ((0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)), "centered", "fft", "legacy", "unity")
    start = perf_counter()
    prepared = prepare_multislice_sample(model, parameters)
    preparation = perf_counter()-start
    t = np.arange(32, dtype=np.float64)
    tilts = np.column_stack((.02*(t-16), .013*((t*7) % 19-9)))
    def repeated():
        return np.stack([simulate_one_tilt(prepared, tilt_h_undilated_fractional=h,
            tilt_k_undilated_fractional=k).amplitudes for h, k in tilts])
    expected = repeated()
    cases = [("repeated one-tilt", repeated),
             ("serial wrapper", lambda: simulate_tilt_series_serial(prepared, tilts).amplitudes)]
    cases += [(f"batched {size}", lambda size=size: simulate_tilt_series_batched(
        prepared, tilts, batch_size=size).amplitudes) for size in (1, 2, 4, 8, 16, 64)]
    print(f"NumPy {np.__version__}; mesh=16x16; tilts=32; stored types=2; occurrences=20; beams=5")
    print(f"Preparation: {preparation:.6f} s (FFT gratings, absorption enabled, centered storage)")
    print("execution | total seconds | milliseconds/tilt | max amplitude residual | max intensity residual")
    for name, function in cases:
        function()
        timings = []
        for _ in range(3):
            start = perf_counter()
            actual = function()
            timings.append(perf_counter()-start)
        elapsed = median(timings)
        amplitude_error = np.max(abs(actual-expected))
        intensity_error = np.max(abs(abs(actual)**2-abs(expected)**2))
        print(f"{name} | {elapsed:.6f} | {elapsed/len(tilts)*1000:.6f} | {amplitude_error:.3e} | {intensity_error:.3e}")


if __name__ == "__main__":
    main()
