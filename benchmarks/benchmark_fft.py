"""Informational 32x32 direct/FFT timing, separate from pytest.

Run with the project's qcbed-python interpreter. Includes validation and
reordering overhead; no timing assertion or legacy equivalence claim.
"""

from statistics import median
from time import perf_counter

import numpy as np

from qcbed.reciprocal import (DirectLattice2D, ReciprocalLattice2D,
                              ReciprocalIndexGrid, ReciprocalVectorGrid)
from qcbed.phase_grating import projected_potential_direct, projected_potential_fft


def main():
    g = ReciprocalVectorGrid(ReciprocalLattice2D(DirectLattice2D(20, 20, 90)),
                             ReciprocalIndexGrid(32, 32))
    rng = np.random.default_rng(20260916)
    c = (rng.normal(size=(32, 32))+1j*rng.normal(size=(32, 32))).astype(np.complex128)
    print(f"NumPy {np.__version__}; centered 32x32 complex128; median of 5 warm calls")
    results = []
    for function in (projected_potential_direct, projected_potential_fft):
        function(c, g)
        durations = []
        for _ in range(5):
            start = perf_counter()
            out = function(c, g)
            durations.append(perf_counter()-start)
        results.append(out)
        print(f"{function.__name__}: {median(durations):.9f} s")
    error = np.abs(results[0]-results[1])
    print(f"maximum absolute residual: {error.max():.9g}")
    print(f"relative L2 residual: {np.linalg.norm(results[0]-results[1])/np.linalg.norm(results[0]):.9g}")


if __name__ == "__main__":
    main()
