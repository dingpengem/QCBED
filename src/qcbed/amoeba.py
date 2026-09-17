"""Float64 translation of readable AMOEBA in legacy_fortran/simpara.o."""

from dataclasses import dataclass
from numbers import Integral

import numpy as np

from .multislice import _snapshot
from .refinement import build_refinement_parameter_layout, _scal, _vector
from .refinement_trial import evaluate_refinement_trial
from .structure_factors import _finite


@dataclass(frozen=True, eq=False)
class LegacyAmoebaResult:
    """Simplex coordinates; best_index identifies the unswapped final best row.

    evaluations includes initialization; search_evaluations is source NUMF.
    iterations counts completed reflection/expansion/contraction/shrink cycles.
    """

    final_simplex: np.ndarray
    final_objectives: np.ndarray
    best_index: int
    evaluations: int
    search_evaluations: int
    iterations: int
    stop_reason: str
    best_evaluated_vector: np.ndarray
    best_evaluated_objective: np.float64

    def __post_init__(self):
        object.__setattr__(self, "final_simplex", _snapshot(self.final_simplex))
        object.__setattr__(self, "final_objectives", _snapshot(self.final_objectives))
        object.__setattr__(self, "best_evaluated_vector", _snapshot(self.best_evaluated_vector))

    @property
    def best_vector(self):
        return self.final_simplex[self.best_index]

    @property
    def best_objective(self):
        return self.final_objectives[self.best_index]

    @property
    def converged(self):
        return self.stop_reason != "max_evaluations"


def _layout(state):
    layout = build_refinement_parameter_layout(state)
    if not 1 <= len(layout.parameters) <= 1000:
        raise ValueError("AMOEBA requires 1..1000 active parameters")
    return layout


def build_legacy_simplex(state, initial_vector, plamda) -> np.ndarray:
    """Driver initialization from physical initial_vector and simplex PLAMDA.

    SCAL stays anchored to state. Rows 1..N displace coordinates 0..N-1.
    Zero or signed PLAMDA is preserved; D1/D2 are not clipping bounds.
    """
    layout = _layout(state)
    n = len(layout.parameters)
    initial = _vector(initial_vector, n)
    steps = _vector(plamda, n)
    a, b = _scal(layout)
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        p = np.tile((initial-b)*a, (n+1, 1))
        for j in range(n):
            p[j+1, j] += steps[j]
    return _snapshot(p)


def _rank(y):
    lo = 0
    hi, nhi = (0, 1) if y[0] > y[1] else (1, 0)
    for i in range(len(y)):
        if y[i] < y[lo]:
            lo = i
        if y[i] > y[hi]:
            nhi, hi = hi, i
        elif y[i] > y[nhi] and i != hi:
            nhi = i
    return lo, hi, nhi


def _run(simplex, evaluate, *, ftol, max_evaluations, fmin):
    """Private scalar-objective seam; source loops retain summation order."""
    p = simplex.copy()
    n = p.shape[1]
    retained = None
    retained_value = np.float64(np.inf)

    def visit(point):
        nonlocal retained, retained_value
        value = np.float64(evaluate(point.copy()))
        if not np.isfinite(value):
            raise ValueError("AMOEBA objective must be finite")
        if value <= retained_value:
            retained, retained_value = point.copy(), value
        return value

    y = np.array([visit(row) for row in p], dtype=np.float64)
    numf = iterations = 0
    while True:
        lo, hi, nhi = _rank(y)
        with np.errstate(divide="ignore", invalid="ignore", over="raise"):
            rtol = 2.0*abs(y[hi]-y[lo])/(abs(y[hi])+abs(y[lo]))
        if rtol < ftol:
            reason = "objective_spread"
            break
        if numf > max_evaluations:
            reason = "max_evaluations"
            break
        if 0.5*(y[hi]+y[lo]) <= fmin:
            reason = "objective_threshold"
            break
        center = np.zeros(n, dtype=np.float64)
        for i in range(n+1):
            if i != hi:
                for j in range(n):
                    center[j] += p[i, j]
        center /= n
        reflected = 2.0*center-p[hi]
        yr = visit(reflected)
        numf += 1
        if yr <= y[lo]:
            expanded = 2.0*reflected-center
            ye = visit(expanded)
            numf += 1
            # Deliberately compare expansion to OLD BEST, not reflection.
            if ye < y[lo]:
                p[hi], y[hi] = expanded, ye
            else:
                p[hi], y[hi] = reflected, yr
        elif yr >= y[nhi]:
            if yr < y[hi]:
                p[hi], y[hi] = reflected, yr
            contracted = 0.5*p[hi]+0.5*center
            yc = visit(contracted)
            numf += 1
            if yc < y[hi]:
                p[hi], y[hi] = contracted, yc
            else:
                for i in range(n+1):
                    if i != lo:
                        p[i] = 0.5*(p[i]+p[lo])
                        y[i] = visit(p[i])
                        numf += 1
        else:
            p[hi], y[hi] = reflected, yr
        iterations += 1
    return LegacyAmoebaResult(p, y, lo, numf+n+1, numf, iterations, reason,
                              retained, retained_value)


def run_legacy_amoeba(state, initial_simplex, experimental, *, ftol,
                      max_evaluations, fmin, phase_grating_method="fft", objective_mode="legacy",
                      evaluation_callback=None, execution=None, before_evaluation=None) -> LegacyAmoebaResult:
    """Run one source AMOEBA stage; initial_simplex is in legacy_simplex coordinates.

    max_evaluations is the soft source MXF limit excluding initial vertices.
    No driver restart or implicit final reevaluation is performed.
    """
    layout = _layout(state)
    _scal(layout)  # Validate the baseline conversion even on immediate return.
    n = len(layout.parameters)
    if not isinstance(initial_simplex, np.ndarray) or initial_simplex.dtype != np.float64:
        raise TypeError("initial_simplex must be a float64 NumPy array")
    if initial_simplex.shape != (n+1, n) or not np.all(np.isfinite(initial_simplex)):
        raise ValueError(f"initial_simplex must be finite with shape ({n+1}, {n})")
    ftol = _finite(ftol, "ftol")
    fmin = _finite(fmin, "fmin")
    if ftol < 0:
        raise ValueError("ftol must be nonnegative")
    if (isinstance(max_evaluations, (bool, np.bool_)) or
            not isinstance(max_evaluations, Integral) or max_evaluations < 0):
        raise ValueError("max_evaluations must be a nonnegative integer excluding bool")

    def evaluate(vector):
        return evaluate_refinement_trial(state, vector, experimental,
            coordinate_mode="legacy_simplex", phase_grating_method=phase_grating_method,
            objective_mode=objective_mode, execution=execution, before_evaluation=before_evaluation, evaluation_callback=evaluation_callback).objective

    with np.errstate(over="raise", invalid="raise", divide="raise"):
        return _run(initial_simplex, evaluate, ftol=ftol,
                    max_evaluations=int(max_evaluations), fmin=fmin)
