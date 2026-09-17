"""Source-derived coordinate parabolic search (atomicscatfact.o:1-194)."""

from dataclasses import dataclass
from numbers import Integral

import numpy as np

from .multislice import _snapshot
from .refinement import (build_refinement_parameter_layout, apply_refinement_parameters,
                         pack_refinement_parameters, _vector)
from .refinement_trial import evaluate_refinement_trial
from .structure_factors import _finite


@dataclass(frozen=True, eq=False)
class LegacyParaminResult:
    """Retained physical point after the final reset, not an all-time archive.

    The last evaluation is always initialization or a cycle-end reevaluation.
    uncertainties is source UNCERT, not a statistical confidence interval.
    """

    best_vector: np.ndarray
    best_objective: np.float64
    last_vector: np.ndarray
    last_objective: np.float64
    bottom_factors: np.ndarray
    top_factors: np.ndarray
    uncertainties: np.ndarray
    average_displacement: np.float64
    cycles: int
    parameter_steps: int
    evaluations: int
    converged: bool

    def __post_init__(self):
        for name in ("best_vector", "last_vector", "bottom_factors", "top_factors", "uncertainties"):
            object.__setattr__(self, name, _snapshot(getattr(self, name)))


def _fraction(c0, c1, c2, bottom, top, cycle, index):
    """Literal AA/BB expression and IEEE fallback, with one-based index parity."""
    b, t = np.float64(bottom), np.float64(top)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        aa = (((c1-c0)/(b**2))+((c2-c0)/((b**2)*(t-1.0)))-
              ((c2-c0)/(b*(t-1.0))))/(1.0-((t+1.0)/b)+(t/(b**2)))
        bb = ((c2-c0)/(t-1.0))-(aa*(t+1.0))
        fraction = -bb/(2*aa) if aa > 0 else 2.0+(bb/(2*aa))
        valid = fraction/fraction == 1
    fallback = 1+((-1)**(cycle+index+1))*(1-b)/2
    if not valid:
        fraction = fallback
    if fraction > 2:
        fraction = np.float64(2)
    if fraction < 0:
        fraction = fallback
    return fraction


def _order(n, random_order, rng):
    if not random_order:
        return list(range(n))
    order, used = [], set()
    while len(order) < n:
        index = int(rng.random()*n)
        if index not in used:
            used.add(index)
            order.append(index)
    return order


def _run(initial, bottom, top, project, evaluate, *, max_cycles, tolerance,
         random_order=False, rng=None):
    bottom, top = bottom.copy(), top.copy()
    n = len(initial)
    uncertainty = abs(1-bottom)
    retained = initial.copy()
    chi = np.float64(1e20)
    evaluations = 0
    last = None
    last_fit = None

    def visit(vector):
        nonlocal retained, chi, evaluations, last, last_fit
        last = project(vector).copy()
        last_fit = np.float64(evaluate(last.copy()))
        evaluations += 1
        if not np.isfinite(last_fit) or last_fit > 1e20:
            raise ValueError("PARAMIN objective must be finite and <= source CHISMIN sentinel 1e20")
        if last_fit <= chi:  # ONEITPARA shared best update, including ties
            retained, chi = last.copy(), last_fit
        return last.copy(), last_fit

    _, c0 = visit(retained)
    best = c0
    cycles = steps = 0
    average = np.float64(1)
    while average > tolerance and cycles <= max_cycles:
        average = np.float64(0)
        cycles += 1
        trial = retained.copy()
        for m in _order(n, random_order, rng):
            steps += 1
            start = trial[m]
            if start == 0 or bottom[m] == 0 or top[m] == 0:
                raise ValueError("PARAMIN unsafe zero start or BOT/TOP divisor")
            trial[m] = bottom[m]*trial[m]
            trial, c1 = visit(trial)
            trial[m] = top[m]*trial[m]/bottom[m]
            trial, c2 = visit(trial)
            fraction = _fraction(c0, c1, c2, bottom[m], top[m], cycles, m)
            trial[m] = fraction*trial[m]/top[m]
            trial, candidate = visit(trial)
            change = abs(start-retained[m])/abs(start)
            trial[m] = retained[m]
            for score in (candidate, c1, c2):
                if score < best:
                    best = c0 = score
                    uncertainty[m] = abs(1-bottom[m])
            displacement = abs(1-bottom[m])
            average += displacement  # old range, before adaptation
            delta = None
            if change >= displacement/2 and change <= .5:
                delta = change
            if change > .5:
                delta = np.float64(.5)
            if change < displacement/2:
                delta = displacement/2  # separate IF, can override the .5 cap
            if delta is None or not np.isfinite(delta):
                raise ValueError("PARAMIN undefined range adaptation")
            bottom[m], top[m] = 1-delta, 1+delta
        average /= n
        chi = np.float64(1e20)
        _, best = visit(retained)
        c0 = best
    return LegacyParaminResult(retained, chi, last, last_fit, bottom, top,
                              uncertainty, average, cycles, steps, evaluations,
                              bool(average <= tolerance))


def run_legacy_paramin(state, initial_vector, experimental, *, max_cycles,
                       tolerance, random_order=False, seed=None,
                       phase_grating_method="fft") -> LegacyParaminResult:
    """Physical-coordinate PARAMIN with BOT=D2 and TOP=D1 from state.ranges.

    max_cycles is literal MAXFPARA: up to max_cycles+1 cycles can execute.
    Random traversal requires an explicit seed; fixed traversal uses no RNG.
    """
    layout = build_refinement_parameter_layout(state)
    n = len(layout.parameters)
    if not 1 <= n <= 1000:
        raise ValueError("PARAMIN requires 1..1000 active parameters")
    initial = _vector(initial_vector, n)
    if any(p.range_factors is None for p in layout.parameters):
        raise ValueError("PARAMIN requires D1,D2 range factors for every active parameter")
    top = np.array([p.range_factors[0] for p in layout.parameters], dtype=np.float64)
    bottom = np.array([p.range_factors[1] for p in layout.parameters], dtype=np.float64)
    if np.any(top == 0) or np.any(bottom == 0):
        raise ValueError("PARAMIN requires nonzero BOT/TOP divisors")
    if isinstance(max_cycles, (bool, np.bool_)) or not isinstance(max_cycles, Integral) or max_cycles < 0:
        raise ValueError("max_cycles must be a nonnegative integer excluding bool")
    tolerance = np.float64(_finite(tolerance, "tolerance"))
    if tolerance < 0:
        raise ValueError("tolerance must be nonnegative")
    if not isinstance(random_order, (bool, np.bool_)):
        raise TypeError("random_order must be bool")
    if seed is not None and (isinstance(seed, (bool, np.bool_)) or not isinstance(seed, Integral) or seed < 0):
        raise ValueError("seed must be a nonnegative integer excluding bool")
    if random_order and seed is None:
        raise ValueError("random_order requires an explicit seed")
    rng = np.random.Generator(np.random.PCG64(seed)) if random_order else None

    def project(vector):
        return pack_refinement_parameters(apply_refinement_parameters(state, vector))

    def evaluate(vector):
        return evaluate_refinement_trial(state, vector, experimental,
            coordinate_mode="physical", phase_grating_method=phase_grating_method).objective

    with np.errstate(divide="raise", invalid="raise", over="raise"):
        return _run(initial, bottom, top, project, evaluate, max_cycles=int(max_cycles),
                    tolerance=tolerance, random_order=random_order, rng=rng)
