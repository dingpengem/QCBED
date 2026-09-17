"""Source-derived anneal.f control flow with explicit float64 and seeded RNG."""

from dataclasses import dataclass
from numbers import Integral

import numpy as np

from .multislice import _snapshot
from .refinement import (build_refinement_parameter_layout, apply_refinement_parameters,
                         pack_refinement_parameters, _vector)
from .refinement_trial import evaluate_refinement_trial
from .structure_factors import _finite


@dataclass(frozen=True, eq=False)
class LegacyAnnealResult:
    """Final PARMIN/CHISMIN, local BESTFIT and last evaluated point.

    best means retained after the source's most recent reset, not a separate
    all-time archive. last_vector can be a rejected proposal or a reevaluation.
    """
    best_vector: np.ndarray
    best_objective: np.float64
    last_vector: np.ndarray
    last_objective: np.float64
    adaptation_objective: np.float64
    iterations: int
    evaluations: int
    temperature: np.float64
    amplitudes: np.ndarray

    def __post_init__(self):
        for name in ("best_vector", "last_vector", "amplitudes"):
            object.__setattr__(self, name, _snapshot(getattr(self, name)))
        for name in ("best_objective", "last_objective", "adaptation_objective", "temperature"):
            object.__setattr__(self, name, np.float64(getattr(self, name)))


def _normal(rng, center, width):
    """oneit.o:26-50; signed SIG is intentional, second normal is discarded."""
    while True:
        first = np.float64(2)*rng.random()-np.float64(1)
        second = np.float64(2)*rng.random()-np.float64(1)
        squared = first**2+second**2
        if squared < 1 and squared != 0:
            scale = np.sqrt(-np.float64(2)*np.log(squared)/squared)
            return width*scale*first+center


def _integer(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer excluding bool")
    return int(value)


def _run(initial, project, evaluate, *, max_iterations, maximum_temperature,
         minimum_temperature, cooling_rate, rng):
    """Private control-flow seam for independent scalar-objective tests."""
    best = initial.copy()
    chi = np.float64(1e20)
    evaluations = 0
    last = best.copy()
    last_fit = np.float64(np.nan)

    def visit(vector):
        nonlocal best, chi, evaluations, last, last_fit
        last = project(vector)
        last_fit = np.float64(evaluate(last))
        evaluations += 1
        # Explicit Python domain policy: never present an undefined/sentinel
        # best point as a successful optimizer result, and never add penalties.
        if not np.isfinite(last_fit) or last_fit > 1e20:
            raise ValueError("ANNEAL objective must be finite and <= the source CHISMIN sentinel 1e20")
        if last_fit <= chi:  # ONEITPARA callback side effect, randnumgen.f:178
            chi = last_fit
            best = last.copy()
        return last_fit

    best_fit = visit(best)  # initial UPDPARS -> ARRANGEPARS -> GOF
    steps = 0
    temperature = np.float64(maximum_temperature)
    amplitudes = np.full(len(best), temperature, dtype=np.float64)
    while steps <= max_iterations and temperature > minimum_temperature:
        steps += 1
        proposal = np.empty_like(best)
        trial_amplitudes = np.empty_like(best)
        with np.errstate(divide="raise", invalid="raise", over="raise"):
            for m in range(len(best)):
                proposal[m] = _normal(rng, best[m], best[m]*amplitudes[m])
                trial_amplitudes[m] = np.float64(1.48258)*abs((proposal[m]-best[m])/best[m])
        fit = visit(proposal)
        change = fit <= best_fit  # exact <=, no tolerance helper in this file
        if change:
            best_fit = fit
            amplitudes = trial_amplitudes
        small = np.count_nonzero(amplitudes < minimum_temperature)  # BEFORE cooling
        if not change:
            amplitudes = amplitudes*cooling_rate
        if small == len(best):
            temperature = temperature*cooling_rate
            amplitudes = np.full(len(best), temperature, dtype=np.float64)
            chi = np.float64(1e20)
            best_fit = visit(best)
        if steps % 10 == 0:  # separate IF: can cause a second evaluation
            chi = np.float64(1e20)
            best_fit = visit(best)
    return LegacyAnnealResult(best, chi, last, last_fit, best_fit, steps,
                              evaluations, temperature, amplitudes)


def run_legacy_anneal(state, initial_vector, experimental, *, max_iterations,
                      seed, maximum_temperature, minimum_temperature, cooling_rate,
                      phase_grating_method="fft", objective_mode="legacy",
                      evaluation_callback=None, execution=None) -> LegacyAnnealResult:
    """Run source-derived ANNEAL in physical coordinates, with no hidden RNG.

    max_iterations preserves MAXITS literally: the source permits MAXITS+1
    proposals. Evaluations additionally include initialization and resets.
    All objectives use evaluate_refinement_trial; unsupported states propagate.
    """
    n = len(build_refinement_parameter_layout(state).parameters)
    if n > 1000:
        raise ValueError("ANNEAL supports at most 1000 active parameters; no truncation")
    initial = _vector(initial_vector, n)
    max_iterations = _integer(max_iterations, "max_iterations")
    seed = _integer(seed, "seed")
    maximum_temperature = np.float64(_finite(maximum_temperature, "maximum_temperature"))
    minimum_temperature = np.float64(_finite(minimum_temperature, "minimum_temperature"))
    cooling_rate = np.float64(_finite(cooling_rate, "cooling_rate"))
    if maximum_temperature <= 0 or minimum_temperature <= 0 or not 0 < cooling_rate < 1:
        raise ValueError("temperatures must be positive and cooling_rate must be in (0,1)")

    def project(vector):
        return pack_refinement_parameters(apply_refinement_parameters(state, vector))

    def evaluate(vector):
        return evaluate_refinement_trial(state, vector, experimental,
            coordinate_mode="physical", phase_grating_method=phase_grating_method,
            objective_mode=objective_mode, execution=execution, evaluation_callback=evaluation_callback).objective

    return _run(initial, project, evaluate, max_iterations=max_iterations,
        maximum_temperature=maximum_temperature, minimum_temperature=minimum_temperature,
        cooling_rate=cooling_rate, rng=np.random.Generator(np.random.PCG64(seed)))
