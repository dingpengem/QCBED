"""Source-derived SIMPARA (atomicscatfact.f), with explicit Python startup/RNG."""

from dataclasses import dataclass
from numbers import Integral

import numpy as np

from .anneal import _normal
from .multislice import _snapshot
from .refinement import build_refinement_parameter_layout, _vector, _entries
from .refinement_trial import evaluate_refinement_trial
from .structure_factors import _finite


@dataclass(frozen=True, eq=False)
class LegacySimparaResult:
    """OUTVAR/CBEST and last callback, with source temperature/direction state.

    best is retained state, not an independent all-time archive. The callback
    may move OUTVAR on ties; periodic reevaluation can change its score.
    """

    best_vector: np.ndarray
    best_objective: np.float64
    last_vector: np.ndarray
    last_objective: np.float64
    central_vector: np.ndarray
    direction: np.ndarray
    perturbation_scales: np.ndarray
    temperature: np.float64
    cycles: int
    evaluations: int
    improvements: int
    stop_reason: str

    def __post_init__(self):
        for name in ("best_vector", "last_vector", "central_vector", "direction", "perturbation_scales"):
            object.__setattr__(self, name, _snapshot(getattr(self, name)))


def _parabola(central, inverse, trial, c0, c1, c2, counter):
    new = np.empty_like(central)
    qualified = 0
    # Source explicitly handles invalid parabola results with a half-step.
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        for i in range(len(new)):
            counter += 1
            a = (((c1-c0)*(trial[i]-central[i])-(c2-c0)*(inverse[i]-central[i])) /
                 ((inverse[i]**2-central[i]**2)*(trial[i]-central[i])-
                  (trial[i]**2-central[i]**2)*(inverse[i]-central[i])))
            fallback = central[i]+(-1.0)**counter*(trial[i]-central[i])/2.0
            if a > 0:
                b = ((c2-c0)-a*(trial[i]**2-central[i]**2))/(trial[i]-central[i])
                new[i] = -b/(2*a)
                if new[i]/new[i] == 1:
                    qualified += 1
                else:
                    new[i] = fallback
                if new[i] > 10000:
                    new[i] = 10000
            else:
                new[i] = fallback
    return new, qualified, counter


def _run(initial, uncertainty, evaluate, *, max_evaluations, high_temperature,
         low_temperature, cooling, random_order, rng):
    """Literal cycle order; evaluate may mutate its vector as ONEITPARA does."""
    out = initial.copy()
    n = len(out)
    temperature = np.float64(high_temperature)
    scales = temperature*uncertainty
    direction = np.zeros(n, dtype=np.float64)  # explicit Python startup policy
    since_change = 0  # source leaves both locals uninitialized
    counter = cycles = evaluations = improvements = 0
    chi = np.float64(1e20)
    last = out.copy()
    last_score = None

    def visit(vector):
        nonlocal chi, evaluations, last, last_score
        score = np.float64(evaluate(vector))
        if not np.isfinite(score) or score > 1e20:
            raise ValueError("SIMPARA objective must be finite and <= source sentinel 1e20")
        if not np.all(np.isfinite(vector)):
            raise ValueError("SIMPARA callback returned non-finite parameters")
        last, last_score = vector.copy(), score
        evaluations += 1
        # Driver passes PARMIN as OUTVAR: emulate the callback's <= write.
        if score <= chi:
            chi = score
            out[:] = vector
        return score

    c0 = best = visit(out)
    central = out.copy()
    change = 0
    while evaluations <= max_evaluations and temperature >= low_temperature*cooling:
        cycles += 1
        if random_order:
            order, used = [], set()
            while len(order) < n:
                counter += 1
                index = int(rng.random()*n)
                if index not in used:
                    used.add(index)
                    order.append(index)
        else:
            order = range(n)
        norm = np.float64(0)
        for j in range(n):
            norm += abs(direction[j])
        perturbation = np.empty(n, dtype=np.float64)
        if change > 0:
            for j in order:
                counter += 1
                direction[j] = (out[j]-central[j])/(abs(out[j]-central[j])+1e-30)
                perturbation[j] = direction[j]*abs(_normal(rng, np.float64(0), scales[j]))
                central[j] = out[j]
        elif norm == 0 or (norm > 0 and since_change % 2 == 1):
            dot = np.float64(1e20)
            while dot > 1.5:
                dot = np.float64(0)
                for j in order:
                    counter += 1
                    perturbation[j] = _normal(rng, np.float64(0), scales[j])
                    dot += direction[j]*perturbation[j]/(abs(perturbation[j])+1e-20)
                dot = abs(dot)
        else:
            for j in order:
                counter += 1
                perturbation[j] = direction[j]*abs(_normal(rng, np.float64(0), scales[j]))
        trial = central+perturbation
        inverse = central-perturbation
        c1 = visit(inverse)
        c2 = visit(trial)
        new, qualified, counter = _parabola(central, inverse, trial, c0, c1, c2, counter)
        c3 = visit(new) if qualified >= 1 else np.float64(1e20)
        change = 0
        for flag, score, vector in ((1,c1,inverse), (2,c2,trial), (3,c3,new)):
            if score < best:
                change, best = flag, score
                out[:] = vector
        if change > 0:
            improvements += 1
            since_change = 0
        else:
            since_change += 1
        if cycles % n == 0:
            best = visit(out)
            for i in range(n):
                central[i] = out[i]
                direction[i] = (out[i]-central[i])/(abs(out[i]-central[i])+1e-30)
            change = 0
        temperature = temperature/cooling if change > 0 else temperature*cooling
        scales = uncertainty*temperature
        c0 = best
    reason = "max_evaluations" if evaluations > max_evaluations else "temperature"
    return LegacySimparaResult(out, best, last, last_score, central, direction,
                              scales, temperature, cycles, evaluations, improvements, reason)


def run_legacy_simpara(state, initial_vector, experimental, *, uncertainties,
                       max_evaluations, high_temperature, low_temperature, cooling,
                       seed, random_order=False, phase_grating_method="fft") -> LegacySimparaResult:
    """One driver-compatible physical-coordinate SIMPARA invocation.

    uncertainties are absolute physical widths, not D1/D2 factors. Startup
    direction and no-improvement counter are explicitly zero in Python.
    max_evaluations is the soft source MAXITS limit, including initialization.
    """
    layout = build_refinement_parameter_layout(state)
    n = len(layout.parameters)
    if not 1 <= n <= 1000:
        raise ValueError("SIMPARA driver callback supports 1..1000 active parameters")
    initial = _vector(initial_vector, n)
    uncertainty = _vector(uncertainties, n)
    if np.any(uncertainty < 0):
        raise ValueError("uncertainties must be nonnegative absolute widths")
    for value, name in ((max_evaluations, "max_evaluations"), (seed, "seed")):
        if isinstance(value, (bool,np.bool_)) or not isinstance(value, Integral) or value < 0:
            raise ValueError(f"{name} must be a nonnegative integer excluding bool")
    high = np.float64(_finite(high_temperature, "high_temperature"))
    low = np.float64(_finite(low_temperature, "low_temperature"))
    cooling = np.float64(_finite(cooling, "cooling"))
    if high <= 0 or low <= 0 or not 0 < cooling < 1:
        raise ValueError("temperatures must be positive and cooling must lie in (0,1)")
    if not isinstance(random_order, (bool,np.bool_)):
        raise TypeError("random_order must be bool")

    def evaluate(vector):
        result = evaluate_refinement_trial(state, vector, experimental,
            coordinate_mode="physical", phase_grating_method=phase_grating_method)
        # UPDPARS writes PT only for these categories. Reuse mapped values;
        # do not repack all coordinates or apply ARRANGEPARS' near-zero floor.
        for i, (category, _, value) in enumerate(_entries(result.physical_state)):
            if category in ("constant_vacancy", "vacancy_coefficient", "absorption", "debye_waller"):
                vector[i] = value
        return result.objective

    with np.errstate(divide="raise", invalid="raise", over="raise"):
        return _run(initial, uncertainty, evaluate, max_evaluations=int(max_evaluations),
                    high_temperature=high, low_temperature=low, cooling=cooling,
                    random_order=random_order, rng=np.random.Generator(np.random.PCG64(seed)))
