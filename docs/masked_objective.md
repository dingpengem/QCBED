# Explicit masked objective policy

“An explicit non-legacy numerical policy for datasets whose zero-weight masked
points contain zero variance.” No historical Fortran executable equivalence is
claimed for this policy. `compute_legacy_objective` and legacy YNORM arithmetic
are unchanged.

`qcbed.masked_objective.compute_selected_objective` defaults to
`objective_mode="legacy"`. Explicit `objective_mode="masked"` calls the separate
`compute_masked_objective`. `evaluate_refinement_trial`, the four existing
optimizer entry points, `run_legacy_refinement`, `validate_refinement` and
`run_local_optimization` also accept `objective_mode`, defaulting to `"legacy"`.
The selected policy is propagated to every evaluation and diagnostic reevaluation.

Let V be the points whose `experimental.weights != 0`. Select V **before** any
variance division. On V, retain the supplied weights w, calculated values c,
experimental values x, and variances v:

* IREF > 0: s = IREF.
* IREF = 0: s = sum(w*c*x/v) / sum(w*c*c/v).
* IREF < 0: s = sum(w*abs(x)) / sum(w*abs(c)).
* Objective = sqrt(sum(w*(s*c-x)*(s*c-x)/v) / sum(w*x*x/v)).

Every sum uses the same V in beam-major accumulation order. Weight is not
replaced by a boolean multiplier or discarded. Negative nonzero weights are
also selected; if cancellation makes the normalizer nonpositive or the squared
error negative, the calculation fails explicitly rather than clipping values.

Every retained variance must be finite and strictly positive, including in
fixed-scale and equal-count modes. Errors identify the first offending tilt,
beam index/H,K, variance and weight. No variance is repaired, replaced, or given
an epsilon. Empty masks, singular/nonfinite YNORM, and undefined normalized
errors also raise. The existing experimental model already rejects nonfinite
input arrays; this policy does not weaken that contract.

The baseline equations are visible in `private predecessor implementation`; they contain no
zero-weight skip. The new selection rule is an explicitly authorized Python
policy, not a recovered missing Fortran branch. Input arrays remain unchanged.

See [local optimization](local_optimization.md) for a runnable Python example,
checkpoint files and the exact progress semantics.
