# Modern physical refinement catalogue

[Documentation index](README.md) · [Current architecture](ARCHITECTURE.md) ·
[Researcher manual](user_manual.md)

The modern active vector is built from complete, committed bindings in the
compiler's source-defined ordering. It is not a fixed two- or eighteen-parameter
vector. Requested depth is optional. YAML/table order does not define numerical
order. Initial values come from the committed physical model, not a second copy
in a binding. All active modern targets require finite explicit lower/upper
rejection bounds; out-of-bound trials are rejected, never clipped.

## Implementation audit

| Target | Canonical source / units | Physical mapping and domain | Modern execution / GUI |
|---|---|---|---|
| Accelerating voltage | `microscope.energy_kv`, kV | Existing energy-dependent sample rebuild; positive | Available only with absorption disabled; legacy voltage/absorption update remains unsupported |
| Requested depth | `crystallography.depth_angstrom`, Å | Existing independent whole-slice expansion; positive; effective depth is derived | Available, optional |
| Region depths | Ordered pre-sliced `sample.regions`, Å | Existing per-region depth layout | Legacy/pre-sliced workflow; modern crystal bridge currently constructs one region |
| Debye–Waller B | Physical site group `debye_waller_b`, Å² | Group ownership → species thermal factor → all recast/stored occurrences; nonnegative (strictly positive with absorption) | Available for every supported group (currently one group per element) |
| Absorption a/b/c | Explicit site-group absorption triple, legacy coefficient units | Group → species coefficient → existing absorptive grating; nonnegative | Available when explicit absorption is enabled; independent C only. Linked C remains supported by the legacy path, not the modern bridge |
| Scan A/B/C first/second | `scan.corner_a/b/c`, reciprocal-basis coefficients | Existing six scan-corner targets; signed, unchanged basis convention | Available by explicit selection; experimental arrays/weights remain fixed |
| Chemical occupancy | Site-group occupancy, fraction | Static occupancy propagates; no dynamic physical occupancy binding exists | Fixed input only; not interchangeable with legacy vacancy |
| Atom x/y | Legacy stored-slice fractional coordinates | Legacy propagation implemented; modern original-site-to-recast coordinate update absent | Unavailable in modern refinement; must not bind simulation-copy indices as physical sites |
| Reflection multipliers | Ordered legacy reflection records, dimensionless | Legacy elastic/absorptive factors implemented, duplicate records retain identity | Modern strict bridge uses unity factors and has no original-hkl multiplier mapping; unavailable |
| Constant/profile vacancy | Legacy region vacancy records | Layout/mapping preserved, nonzero vacancy evaluator unsupported | Unavailable |
| DILATION | Legacy vacancy-coupling coefficient, dimensionless | Existing folding semantics; not lattice dilation/strain | Deferred: nonzero vacancy propagation unsupported; zero-vacancy scale is unity |

All four existing numerical optimizers accept dynamically sized vectors. Modern
execution reuses those engines through a physical-vector evaluation adapter;
no optimizer equations are replaced. PARAMIN/SIMPARA now use the same durable
journal and safe cancellation boundaries as ANNEAL/AMOEBA. The legacy two-target
B7 diagnostic adapter remains separately supported.

## D1, D2, PLAMDA: exact meanings

For the legacy packed scaling reference `x0`, simplex scaling is
`a = 1 / ((D2 - D1) * x0)`, `b = D1 * x0`, `u = a * (x - b)`.
D1 and D2 are dimensionless factors, **not** hard bounds. AMOEBA adds PLAMDA
to a simplex coordinate, so its physical displacement is
`PLAMDA * (D2 - D1) * x0`. A zero starting value or equal factors makes this
scaling singular and is rejected in modern readiness; signed PLAMDA retains existing semantics.
Legacy packing uses a `1e-10` near-zero floor for its scaling reference. Modern
run initialization reads the raw committed catalogue value (including zero),
not that packing floor. The legacy compiler and packing API are unchanged.

PARAMIN uses `BOT = D2`, `TOP = D1` as multiplicative coordinate-search factors;
its starting values and BOT/TOP divisors must be nonzero. PLAMDA is not used.
SIMPARA uses PLAMDA as a nonnegative **absolute physical uncertainty width**,
multiplied by its current temperature. ANNEAL uses its existing temperature /
adaptive-amplitude proposal rule, not PLAMDA. None of these fields substitutes
for explicit modern rejection bounds. Modern proposals are checked before
physical folding; legacy callback behavior remains on the legacy path.

## Journals, results and compatibility

New modern runs write version-2 variable-vector stage and parent journals, with
physical identity, initial value, bounds and each evaluated vector. Results show
initial, displayed, latest, best and retained values separately, plus per-target
history. Viewing results never recomputes the objective. Version-1 two-target
journals keep their original replay branch; checkpoint ownership/data checks
remain mandatory. New aliases enter configuration identity when active; merely
making additional targets available does not invalidate older two-target IDs.

Unsupported families remain unavailable with reasons. In particular this change
does not implement symmetry-aware coordinate refinement, independent same-element
thermal species, vacancy physics, lattice dilation, or modern original-hkl
reflection-factor ownership.

## Recommended frontend metadata (stable defaults)

The GUI and compact YAML loader share `optimizer_defaults.recommended_binding`.
Explicit values always take precedence. Profiles are independent of physical hard
bounds and never adjust a physical starting value:

| Target category | D1 | D2 | PLAMDA | Existing source |
|---|---:|---:|---:|---|
| Requested depth | 1.05 | 0.95 | 0.5 | `depth_refinement.enable_requested_depth`, reference fixture |
| Debye–Waller B | 1.05 | 0.95 | 0.0002 | `examples/b7a_refinement_layout.py`, reference fixture |
| Absorption coefficients | 1.05 | 0.95 | 0.005 | Reference absorption search profile (b/c); shared coefficient-category frontend profile |
| Scan coordinates | 1.05 | 0.95 | 0.0005 | `tests/data/reference_config.yml` |
| Voltage | 0.9 | 1.1 | 0.01 | Existing `test_modern_dynamic_refinement.selection` voltage regression profile |

These are reproducible initialization profiles, not newly optimized step sizes or
universal convergence guarantees. No hard bounds are inferred. Unsupported future
categories require an audited profile or explicit advanced metadata. The same
stored PLAMDA retains each optimizer's interpretation described above: choosing a
workflow does not rescale it. Zero-start restrictions and all domain/readiness
checks remain unchanged. High-level defaults remain `modern_yaml.OPTIMIZER_DEFAULTS`;
loaded/custom budgets, temperatures and tolerances are not reset by this UI change.

## C-stage editor sequence

C has four sub-tabs: **Simulation Setup**, **Refinement Targets**, **Optimizer**,
and **Execution / Output**. Validate A and B, then Validate Simulation Setup to
generate targets from committed physical inputs. Select targets and hard bounds;
recommended metadata constructs the existing binding objects automatically.
Changing setup requires re-validation and preserves compatible committed bindings.

The normal target table never displays D1/D2/PLAMDA. The collapsed **Expert / Legacy settings — per-target search metadata** section in Optimizer shows them, with optimizer-specific
usage and editable overrides when recommendations are unchecked. Selecting
recommendations explicitly replaces per-target overrides; merely opening a file
or changing workflow never does. Inactive fields need no manual input.

**Validate C / Run Readiness** commits targets and validates execution settings
and external paths. Existing-output replacement remains an explicit confirmation
condition, separate from physical validity. No output filename input is required.
No optimizer evaluations occur during these validations.

High-level search budgets and tolerances are in **Advanced optimizer settings**;
per-target D1/D2/PLAMDA are explicitly **Expert / Legacy**. Fresh modern drafts use
the existing rejecting-bound policy and batched default so these implementation
choices do not block the basic physical-target workflow. Loaded policies and all
checkpoint identity checks remain unchanged.

## Hard-bound audit and immediate target display

Physical hard bounds have no universal defaults in `BoundsConfig` or
`RefinementBoundRules`. Legacy D1/D2 and PLAMDA tables describe search/scaling,
not hard rejection intervals. The GUI now resolves **configured** bounds for
unselected as well as selected targets using the existing backend resolver:
specific target > category absolute > category relative > no bound.

| Supported family | Existing evidence | Default applicability |
|---|---|---|
| Requested depth | `examples/reference_configuration.py`, `depth_refinement.py`: 600–1200 Å; small tests use 1–20 Å | Specimen/test-specific; no universal finite interval |
| Group Debye–Waller B | Reference category 0.5–1.5 Å²; `b7a_refinement_layout.py` explicitly calls these project bring-up limits | Display when configured; not a universal material bound |
| Absorption a/b/c | Reference category 0.5–1.5 in legacy coefficient units | Display when configured; other valid starts need different bounds |
| Voltage | Dynamic regression uses 100–300 kV | Test-specific; positivity alone supplies no finite upper bound |
| Scan A/B/C components | Reference relative fraction 0.2, absolute fallback 0.05 | Explicit user-selected rule; not universal geometry limits |

When a rule is present its lower/upper values appear immediately, before Refine
is checked. Selection fills only the audited optimizer metadata; validation makes
a complete supported selection Active. A row lacking a canonical/configured bound
says **no canonical hard-bound default** and requires explicit specimen-appropriate
input. No reference-campaign interval is silently injected into fresh projects.
Unsupported families remain disabled.

Unchanged inherited category/relative bounds stay as rules in compact YAML;
the GUI does not duplicate them into a specific bound for every row. Edited bounds
are saved as target-specific overrides. Existing specific bounds retain priority.
