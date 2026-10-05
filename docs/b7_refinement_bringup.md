# B7a — physical site groups and refinement binding bring-up

B7a defines a nonempty layout, validates bounds, reapplies the initial vector and
performs one deterministic transient perturbation. **No optimizer runs.** The
initial mode is conventional because it is the simpler validated signal path;
its smaller B6b objective is not the justification. Gradient remains available.

## Audit and architecture

The existing `ParameterTarget`/`ParameterBinding`, compiler target map,
`RefinementParameterLayout`, masks, packing, D1/D2, PLAMDA and bounds remain the
numerical parameter system. `available_refinement_targets` and `set_refinement`
now understand the opt-in physical group aliases. Optimizers still see dynamic
numeric vectors, and none is changed or invoked here. Existing non-group
configuration paths retain their previous behavior.

| Family | B7a status | Physical ownership | Existing binding and remaining work |
|---|---|---|---|
| Debye–Waller B | READY, enabled | Physical site group | Existing region/species B resolved from a stable physical group; original/recast/slab metadata and scattering rebuild tested |
| Global lattice/dilation | DEFERRED | Original conventional lattice | Existing legacy dilation folds values and changes derived state; a physical original-cell rebuild and reviewed bounds are not implemented here |
| Requested sample depth | DEFERRED, not selected | Sample/region | Existing strict-< cyclic sequence until requested depth is reached/exceeded; reference bounds600–1200 Å exist. Changes can leave sequence unchanged until a boundary; this is not continuous effective depth |
| Occupancy/vacancy | DEFERRED | Physical site group | Fixed chemical occupancy is explicit. Legacy vacancy fractions/profile coefficients have distinct semantics and nonzero trial hazards; do not substitute projected weights |
| Absorption a/b/c | DEFERRED as a target | Group → shared species table | Existing independent coefficients are physically explicit. Group metadata carries (0,1,1), but B7a application enables only B. Linked C unsupported in crystal bridge; no new absorption model |
| Accelerating voltage | DEFERRED | Global microscope | Wavelength binding exists; changing voltage with absorption is rejected by legacy trial semantics. No reviewed voltage bounds |
| Scan A/B/C | DEFERRED | Global scan | Existing six dynamic targets work; retain validated A10 geometry for first bring-up |
| Independent atom x/y | NOT APPLICABLE to symmetric FCC bring-up | Future physical symmetry modes | Legacy stored-atom targets exist, but four independent site edits would break intended symmetry |
| Reflection multipliers | NOT APPLICABLE to first real layout | Original reflection records | Dynamic legacy records supported; real configuration has none, with unity multipliers. Do not add them to absorb model disagreement |

The source-defined numerical ordering remains authoritative. Group or binding
YAML declaration order cannot reorder parameters. A two-element/two-group test
verifies this and verifies that perturbing one group leaves the other unchanged.
The current scattering backend allows one shared B per element. Independently
parameterized same-element groups are **explicitly rejected**, even if their
starting B values happen to match; this avoids conflating physical owners.

## Canonical site-group domain model

`qcbed.configuration.AtomicSiteGroup` is immutable and has:

- `id`: stable physical group ID;
- `element`: canonical element symbol;
- `site_indices`: ordered indices in the original immutable structure;
- `debye_waller_b`: isotropic B in Å², not RMS displacement;
- `occupancy`: chemical fraction0<occupancy≤1;
- `absorption`: optional existing independent (a,b,c) coefficients.

The structure file/ASE owns cell, chemical labels and unwrapped original atomic
coordinates. Explicit QCBED groups own B, occupancy, absorption and refinement
ownership, and may supply/override structure-file scattering metadata. If groups
are supplied, old `debye_waller` species metadata must be empty and the old
absorption settings disabled/empty: two simultaneous authorities are rejected.
Without groups the earlier strict species/structure conflict checks are unchanged.

Group validation rejects duplicate/empty IDs, duplicate/noninteger/out-of-range
member indices, overlap, missing coverage, inconsistent elements, invalid B or
occupancy, nonfinite coefficients and unsupported independent same-element groups.
B is optional in the domain model but mandatory when compiling a scattering sample.
Absorption must be supplied for all groups or absent for all; signed finite
coefficients retain existing input semantics. There is no invented fallback.
Grouping asserts user-provided equivalence; it does not run symmetry discovery.
Reordering the source file requires updating membership explicitly.

The real group is `Al_fcc`, elementAl, members(0,1,2,3), B0.86 Å², occupancy1,
absorption(0,1,1). All four sites share these values. Legacy projected half-weights
are not chemical occupancy. `AtomicSite` adds optional `physical_group_id` and
`original_site_index`; recasting preserves them and depth-slab construction copies
them without changing positions or numerical arithmetic. Simulation provenance is
compact: original group membership + each stored slab's original indices/group
IDs + the existing occurrence sequence. Each of830 atom occurrences is traceable;
there is no binding to a temporary simulation atom number.

## Existing-layout binding and bounds

The canonical target is:

```python
ParameterTarget(ParameterCategory.DEBYE_WALLER, ('site_group', 'Al_fcc'))
```

Compiler alias resolution maps it to the existing region/species target
`('crystal','Al')`, then the existing numerical `(region0,species0)` slot. The
catalogue exposes the physical target, not that numerical slot. B7a grouped
bindings reject non-B categories; their catalogue entries show deferred status.
No second layout, optimizer vector or positional array is introduced.

| Target ID | Scope | Initial | Bounds | Units |
|---|---|---:|---|---|
| site_group.Al_fcc.debye_waller_b | site-group | 0.86 | [0.5,1.5] inclusive | Å² |

These are existing **project bring-up bounds**, explicitly adopted from
`examples/reference_configuration.py`'s Debye–Waller category absolute bounds.
They are not claimed as universal Al limits or experimental confidence intervals.
D1/D2=(1.05,0.95), PLAMDA=0.0002 also come from that reference. The starting
physical parameters remain the B6b inputs. Depth and lattice are not enabled.

`apply_crystallographic_parameters(compiled, vector)` validates finite float64,
dynamic length, supported bindings and compiled rejection bounds **before**
legacy scalar application. Negative inputs are rejected, not absolute-value
folded. Inclusive endpoints work; below/above/NaN/Inf are rejected before rebuild.
The existing physical mapping is used for accepted values; resulting B updates
immutable original group metadata. Pure `compile_config(..., structure=...)`
rebuilds recast cells, strict sample and scattering from that new physical state.
It performs no file I/O and never mutates the baseline.

`physical_target_layout` reports physical identity, scope, initial value, units,
bounds and internal binding in existing deterministic order.
`site_occurrence_provenance` exposes compact ownership/sequence metadata.

## Current YAML and GUI status

The additive `crystallography.site_groups` field is implemented in the typed
model and automatically supported by the existing strict YAML loader/serializer:

```yaml
crystallography:
  # Existing structure, zone, basis, slicing, beams, etc. are still required.
  debye_waller: []
  absorption:
    enabled: false
    coefficients: []
  site_groups:
    - id: Al_fcc
      element: Al
      site_indices: [0, 1, 2, 3]
      debye_waller_b: 0.86
      occupancy: 1.0
      absorption: [0.0, 1.0, 1.0]
```

Refinement uses the existing binding/bounds schema with
`category: debye_waller`, `identifiers: [site_group, Al_fcc]`. Membership and source
order survive Save/reload exactly. Future GUI controls can present Group,
Element, Sites, B, Occupancy, Absorption, Refine, Bounds using these objects and
the existing catalogue. GUI editing, same-element independent scattering groups,
linked absorption and symmetry-breaking displacement modes remain deferred.

## Snapshot/checkpoint boundary

The immutable compiled snapshot contains the original applied sites, group
metadata, aliases, source configuration and numerical layout. Diagnostics store
full original coordinates/lattice plus compact occurrence provenance. YAML and
JSON snapshots reconstruct the baseline (given the referenced structure file);
tests also check identical packed-vector reapplication and state equality.

Existing version1 checkpoint identities are byte-for-byte unchanged for configs
without groups. Group configs opt into version2: the digest additionally records
original cell/sites, zone, original hkl basis, group ownership and target aliases.
Changing only a group ID invalidates compatibility, even if numerical arrays match.
Old identities are never silently reinterpreted. Binding declaration order does
not affect identity. This establishes **identity support**, not optimizer resume:
canonical Run/Resume explicitly rejects site-group configs as diagnostic-only
until physical-group optimizer integration is separately validated in B7b. No
checkpoint is fabricated and existing non-group runner behavior is unchanged.

## Measured real bring-up

```bash
python examples/b7a_refinement_layout.py --output tmp/b7a_refinement_layout
```

The output must be new. It contains `configuration.yml`, `layout_summary.json`
and `layout_summary.md`; baseline source YAML and reference inputs are untouched.

- Active parameter count1; ordered vector[0.86].
- Baseline YNORM14870390.0181119; objective0.939675432747781 — exact B6b match.
- Identical-vector application: exact physical state, amplitudes, intensities,
  YNORM and objective; no drift.
- Transient B0.8602: YNORM14872521.395027455; objective0.9396971244307921.
- Maximum phase-grating change1.1442451797009533e-5; amplitude change
  3.888962983081606e-5. Four original sites, four stored sites and830 simulation
  atom occurrences share the new B. Occupancy, absorption, original positions,
  lattice, depth sequence, beam identities, signal/variance/support and A10 remain
  unchanged. The worsening objective is reported without tuning.

The diagnostic performs exactly three fixed evaluations: baseline, identical
vector reapplication, and one perturbation. B6b's default repeat behavior remains
unchanged; B7a explicitly disables that extra repeat and reports it unmeasured.

Recommended B7b candidate: **shared Al_fcc B only**, with the above project bounds,
provided the physical-group application path is integrated with a controlled
optimizer and checked before enabling Run/Resume. No optimizer, GUI or app build
is started in B7a.

## Validation results

**266 focused tests passed; 2759 full-suite tests passed in159.99 seconds.**
`git diff --check` is clean. Coverage includes membership/domain rejection,
physical group propagation, original/recast/slab/occurrence provenance, dynamic
two-group ordering, unrelated-group isolation, inclusive bounds and nonfinite
rejection before rebuilding, catalogue/edit helper integration, YAML snapshots,
version1 compatibility/version2 ownership identity, exact real baseline and
one-target sensitivity. No optimizer was invoked.

# B7b — controlled single-target PARAMIN

B7b enables a dedicated backend Run/Resume/Cancel integration for exactly
`site_group.Al_fcc.debye_waller_b`, conventional mode, initial0.86 Å²,
bounds[0.5,1.5] Å². No other target or physical input changes. General
`run_config`, GUI group editing and app rebuilding remain deferred; the B7a
restriction on that general runner remains intentional.

## Binary mask contract

The real conventional and gradient prepared datasets both contain only weights0
and1. Positive unique weights are exactly[1.0], with7974 conventional and7281
gradient active points. The current extraction already thresholds integrated
support with `support_fraction > 0.875` (`threshold_support`), so8/9 or0.9
coverage becomes1, not a fractional chi-square multiplier. Exactly0.875 is
excluded. No extraction/scientific behavior was changed.

`audit_binary_masks` makes this a hard B7b input invariant, rechecked alongside
all signal/variance/support hashes for every physical optimizer evaluation.
**Mask/support is an inclusion gate; variance is the statistical uncertainty
denominator.** Included error terms are residual²/variance, with the existing
global YNORM and root-normalized error unchanged. Excluded terms remain exactly
zero, regardless of placeholder observations/variance. Generic legacy APIs retain
their historical ability to accept explicit nonbinary weights; this patch does
not silently reinterpret historical inputs. Gradient is audited but not optimized.

## Optimizer and proposal application

The existing source-derived `paramin._run` engine is used unchanged. Its public
legacy wrapper normally supplies a legacy-state objective; this integration supplies
the B7a physical-group callback instead. There is no new search algorithm or SciPy
optimizer. Numerical proposals remain one-dimensional arrays with deterministic
traversal, BOT0.95/TOP1.05 from the compiled D2/D1 records. No random seed is needed.

The budget is literal MAXFPARA1, tolerance0: up to2 source cycles, at most9
optimizer proposals. Proposals outside[0.5,1.5] are recorded as rejected and return
positive infinity through PARAMIN's existing `allow_invalid` policy, before any
physical rebuild. They have no fabricated YNORM or calculated arrays. Accepted
proposals call `apply_crystallographic_parameters` and the existing B6b forward/
objective diagnostic path exactly once, with repeat disabled. All arrays/objectives
must be finite. Unexpected failures remain failures, not arbitrary finite penalties.

Each accepted proposal validates all four original sites, recast/slab copies and
830 atom occurrences sharing B. Restoring only B in the source configuration must
reproduce the baseline configuration exactly. Original site metadata/positions,
lattice, slice sequence, indexing, scan and experimental manifests remain fixed.
B6b decomposition checks run at every accepted point, including baseline/best/final.

## Measured local checks and result

| Diagnostic B (Å²) | YNORM | Objective |
|---|---:|---:|
| 0.8598 | 14868258.707273409 | 0.9396537400324485 |
| 0.86 | 14870390.0181119 | 0.939675432747781 |
| 0.8602 | 14872521.395027455 | 0.9396971244307921 |

Four pre-optimizer evaluations comprise baseline, identical baseline repeat,
and the two fixed PLAMDA-sized perturbations. Baseline amplitudes, raw/objective
intensities, YNORM and objective repeat exactly. The B6b baseline matches exactly.
Local sensitivity neither adjusts bounds nor hard-codes a search direction.

PARAMIN proposed, in order:
0.86, 0.817, 0.903, 0.8814999999999998, 0.817,
0.7761499999999999, 0.8578499999999999, 0.7965749999999998,
0.7761499999999999. All9 were within bounds and physically evaluated.
Synthetic callback tests separately force rejected proposals and verify zero
physical calls for them.

Best, retained final and last evaluated happen to coincide in this run:
B=0.7761499999999999 Å², YNORM13982329.708364274,
objective0.9304373947054256. They are tracked separately in the journal; the
initial objective is0.939675432747781. PARAMIN reports **not converged** after2
cycles. No global optimum or experimental best-fit claim is made.

## Checkpoint/replay and cancellation

`run_single_target_refinement(compiled, datasets, output_directory, resume=False,
stop_after=None, cancel_requested=None)` is the dedicated B7b backend API.
Fresh runs require a new directory. Resume validates configuration/data/layout
identity before execution. `load_single_target_checkpoint` validates the explicit
`qcbed-b7b-paramin-replay`, version2 JSON schema/checksum. It contains the current
ownership-aware version2 canonical identity, signal/data fingerprint, fixed
optimizer settings, physical target/bounds, baseline diagnostics/provenance,
ordered history, last completed physical evaluation and all-time best.

PARAMIN local variables are reconstructed by replay: restart the same algorithm,
verify every proposed B against the saved prefix, and return the saved objective
without propagating those points again. Then evaluate new proposals normally.
This is a deterministic journal checkpoint, not a claim of serialized native
PARAMIN stack state. Code/format changes that alter proposal semantics must not
reuse this format/version. It is intentionally distinct from general runner
checkpoints; existing version1 readers and identities are unchanged.

Each completed evaluation or rejection is written through the existing atomic,
fsync-backed JSON writer. Cooperative cancellation is checked before the next new
evaluation, leaving the last completed result, best, history and checkpoint intact.
There is no force-kill or mid-propagation cancellation. A failure is separately
recorded as failed and preserves the prior completed journal. The hook can use a
backend cancellation flag; GUI integration remains deferred.

A real run stopped after3 completed proposals and resumed through6 more matched
uninterrupted proposal/history/best/final records **exactly**. Completed-prefix
replay does not duplicate history or physics. Tests cover cancellation after2
proposals, resumption, pre-evaluation cancellation, failed evaluations, missing/
corrupt checkpoints, incompatible identity and fresh-directory overwrite protection.

## Command and output

```bash
python examples/b7b_single_target_refinement.py --output tmp/b7b_single_target_refinement
```

The command performs one main run and a separate interrupted/resumed validation
run under `resume_validation/`. Each logical run uses4 diagnostic evaluations and9
optimizer evaluations (26 physical evaluations total across both validation runs).
It writes `refinement_summary.md/.json`, `history.csv`, `checkpoint_latest.json`
and a separate `refined_snapshot.yml` labeled B7b diagnostic refinement output.
Source configuration remains B0.86. No production/reference file is overwritten.

The next candidate remains controlled shared-B work if a longer optimization is
wanted; do not infer convergence from this pilot. For a separately reviewed B7c,
requested depth is the clearest additional candidate with existing bounds and
explicit whole-slice discretization. Global lattice dilation still needs a proper
original-cell rebuild. No additional target is enabled here.

## B7b validation results

Six focused B7b tests plus the retained zero-weight placeholder regression passed
(**7 focused tests**). The full suite passed **2765 tests in293.49 seconds**.
`git diff --check` is clean. Real uninterrupted versus partial/resumed execution
was compared exactly, in addition to callback-level rejection/cancel/failure tests.
No existing test, objective formula or optimizer algorithm was weakened or changed.

## B7c — requested depth audit and controlled B/depth profile

B7c supersedes the earlier B7a/B7b depth deferral **only for the controlled
backend profile described here**. Lattice dilation and all other target families,
GUI controls, and app packaging remain deferred. The canonical real input stays
at B = 0.86 Å² and requested depth = 840 Å.

### Depth semantics and activation gate

`qcbed.depth_geometry.requested_depth_to_sample_geometry(request, thicknesses)`
is now the shared mapping used by strict B5 construction and legacy refinement
sample preparation. Starting at zero, append stored slices in cyclic order while
`accumulated_depth < requested_depth`, adding each thickness sequentially.
The effective depth is that accumulated sum. Stored slice count describes the
pattern; occurrence count describes its expansion. Neither is a fixed capacity.
No fractional slice, interpolation, `ceil(request / thickness)`, or multiplied
sum replaces the source floating-point addition order.

An expanded sequence owns the interval `(previous_sum, effective_depth]`.
For the real FCC Al cell there are two stored slices, each 2.0245 Å thick:

| Requested-depth interval (Å) | Occurrences | Effective depth (Å) |
|---|---:|---:|
| (836.1184999999964, 838.1429999999964] | 414 | 838.1429999999964 |
| (838.1429999999964, 840.1674999999964] | 415 | 840.1674999999964 |
| (840.1674999999964, 842.1919999999964] | 416 | 842.1919999999964 |

The audit evaluates 838, 839, 840, 841, 842 Å and exact/nextafter boundary
requests, all at B = 0.86. Equal sequences reproduce YNORM and objective exactly:

| Occurrences | YNORM | Objective |
|---:|---:|---:|
| 414 | 14839764.173143458 | 0.9409387474537082 |
| 415 | 14870390.0181119 | 0.939675432747781 |
| 416 | 14900910.798010143 | 0.9384327629151047 |

Requested depth is READY for explicit candidate search: semantics, deterministic
mapping, finite response, bounds, and checkpoint/snapshot provenance are verified.
It is **not** certified as a smooth PARAMIN coordinate. PARAMIN's parabolic
interpolation and SIMPARA's continuous simplex operations do not resolve the
physical identifiability of a plateau. B7c therefore enumerates three distinct
sample sequences and uses unchanged one-dimensional PARAMIN for B at each.

Depth bounds [600, 1200] Å come from the existing Python reference configuration
(`examples/reference_configuration.py`) and full-refinement runner constants
(`examples/run_full_refinement.py`). These constrain requested region depth with
the same cyclic whole-slice semantics; they are not newly inferred Fortran hard
limits. B bounds remain [0.5, 1.5] Å².

### Layout, search, and invariants

The dynamic source-defined numerical order remains **requested depth, then B**:
`sample.requested_depth` (Å), `site_group.Al_fcc.debye_waller_b` (Å²).
There are exactly two active parameters, initially `[840, 0.86]`. Physical target
aliases map into the existing region-depth/Debye–Waller layout; declaration order
does not reorder that vector. Parameter application replaces the canonical depth
and site-group B, then rebuilds derived crystallography and the strict sample.

The candidates are 838.1429999999964, 840, and 840.1674999999965 Å, deliberately
covering 414/415/416 occurrences. The pre-search grid crosses them with B values
0.84, 0.86, and 0.88. Search starts from canonical B = 0.86 independently at each
depth. PARAMIN MAXFPARA = 0 still executes one source cycle: five proposals per
candidate, 15 total. This is a local bring-up, not an exhaustive search or a
convergence claim. Requested values are recorded even when they share geometry;
no evaluation cache was introduced.

Every physical evaluation checks the 7974 conventional active points and binary
weights, prepared-data manifest, unchanged variance/support, original lattice,
site properties other than B, stored-slice geometry, and unrelated configuration
fields. Depth changes only the occurrence sequence; B changes scattering without
changing that sequence. Objectives and normalization use existing routines.

### Execution, checkpoint, and snapshots

```bash
python examples/b7c_two_target_refinement.py --output tmp/b7c_two_target_refinement
```

The example runs the depth audit before adding the depth binding, then a baseline,
nine-point grid, and 15-proposal profile. A second profile is interrupted after
seven proposals and resumed for eight more. In total this validation performs
59 physical evaluations (9 audit + 25 main + 25 split/resumed). Replay does not
repeat completed physical evaluations.

`qcbed.depth_refinement` exposes `enable_requested_depth`, `sample_extent`,
`evaluate_depth_point`, `run_depth_profile`, and `load_depth_checkpoint` for this
scoped backend workflow. Generic canonical/GUI execution remains unchanged.
The outer checkpoint format is `qcbed-b7c-depth-profile`, version 1, containing
canonical ownership identity version 2, the exact layout, data fingerprint,
audited candidates, initial/grid results, and aggregated history. Keep the
`depth_0/`, `depth_1/`, and `depth_2/` child PARAMIN replay journals with it.
Child identities include exact requested/effective geometry and the parent
identity. Existing B-only canonical identity and B7b journal defaults remain
compatible; old B7b checkpoints are not reinterpreted as B7c checkpoints.

Cooperative cancellation is checked before new optimizer evaluations, after the
previous complete child journal is durable. The baseline/grid diagnostics precede
these optimizer cancellation boundaries. Cancelled history and best-so-far are
preserved; resume replays completed proposals through the same optimizer and
checks compatibility. A child journal is authoritative for an evaluation completed
before an unexpected outer-process failure; retain the whole output directory.

History records B, requested/effective depths, occurrence/stored-slice counts,
sequence identity for evaluated proposals, YNORM, objective, and bound status.
Initial, best, last evaluated, and retained final are separate. The explicitly
selected retained final is the best sampled point across the finite profile;
`stage_retained_final` also retains the last inner optimizer's own result.
`refined_snapshot.yml` is a separate diagnostic configuration; derived effective
depth/count/sequence identity live in `refined_geometry.json`. No baseline or
reference input is overwritten. CSV and Markdown/JSON reports contain the full
audit, grid, history, and exact resume comparison.

### Measured B7c result

The real validation completed with an exact B6b baseline and exact uninterrupted
versus interrupted/resumed equality for history, best, last evaluated, and
retained final. The nine-point grid was finite:

| B (Å²) | Objective, 414 occurrences | Objective, 415 | Objective, 416 |
|---:|---:|---:|---:|
| 0.84 | 0.938896580947465 | 0.937500404169553 | 0.936127031762307 |
| 0.86 | 0.9409387474537082 | 0.939675432747781 | 0.9384327629151047 |
| 0.88 | 0.9429767976498928 | 0.9418400976552007 | 0.9407220452015184 |

Best and retained final: B = 0.817 Å², requested depth = 840.1674999999965 Å,
effective depth = 842.1919999999964 Å, 416 occurrences / 2 stored slices,
YNORM = 14442532.573154794, objective = 0.9334495059717921. The last evaluation
has the same physical values but a distinct history index (14 versus best's 11,
zero-based). This shorter B search does **not** outperform B7b's two-cycle result
0.9304373947054256; it establishes mixed-parameter plumbing, not an improved fit
or convergence. The best depth is at the edge of the tiny candidate set.

The canonical input remains unchanged. All evaluations retained 7974 active
weight-one points, exact zero contribution from excluded points, and unchanged
experimental variance/support. Artifacts are under
`tmp/b7c_two_target_refinement/` (generated, not source inputs).

Recommended follow-up, subject to a separate request: examine a justified wider
requested-depth profile and B convergence before adding target families. Global
lattice dilation needs its own original-cell → reciprocal → recast → strict B5
sample → HK-consistency audit. Occupancy/absorption should remain deferred pending
identifiability and physical bounds; none is enabled by B7c.

The focused mapping/application suite passed 59 tests; all seven existing reuse
profiler tests and all six B7b tests also passed. The profiler's source-line marker
was updated to recognize the extracted shared sequence mapper; its timing phases
and disabled behavior are unchanged. Existing B7b saved canonical identity was
compared directly and matched exactly. Cancellation tests cover zero, two, and
five completed optimizer evaluations (including a depth-stage transition), then
resume to all 15 proposals with preserved history using deterministic synthetic
objectives. The real objective was used for the separate interrupted/resumed
validation above.

B7c source additions: `src/qcbed/depth_geometry.py`,
`src/qcbed/depth_refinement.py`, `examples/b7c_depth_audit.py`,
`examples/b7c_two_target_refinement.py`, `tests/test_depth_geometry.py`, and
`tests/test_depth_refinement.py`. Integration changes are confined to
`crystallography/simulation_sample.py`, `crystallography/site_groups.py`,
`crystallography_refinement.py`, `refinement_trial.py`, `config_targets.py`,
`canonical_run.py`, and the existing journal helper in
`single_target_refinement.py`, plus the profiler marker and the three status
Markdown files. No optimizer implementation or objective/propagation kernel changed.

Final B7c full-suite regression: **2792 passed in 288.64 seconds**. The first
full run exposed six stale profiler-marker failures; after updating the marker,
all seven profiler tests and the entire suite passed without weakening tests.
`git diff --check` is clean. **B7c PASS**; stop here, with no automatic B7d work.

## B7d — existing QCBED ANNEAL → local AMOEBA

B7d uses the existing source-derived `qcbed.anneal._run` and
`qcbed.amoeba._run` execution seams, unchanged. These are the engines behind the
canonical `ANNEAL` and `AMOEBA` optimizer options (`ANNE` and `SIMP`, or the
historical combined `ANSI` workflow). No new annealing algorithm or third-party
optimizer is introduced. The adapter is scoped to exactly the B7c source-ordered
physical targets `[sample.requested_depth, site_group.Al_fcc.debye_waller_b]`.

Thickness-dependent dynamical diffraction can produce multiple objective minima;
neither monotonicity nor local convexity is assumed. ANNEAL supplies the global
exploration stage before an existing local method. This label does not guarantee
that a finite stochastic search finds a global optimum.

### Source semantics and explicit case settings

`private predecessor implementation` and the existing Python translation propose each
coordinate from a normal distribution centered on retained PARMIN, with signed
width `PARMIN[m] * TEMP[m]`. On `FIT <= BESTFIT`, widths become
`1.48258 * abs((proposal - best) / best)`; otherwise widths are multiplied by
COOLRATE. The test for all widths below MINTEMP occurs **before** cooling. A
small-width reset cools the overall temperature, restores coordinate widths,
and reevaluates PARMIN. Every tenth iteration also reevaluates PARMIN. These
independent resets can cause two evaluations in one iteration. MAXITS permits
MAXITS+1 proposal iterations. The source has improvement-based acceptance, **not
Metropolis uphill acceptance**; no acceptance law is substituted by B7d.

The canonical real input remains requested depth840 Å and B0.86 Å². Because 840
is outside this run's informed window, the ANNEAL initial state is separately
specified as `[1050 Å, 0.86 Å²]`, never written over the baseline input.

| Setting | B7d value / interpretation |
|---|---|
| General legal requested-depth bounds | [600,1200] Å, unchanged |
| Current FCC Al search window | [1000,1100] Å for both stages |
| B bounds | [0.5,1.5] Å², unchanged |
| Independent seeds | 2026,2027,2028, NumPy PCG64 as in existing Python ANNEAL |
| MAXTEMP | 0.05 relative width: initial depth sigma52.5 Å; B sigma0.043 Å² |
| MINTEMP | 0.0001 relative width |
| COOLRATE | 0.9, retained from reference config |
| MAXITS / maxfsimp | 59: up to60 proposals per seed, plus initial/reset evaluations |
| D1/D2 | (1.05,0.95), unchanged for both targets; not used by ANNEAL |
| PLAMDA | depth0.5, B0.0002, unchanged; not used by ANNEAL |
| Local optimizer | AMOEBA, both targets free, MXF40 soft search-evaluation limit |
| Local tolerances | reference FTOL1e-6 and CHIMIN1e-6 |

`case_settings(config)` returns an immutable settings replacement; canonical
settings/files remain untouched. The reference temperatures100 and1 are literal
fractional widths in this implementation, unsuitable for a100 Å-wide interval
near1050 Å. The explicitly declared .05/.0001 overrides enable meaningful
exploration and adaptation inside that window, without changing the optimizer's
schedule equations. Seeds are fixed before observing results, not tuned.

AMOEBA is chosen over PARAMIN because requested depth produces plateaus;
AMOEBA does not impose a local parabolic fit or require derivatives. Its source
simplex initialization, ranking, expansion/contraction/shrink, tolerances and
soft evaluation budget are unchanged. Simplex coordinates use the existing SCAL
conversion anchored to the compiled baseline and original D1/D2/PLAMDA. A
simplex origin round trip can differ from the handoff by floating-point conversion
roundoff; no input rounding is introduced. Both physical parameters stay free.

### Domain, provenance, and diagnostics

The adapter rejects out-of-window/out-of-B-bounds proposals **before** physical
rebuilding, returning the existing invalid-trial infinity score to the optimizer.
Rejected rows have null objective/YNORM, never fabricated physical outputs.
Positive rejected requested depths may still have a derived geometric extent;
nonpositive invalid requests have null geometry. Physical proposals are never
rounded to integer Å, clipped, abs-folded or converted into fractional slices.
The unchanged canonical mapper expands whole occurrences.

A sparse five-point depth profile (1000,1025,1050,1075,1100 Å at B0.86) is diagnostic
only and does not select the handoff. The lowest finite valid evaluated ANNEAL
objective across completed seeds selects the global best. Initial, best,
last evaluated, last proposed, and optimizer-retained final are distinguished.
The boundary detector flags endpoint plateaus and their immediately adjacent
slice interval (using the local slice thickness, not an arbitrary Å tolerance). A boundary result stops local refinement and calls for
a separately reviewed wider window; it is not reported as a thickness optimum.

Every real evaluation reuses B7c canonical target application and invariants:
7974 conventional included points with weights1, unchanged experimental signal,
variance/support, beam order, A10 corners, lattice/zone/basis/occupancy/absorption,
and shared B inherited through physical/recast/stored-slice ownership. Requested
and effective depths, occurrence count, sequence identity, B, YNORM and objective
are recorded separately. No dense thickness grid or evaluation cache is added.

Engine callbacks expose evaluated values, not each internal acceptance decision
or temperature-stage transition. History therefore distinguishes hard-bound
acceptance from best-so-far updates; unavailable optimizer acceptance is null,
not inferred or mislabeled. Final temperature/amplitudes and iteration/evaluation
counts are recorded from ANNEAL's existing result. Figures mark best-so-far states
and all valid evaluated states rather than inventing an accepted-state trajectory.

### Checkpoint, replay, and cancellation

`qcbed.anneal_local_refinement.run_optimizer_stage()` accepts compiled physical
configuration and prepared datasets. `load_stage_checkpoint()`,
`select_global_best()`, and `boundary_hit()` expose supporting diagnostics.
The distinct `qcbed-b7d-replay` checkpoint format/version1 stores canonical/data
identity, ordered target layout, search window, settings, initial point, seed,
all completed evaluations, RNG state at each callback, current recorded result,
and best-so-far. Completed ANNEAL additionally stores temperature/amplitudes;
completed AMOEBA stores its final simplex/objectives and stop reason.

Resume reruns the **unchanged algorithm** from the original initial point/seed or
simplex, returning saved completed answers and verifying every proposal and RNG
state. This reconstructs internal adaptive/simplex state exactly, without
recomputing completed physical evaluations. It does not restart from best-so-far.
Cooperative cancellation and controlled interruption occur before the next new
evaluation; durable completed answers and their geometry remain available for
replay. Failures are recorded distinctly. Old B7b/B7c checkpoint APIs and formats
are unchanged and are never interpreted as B7d checkpoints.

```bash
python examples/b7d_anneal_local_refinement.py --output tmp/b7d_anneal_local_refinement
```

The fresh-directory command performs baseline/profile checks, three ANNEAL runs,
local refinement only when the global best is interior, then independent real
interrupted/resumed comparisons. It writes separate ANNEAL/local diagnostic YAML
snapshots with state/settings sidecars, full histories/checkpoints, numeric
per-beam squared-objective contributions and landscape figures. Per-beam
contributions sum to objective squared; no per-beam fitted weights are added.
Canonical FCC Al input and structure files are never overwritten.

Basin reproducibility, optimizer stop criteria, and physical convergence are
separate questions. A finite interior result or AMOEBA objective-spread stopping
criterion alone does not establish physical convergence. Lattice dilation,
occupancy, absorption, scan, voltage, reflection/coordinate targets, GUI and app
packaging remain deferred.

The measured pre-ANNEAL sparse profile at B0.86 (diagnostic, not handoff selection):

| Requested Å | Effective Å | Occurrences | YNORM | Objective |
|---:|---:|---:|---:|---:|
| 1000 | 1000.1029999999955 | 494 | 17947145.100170255 | 0.8729005921762424 |
| 1025 | 1026.4214999999954 | 507 | 18560708.485801358 | 0.9119900364458351 |
| 1050 | 1050.7154999999952 | 519 | 19130365.207324605 | 0.9516376833068675 |
| 1075 | 1075.0094999999951 | 531 | 19697158.29052972 | 0.9909735773431866 |
| 1100 | 1101.327999999995 | 544 | 20301364.179723773 | 1.0292577775655067 |

Effective depth can exceed1100 Å for a legal requested1100 Å proposal, because
whole-slice mapping is unchanged. Such endpoint geometry is explicitly a boundary
hit; the window constrains the requested physical parameter, not a clipped
simulated extent. This sparse fixed-B profile alone cannot establish landscape
monotonicity or the best joint B/depth basin.

The example also supports process-level continuation with `--resume --output
<existing-directory>`. It requires a first-stage checkpoint, verifies compatible
stage identities, replays completed stages without physical reevaluation, and
rechecks the boundary policy before any local handoff. It never restarts ANNEAL
from the best-so-far vector. Stage-specific backend interruption/cancellation is
available through `run_optimizer_stage(..., stop_after=..., cancel_requested=...)`.

Boundary adjacency is conservative: the plateau containing an endpoint and the
immediately neighboring slice interval are flagged. The interval width comes
from `effective_depth - plateau_lower_exclusive`; no requested-depth rounding or
fixed arbitrary Angstrom tolerance is used. This policy only governs reporting
and whether local refinement may start; it does not modify ANNEAL proposals.

### Real ANNEAL/local results

| Seed | Requested Å | Effective Å | Occurrences | B Å² | YNORM | Objective | Boundary-adjacent |
|---:|---:|---:|---:|---:|---:|---:|---|
| 2026 | 1004.8311315846368 | 1006.1764999999955 | 497 | 0.5437869170250926 | 13939370.73628736 | 0.8588311454653902 | No |
| 2027 | 1001.0926255835353 | 1002.1274999999955 | 495 | 0.8204646978746574 | 17461256.792640768 | 0.8707128065526254 | Yes |
| 2028 | 1002.0329740580512 | 1002.1274999999955 | 495 | 0.807407856773129 | 17286634.024007335 | 0.8691507634931115 | Yes |

Each ran60 proposal iterations. Total evaluations were67/68/68, including source
initialization and reset reevaluations. Valid physical evaluations were61/47/56;
bound rejections were6/21/12. Distinct valid sequence geometries were11/16/11.
Final overall temperatures were.05/.045/.045. Histories retain every requested
proposal; no repeated-geometry collapse is used. Seed2026's depth width adapted
to about1.36e-8 relative width under the unchanged source rule. This very fine
late width is an observed algorithm effect, not an imposed dense grid or rounding.

Seed2026's lower finite objective selects the handoff. Its best, last evaluated
and retained final happen to coincide; they remain separate records. The
other seeds approach the lower boundary and substantially different B values.
These are different retained states, **not proof of reproducible basin convergence
or proven distinct local minima**. A wider informed window should be reviewed
before interpreting a physical thickness optimum; B bounds were not changed.

AMOEBA kept both targets free and completed7 iterations /16 evaluations (10
valid physical evaluations,6 bound rejections), stopping on `objective_spread`.
Best/retained final: requested1002.2061315846368 Å,
effective1004.1519999999955 Å (496 occurrences), B0.5437788545250926 Å²,
YNORM13906515.486930933, objective0.8541807283629408. This result does not hit the
conservative depth-boundary guard. Last evaluated is distinct:
requested1003.5186315846368 Å, same effective depth, B0.5437871857750927 Å²,
YNORM13906621.548532056, objective0.8541807834310685.

Local B moved only about8.06e-6 Å²: the inherited B PLAMDA is small. Depth moved
one effective occurrence from handoff. Thus AMOEBA's spread criterion is **not
physical convergence**; a plateau and small simplex width can satisfy it. The
seed differences and proximity of other runs to the search boundary reinforce
that limitation. No local jump to an unrelated basin is demonstrated.

| Comparison | Requested Å | Effective Å | B Å² | YNORM | Objective |
|---|---:|---:|---:|---:|---:|
| Canonical B6b baseline | 840 | 840.1674999999964 | 0.86 | 14870390.0181119 | 0.939675432747781 |
| Historical B7b best | 840 | 840.1674999999964 | 0.77615 | 13982329.708364274 | 0.9304373947054256 |
| B7d ANNEAL best | 1004.8311315846368 | 1006.1764999999955 | 0.5437869170250926 | 13939370.73628736 | 0.8588311454653902 |
| B7d local best | 1002.2061315846368 | 1004.1519999999955 | 0.5437788545250926 | 13906515.486930933 | 0.8541807283629408 |

B7b and B7d use different thicknesses; they are not equivalent B-only fits.
Canonical input remains840/.86. No search beyond1000–1100 Å or outside B bounds
was performed: rejected proposals did not run physical simulations.

B7d source additions: `src/qcbed/anneal_local_refinement.py`,
`examples/b7d_anneal_local_refinement.py`, and
`tests/test_anneal_local_refinement.py`. Only this documentation, README and
PROJECT_SPEC were updated alongside them; no optimizer, physics, objective,
compiler, B7b/B7c checkpoint implementation, GUI or packaging module was changed.

Focused validation: **81 passed in4.09 seconds**, including existing ANNEAL and
AMOEBA numerical tests, B7c mapping/invariants, new stage replay/cancel/handoff,
bounds, seed reproducibility, boundary adjacency and production control validation.
Full-suite validation: **2816 passed in294.46 seconds**. `git diff --check` is clean.
Cancellation tests use deterministic synthetic objectives for both engines,
before the first physical evaluation and after two completed physical evaluations;
resumption preserves the completed prefix and finishes successfully. Exceptions
remain failures, not cancellations. Real interrupted/resumed comparisons are
reported separately below.

### Per-beam comparison

These contributions add to objective squared, not objective. Beam order is the original prepared-data order.

| Beam ID | Original HKL | Baseline | ANNEAL best | Local best |
|---|---|---:|---:|---:|
| lattice_0_0 | (0, 0, 0) | 0.26883262497 | 0.0339689424965 | 0.0340089380395 |
| lattice_-1_0 | (-2, 0, 0) | 0.204378462272 | 0.248819292928 | 0.245383131928 |
| lattice_1_0 | (2, 0, 0) | 0.0358972252678 | 0.0492555432498 | 0.0483043113317 |
| lattice_0_-1 | (0, -2, 0) | 0.0865945262484 | 0.104841622624 | 0.103259868979 |
| lattice_0_1 | (0, 2, 0) | 0.08331367794 | 0.10810133078 | 0.10645357278 |
| lattice_-1_-1 | (-2, -2, 0) | 0.0285061259441 | 0.0278879428027 | 0.0278120427968 |
| lattice_1_1 | (2, 2, 0) | 0.126024162853 | 0.118304838713 | 0.118179770702 |
| lattice_-1_1 | (-2, 2, 0) | 0.0241222019497 | 0.0218907120517 | 0.0217854019043 |
| lattice_1_-1 | (2, -2, 0) | 0.0253209114643 | 0.0245207107764 | 0.0244376782451 |

Most net improvement comes from the direct beam: its squared-objective contribution drops
from0.268832625 to0.034008938. Four of the nine beams worsen relative to baseline,
while five improve. Thus the aggregate improvement is not uniform across beams.
Local refinement improves eight contributions relative to ANNEAL best; the direct beam worsens slightly.
No per-beam weighting or normalization was fitted to produce this comparison.

Both real replay validations passed **exactly** after interruption at12 scoring
evaluations: ANNEAL and AMOEBA matched uninterrupted proposal history, RNG records,
requested/effective geometry, objectives, best/last/retained state and final
optimizer metadata. The validation command also resumed completed ANNEAL journals
without reevaluating their physical prefix before the reviewed local handoff.
No B7b/B7c checkpoint was changed or reinterpreted.

**B7d PASS for integration**, not physical convergence. Mask/variance invariants
held through every valid objective; excluded points retain zero contribution.
Recommended follow-up is a separately reviewed thickness-window and local-step
adequacy study with the same two targets. A later parameter-family task may audit
global lattice dilation from the original cell through reciprocal/recast/strict
sample/HK consistency; it is not enabled here. Occupancy, absorption, GUI and app
rebuilding remain deferred.

## B7d informed-window rerun: 900–1100 Å

The existing runner now accepts an explicit case window; the original default
1000–1100 Å and its checkpoint identity remain compatible. General legal bounds
stay600–1200 Å, with B bounds0.5–1.5 Å². The new case window is included in stage
identity, so checkpoints cannot silently resume against a different window.

```bash
python examples/b7d_anneal_local_refinement.py \
  --depth-window 900 1100 --output tmp/b7d_anneal_local_refinement_900_1100
```

This rerun retains seeds2026/2027/2028,60 ANNEAL proposals per seed, relative
temperatures.05/.0001, cooling.9, and AMOEBA settings/D1/D2/PLAMDA. The deterministic
window-midpoint initialization is now1000 Å/B0.86; canonical baseline840/.86 is
unchanged. Sparse diagnostics are900/950/1000/1050/1100 Å at B0.86.

Mask-filtered figures use existing evaluated experimental, globally scaled
calculated, and residual arrays. Each state has one9-beam ×3-column figure,
shared intensity/residual color scales, and gray excluded support. The saved NPZ
contains beam/row/column arrays, explicit Boolean mask and global YNORM; excluded
entries are NaN, not displayed zeros. No per-beam normalization is introduced.
Files are `baseline_masked_pattern_match`, `anneal_best_masked_pattern_match`,
and `local_best_masked_pattern_match` with PNG/NPZ suffixes under `figures/`.
`per_beam_residuals.csv` preserves ordered HKL identities and contributions to
objective squared, including differences from canonical baseline.

At B0.86 the wider sparse profile gives objectives0.8698570679639174,
0.8353568122965725,0.8729005921762424,0.9516376833068675,
1.0292577775655067 at requested900/950/1000/1050/1100 Å respectively. This
reveals lower-depth structure unavailable in the old search window, without
selecting an optimizer result from the sparse profile.

Changed-setting/figure tests:26 passed in2.06 seconds. Full regression:2818
passed in297.96 seconds. Masked-array tests verify active values exactly,
column-fast scan reshaping, excluded NaNs, saved mask/YNORM, duplicate labels,
unchanged inputs and PNG/NPZ generation. Window tests verify900–1100 constraints
and incompatible-window resume rejection while preserving general legal bounds.

The 900–1100 Å controlled search completed with these ANNEAL best states:

| Seed | Requested Å | Effective Å | B Å² | Objective |
|---:|---:|---:|---:|---:|
| 2026 | 943.3711350403735 | 943.4169999999958 | 0.7514837950642871 | 0.8069375536245184 |
| 2027 | 945.2443952812022 | 945.4414999999958 | 0.7412699071104649 | 0.8044060827967126 |
| 2028 | 948.7454144090252 | 949.4904999999958 | 0.8220501585908963 | 0.8257054413081869 |

All are interior, neither near900 nor near1100 Å. Seed2027 supplies the lowest
valid objective and the local handoff. All three ran60 proposal iterations /67
scoring evaluations, preserving the previous budget and source reset semantics.

AMOEBA best/retained final: requested941.9631452812022 Å,
effective943.4169999999958 Å, B0.7403503789854657 Å²,
YNORM15236619.368312828, objective0.8039437253022698. It completed21 iterations
and45 scoring evaluations, stopping at the source's soft MXF40 limit, **not
convergence**. Both requested and effective depth remain interior. This improves
on the previous window's local objective0.8541807283629408: the old interval
excluded better lower-depth states. It does not prove a global physical optimum.
The window is not widened further and no new target families are enabled.

The wider run's ANNEAL and local replay comparisons both matched exactly.
All9 beams improve versus canonical baseline for both selected states; none
worsen. The masked NPZ exports have identical experimental values and identical
7974-point masks across all three states, exact calculated-minus-experimental
residuals, and NaN outside support. PNGs were visually inspected. Complete
per-beam numeric differences are in the new directory's `per_beam_residuals.csv`.
The requested/effective best depths are interior; the range was not widened again.
`git diff --check` is clean. No new target families, GUI or app rebuild was started.

## B7e: gradient-mode dilation audit — deferred

Gradient is the primary mode for B7e. The canonical fixed baseline is 7281
active points, binary active weights `[1.0]`, YNORM 13447206.106343817 and
objective 1.1896911365294378. Conventional B7d scores are historical results,
not equivalent gradient scores. Support fraction > 0.875 determines inclusion;
variance remains the statistical denominator.

The scientific gate fails before optimization. In `private predecessor implementation`,
DILATION is a dimensionless **vacancy-coupling coefficient**, with the source
physical domain 0–1. Selected values are mapped by absolute value and reciprocal
above 1, not by a hard-bound rejection. For S records, APPLIEDDILATION starts at 1
and becomes `(1 - VAC*DILATION)**(1/3)` only when `0 < VAC*DILATION <= 0.5`.
`private predecessor implementation` corroborates the initial-state equation;
its lines 2156–2165 corroborate packing. The source uses `DBLE(1.0/3.0)`;
the formula here describes semantics, not a replacement numerical implementation.
The reference coefficient and compiled real-case coefficient are 0; the applied
scale is 1. Packing a zero coefficient uses the existing 1e-10 near-zero floor.

This coefficient is neither an absolute lattice parameter in angstroms, a strain,
nor an independent multiplicative lattice scale. With zero vacancies it has no
lattice sensitivity. Its 0–1 domain must not be transferred to a new scale target.
No validated independent original-cell lattice-scale bounds were found in the
canonical configuration, reference bounds, or source. D1/D2 and PLAMDA are
optimizer coordinates/steps, not defensible physical lattice bounds.

The current crystallographic application rebuilds group B and requested depth
using the original structure; it does not implement original-cell dilation.
B5 explicitly uses the no-vacancy/unit-applied-dilation subset. Accordingly no
third target, perturbation sweep, ANNEAL/AMOEBA run, refined snapshot or new
checkpoint is created. Existing numerical code and checkpoint formats remain
unchanged. No original-cell/reciprocal/recast rebuild equivalence is claimed.

Reproduce the audit and gradient baseline (new output directory required):

```sh
python -m examples.b7e_gradient_audit --output tmp/b7e_gradient_three_target_refinement
```

The audit writes `audit_summary.md` and `audit_summary.json`. A future separate
original-cell lattice-scale or strain target requires an explicit convention,
baseline, and scientifically justified bounds first. Then validate the full
original direct/reciprocal → ZoneAxisCell indexing/recast → strict physical
MultisliceSimulationCell → original-hkl/internal-HK/sample rebuild before
optimizing. Preserve fractional sites, site-group ownership and requested depth.
Do not substitute a derived-cell scale or a full-period projection.

B7e validation: focused tests **37 passed** (20.02 s); full suite **2820 passed** (291.38 s); `git diff --check` clean. Gradient baseline absolute errors were exactly zero for both YNORM and objective.

## B7e continuation: gradient depth/B refinement

After the dilation audit, the primary refinement mode for this stage is gradient
QCBED with **exactly** `sample.requested_depth`, then
`site_group.Al_fcc.debye_waller_b`, in compiler order. Legacy DILATION remains
vacancy coupling, deferred; no independent lattice scale is introduced.

```sh
python -m examples.b7e_gradient_depth_b_refinement --output tmp/b7e_gradient_depth_b_refinement
# Compatible existing stage journals only:
python -m examples.b7e_gradient_depth_b_refinement --output tmp/b7e_gradient_depth_b_refinement --resume
```

The source FCC Al configuration remains unchanged. The example writes a separate
17-digit gradient experimental body and diagnostic YAML snapshots with gradient
metadata. The canonical baseline is depth840 Å/B0.86 Å², YNORM
13447206.106343817, objective1.1896911365294378. The adapter selects only the
prepared gradient dataset. Its 7281 active weights are exactly1: support fraction
>0.875 determines binary inclusion, and variance is the statistical denominator.
Signal/variance/weights/beam order/scan geometry are fingerprinted and checked.
No conventional optimization or cross-mode ranking is performed.

The fixed-B sparse profile uses900/950/1000/1050/1100 Å. ANNEAL starts at
[1000,.86], seeds2026/2027/2028, MAXTEMP.05, MINTEMP.0001, cooling.9,
MAXITS59 (60 proposal iterations plus source initial/reset evaluations).
The case window is900–1100 Å, separate from legal600–1200 Å; B bounds.5–1.5 Å².
All finite bound-valid trials use the existing physical rebuild and gradient
objective. Out-of-bound proposals are recorded with null objective, not fabricated
finite scores. Bound acceptance is distinct from unexposed optimizer acceptance.

AMOEBA receives the lowest finite ANNEAL result. Existing SCAL/D1/D2(1.05,.95)
remain anchored at canonical [840,.86]. Explicit local PLAMDA magnitudes .05/.25
produce nominal physical steps4.2 Å/.0215 Å²; signs direct each step toward its
bound midpoint. All three physical simplex vertices are measured before handoff,
requiring a distinct depth occurrence count and a changed gradient objective for
B. No clipping/rounding is introduced. The legacy inverse coordinate map can
round an exact endpoint by an ulp; normal strict rejection remains in force.
MXF100 is a source soft budget, not a strict count including initialization.
Optimizer stopping and physical convergence are reported separately. Depth
nearness uses the existing endpoint/neighbor slice-interval rule; B nearness is
explicitly one nominal local B step (.0215 Å²). No window is widened.

Gradient checkpoints reuse the existing engine and completed-prefix replay.
Their problem identity explicitly adds signal mode, PLAMDA, range factors and
site-occurrence provenance. Historical conventional identities remain unchanged.
The example verifies real interrupted/resumed ANNEAL and AMOEBA histories,
best/last/retained states and optimizer metadata exactly. Focused tests also cover
cooperative cancellation before/between evaluations and incompatible simplex
settings on resume. Replay does not reevaluate the completed prefix.

Figures/NPZ use actual globally normalized gradient calculated signals, with
excluded pixels NaN and common display scales across beams. NPZ records ordered
beam IDs, lattice coordinates and original conventional HKL. Per-beam CSV retains
active counts, squared-objective contributions, fractions and changes against
the gradient baseline. There are no fitted per-beam weights/scales. All outputs
are diagnostic artifacts; this does not enable a new GUI or general runner path.

### Measured B7e gradient result

Workflow **PASS**, with a **B-boundary-limited** local result; not physical convergence.
Canonical gradient baseline reproduced exactly (zero YNORM/objective difference).
All three ANNEAL seeds used68 scoring evaluations; seed2026 supplied the lowest finite state.

| State | Requested Å | Effective Å | Occurrences | B Å² | YNORM | Gradient objective |
|---|---:|---:|---:|---:|---:|---:|
| ANNEAL 2026 | 930.364980613540 | 931.269999999996 | 460 | 0.752891188781 | 11741665.264922743663 | 0.943139477890898 |
| ANNEAL 2027 | 943.356755683338 | 943.416999999996 | 466 | 0.743498424017 | 11755258.690354241058 | 0.948602506012085 |
| ANNEAL 2028 | 936.360434477903 | 937.343499999996 | 463 | 0.868394715893 | 12918110.526036154479 | 0.958043853114421 |
| AMOEBA best = retained final | 922.084336397803 | 923.171999999996 | 456 | 0.500000094220 | 9333875.106310537085 | 0.912195417120245 |

Actual local PLAMDA was(-.05,-.25): initial physical steps +4.2 Å/+0.0215 Å². The depth vertex changed460→462 occurrences; the B vertex changed gradient objective0.9431394778908981→0.9460706307848331.
AMOEBA stopped on objective spread after41 iterations/82 scoring calls, including28 bound rejections (54 physical evaluations). Requested depth is interior to900–1100 Å, while B=0.5000000942198216 Å² is adjacent to its0.5 lower bound. No bound/window was widened.
All9 per-beam squared-objective contributions improved from the gradient baseline for both ANNEAL best and local best;0 worsened. The sparse fixed-B profile minimum was950 Å, objective0.9682704100351635; final states came from optimization, not profile selection.
Both real interrupted/resumed stages produced checkpoint dictionaries exactly equal to their uninterrupted counterparts. Focused cancellation tests preserve completed history before/after interruption. Both refined YAML snapshots and the exported gradient body reload with exact parameters and signal/variance/weights.
All three PNG/NPZ comparisons are under `tmp/b7e_gradient_depth_b_refinement/figures/` as `gradient_baseline_masked_pattern_match`, `gradient_anneal_best_masked_pattern_match`, and `gradient_local_best_masked_pattern_match`. All masks contain7281 points; outside-mask arrays are NaN, ordered beam IDs and active residuals are exact.
Validation: **37 focused tests passed** (26.62 s), **2829 full-suite tests passed** (322.66 s), `git diff --check` clean.
Next scientific review: inspect residuals and justify the B lower bound before extending searches or enabling another parameter family. Independent lattice scale remains deferred; occupancy/absorption are not enabled. Historical conventional scores are not directly comparable to gradient scores.
