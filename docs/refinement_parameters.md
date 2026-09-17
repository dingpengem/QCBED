# Refinement parameter mapping and trial evaluation

## Evidence and scope

This implementation translates scalar parameter packing/unpacking and an
explicitly supported trial-evaluation path. It implements no optimizer and
claims no historical executable equivalence. Readable executable statements
take precedence over names, comments and the related-version private predecessor implementation.

| Source | Evidence |
| --- | --- |
| private predecessor implementation| ARRANGEPARS selection, ordering and near-zero replacement |
| private predecessor implementation| UPDPARS assignments, linked parameters and dependent updates |
| private predecessor implementation| simplex SCAL setup, conversion and vertex calls |
| private predecessor implementation| ONEIT inverse conversion, three update/sequence passes, vector restoration |
| private predecessor implementation| ONEITPARA physical-vector update/sequence passes |
| private predecessor implementationprivate predecessor implementation| propagation conventions, optional GRADIENT and objective return |
| private predecessor implementation| SFMAP2 exact/reflection-symmetry factor matching |
| private predecessor implementation| readable GRADIENT body despite .o extension |
| private predecessor implementation| masks, reflection records, IREF, PLAMDA and D1/D2 input |

All scalar categories are represented and mapped. Trial evaluation supports
S atoms, constant zero vacancy and no free-surface correction. As explicitly
agreed for this task, unresolved derived-state cases raise errors rather than
silently choosing an interpretation. The existing forward/scaling/objective
modules and parser semantics are unchanged.

## Active vector order and targets

NPAR starts at zero. Each selected scalar increments it once; YY appends two
in order. Inactive values are omitted entirely, not represented by zero slots.
The Python vector and all physical target indices are zero-based. Species
identifiers retained on AtomSite remain one-based, as in the existing API.

| Order | Category / target | Selection | Packing | Applying a trial |
| --- | --- | --- | --- | --- |
| 1 | voltage KV | KVOPT=1 | KV | abs(value) |
| 2 | depth DEPTH(region) | TMASK=1 | depth | abs(value) |
| 3 | vacancy VA,VB per region | CONCFLAG>0, respective mask=1 | abs(value), VA before VB | abs(value) |
| 3 alternative | CONSTVACCONC(region) | CONCFLAG=0, mask=1 | value | cap >=1 to 1, then abs if <=0 |
| 4 | DILATION | MASKDILATION=1 | abs, reciprocal if >1 | same folding |
| 5 | INELA,INELB,INELC per species | raw absorption starts YE, mask=1 | abs(value) | abs(value) |
| 6 | DW(region,species) | MASKDW=1 | value (no abs) | abs(value) |
| 7 | XPS(slice,atom,x then y) | coordinate mask=1 | signed fractional coordinate | signed assignment |
| 8 | scan A,B,C, component 1 then 2 | YY both, YN first, NY second, NN neither | signed scan-basis coefficient | signed assignment |
| 9 | SFFAC,ABSFAC per region/reflection | same character-mask convention | signed multiplier | signed assignment |

Outer ordering within each category is exactly the row ordering shown above.
After packing **every** active value with abs(value)<1e-10 becomes positive
1e-10; the strict boundary is tested. Applying a physical vector does not
perform that packing step again. Round trips are therefore not identities
for sign-folded or near-zero values. Constant vacancy -2 becomes +2 after
the source's ordered operations; it is not a conventional [0,1] clamp of the
parameter. The separate derived VAC array is clipped later in Fortran.

INELC mask=-1 omits C from the vector and copies the current B into C during
unpacking, even though C is not independently active. Other inactive scalar
values remain unchanged. Atom starting occupancy is not directly packed;
vacancy parameters control derived occupancy. Energy is kV, depths angstrom,
DW angstrom squared; atom positions are direct-cell fractions, scan values
are dimensionless scan-basis coefficients, and SF factors are dimensionless.
There are no direct lattice-a/b/gamma, slice-shift, background or independent
YNORM parameters in ARRANGEPARS. IREF controls the existing YNORM calculation.

SF entries retain record identity and duplicate H/K values. They are elastic
and absorptive multipliers, not independent real/imaginary components:
SFMAP2 multiplies both SFRE/SFIM by SFFAC and AFRE/AFIM by ABSFAC. UPDPARS
does not copy factors among reflection records. Symmetry acts later during
grid construction. All exact H/K matches are searched before symmetry matches;
the first matching record wins, so duplicate slots need not affect the output.
NONE, INVE, RECT, SQUA and HEXA follow the literal source transforms. Tests
distinguish exact-match priority from an earlier symmetry-related record.

## SCAL and range semantics

SCAL is an array, not a called routine. In simplex branches it is constructed
once from the packed starting values p0 and input factors D1,D2:

```text
a = 1.0 / (D2-D1) / p0
b = D1 * p0
q = (physical_best - b) * a
physical_trial = q / a + b
```

ONEIT receives q, saves it, inverse-transforms it for UPDPARS, then restores
the optimizer vector. ONEITPARA receives physical values directly. The new
API names these modes `legacy_simplex` and `physical`. Both use a stable
baseline; repeated evaluations must pass that baseline, not the preceding
trial's state. No SCAL conversion is applied to physical-mode trials.

D1,D2 are stored as `range_factors` metadata and used only for the affine
conversion, not clipped bounds. The reference pairs are (1.05,0.95), with
the apparently reversed order preserved. ARRANGEPARS itself does not apply
them. New Python validation rejects singular SCAL rather than allowing an
undefined coordinate conversion. PLAMDA and optimizer controls are consumed
by the refinement reader only to reach the D1/D2 records; no optimizer step
or convergence behavior is implemented. Missing/uncertain optimizer linkage
does not affect the readable affine formula.

## Immutable Python data flow

```python
state = parse_refinement_state("reference_data/FOR007.DAT")
layout = build_refinement_parameter_layout(state)
physical_vector = pack_refinement_parameters(state)
mapped = apply_refinement_parameters(state, physical_vector)
result = evaluate_refinement_trial(state, physical_vector, experimental)
```

RefinementState holds the existing model/scan plus immutable region, absorption,
reflection and range data. RefinementParameterLayout contains frozen specs:
vector_index, category, target, initial packed value, and optional range_factors.
Arrays use parsed counts, with no legacy 1000/2000 capacity. The dedicated
refinement reader retains validated prefix fields through an observer around
the existing READ cursors, then reads the later masks/ranges. It does not alter
the first-stage FOR007 or scan parser APIs. Unknown masks/symmetries are rejected
as Python validation policy. No optimizer is selected or executed by parsing.

Packing returns an immutable float64 vector. Applying requires a finite
float64 NumPy vector of exactly the layout length, and returns a new state.
Zero active parameters yield an empty vector: scalar mapping/evaluation does
not divide by NPAR. This is not a promise that legacy optimizers accept NPAR=0.
Packing's Fortran mutation of shared DILATION and UPDPARS' mutation of PT are
represented by values in the returned vector/state, never by baseline mutation.
There is no CHISMIN/PARMIN best-state bookkeeping, UPDATE output, or GD lookup.

The evaluator maps the vector, validates supported derived physics, reconstructs
the slice sequence, rebuilds gratings, generates scan tilts, and calls the
existing serial slice simulator with existing per-tilt propagators. Unity
gratings reuse preparation; nonunit reflection factors are applied to the
existing elastic/absorptive kernel outputs before the existing grating builder.
The base prepared sample is used for geometry/propagators only when factors
are nonunit; its unity metadata is not relabeled or returned as a nonunit sample.
Complex128 amplitudes and float64 intensities retain `(tilt,beam)` order.

For the supported subset, VAC=0 and APPLIEDDILATION=1. Each region cycles its
stored slices until accumulated thickness is >=DEPTH, using repeated additions
and source order, not a rounded slice count. The three Fortran UPDPARS/MSEQ
passes give identical derived state here; one construction suffices. Wavelength,
interaction constant, grids, DW-weighted factors and gratings are recalculated
from the mapped state using existing numerical APIs. No stale cache is reused.

The discovered readable GRADIENT applies an eight-neighbor mean minus the
center with zero padding, to calculated intensities only when GRADIENTOPT=1.
The wrapper retains both raw calculated and objective-input intensities and
preserves source neighbor-addition order. It uses dynamic storage instead of
the old 250x250 capacity. Other integer flag values follow the source's no-call
branch. This preprocessing does not modify the experimental dataset.

Finally compute_legacy_objective computes/applies YNORM through the existing
scaling APIs and returns the scalar objective. The wrapper does not fit/apply
YNORM twice. RefinementTrialResult exposes objective, scaling, physical_state,
amplitudes and both intensity arrays; all array snapshots are immutable. Same
baseline/vector evaluations are deterministic; failure leaves no partial state.

## Reference layout and validation limits

Both reference FOR007 files have this exact 18-entry layout (indices zero-based):

| Vector indices | Physical target |
| --- | --- |
| 0 | region 0 DEPTH |
| 1,2 | species 0 INELB,INELC |
| 3 | DW(region 0,species 0) |
| 4,5 | A(1),A(2) |
| 6,7 | B(1),B(2) |
| 8,9 | C(1),C(2) |
| 10-17 | region 0 SFFAC records 0-7 |

Reflection records are (1,1),(1,-1),(-1,1),(-1,-1),(2,0),(0,-2),(-2,0),(2,0).
The last duplicate is preserved, not deduplicated. FOR007 starting values are
840,1,1,.86,-.508156,.407641,.407641,.508156,-.407641,-.508156, then eight ones.
Tests verify these independently recorded values and ordering. They do not
claim a successful reference fit: zero variances can still produce NaN in
the unchanged objective, and no historical executable comparison was run.

## Explicit unresolved/rejected trial states

- Nonzero constant vacancy, vacancy profiles and free-surface correction are
  rejected by the evaluator, while scalar mapping remains available. UPDPARS
  uses a region-loop I in RTYPE(I) within a slice loop and reads SLNO on paths
  without an assignment (private predecessor implementation). Neighbor-region vacancy formulas
  and previous applied dilation also require state-order decisions. These
  must be resolved before enabling these derived updates. No correction to
  these equations is silently introduced.
- G records are rejected: legacy RNG/gas behavior is not deterministic from
  the supplied vector alone. No randomization is implemented.
- Voltage changes with absorption enabled are rejected: UPDPARS updates KV,
  WAVEL and SIGMA but not VOLTS, whereas SFMAP2 absorption uses VOLTS. The
  existing forward API uses one coherent voltage. Choosing one interpretation
  here would silently change the legacy state relationship.
- Nonpositive mapped region depths and unsupported existing forward domains
  are rejected before simulation. These are explicit Python domain limits,
  not invented optimization bounds. Singleton scans remain unsupported by
  the existing source-backed scan generator.
- Exact Fortran default-literal rounding, compiler options, shared-state
  side effects in optimization and executable equivalence remain unverified.

Tests cover literal ordering/transformations, multiple regions/species/slices,
near-zero packing, linked/inactive values, exact affine conversion, duplicate
reflections/symmetry, dynamic counts, invalid vectors, deterministic immutable
trials, a manually assembled complex-wave oracle, nonunit origin factors,
GRADIENT's hand stencil and all YNORM branches. No optimizer tests are fabricated.

Validation used exactly
`C:\Users\peng863\AppData\Local\miniconda3\envs\qcbed-python\python.exe`:
**34 focused tests passed in 0.89 s; 1391 complete-suite tests passed in
26.66 s.** The manually assembled complex-wave comparisons match exactly;
the independent nonunit-origin phase check passes at absolute tolerance 1e-15.
The first expanded focused run hit pytest temporary-directory permissions;
the run with required access passed. Python/documentation diff checks and
new-file trailing-whitespace checks passed. Full-tree `git diff --check`
reports only the user's pre-existing private predecessor implementation whitespace at lines 4,101,105.
Protected source/reference data were not modified; no files were staged or
committed. The pre-existing optimizer-architecture document diff was untouched.
