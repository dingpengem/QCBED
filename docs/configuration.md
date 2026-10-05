# Canonical Python configuration — Stages 1–5

[Documentation index](README.md) · [Current architecture](ARCHITECTURE.md) ·
[Researcher manual](user_manual.md)

`qcbed.configuration` defines immutable, typed user configuration objects.
`qcbed.compile_config(config)` translates them into existing internal objects
without reading files, loading experimental arrays, constructing physics, or
running an optimizer. Experimental I/O is a separate explicit operation. No
GUI dependency, exporter, cache, or replacement legacy entrypoint is introduced.
Stage 3 adds an optional [strict YAML frontend](yaml_configuration.md), while
these typed objects and the compiler remain format-independent.

```python
from qcbed import compile_config
from reference_configuration import reference_configuration  # examples directory

raw = reference_configuration()
resolved = compile_config(raw)
state = resolved.state
controls = resolved.optimizer
vector = resolved.packed_starting_vector
```

`examples/reference_configuration.py` independently constructs the current Al
full-refinement configuration with Python objects. It does not parse FOR007,
inspect the experimental header, or import the editable runner. It includes the
original 18 bindings, the current depth/scan category selection, ANSI settings,
batch 64 and the explicit production bounds. The reference's duplicate reflection
is preserved. Paths remain unresolved strings during compilation. The default
experimental specification selects a headerless body under `tmp/`;
`reference_configuration(legacy_header=True)` retains the legacy file metadata.
This example does not change `examples/run_full_refinement.py` or its resume path.

## Raw and resolved objects

`QcbedConfig` holds schema version 1, microscope/sample/scan/experimental settings,
refinement controls, execution settings, run settings, and optional legacy-export
metadata. Nested configuration uses frozen dataclasses and enums. Collection
inputs are copied to tuples; finite real values and strict integer types are
checked without accepting booleans as numbers. Cross-references are checked by
the compiler. Existing `MeshConfig`, `SliceGeometry` and `ForwardExecutionConfig`
are reused. The canonical `configuration.MicroscopeConfig` adds the zone axis;
the existing `qcbed.MicroscopeConfig` is unchanged.

`CompiledConfig` is a distinct frozen object containing:

- The original raw configuration and external experimental-data specification.
- `RefinementState`, source-ordered `RefinementParameterLayout`, and
  `LegacyRefinementConfig` (including existing execution settings).
- Resolved bounds, bound rules and `FullRefinementConfig` runner settings.
- `range_factors` and `plamda` tuple properties in source parameter order.
- `packed_starting_vector`, obtained through the existing immutable packing API.

The output directory and optional export metadata remain in `configuration.run`
and `configuration.legacy_export`: the compiler never creates a directory or
exports a file. Experimental format/header fields are specifications only.

## Stable references and ordering

Species, regions, globally named stored slices, and per-slice atoms have stable
IDs. Reflection IDs are local to a region. Duplicate IDs in their scope are
errors; duplicate reflection H/K and beam H/K entries are allowed. Regions are
flattened into stored slices in declared region/slice order. Species, atoms,
beams and reflection records preserve their explicitly ordered sequences.
Debye–Waller and absorption entries refer to species IDs, so their declaration
order does not affect species association.

`ParameterTarget(category, identifiers)` resolves these stable identifiers:

| Category enum value | Identifiers |
| --- | --- |
| `voltage`, `dilation` | `()` |
| `depth`, `constant_vacancy` | `(region_id,)` |
| `vacancy_coefficient` | `(region_id, "a" or "b")` |
| `absorption` | `(species_id, "a" or "b" or "c")` |
| `debye_waller` | `(region_id, species_id)` |
| `atom_coordinate` | `(slice_id, atom_id, "x" or "y")` |
| `scan_coordinate` | `("A" or "B" or "C", "first" or "second")` |
| `reflection_multiplier` | `(region_id, reflection_id, "elastic" or "absorptive")` |

Use `ParameterCategory` enum members, not arbitrary strings. Each
`ParameterBinding` supplies a target, a `(D1,D2)` pair and PLAMDA. Duplicate or
unknown targets are rejected. Bindings set the existing masks; the existing
layout builder determines vector order. Reordering bindings therefore leaves
layout, packed vector, D1/D2 and PLAMDA unchanged. Optional category selection
uses the existing subset function after this full source-order assembly.

Changing the explicitly ordered physical species/region/slice/atom/reflection
lists is different from changing binding declaration order: those lists define
the scientific source order and are deliberately not sorted by ID or H/K.

## Numerical and policy boundaries

Configuration values are starting inputs, not mapped or packed values. Negative
or near-zero values are not folded, clipped, wrapped or replaced in the raw
configuration or compiled baseline. Calling `packed_starting_vector` invokes
the existing ARRANGEPARS projection. Evaluating a trial invokes the existing
UPDPARS mapping. Neither operation is duplicated in this layer.

Linked absorption C retains its configured input C value and compiles to the
existing `-1` C mask. Its B-to-C assignment still happens in the existing mapper.
A binding that independently refines linked C is rejected. S/G record semantics
remain explicit, and G atoms cannot acquire coordinate-refinement masks.

Vacancy profiles, free-surface settings and gas records remain representable.
Compiling them does not claim that the current evaluator supports them. Existing
runtime rejection is unchanged. Singleton scans remain representable but are
still rejected by scan generation. There is no new physical-domain workaround.

Units remain kV, angstroms, angstrom squared and lattice degrees. Atom positions
and slice translations are fractional. Scan corner components are coefficients
in the supplied scan basis, not angular tilts. Third basis components are retained.
`gradient_option` remains the original integer, not a Boolean conversion.

IREF and legacy optimizer control names/semantics are retained, including shared
MAXFSIMP. Range factors remain search/scaling inputs, never hard bounds. Optional
hard bounds, masked objective, numerical rejection, seed, execution policy and
runner budgets remain separately selected Python behavior. Bound resolution
delegates to the existing mapped-start resolver on a separate state; specific,
category-absolute, category-relative and unbounded precedence is unchanged.

Legacy raw strings can be retained as metadata. Absorption's raw flag must agree
with the explicit enabled flag because the existing model uses its literal `YE`
prefix. Raw title/filename metadata must agree with the corresponding configured
values after the existing trailing-space convention. DSMAX is retained for a
future exporter and has no newly introduced runtime effect.

## Validation

`tests/test_configuration.py` compares the direct Al configuration against the
unchanged parser/subset/current-example path, including complete physical state
and raw records, all 18 parameters and the selected seven, scan arrays, ordered
duplicate reflections/beams, optimizer and runner controls, bounds and packing.
Additional tests cover multiple species/regions/slices, every parameter category,
linked C, unmodified starting signs, invalid references/types, empty layouts,
deep immutability, filesystem-I/O prohibition and existing unsupported-state
rejection. A small two-slice, four-tilt forward comparison checks exact complex
amplitudes, intensities, objective and YNORM against the legacy-derived state.

This is configuration-to-existing-Python equivalence, not new compiled-Fortran
equivalence. Legacy export and runner migration remain later stages.

Validation on the Apple Silicon `qcbed-python` environment: **20 focused tests
passed in 0.39 s; all 1772 tests passed in 20.72 s**. The small forward comparison
has zero complex-amplitude/intensity difference and exact objective/YNORM equality.

## Dynamic-count and extensibility audit

`CompiledConfig.parameter_count` is exactly `len(layout.parameters)`. There is
no 18-entry vector or eight-reflection capacity in the compiler or execution
pipeline. Those numbers describe the Al fixture only. The reference example now
generates reflection bindings using the actual reflection-list length. Generic
reflection records may contain any number of independently enabled elastic and
absorptive targets, including distinct records sharing the same H/K.

Three concepts remain separate:

1. **Category:** `ParameterCategory` and immutable `PARAMETER_SEMANTICS` in
   `refinement_categories.py` define target arity, existing source-order group
   and selective mapped-value writeback behavior. Registry size is not vector size.
2. **Instance:** canonical `ParameterTarget` identifies concrete stable IDs and
   components; compiled `RefinementParameterSpec.target` contains source indices.
   Debye–Waller instances are `(region, species)` pairs, not reserved vector slots.
3. **Ordering:** the existing `_entries` traversal generates the active layout.
   It remains authoritative; no new sorting changes summation or mapping order.
   Vacancy variants interleave by region within their common order group.

The shared taxonomy removes independent target-arity dictionaries in layout and
bounds validation, and the duplicate category whitelist in subset selection.
Selective UPDPARS writeback now lives in the parameter-mapping layer; driver and
SIMPARA adapters call it without maintaining physical-category lists. Its original
operation order and in-place update semantics are unchanged. Numerical optimizer
search loops are unchanged and operate on vectors/callbacks sized from N.

| Path | Count/order source |
| --- | --- |
| Compiler masks | Actual regions/species/slices/atoms/reflections and resolved bindings |
| D1/D2 and PLAMDA | Binding lookup in the generated layout's order, then existing subset selection |
| Packing and physical mapping | Existing source `_entries` / layout traversal |
| Bounds | Concrete targets, category rules expanded over the active layout |
| Simplex and optimizer arrays | N from layout/vector length; simplex shape `(N+1,N)` |
| Callback writebacks | Source-ordered mapped entries and shared category semantics |
| Checkpoints/resume | Serialized complete layout, vector and matching labels; mismatched layouts rejected |
| Reports/trajectories | Layout-derived labels and vector components |

An existing limitation is deliberately isolated rather than silently removed:
`LEGACY_OPTIMIZER_PARAMETER_LIMIT = 1000` preserves the current adapter/driver
validation contract. It is not an allocation size and does not limit canonical
compilation; a test compiles 1003 active reflection instances. Thus the compiler
supports arbitrary record counts, while current legacy optimization APIs are not
advertised as unlimited. Changing this guard needs separate validation. The
runner's default 1000-evaluation budget and SIMPARA's numerical 10000 cap have
different meanings and were not altered.

Zero active parameters are supported by the existing evaluation-only route
(`loops <= 0`, no bindings, no category filter). Positive-loop driver refinement
still requires at least one parameter. Explicit category filtering that selects
no parameters retains its existing error instead of inventing a new policy.

### Adding a future physical category

Define its input representation, stable target semantics and supported numerical
domain. Add category metadata and compiler target/mask resolution, including
subset clearing and valid bound targets. Implement source-ordered enumeration,
packing and physical application in the parameter layer; choose selective
writeback explicitly. Existing categories must retain their ordering. Add tests
for physical mapping, ordering, invalid targets and dynamic selection/counts.

These are domain-layer changes, not new optimizer dimensions or new optimizer
branches. No optimizer search loop, checkpoint vector format or report column
count needs rewriting. This is a source-maintained extension boundary, not a
runtime plugin mechanism: merely registering a category does not implement its
physics. Occupancy or additional coordinates are not enabled by this audit.

New dynamic tests exercise 0/1/3/11 reflection records (0/2/6/22 active components),
1003 compiled reflection targets, 2/12/6 active thermal targets across variable
region/species counts, both reference layouts, and all four optimizers with 1 and
5 active instances. Evaluation-only runner tests verify zero- and five-element
checkpoint layouts/vectors, report changes and resume. Existing mapping/optimizer
equivalence tests remain the numerical regression baseline.

Audit validation on Apple Silicon: **228 focused tests passed in 1.58 s; all
1791 tests passed in 22.14 s**. `git diff --check` passed. Existing exact numerical
comparisons, including the configuration-to-legacy forward comparison, pass;
no physics formulas or optimizer search algorithms changed.

## Stage 2: configuration-driven experimental data

```python
from pathlib import Path
from qcbed import compile_config, load_experimental_data

resolved = compile_config(reference_configuration())
data = load_experimental_data(resolved, base_directory=Path.cwd())
```

`ScanConfig.scan_columns` / `scan_rows` and
`ExperimentalDataConfig.gradient_option` are authoritative. The explicit
`ExperimentalDataConfig.header` mode selects:

- `"absent"`: the external file contains only intensities, variances and weights.
- `"validate"`: consume the legacy logical NUMP/NUML/GRADIENTOPT header and require
  all three integers to equal configuration. A mismatch raises a contextual
  `ValueError`; the header never overrides configuration.

There is no header autodetection. Relative paths resolve against the required
`base_directory`; absolute paths are used directly. Compilation remains free of
filesystem I/O. `parse_experimental_intensities` retains its legacy signature and
header-driven behavior. Both APIs share the original body parser, including
logical-record continuation and discarded surplus fields on the final physical
record of each READ. Beam-major file order becomes public `(tilt, beam)` arrays;
duplicate beams, zeros, signed weights, scan order and GRADIENT semantics remain
unchanged. Reading does not apply GRADIENT or normalize data.

Run `python examples/reference_configuration.py` once from the development
environment to create `tmp/configuration_reference/in_body.dat` and load it via
the new API. This explicit demonstration validates and removes the complete
logical header from a separate copy, leaves reference data untouched, and refuses
to overwrite an existing destination. Merely importing or constructing the example
performs no I/O. The production legacy runner and FOR007 inputs remain supported.

`tests/test_experimental_config.py` compares all three modes, including continued
headers/body records, mismatch errors, no autodetection, duplicate beams, zeros,
non-Boolean gradient options and immutable arrays. Reference-data tests compare
all 1521 tilts and nine beams, compiled refinement setup, and exact amplitudes,
intensities, YNORM and objective using a shortened depth with the full scan.
These establish equivalence to the existing Python path, not a new Fortran run.

Stage 2 validation on Apple Silicon using `qcbed-python`: **94 focused tests
passed in 2.23 s; all 1807 tests passed in 21.40 s**. Exact comparisons passed
with no amplitude/intensity differences and equal YNORM/objectives. The executable
reference demonstration loaded `(1521, 9)` arrays with GRADIENTOPT=1.
`git diff --check` passed. No YAML Stage 3 work was started.

## Multi-region audit before YAML

Each `RegionConfig.depth_angstrom` compiles to its own `RefinementRegion` depth.
`SampleConfig` has no shared thickness field. Regions and their stored slices
retain declaration order; slice IDs are globally unique. Regions may have
unequal stored-slice counts and unequal thickness patterns.

`refinement_trial._prepare_trial` resets accumulated thickness and pattern
position for each region. It repeats that region's stored pattern until its
own depth is reached or exceeded, then appends the next region's sequence.
The last slice remains a whole occurrence: it is not shortened to fit depth.
An exact boundary includes no extra occurrence; moving just above it adds the
next occurrence. This is existing numerical behavior, not a new depth policy.

A `ParameterTarget(DEPTH, (region_id,))` resolves to exactly one region index.
Bindings may select any subset; vector ordering follows configured region order,
not binding declaration order. Independent target bounds resolve through the
same IDs. Existing packing, mapped-value handling, D1/D2 and PLAMDA retain their
semantics. Fixed region depths remain part of the full physical state and resume
fingerprint even though they have no active vector entry. Callback states,
parameter snapshots and checkpoints preserve region identity through ordered
indices; reports label active depths with those distinct targets.

`tests/test_multiregion_configuration.py` covers one region, two fixed regions,
one/both independently refined depths, three unequal patterns, independent bounds,
exact adjacent-float sequence boundaries, and a 23-region layout. Runner tests
check callbacks, saved physical depths, report labels, completed resume and an
interrupted stage-boundary resume against an uninterrupted run. Changing a fixed
region depth is rejected on resume by the problem fingerprint.

No production architecture or numerical code needed changes. Existing supported
trial restrictions remain: positive region depths, advancing positive slice
thicknesses, and no vacancy/free-surface derived-state hazards or G-atom
randomization. These tests do not expand that physics support or establish a new
compiled-Fortran comparison. YAML Stage 3 remains unstarted.

Multi-region audit validation on Apple Silicon (`qcbed-python`): **73 focused
 tests passed in 2.85 s; all 1820 tests passed in 21.55 s**. Thirteen new regression
cases passed, including exact resumed/uninterrupted state and vector comparisons.
`git diff --check` passed.

## Stage 3: YAML frontend

`load_yaml_config(path)` now constructs these same immutable objects from strict
YAML, validates through the existing compiler, and resolves configured paths
relative to the YAML source directory. See the [schema and GUI mapping reference](yaml_configuration.md)
and [complete example](../examples/qcbed_config.yml). The canonical model, compiler,
dynamic layout and numerical behavior are unchanged. Validation: **107 focused
and all 1859 tests passed**; no main-runner migration or FOR007 export was added.

## Stage 4: canonical production entry point

`python -m qcbed.run --config examples/qcbed_config.yml` runs the compiled
configuration through the existing production runner with no FOR007 input.
`--resume` uses the configured output directory and explicitly marked canonical
checkpoint identity. See the [recommended workflow](yaml_runs.md) for path rules,
outputs, resume guarantees and the retained legacy compatibility workflow.

## Stage 5: optional compatibility tools

`qcbed.legacy_config` provides explicit one-time import, YAML writing and optional
FOR007 export. `python -m qcbed.legacy import ...` requires explicit Python policy
choices; `python -m qcbed.legacy export --config ... --destination ... --bundle`
creates separate legacy input/data files and an omissions report. See
[legacy migration and export](legacy_migration.md). Canonical configuration,
compilation, YAML loading and production execution remain unchanged.

## Crystallographic source compilation

`compile_config(config, structure=AtomicStructureModel(...))` is the pure B4b
entry point. `compile_config_from_sources(config)` is the explicit file-reading
frontend for document-resolved paths. The result includes a crystallography bundle,
resolved pre-sliced execution configuration and original source configuration.
See [B4b schema and restrictions](crystallography.md#b4b-production-configuration-integration).
No legacy field is reinterpreted as an original-cell hkl basis.

## B7a physical site groups

Opt-in `crystallography.site_groups` owns original physical sites and their B,
occupancy and independent absorption metadata. The existing target catalogue and
compiler resolve group B aliases into the existing dynamic layout; no simulation
atom index is user input. Recast/slab copies retain original index and group ID.
See [B7 bring-up](b7_refinement_bringup.md) for the implemented schema, one-group-per-
element limitation, explicit bounds, exact baseline tests and version2 identity.
Existing non-group semantics and checkpoint identities remain unchanged.
Group optimizer/resume execution and GUI editing are deferred.

Modern editable documents use the [A/B/C frontend](modern_yaml.md). This is a
serialization boundary over the existing typed model, not a second configuration
model. Compilers and immutable run snapshots retain derived runtime state; GUI
Save writes only user inputs and a preparation-descriptor reference. Equivalent
verbose and compact documents preserve checkpoint identity.
