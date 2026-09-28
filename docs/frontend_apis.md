# Stage 6A: frontend-independent configuration APIs

These APIs use the existing `QcbedConfig` dataclasses. They add no GUI toolkit,
widget state, alternate scientific schema, background worker or cancellation API.
The compiler, numerical driver and legacy workflows remain authoritative.

## Structured validation

```python
from qcbed import validate_config, validate_fields, ValidationPhase
from qcbed.configuration import MicroscopeConfig

field_issues = validate_fields(
    MicroscopeConfig, {'energy_kv': 200.0, 'zone_axis': (0, 0, 1)},
    section='microscope', object_path=('microscope',),
)
issues = validate_config(config)
file_issues = validate_config(
    config, phases=(ValidationPhase.FILES,), base_directory=project_directory,
    resume=False,
)
export_issues = validate_config(config, phases=(ValidationPhase.EXPORT,))
```

`validate_fields` accepts constructor keywords with already-typed nested objects;
it checks declared scalar/sequence types, unknown fields and missing required
fields before accepting immutable objects, then reuses the constructor for local
domain checks. Multi-field constructor failures have an object-level location. It
does not coerce strings or
introduce a dictionary-based configuration frontend. Existing constructors retain
their domain validation and exception behavior.

`validate_config` returns a tuple of frozen `ValidationIssue` objects with `code`,
`severity`, `phase`, `location`, `message` and `related_locations`. A location has
`section`, `object_path`, `field` and optional `occurrence`. Paths use field names,
`Id` selectors and, for refinement entries, `ParameterTarget` identities. Occurrence
is only a diagnostic disambiguator for duplicate records, never a refinement index.
`section`, `object_path` and `field` are also available directly on an issue.

The phases are:

| Phase | Purpose |
| --- | --- |
| `FIELD` | Type/shape checks; also available for uncommitted constructor input |
| `COMPILE` | Stable IDs, species references, target references, linked masks, metadata consistency and existing compiler checks |
| `EXECUTION` | Mapped starting-state support, singleton scans, unsupported gas/vacancy/free-surface states, optimizer counts and singular AMOEBA range factors |
| `FILES` | Explicit file loading, output readiness and canonical resume identity checks; requires an absolute base directory |
| `EXPORT` | Prospective legacy representability and explicit omission warnings; defaults to a separate bundle |

Default validation includes FIELD, COMPILE and EXECUTION and does no filesystem I/O.
Dependent phases first perform field and compilation checks and return blockers if
those prerequisites fail. Field-only checks do not compile. File checks read data
but never create directories or outputs. Export checks render records in memory;
`export_bundle=False` requires verified `filename_limits` and checks the configured
standalone filenames.

Duplicate-ID diagnostics include both stable paths and physical occurrences.
Unknown atom species and refinement targets identify the offending object; linked-C
errors identify the binding. Compiler domain errors without an explicit location
retain a structured `compile_error` at configuration scope with the original message.
Consumers should use codes/locations rather than parse exception text. This is a
preflight, not proof that every optimizer trial will succeed, nor a reservation of
filesystem destinations. Existing compiler/evaluator/runner exceptions remain intact.

## Available physical targets

```python
from qcbed import available_refinement_targets

for entry in available_refinement_targets(config):
    print(entry.target, entry.starting_value, entry.selected, entry.active)
```

Accepts `QcbedConfig` or `CompiledConfig` and returns frozen `RefinementTargetInfo`
entries. Each includes:

- `target` (`ParameterTarget`) and `category`;
- `references`: identifier-position/stable-object-path pairs;
- `stable_ids`: referenced object IDs, excluding axis/component tokens;
- a label, units and the **raw configured starting value**;
- `selected`: a binding exists;
- `active`: that binding survives category filtering;
- `available` and a reason when independent selection is unavailable/unsupported.

Enumeration uses the compiler's extracted `_target_map`; no frontend maintains a
second enumeration rule set. Display ordering follows the existing source category
order and concrete target indices, without changing compilation or packing. Starting
values are not folded, wrapped or replaced with packed near-zero values. Linked
absorption C, G-atom coordinates, unsupported vacancy refinement and voltage with
absorption remain visible with reasons. A target's availability does not certify
whole-experiment readiness.

The catalogue requires a compiler-valid configuration. Incomplete edit snapshots
should first be diagnosed; the UI may retain its last valid catalogue while editing.
Counts remain dynamic, including configurations above 1,000 active targets. The
existing driver limit is reported by execution readiness, not imposed on editing or
compilation. No fixed reflection, atom, region or thermal-target slots are added.

## Pure stable-ID edits

```python
from qcbed import Id, update_object, reorder_objects, duplicate_object

region = ('sample', 'regions', Id('film'))
changed = update_object(config, region, depth_angstrom=120.0)
reordered = reorder_objects(changed, ('sample', 'regions'), ('substrate', 'film'))
copy = duplicate_object(changed, region, copy_refinement=True)
```

`get_object` resolves paths of field names and `Id(value)` selectors; collection
indices are not edit identities. `Id(value, 'species')` addresses thermal/absorption
entries associated with a species.

| Helper | Behavior |
| --- | --- |
| `add_object` | Append an ID-bearing object, or insert at explicit `before_id`; reject duplicate IDs and globally repeated stored-slice IDs |
| `replace_object` / `update_object` | Rebuild nested immutable snapshots; retain the selected object's ID |
| `reorder_objects` | Require an exact permutation of IDs; explicitly changes scientific order |
| `deletion_dependencies` | Report bindings, specific bounds, atom assignments and thermal/absorption references outside the selected subtree |
| `delete_object` | Raise `DependencyError` carrying those dependencies until the caller explicitly resolves them |
| `duplicate_object` | Append a copy with deterministic fresh `-copy-N` IDs throughout the owned subtree |
| `set_refinement` | Enable/disable or update a binding by `ParameterTarget`, never by vector position |

Duplication leaves shared species references shared. By default it copies neither
bindings nor bounds. Explicit `copy_refinement=True` remaps affected bindings and
specific bounds to the new IDs; unrelated targets remain untouched. Disabling a
binding does not remove its specific bound or change `active_categories`. New bindings
require explicit D1/D2 and PLAMDA. Existing signed/zero values remain preserved.
Independent selection of catalogue-unavailable targets is refused.

These operations preserve dataclass field types, not necessarily all cross-object
completeness after each edit. For example, adding a species requires explicit thermal
values for every region and absorption coefficients when enabled; the helper does
not guess them. `validate_config` identifies the incomplete relationships before
compilation/execution. No hidden deletion cascade, species reassignment, ID rename,
parameter defaulting or scientific sorting occurs. ID renaming is outside this stage.
Undo/redo may retain old snapshots; no UI state is implemented here.

## Paths, Save and Save As

```python
from qcbed import resolve_config_paths, save_yaml_config, load_yaml_config

resolved = resolve_config_paths(config, base_directory=project_directory)
save_yaml_config(resolved, 'experiment.yml', base_directory=project_directory)
save_yaml_config(edited, 'experiment.yml', base_directory=project_directory,
                 overwrite=True)
```

`base_directory` is mandatory and must be absolute. Both relative configured paths
and a relative save destination resolve against that base, never incidental CWD.
The existing YAML loader delegates to the same resolver using its source directory.
`run.output_filename` remains independently relative to that directory, not to the
run output directory. Legacy raw filename padding is retained consistently.

Saving emits resolved absolute paths and round-trip numeric precision, preserving
IDs, enum values, ordered records, duplicates and settings. It returns the resolved
snapshot. Thus an already-resolved configuration satisfies exact
`load_yaml_config(saved) == config`. For a relative-path input, equality is to
`resolve_config_paths(config, base_directory=...)`, because the unchanged loader
returns absolute paths. Frontends should normalize a document's paths once when
establishing its base.

Save As uses exclusive atomic publication and refuses an existing destination.
`overwrite=True` stages a file beside the destination, flushes it, and uses atomic
replacement; a failure before publication leaves the prior file intact. Staging files
are cleaned up. Protected source/reference files and the configured experimental
file cannot be overwritten. Save does not create parent directories. No preservation
of YAML comments, key formatting, or original relative path spelling is promised.
The older migration-only `write_yaml_config` remains a compatibility API; new
frontends should use `save_yaml_config` with an explicit base.

## Canonical execution and immutable snapshots

```python
from qcbed import prepare_run_snapshot, run_snapshot, run_config

snapshot = prepare_run_snapshot(config, base_directory=project_directory)
# Later edits produce a different config; they cannot change snapshot.
result = run_snapshot(snapshot)
# Convenience equivalent:
result = run_config(config, base_directory=project_directory, resume=False)
```

`RunSnapshot` is frozen and contains a compiled immutable configuration plus the
existing read-only experimental arrays. Preparing it resolves paths, compiles and
loads data but does not run the objective. Editing either configuration objects or
the source experimental file afterward does not alter that snapshot. Running it
uses the captured data, not a fresh file read.

The YAML runner is a thin loader over `run_config`. Both entry points use the exact
same existing canonical identity format/version/digest fields and call the same
`run_full_refinement` with identical output, callback and resume behavior. Tests
compare the extracted identity against the pre-extraction formula and exercise
canonical → YAML → canonical resume. No temporary YAML is used by canonical runs,
and PyYAML is not required to prepare/run a Python configuration.

Resume with a newly prepared configuration/data snapshot still checks all existing
fingerprints. A later edit neither mutates an active run nor grants compatibility
with its checkpoint. Existing legacy/canonical checkpoint separation is unchanged.
Existing progress hooks are retained; background scheduling and cancellation are
explicitly deferred. No optimizer callbacks, evaluation counts or numerical kernels
were changed by this stage.

## Stage 6A validation

On Apple Silicon using `qcbed-python`, **138 focused tests passed in 7.69 s** and
**all 1906 tests passed in 28.16 s**. Twenty-seven new regression cases cover
structured field/duplicate/reference diagnostics, dynamic target counts and the
existing 1,000-parameter execution guard, scoped editing/copy/deletion, atomic save
failures, path-base independence, YAML round trips, snapshot isolation, and exact
canonical/YAML resume equivalence. Existing configuration, migration, numerical and
runner regressions pass. The production CLI help and `git diff --check` also passed.
