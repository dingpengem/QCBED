# Modern editable YAML: A/B/C inputs

`schema_version: 2` is a **document format**, decoded into the same immutable
`QcbedConfig` used by Python and the GUI. It does not change the typed model,
compiler, canonical run identity, numerical ordering, or optimizer behavior.
Legacy/pre-sliced YAML and verbose modern version-1 documents still load through
the existing strict decoder. Saving a modern object writes the new format; saving
a legacy object preserves its existing format. No conversion guesses missing
physical inputs.

The [FCC Al teaching project](../examples/fcc_al.yml) starts from the real raw
pattern and CIF. Open it in the GUI and complete A before running or using the CLI.
It uses 1000 Å, one stored slice, AMOEBA and an explicit two-evaluation teaching
limit. See [the study](EXAMPLE_FCC_AL.md) for all example choices and bounds.

A raw document uses `data_prep.pattern` (plus optional `aperture` and `disk_profile`).
After preparation, Save writes a `data_prep.prepared_dataset` descriptor reference.
The verbose `fcc_al_grouped_gui.yml` remains a historical compatibility/replay input:
it uses a previously generated gradient body, 840 Å and two slices. It is not the
same physical/tutorial configuration as the current minimal raw example.

## A: prepared input reference

`data_prep.prepared_dataset` references a small `qcbed-prepared-input-v1` JSON
input descriptor containing:

- the existing experimental body path, header mode and signal mode;
- recorded scan geometry/dimensions and ordered stable beam IDs/(n1,n2);
- the explicit experimental-to-scan basis convention;
- when present, the companion preparation-metadata path and its SHA256.

The descriptor contains no intensity/variance/weight arrays or compiled structure.
Its metadata link retains the existing Patch A source pattern/manifest, original
support source, manual exclusions, and manual/confirmed zone-axis result. It uses
the existing `.dat.prep.json` provenance; no second mask or marker authority is
introduced. Sigma remains the uncertainty in the referenced experimental body;
unknown raw sigma/noise sources are not invented. Existing overview manifests
retain raw-source provenance where it is available.

Save publishes an immutable content-addressed descriptor in `.qcbed-data` beside
the destination YAML. A/scan edits produce a new descriptor. Masking still produces
the existing immutable weight-only body/provenance snapshot. Save As retains exact
input references and canonical values. Keep these assets with the project; YAML
alone is not an export bundle. Generated assets remain ignored by Git unless
explicitly selected for distribution. The teaching example generates its descriptor during preparation; no machine-local
prepared descriptor is shipped as a portable input.

YAML paths resolve against the YAML directory; descriptor paths resolve against
the descriptor directory. GUI saves resolved absolute references, as before.
There is no environment substitution or incidental-CWD interpretation. A changed
provenance checksum or mismatched signal mode is rejected. User validation booleans,
hover state and current widget selections are not saved.

## B: physical model

`atomic_model` contains `structure` (path string or path/format mapping), original
`zone_axis`, conventional `cbed_basis`, and ordered `site_groups`. `sites` means
**original structure site indices**, never recast occurrences. Group IDs, B values,
occupancies and explicit absorption coefficients survive unchanged. Explicit site
groups require `occupancy` as well as ID, element, sites and B. Optional
`thermal_factors` and `absorption` retain pre-group modern input compatibility.

No internal HK, zone transform, recast lattice, simulation-cell atom IDs, effective
depth or generated slice counts are written to editable YAML.

## C and D: simulation, refinement and execution inputs

`qcbed` contains the title, voltage (kV), signal mode, requested depth (Å), explicit
`numerics`, `workflow`, optional optimizer overrides, dynamic refinement controls,
execution settings and output settings. Output requires only `directory`;
`final_out.dat` is derived within it. An explicit legacy/custom `filename` remains readable for compatibility; the modern
GUI exposes only the directory.
Snapshot/checkpoint intervals and execution budgets remain optional user settings.

`refinement.targets` is an ordered mapping from existing catalogue identities to
`range_factors: [D1, D2]` and `step_size` (unchanged PLAMDA). Presence selects the
target; omission disables it. `refinement.bounds` is a separately ordered mapping
from target identity to `[lower, upper]`, preserving inactive bounds and arbitrary
binding/bound declaration order. Numerical order still comes from the compiler.

Examples: `sample.requested_depth` and
`site_group.<actual stable group ID>.debye_waller_b`. Group names/counts are never
hard-coded. Names resolve through the public target catalogue; unknown IDs are
rejected. Compatibility targets without a physical alias use the catalogue's
category plus stable identifiers. This is not a promise of new runtime capabilities.
Category absolute/relative bounds and category filters remain optional advanced
inputs. Bounds, step sizes and range factors are not normalized or folded.

Version-2 optimizer defaults (omitted on Save when equal):

| YAML control | Existing control | Default |
|---|---|---|
| cycles | loops | 1 |
| stage_budget | maxfsimp | 20 |
| coordinate_search_budget | maxfpara | 3 |
| objective_tolerance | ftol | 1e-6 |
| stop_objective | chimin | 1e-6 |
| maximum_temperature | maxtemp | 100 |
| minimum_temperature | mintemp | 1 |
| cooling_rate | coolrate | 0.9 |
| seed | seed | 2026 |
| coordinate_search_option | paraopt | 0 |
| verbosity | loud | 0 |

`stage_budget` retains the existing algorithm-specific iteration/soft-budget
semantics; it is **not** a new strict objective-evaluation limit. Compatibility
coordinate-search controls do not enable unsupported modern optimizers.
Workflow maps the existing ANNEAL/AMOEBA or compatibility method names without
changing algorithms. Other defaults: title `QCBED`, serial execution, masked
objective, reject invalid trials/bounds, IREF -1, no active targets, and existing
RunConfig defaults (snapshot 100, checkpoint 25, 1000 evaluations, 3600 seconds).
These are format defaults, not inferred FCC Al parameters. Changing their meaning
requires a schema-version decision. Non-default values are always retained.

## Compile, snapshots and resume

User YAML = persistent inputs. The unchanged compiler produces the derived
physical state, internal beam mapping, effective depth, target layout and bounds.
The unchanged immutable `RunSnapshot` carries compiled configuration plus copied
experimental arrays; checkpoint/results retain execution state and provenance.
The descriptor's filename and frontend schema version do not enter the typed
configuration or checkpoint identity. Equivalent old/new documents therefore
have identical canonical/grouped checkpoint identities. A real mask/geometry/model
change still changes the relevant identity and invalidates incompatible Resume.

The GUI continues editing `QcbedConfig` directly; only its load/save frontend
changes. A → B → C validation remains transient and ordered. D validates output separately. Fresh modern GUI runs use unique managed child
directories; Resume requires a separate compatible checkpoint. The subsequent [packaging closure](../packaging/macos/VALIDATION.md) verifies
this frontend in the rebuilt macOS app without a refinement campaign.

## Historical schema-migration validation

79 focused schema/configuration tests and 98 GUI tests passed; the normal full
suite passed **2860 tests** (340.26 seconds). Fixed small-fixture evaluations
preserved compiled state, phase gratings, intensities, objective/YNORM and checkpoint
identity exactly. No optimization campaign or app rebuild was performed.
`git diff --check` passed.
