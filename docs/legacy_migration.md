# Optional legacy compatibility and migration

Canonical YAML remains the production input. FOR007 import/export is a separate,
one-time tool; `qcbed.run` never imports these tools or reads FOR007.

```sh
python -m qcbed.run --config experiment.yml
python -m qcbed.run --config experiment.yml --resume
```

## Import once, then run canonically

Legacy files do not define Python execution, seeds, hard bounds, objective policy,
rejection behavior or run budgets. The import command requires explicit choices:

```sh
python -m qcbed.legacy import \
  --for007 reference_data/FOR007.DAT --experimental reference_data/in.dat \
  --yaml migrated.yml --seed 2026 --execution batched --batch-size 64 \
  --objective masked --invalid-trial-policy reject --bounds none \
  --output-directory tmp/migrated_run --output-filename tmp/migrated_out.dat \
  --snapshot-every 25 --checkpoint-every 25 \
  --max-evaluations 150 --max-seconds 1800
python -m qcbed.run --config migrated.yml
```

These example choices are explicit Python policies, not values inferred from
Fortran. Use `--objective legacy --invalid-trial-policy raise` for those legacy
policies. Zero snapshot cadence disables snapshots; zero evaluation/time budgets
mean unlimited. Checkpoint cadence must remain positive. Import preserves the
legacy optimizer method, loop counts, controls, PLAMDA and D1/D2. The seed is an
explicit Python reproducibility setting; it does not promise historical Fortran
random-number sequence equivalence. Hard bounds have no legacy input source;
the CLI requires `--bounds none`, after which bounds can be edited in YAML.

The data path is mandatory so the caller chooses the correct experimental file.
Its logical header supplies NUMP/NUML/GRADIENTOPT during migration. The generated
canonical specification uses `header: validate`, absolute file paths and the
existing external data body. Canonical runs subsequently require only that YAML
and experimental file, not the original FOR007. Output destinations are explicitly
relocated by the supplied run settings. The migrated YAML is created exclusively;
existing YAML, experimental inputs, and reference files are never overwritten.

Python API:

```python
from qcbed.legacy_config import MigrationPolicy, import_legacy_config, write_yaml_config
# policy must explicitly supply seed, execution, RunConfig, objective_mode,
# invalid_trial_policy and BoundsConfig.
raw = import_legacy_config('FOR007.dat', experimental_path='in.dat', policy=policy)
write_yaml_config(raw, 'migrated.yml')
```

Stable IDs are assigned deterministically in source order: `region-0`, `slice-0`,
`species-0`, per-slice `atom-0`, and per-region `reflection-0`. No records are
sorted or deduplicated. Reflection identity is record identity, not H/K identity.
Depth targets and thermal targets resolve to the corresponding region/species.
Linked absorption C retains its input value and link flag; import does not apply
UPDPARS. DSMAX and empty-slice S/G types are recovered by capturing validated READs,
even where intermediate models do not retain them. Unknown trailing records are
rejected rather than discarded. Unsupported physical states remain representable,
but retain existing evaluator rejection behavior.

## Optional export bundle

```sh
python -m qcbed.legacy export --config experiment.yml \
  --destination tmp/legacy_bundle --bundle
```

The destination must not exist. This creates `FOR007.dat`, a separate header-bearing
`in.dat`, and `compatibility.json`. FOR007 uses short reference-style `in.dat` and
`out.dat` names, intended to run with the bundle as working directory. Experimental
values are read via the existing configured reader and written beam-major, preserving
all variances, weights, zeros, signs and duplicate beams. GRADIENT is not applied.
Canonical files remain untouched. The export command does not launch a Fortran
executable and does not create FOR007 runtime dependencies for Python.

```python
from qcbed.legacy_config import export_for007
report = export_for007(raw_or_compiled_config, 'new_bundle', bundle=True)
```

Standalone export requires `header: validate` and explicit verified filename widths
for the intended target executable. For example, **only if those widths have been
verified for your build**:

```sh
python -m qcbed.legacy export --config experiment.yml --destination exported.dat \
  --experimental-filename in.dat --output-filename out.dat --filename-limits 9 10
```

The corresponding API keywords are `experimental_filename`, `output_filename`,
and `filename_limits=(EXPINT_width, BMNAME_width)`. The caller supplies/relocates
the matching header-bearing data file; standalone export does not copy it.
Names exceeding the explicit widths are rejected, never truncated. Bundle names
are proven by reference inputs, but target executable array capacities still need
checking for large user configurations. Related-version private predecessor implementation declarations
are not imposed as Python capacities or claimed as universal binary limits.

## Representation and omissions

The exporter follows readable `private predecessor implementation` input READs and `private predecessor implementation` record order.
Numeric records use round-trip binary64 decimal representations rather than UPDATE's
short display precision. Species labels retain A2 semantics; non-ASCII text,
embedded newlines, titles longer than 80 characters, nonpositive depths/voltage,
missing DSMAX, non-four-character METHOD selectors and fractional scan-basis
components are rejected. Integer scan-basis output is corroborated by UPDATE's
integer format; coordinates and shifts are never rounded or wrapped.

Export uses **compiled active masks and layout**. Source ordering, region depths,
unequal stored-slice patterns, species/atoms, duplicate reflections, linked masks,
D1/D2 and PLAMDA are retained. If category filtering disabled bindings, those inactive
bindings and the category-selection mechanism cannot be recovered from FOR007.
Stable ID spellings likewise have no FOR007 representation; reimport generates IDs.

Every export returns (and the CLI prints) a report containing the omitted execution
backend/batch size, seed, bounds, objective/rejection policies, run directory,
checkpoint/snapshot settings and time/evaluation budgets. Bundles persist that report
as `compatibility.json`. These settings are not smuggled into legacy records or
silently treated as equivalent behavior. FOR007 represents legacy settings only;
retain canonical YAML as the authoritative configuration.

## Validation scope

Tests compare reference legacy → canonical → export → legacy parse, including
all physical state and masks, layout, ranges, PLAMDA and optimizer controls (filename
relocations are compared separately). They also cover variable reflection counts,
multiple regions with different patterns, variable thermal targets, linked absorption,
empty G slices, duplicate beams, full parameter categories, high-precision values,
no-overwrite checks and CLI migration. Representative forward comparisons include
complex amplitudes, intensities, YNORM and objective, using a shortened reference
depth and small multi-region cases. This validates the existing Python parser and
numerical path, not a new compiled-Fortran executable comparison or binary capacity
guarantee. No physics, propagation, objectives or optimizers were modified.

Stage 5 validation on Apple Silicon (`qcbed-python`): **72 focused tests passed
in 5.26 s; all 1879 tests passed in 27.20 s**. Twelve new regression cases passed.
Round-trip amplitude/intensity arrays matched exactly, as did YNORM/objective in
representative comparisons. The CLI help check and `git diff --check` passed.
Canonical configuration/compiler/runner and all numerical modules are unchanged
by this stage; the new compatibility modules are outside their import path.
