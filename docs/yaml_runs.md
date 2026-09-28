# Production refinement from YAML

The recommended entry point is now:

```sh
python -m qcbed.run --config examples/qcbed_config.yml
python -m qcbed.run --config examples/qcbed_config.yml --resume
```

Use the `qcbed-python` environment, or install the package with
`pip install -e '.[yaml,visualization]'`. See the [schema reference](yaml_configuration.md)
for the complete typed schema. The supplied YAML example expects a headerless
experimental body; `python examples/reference_configuration.py` creates that
separate copy once without changing reference data. For a new experiment, supply
its intensity/variance/weight body yourself, or explicitly select `header: validate`
for a file with a matching legacy header. No FOR007 file is needed by this workflow.

Before running, choose a fresh `run.output_directory` and `run.output_filename`.
Every path is resolved relative to the YAML file, including `output_filename`
(which is not relative to `output_directory`). The runner uses the configured
output directory exactly; it does not append a timestamp or apply quick/example
overrides. Set budgets, snapshots, checkpoints, execution, seed, refinement targets
and bounds in YAML. A fresh run refuses an existing run directory or output file.
The shipped reference settings select the established masked objective and bounds;
these remain explicit user-selected policies, not changes to legacy semantics.

Programmatic callers use:

```python
from qcbed.run import run_yaml_refinement
result = run_yaml_refinement('experiment.yml', resume=False)
```

`progress_hook` is the existing completed-record observer and may be supplied to
this API. It does not change optimizer callbacks or evaluation counts. The CLI
prints the final summary and exits zero for completed runs, one for other returned
statuses. Configuration and input errors propagate before refinement begins.

## Shared execution and outputs

The adapter calls `load_yaml_config`, `compile_config`, `load_experimental_data`,
then `run_full_refinement`. The latter remains the common production execution
path used by the legacy example. No physics, objective or optimizer implementation
is duplicated. All sample regions, selected parameters, bounds and scheduling
settings come from the compiled configuration.

Existing run artifacts remain: progress JSONL, parameter snapshots, active vectors,
intensities, retained-stage checkpoints, summaries, previews, plots, report and
`final_out.dat`. Final simulated output is also written to `run.output_filename`
using the same existing output writer. This is simulated-intensity output, not
FOR007 export. Protected source/reference directories remain forbidden output
locations. No GUI or FOR007 exporter is introduced.

Canonical Python callers can also use `qcbed.run_config(config, base_directory=...)`
or prepare an immutable `RunSnapshot`; see [frontend APIs](frontend_apis.md). The
YAML entry point delegates to this same canonical path.

## Identity and resume

Canonical checkpoints carry `checkpoint_identity.format: qcbed-canonical-run`
and identity version 1. Its digest covers the resolved compiled state, layout,
runner settings, experimental specification and run settings. The shared runner
also fingerprints the parsed experimental metadata and all intensity, variance and
weight arrays. It does not hash YAML bytes or require FOR007. Comments and mapping
key order may change without invalidating resume; changing data, paths, physical
values, optimizer/execution settings or run settings is rejected. Ordered physical
records remain meaningful and must not be reordered to change the problem.

Resume uses the existing deterministic retained-stage restart policy; it is not
an arbitrary instruction-level optimizer restart. Use the same config and output
directory with `--resume`. Existing concurrency locks and interruption handling
remain in force. Canonical checkpoints cannot be resumed through the legacy entry
point, and unmarked legacy checkpoints cannot be resumed through the canonical
entry point. There is no implicit checkpoint conversion.

## Legacy compatibility workflow

`python examples/run_full_refinement.py` retains its FOR007/experimental-input
path, explicit example options, output-directory handling and resume behavior.
Use its `--help` for the existing options. Unmarked version-1 legacy checkpoints
remain supported by that workflow. The new adapter never calls a FOR007 parser;
legacy record classes are reused internally by the validated compiler.

## Integration coverage

Tests exercise YAML-only fresh runs without a FOR007 file and with FOR007 parser
calls made fatal, 1/2/3/5 regions and independent parameter subsets, exact-bound
rejection, snapshots/checkpoints/output files, and interrupted resume matching an
uninterrupted run. Reference-path comparison uses the parsed legacy Al state with
an explicitly shortened depth and evaluation-only controls, checking identical
final intensities, YNORM and objective. Existing numerical regression suites remain
the broader equivalence baseline. Test runs do not claim a new compiled-Fortran
comparison or a completed full-depth reference optimization.

Stage 4 validation on Apple Silicon (`qcbed-python`): **81 focused tests passed
in 5.44 s; all 1867 tests passed in 24.98 s**. The eight new integration cases
passed, including exact resumed/uninterrupted vector and physical-state comparisons.
The module CLI help command and `git diff --check` also passed.
