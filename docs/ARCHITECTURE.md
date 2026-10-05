# Developer architecture

[Documentation index](README.md) · [Specification](../PROJECT_SPEC.md)

## Input and scientific boundaries

`configuration.py` defines immutable typed configuration. `yaml_config.py`,
`modern_yaml.py`, `config_serialization.py` and `config_paths.py` decode/encode
references against explicit document bases. The legacy FOR007 frontend remains
separate. `config_compiler.py` resolves physical inputs and target layout;
`config_validation.py` validates the typed contract.

`data_preparation/` and `detector/` own measured-pattern geometry, calibration,
uncertainty and weight preparation. `crystallography/` owns original-site identity,
zone recasting and strict simulation-sample construction. The multislice kernels
and the existing refinement engines retain their tested units, ordering and FFT
conventions. Consult [Fortran architecture](fortran_architecture.md) before changing
numerical code; consult [source provenance](source_provenance.md) before interpreting
related-version declarations.

## GUI ownership

`gui/window.py` coordinates A/B/C/D. `preparation_view.py`/`preparation_model.py`
handle interactive raw preparation; `data_view.py`/`data_model.py` inspect prepared
signals and masks. `modern_widgets.py` owns atomic-model/site-group presentation.
`document.py` supplies editable fields and committed state. `workflow.py` owns
transient validation receipts. `optimizer_panel.py` supplies optimizer controls.

C contains only Simulation Setup, Refinement Targets and Optimizer. D contains
`run_dashboard.py`, output validation via `output_paths.py`, and Advanced execution.
Moving a widget between scroll areas must detach it from its previous owner before
that owner is deleted; visibility tests exercise this Qt ownership requirement.
Do not duplicate execution/output state in C or serialize readiness receipts.

## Execution and results

`prepare_run_snapshot` → `run_snapshot` in `canonical_run.py` is the shared execution
entry point. GUI `run_controller.py` launches the immutable snapshot in a child
process. `grouped_run.py` adapts modern physical targets to the existing optimizers
and durable per-evaluation journals. `gui/readiness.py` checks capability and Resume
identity separately from writable-directory validation. `output_paths.py` maps a
GUI output container to a unique fresh child or an existing checkpoint directory.

Batched scheduling controls memory/work partitioning; it does not define physics.
The canonical modern backend requires batched execution. Standalone tilt-process
APIs do not imply a worker setting in canonical configuration. Grouped execution
checkpoints each completed evaluation; legacy interval fields remain compatible
but do not change this modern cadence. Evaluation/time budgets apply per invocation.

Results readers consume recorded histories, checkpoints and comparison arrays;
viewing results does not recompute simulations. Retained, best and latest states
are distinct. A pre-evaluation cancellation can validate lifecycle wiring but cannot
supply a fitted result or demonstrate convergence.

## Distribution and tests

`packaging/macos/build.py` drives the PyInstaller spec and records source/build
provenance. `app_smoke.py` dispatches explicit external-input acceptance manifests;
`usability_smoke.py` exercises widget visibility and A/B/C/D without a campaign.
Tests cover numerical references, preparation, crystallography, dynamic targets,
checkpoint identity and GUI state. See [repository audit](REPOSITORY_AUDIT.md) and
[release report](RELEASE_REPORT.md). Keep historical scientific fixtures even when
new tutorials use a simpler configuration.
