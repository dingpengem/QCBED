# YAML configuration frontend (schema version 1)

[Documentation index](README.md) · [Current architecture](ARCHITECTURE.md) ·
[Researcher manual](user_manual.md)

Install the optional frontend with `pip install -e '.[yaml]'`. The development
conda environment already declares `pyyaml`. Numerical modules and the Python
configuration API do not import or require PyYAML.

```python
from qcbed import load_yaml_config, compile_config, load_experimental_data

raw = load_yaml_config('examples/qcbed_config.yml')
compiled = compile_config(raw)
data = load_experimental_data(compiled, base_directory='.')
```

The [complete example](../examples/qcbed_config.yml) represents the independently
constructed [Python reference](../examples/reference_configuration.py), with paths
relative to the YAML file. Run that Python example once to create its separate
headerless data body. Loading YAML itself does not create files, load experimental
arrays, or run a simulation. The [Stage 4 production entry point](yaml_runs.md) consumes this configuration.

## Mapping to typed objects and GUI controls

Every mapping is one existing immutable dataclass; its keys are exactly the field
names below. Sequences become immutable tuples. No parallel YAML configuration
model exists. A GUI can construct the same classes directly, using the same
required fields, defaults, enum choices, ID references, and compiler validation.
`load_yaml_config` validates through the existing pure compiler and returns only
`QcbedConfig`; the caller compiles it for execution. The compiler stays independent
of YAML and path resolution.

Required fields must be present. Optional fields use their dataclass defaults;
`null` is accepted only for nullable types. Integer fields reject floats, booleans,
and quoted numbers. Real fields accept finite integer or floating numeric scalars,
never strings or booleans. Use scientific notation such as `1.0e-6`. Unknown
fields, duplicate mapping keys, invalid enums, duplicate IDs, missing references,
invalid target arity/components and invalid dimensions are errors. Error messages
from frontend type checks include a field path; compiler errors describe the
invalid physical reference or setting.

Parameter category enum values are the lowercase strings in the
[canonical target table](configuration.md#stable-references-and-ordering).
Vacancy and surface enum values are integers `0`, `1`, `2`, not enum names.
Literal choices in the field table are exact and case-sensitive. Optimizer method
strings are validated by the existing compiler/optimizer configuration.

Use ordered list controls for regions, slices, species, atoms, reflections and
beams. Their YAML sequence order is retained exactly. Duplicate beam and reflection
H/K values remain separate records; reflection IDs must remain distinct. Bindings
and bounds select stable IDs, while compiled vector order follows the existing
source-defined layout, irrespective of mapping-key or binding order. Depth and
Debye–Waller targets belong to specific regions; no fixed region or parameter
count is introduced. See the canonical documentation for ID scopes and target
components.

Units: depths, slice thicknesses and cell lengths are angstroms; thermal factors
are angstrom squared; energy is kV; lattice angles are degrees. Atom coordinates,
slice shifts and scan basis/corner coefficients are dimensionless. Values are not
wrapped, clipped, folded, normalized or deduplicated by the frontend. Existing
unsupported physical states remain subject to the same evaluator restrictions.

## Filesystem paths

`experiment.path`, `run.output_directory`, `run.output_filename` and optional
`legacy_export.filename` resolve to absolute paths relative to the YAML file's
parent directory (including normal `..` and symlink resolution). Absolute paths
remain absolute. `run.output_filename` is independently relative to that directory,
not relative to `output_directory`. No existence check or directory creation is
performed for these configured paths. The source YAML file must exist.

Optional raw experimental/output filename metadata must first agree with its
unresolved configured filename. Its path portion is then resolved identically,
retaining trailing whitespace, to keep the compiler's legacy consistency checks.
Other legacy metadata is untouched. Equivalence to direct Python configurations
therefore compares the same explicitly resolved paths.

No environment substitution, tilde expansion, Python expression evaluation or
executable tags are supported. A string such as `${HOME}/in.dat` remains a literal
path component. Only one YAML document is allowed. SafeLoader is used with duplicate
key rejection; aliases and merge keys are rejected. Object-construction tags are
never enabled. Quote identifiers that YAML might otherwise interpret as booleans,
numbers or timestamps.

## Field reference

The following tables mirror the existing dataclasses. `required` means no default;
all other entries show the Python default. Nested dataclass names denote mappings,
`tuple[T, ...]` denotes an ordered list, fixed tuples require that exact length,
and `| None` permits YAML `null`.

### MeshConfig

| Field | Type | Default |
| --- | --- | --- |
| `mesh_x` | `int` | `required` |
| `mesh_y` | `int` | `required` |

### SliceGeometry

| Field | Type | Default |
| --- | --- | --- |
| `thickness_angstrom` | `float` | `required` |
| `x_shift_fractional` | `float` | `0.0` |
| `y_shift_fractional` | `float` | `0.0` |
| `cell_a_angstrom` | `float \| None` | `None` |
| `cell_b_angstrom` | `float \| None` | `None` |
| `lattice_gamma_degrees` | `float \| None` | `None` |

### ForwardExecutionConfig

| Field | Type | Default |
| --- | --- | --- |
| `execution_backend` | `typing.Literal['serial', 'batched']` | `'serial'` |
| `batch_size` | `int \| None` | `None` |

### MicroscopeConfig

| Field | Type | Default |
| --- | --- | --- |
| `energy_kv` | `float` | `required` |
| `zone_axis` | `tuple[int, int, int]` | `required` |

### SpeciesConfig

| Field | Type | Default |
| --- | --- | --- |
| `id` | `str` | `required` |
| `label` | `str` | `required` |
| `atomic_number` | `int` | `required` |
| `label_raw` | `str` | `''` |

### AtomConfig

| Field | Type | Default |
| --- | --- | --- |
| `id` | `str` | `required` |
| `species` | `str` | `required` |
| `x_fractional` | `float` | `required` |
| `y_fractional` | `float` | `required` |
| `starting_occupancy` | `float` | `required` |

### StoredSliceConfig

| Field | Type | Default |
| --- | --- | --- |
| `id` | `str` | `required` |
| `geometry` | `SliceGeometry` | `required` |
| `atoms` | `tuple[AtomConfig, ...]` | `required` |
| `record_type` | `Literal['S', 'G']` | `'S'` |
| `record_type_raw` | `str` | `''` |

### VacancyConfig

| Field | Type | Default |
| --- | --- | --- |
| `kind` | `VacancyKind` | `<VacancyKind.CONSTANT: 0>` |
| `values` | `tuple[float, ...]` | `(0.0,)` |

### SurfaceConfig

| Field | Type | Default |
| --- | --- | --- |
| `kind` | `SurfaceKind` | `<SurfaceKind.NONE: 0>` |
| `coefficients` | `tuple[float, float, float]` | `(0.0, 0.0, 1.0)` |

### ThermalFactorConfig

| Field | Type | Default |
| --- | --- | --- |
| `species` | `str` | `required` |
| `value_angstrom_squared` | `float` | `required` |

### ReflectionConfig

| Field | Type | Default |
| --- | --- | --- |
| `id` | `str` | `required` |
| `h` | `int` | `required` |
| `k` | `int` | `required` |
| `elastic` | `float` | `required` |
| `absorptive` | `float` | `required` |

### RegionConfig

| Field | Type | Default |
| --- | --- | --- |
| `id` | `str` | `required` |
| `depth_angstrom` | `float` | `required` |
| `slices` | `tuple[StoredSliceConfig, ...]` | `required` |
| `debye_waller` | `tuple[ThermalFactorConfig, ...]` | `required` |
| `reflections` | `tuple[ReflectionConfig, ...]` | `()` |
| `symmetry` | `Literal['NONE', 'INVE', 'SQUA', 'RECT', 'HEXA']` | `'NONE'` |
| `vacancy` | `VacancyConfig` | `VacancyConfig(kind=<VacancyKind.CONSTANT: 0>, values=(0.0,))` |
| `free_surface` | `SurfaceConfig` | `SurfaceConfig(kind=<SurfaceKind.NONE: 0>, coefficients=(0.0, 0.0, 1.0))` |

### AbsorptionSpeciesConfig

| Field | Type | Default |
| --- | --- | --- |
| `species` | `str` | `required` |
| `a` | `float` | `required` |
| `b` | `float` | `required` |
| `c` | `float` | `required` |
| `c_mode` | `Literal['independent', 'linked_to_b']` | `'independent'` |

### AbsorptionConfig

| Field | Type | Default |
| --- | --- | --- |
| `enabled` | `bool` | `False` |
| `coefficients` | `tuple[AbsorptionSpeciesConfig, ...]` | `()` |
| `legacy_raw` | `str \| None` | `None` |

### SampleConfig

| Field | Type | Default |
| --- | --- | --- |
| `mesh` | `MeshConfig` | `required` |
| `species` | `tuple[SpeciesConfig, ...]` | `required` |
| `regions` | `tuple[RegionConfig, ...]` | `required` |
| `dsap_inverse_angstrom` | `float` | `required` |
| `dilation` | `float` | `0.0` |
| `absorption` | `AbsorptionConfig` | `AbsorptionConfig(enabled=False, coefficients=(), legacy_raw=None)` |

### ScanConfig

| Field | Type | Default |
| --- | --- | --- |
| `scan_columns` | `int` | `required` |
| `scan_rows` | `int` | `required` |
| `basis_1` | `tuple[float, float, float]` | `required` |
| `basis_2` | `tuple[float, float, float]` | `required` |
| `corner_a` | `tuple[float, float]` | `required` |
| `corner_b` | `tuple[float, float]` | `required` |
| `corner_c` | `tuple[float, float]` | `required` |

### ExperimentalDataConfig

| Field | Type | Default |
| --- | --- | --- |
| `path` | `str` | `required` |
| `beam_indices_hk` | `tuple[tuple[int, int], ...]` | `required` |
| `gradient_option` | `int` | `required` |
| `header` | `Literal['absent', 'validate']` | `'absent'` |
| `format` | `Literal['legacy_intensity_blocks']` | `'legacy_intensity_blocks'` |

### ParameterTarget

| Field | Type | Default |
| --- | --- | --- |
| `category` | `ParameterCategory` | `required` |
| `identifiers` | `tuple[str, ...]` | `()` |

### ParameterBinding

| Field | Type | Default |
| --- | --- | --- |
| `target` | `ParameterTarget` | `required` |
| `range_factors` | `tuple[float, float]` | `required` |
| `plamda` | `float` | `required` |

### OptimizerConfig

| Field | Type | Default |
| --- | --- | --- |
| `method` | `str` | `required` |
| `loops` | `int` | `required` |
| `maxfsimp` | `int` | `required` |
| `maxfpara` | `int` | `required` |
| `ftol` | `float` | `required` |
| `chimin` | `float` | `required` |
| `maxtemp` | `float` | `required` |
| `mintemp` | `float` | `required` |
| `coolrate` | `float` | `required` |
| `seed` | `int` | `required` |
| `paraopt` | `int` | `0` |
| `loud` | `int` | `0` |

### TargetBound

| Field | Type | Default |
| --- | --- | --- |
| `target` | `ParameterTarget` | `required` |
| `lower` | `float \| None` | `None` |
| `upper` | `float \| None` | `None` |

### CategoryBound

| Field | Type | Default |
| --- | --- | --- |
| `category` | `ParameterCategory` | `required` |
| `lower` | `float \| None` | `None` |
| `upper` | `float \| None` | `None` |

### RelativeBound

| Field | Type | Default |
| --- | --- | --- |
| `category` | `ParameterCategory` | `required` |
| `fraction` | `float` | `required` |
| `absolute_fallback` | `float` | `required` |

### BoundsConfig

| Field | Type | Default |
| --- | --- | --- |
| `policy` | `Literal['none', 'reject']` | `'none'` |
| `specific` | `tuple[TargetBound, ...]` | `()` |
| `absolute` | `tuple[CategoryBound, ...]` | `()` |
| `relative` | `tuple[RelativeBound, ...]` | `()` |

### RefinementConfig

| Field | Type | Default |
| --- | --- | --- |
| `optimizer` | `OptimizerConfig` | `required` |
| `parameters` | `tuple[ParameterBinding, ...]` | `required` |
| `iref` | `float` | `required` |
| `objective_mode` | `Literal['legacy', 'masked']` | `'legacy'` |
| `invalid_trial_policy` | `Literal['raise', 'reject']` | `'raise'` |
| `active_categories` | `tuple[ParameterCategory, ...] \| None` | `None` |
| `bounds` | `BoundsConfig` | `BoundsConfig(policy='none', specific=(), absolute=(), relative=())` |

### RunConfig

| Field | Type | Default |
| --- | --- | --- |
| `output_directory` | `str` | `required` |
| `output_filename` | `str` | `'out.dat'` |
| `snapshot_every` | `int \| None` | `100` |
| `checkpoint_every` | `int` | `25` |
| `max_total_evaluations` | `int \| None` | `1000` |
| `max_wall_time_seconds` | `float \| None` | `3600.0` |

### LegacyExportConfig

| Field | Type | Default |
| --- | --- | --- |
| `dsmax_inverse_angstrom` | `float \| None` | `None` |
| `filename` | `str \| None` | `None` |
| `title_raw` | `str \| None` | `None` |
| `output_filename_raw` | `str \| None` | `None` |
| `experimental_filename_raw` | `str \| None` | `None` |

### QcbedConfig

| Field | Type | Default |
| --- | --- | --- |
| `title` | `str` | `required` |
| `microscope` | `MicroscopeConfig` | `required` |
| `sample` | `SampleConfig` | `required` |
| `scan` | `ScanConfig` | `required` |
| `experiment` | `ExperimentalDataConfig` | `required` |
| `refinement` | `RefinementConfig` | `required` |
| `execution` | `ForwardExecutionConfig` | `required` |
| `run` | `RunConfig` | `required` |
| `legacy_export` | `LegacyExportConfig` | `LegacyExportConfig(dsmax_inverse_angstrom=None, filename=None, title_raw=None, output_filename_raw=None, experimental_filename_raw=None)` |
| `schema_version` | `int` | `1` |

## Stage 3 validation

On Apple Silicon using `qcbed-python`, **107 focused tests passed in 3.10 s**
and **all 1859 tests passed in 23.01 s**. The 39 new frontend cases cover strict
rejections, default handling, path and raw-metadata resolution, dynamic reflection
and thermal targets, and ordered multi-region configurations. Complete typed and
compiled reference equivalence uses explicitly identical resolved paths. Forward
comparisons use the full reference scan with shortened depth and match amplitudes,
intensities, YNORM and objective exactly. This is equivalence to the existing Python
path, not a new compiled-Fortran comparison. `git diff --check` passed.

No numerical kernels, canonical dataclasses, compiler, dynamic layout, experimental
reader, legacy parser or production runner were changed for this frontend. FOR007
retirement/export and GUI implementation remain deferred.

## Crystallography-derived mode (B4b)

The optional `crystallography` section selects original-structure input with
`sample: null`. In this mode, omit scan `basis_1/basis_2` and experimental
`beam_indices_hk`; these are derived. Required original hkl basis, beam coordinates,
explicit physical slicing, thermal factors and source-loading rules are documented
in [the crystallography schema](crystallography.md#b4b-production-configuration-integration).
See `examples/crystallography_config.yml`. Legacy YAML remains supported unchanged.
The GUI editor does not yet support this mode (B4c).

## Modern A/B/C document frontend

Modern Save/Save As now uses `schema_version: 2` with `data_prep`, `atomic_model`
and `qcbed`; see [the schema and versioned defaults](modern_yaml.md). This frontend
constructs the same `QcbedConfig`; the version-1 schema documented above remains
the legacy/verbose compatibility input. Prepared arrays and interactive provenance
are referenced externally, not embedded in user YAML.

### Compact modern search defaults

A selected entry in `qcbed.refinement.targets` may be `{}` to use the documented
[reference search profile](modern_refinement_catalogue.md). `range_factors` and
`step_size` are optional overrides, resolved into the same typed ParameterBinding.
Compact Save omits each value only when exactly equal to its canonical default;
old explicit compact and verbose configurations still load unchanged. Bounds
remain separate in `qcbed.refinement.bounds`; no physical bounds are inferred.
Defaults are stable frontend metadata, not changes to optimizer algorithms.
