# Canonical regression data

`reference_config.yml` is the fixed YAML counterpart of
`examples/reference_configuration.py::reference_configuration`. Keep complete
canonical equality here, including presentation metadata. Intentional reference
changes require reviewing both representations and their numerical comparisons.
Relative paths are resolved against this fixture directory; forward tests supply
a temporary experimental body derived from the read-only reference data.

`examples/qcbed_config.yml` is a user-editable configuration, not the immutable
regression baseline. Editor tests may open it without requiring its title to match
this fixture. The title-only regression separately verifies exact compiled physics,
forward amplitudes/intensities, YNORM, and objective invariance.
