# Scientific provenance and private reference material

QCBED's Python implementation was developed against an earlier multislice program
and recorded scientific reference data. The shared repository retains the Python
implementation, numerical/functional regressions and scientific reference inputs.
Legacy implementation sources, source listings and detailed inventories are held
privately by the maintainer; they are not distributed in this repository.

Historical validation statements describe the inputs, software and private assets
available at the time. Related implementation versions must not be treated as
interchangeable evidence. No claim of executable equivalence extends beyond the
comparisons actually performed. Preserve documented units, axes, array ordering,
Fourier signs, normalization, boundary conditions and parameter semantics.

## Validation without private sources

Normal Python/app use does not require the private implementation or a Fortran
compiler. Most numerical, GUI and configuration tests use the shared Python code
and recorded data. Five direct source-integrity/comparison checks retain their
original scientific assertions. They require maintainer
access to the private sources; the inventory comparison also needs the private
original inventory. A shared checkout alone cannot reproduce the earlier complete
2,965-test acceptance result. That result remains a historical measurement.

The affected checks are the three source-dependent checks in
`tests/test_legacy_source_inventory.py`, the source-table provenance comparison in
`tests/test_scattering.py`, and the compiled predecessor-sequence comparison in
`tests/test_one_slice_simulation.py`. Do not fabricate private inputs or weaken
these assertions to obtain a green test count.

The public JSON inventory is a provenance notice, not an embedded copy of private
source statements. Maintainers retain its original privately for integrity work.
See [current architecture](ARCHITECTURE.md), [configuration](configuration.md), and
[the refinement catalogue](modern_refinement_catalogue.md) for shared contracts.

Maintainers may set `QCBED_PRIVATE_REFERENCE_ROOT` to an external evidence directory
when running pytest. Private evidence is never copied into the shared tree or app.
Without it, five optional source checks are explicitly skipped; public numerical
regressions still run. A configured but invalid directory fails rather than skipping.
