# Shared-parameter layout audit

The audited tree already reduces shared coefficient groups to optimizer coordinates.
The reported unchanged count was not reproduced as duplicate backend coordinates:
Real-only and Magnitude-only pairs contribute one coordinate, and Complex pairs
contribute two. The GUI still used ambiguous **Active parameters** labels and did
not expose the resolved group layout. Those labels are now split into **Active
physical targets** and **Optimizer parameters**, both derived from compilation.
No particular user run or saved project was supplied for this report, so the
precise state behind that observation cannot be reconstructed.

A separate functional defect was reproduced: explicit coefficient-only category
selection was rejected by the legacy preview's nonempty ordinary-vector check.
The modern bridge now permits a zero-coordinate ordinary preview before appending
the selected elastic coordinates. This also permits YAML target discovery before
ordinary target bindings have been restored. Legacy standalone subset validation
and numerical propagation are unchanged.

## Trace through the implementation

1. `gui/projected_elastic.py` commits one `equality_group` identifier per original
   hkl/slice record. Duplicate hkl/slice records are rejected; a record cannot have
   two group identifiers.
2. `projected_elastic.parameter_groups` resolves groups, verifies mode/slice/start
   compatibility and intersects hard bounds. Multiple groups remain distinct.
3. `attach_coefficients` appends one specification and step size per shared
   component, plus its bound entry. It does not append one specification per member.
4. `compiled.packed_starting_vector` follows that actual layout. The GUI's new
   `resolved_layout` presentation reads the same compiled component metadata.
5. `trial_records` reads each shared coordinate once and reconstructs each physical
   member. Real-only retains each fixed imaginary part; Magnitude-only retains each
   IAM phase; Complex shares both components. `apply_grids` preserves conjugate
   Friedel partners separately.
6. YAML retains `equality_group`, original hkl, slice and mode. Checkpoints retain
   the existing layout identity. Results presents group coordinates and individual
   physical members without modifying stored checkpoint labels or beam identities.

## Regression examples

- Two coefficient-only Real or Magnitude targets: 2 independent → 1 shared.
- Add depth: 3 independent → 2 shared; Make Independent restores 3.
- Two coefficient-only Complex targets: 4 independent → 2 shared.
- Depth plus four Real targets in two groups: 5 physical targets, 2 groups,
  1 independent target, 4 mandatory Friedel relationships, 3 optimizer coordinates.
- Three Complex members in one group: 3 physical targets, 2 coordinates.

Focused checks cover validation, starts and hard-bound lengths, independent group
writeback, YAML round-trip, checkpoint compatibility and Results member identities.
No refinement campaign is involved. See the release record for measured test and
packaging results.


## Follow-up: actual validation-button defect

The later real-widget/frozen audit found a missed path in the earlier tests:
`validate_c_part('refinement')` bypassed `elastic_panel.apply`, unlike the general
Apply/Validate route. Starting with unchecked/uncommitted coefficients, pressing
**Validate Refinement Targets** could report success while omitting both visible
coefficient drafts, then refresh away their values. Earlier tests began with
already-committed coefficient records or explicitly applied the panel and therefore
missed this defect. The button now commits the coefficient editor before ordinary
refinement validation. See `docs/sharing_gui_state_audit.md` for the before/after
state-layer evidence. This supersedes any implication that backend-only checks
established correctness of every GUI validation path.
