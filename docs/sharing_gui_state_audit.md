# Actual sharing-button and validation audit

## Confirmed defect and scope

The real validation button bypassed coefficient Apply/commit. In the instrumented pre-fix frozen app, checking depth and both Real-only coefficients then pressing Validate reported success but committed only depth (vector length 1), discarding the coefficient drafts during refresh. This is a GUI draft/validation synchronization defect, not duplicate allocation in the shared-coordinate compiler. The exact reported 3 → 3 observation was not independently reproduced; the app path used for that observation has not been identified.

Earlier tests started with committed coefficient records or applied the editor explicitly. The new regression starts with no coefficient records, edits real widgets, clicks Validate, selects physical rows with mouse events, clicks Share, accepts the actual dialogs, revalidates, clicks Make Independent, and revalidates again. The frozen harness performs the same route after real A/B/Simulation Setup validation.

## Fix

Validate Refinement Targets now commits coefficient drafts before validating the complete layout. Share/Make Independent commit immediately, invalidate the refinement receipt and C readiness, disable Run/Resume, refresh, and display the compiled dimension transition. No hidden second Apply is required. The coefficient panel, D, and resolved summary use compiled dimensions; the panel no longer labels hypothetical checked-row counts as active optimizer parameters. Physics and stored group identities are unchanged.

## Frozen Real-only trace

All rows below retain three active physical targets. Bounds length counts one lower/upper bound entry per optimizer coordinate.

| State | Starting-vector length | Bound entries | Displayed optimizer count | Refinement receipt |
|---|---:|---:|---|---|
| independent_validated | 3 | 3 | 3 committed | Present |
| made_independent_validated | 3 | 3 | 3 committed | Present |
| selected_before_share | 3 | 3 | 3 committed | Present |
| shared_before_validation | 2 | 2 | 2 committed | Invalidated |
| shared_validated | 2 | 2 | 2 committed | Present |
| independent_before_validation | 3 | 3 | 3 committed | Invalidated |
| independent_revalidated | 3 | 3 | 3 committed | Present |

## Exact compiled coordinate order

### independent_validated

1. `sample.requested_depth` — start `1000.0`, bounds `[500.0, 1500.0]`.
2. `Projected elastic coefficient E_g (0, 2, 0), slice [0, 1] / real` — start `19.5`, bounds `[0.0, 40.0]`.
3. `Projected elastic coefficient E_g (2, 0, 0), slice [0, 1] / real` — start `19.5`, bounds `[0.0, 40.0]`.
### shared_validated

1. `sample.requested_depth` — start `1000.0`, bounds `[500.0, 1500.0]`.
2. `Shared refinement parameter group Audit pair {(0, 2, 0), (2, 0, 0)}, slice (0.0, 1.0) / real` — start `19.5`, bounds `[0.0, 40.0]`.

## Physical selection, persistence and writeback

Both selected row identities are original `(0,2,0)` and `(2,0,0)`, slice `[0,1]`. Both Refine checkboxes remain checked after validation. Independent rows display Independent and store `equality_group: null`; shared rows display Shared: Audit pair and store `equality_group: Audit pair`. Compiled group metadata contains both original hkl/slice identities. The real group contributes exactly one coordinate. A small in-memory coordinate perturbation writes the same changed real component to both owners while preserving their own fixed imaginary parts and mandatory Friedel partners. No objective evaluation is invoked.

JSON records retain every selected row, checkbox, visible sharing state, committed record, receipt, compiled label/identity, vector and bound lengths, displayed counts, readiness flag and writeback value. Local records: `tmp/sharing_actual/after/sharing-real.json`, `sharing-complex.json`, and `sharing-magnitude.json`. The pre-fix frozen record is `tmp/sharing_actual/before2/sharing-real.json`; that diagnostic run stopped on an invalid-bound writeback after the omitted drafts. It is not a passing smoke.

## Validation

- Final focused suite: **83 passed**, 176 Qt deprecation warnings, 43.64 seconds.
- Final frozen actual-widget smoke: **PASS**, zero optimizer evaluations.
- Real/Magnitude with depth: **3 → 2 → 3**; Complex with depth: **5 → 3 → 5**.
- App and ZIP rebuilt; **331 arm64 Mach-O binaries** (includes QtTest for the opt-in widget smoke), ad-hoc signing, deep/strict verification and ZIP integrity passed.
- Full regression was not rerun for this focused fix. No refinement campaign. Apple Silicon only; not notarized.
