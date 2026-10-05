# Atomic Model and mask inspection

Loading a different structure clears previous bindings and thermal metadata.
Initial groups collect sites with the same element, source occupancy and source
Debye–Waller B. This is an editable starting partition, not a symmetry-equivalence claim.
The Al FCC CIF produces one group with four sites. Debye–Waller B remains blank if the source
does not supply it. Enter Debye–Waller B (Å²), select Apply Site-Group Parameters, then validate B. Atomic Model.
Group ID, element, source member indices, occupancy and absorption coefficients
are editable there. Blank absorption means disabled. Group membership must cover
every original site exactly once; bindings must be removed before renaming IDs.
Independent same-element groups remain subject to the existing backend limitation.

The selected-beam mask selector exposes:

- Geometric Mask: existing extraction ellipse/radial support intersected with
  detector support, integrated and thresholded onto the output pixel grid.
- Processing Validity: additional eligibility inside geometric support.
  Outside geometric support this diagnostic is neutral (one), not a claim of
  measurable data.
- Manual Mask: explicit user exclusions, initially all ones.
- Final Weight: the exact binary objective input, equal to the AND of these layers.

The existing radius factor (0.9) and coverage threshold (0.875) are unchanged.
Gradient mode propagates numerical validity through its neighbor stencil.
Conventional mode also retains positive-intensity/noise-validity checks: its final
base weight is not guaranteed to equal geometric support for defective inputs.
For positive valid conventional input, the base weight equals geometric support.

New live preparations record geometric/processing diagnostics in their existing
provenance. Older preparations without these records explicitly show the layers
as unavailable; they are not reconstructed from guessed geometry. Manual mask
snapshots retain the original diagnostic provenance. No scientific arrays are
modified to make the mask display look complete.

Slice boundaries are editable under C → Numerical / Simulation cell. They are
fractional depth coordinates, strictly increasing from 0 to 1, with at least
three entries (two slabs). New generic drafts leave them unset; loading a CIF
does not infer a slicing scheme. The established FCC Al reference uses
0, 0.5, 1. New preparation retains an existing B/C slicing choice, but invalidates
validation receipts. Enter and apply explicit boundaries before validating C. Main QCBED Parameters.

Site-group table edits are drafts until Apply Site-Group Parameters commits all rows. Validate and Save never implicitly apply these drafts. For the conventional FCC Al reference, the initial Al_fcc group needs a user-supplied Debye–Waller B of 0.86 Å². No thermal value is invented. Applied values survive Save / Save As and reload.

## Stage-local validation

Validate A checks preparation only and records its snapshot. Validate B checks
the loaded structure, beam-frame orientation, original/conventional reciprocal
basis, site ownership and explicit thermal/occupancy/absorption metadata.
It does not build a simulation sample or require slice boundaries, depth,
optimizer, execution or output settings. No placeholder numerical settings are
used. Successful validation displays “B. Atomic Model: Valid”.

C owns slice boundaries and the remaining simulation/run controls. Invalid C
does not erase a valid B receipt. Changing A invalidates A/B/C; changing atomic
model fields invalidates B/C; changing numerical/run fields invalidates C only.
Run remains disabled until all receipts are valid and retains the canonical
full-configuration compiler and compatibility checks.

## Modern output directory

C exposes only Output directory. The final numerical export is managed as
final_out.dat inside that directory; histories, checkpoints, reports and
diagnostics retain the runner's existing internal filenames. Changing the
directory also changes the derived final export path. Modern Save/Save As writes
the directory and any run-policy settings, without a separate output filename.
Older compact filename/output_file keys remain readable and are normalized to
the managed path. Legacy-mode configuration retains its explicit output filename.
Existing checkpoint compatibility checks remain authoritative; a migrated output
path is not a promise of compatibility with a checkpoint using different paths.

## One project YAML

Normal launch loads no project. File → Load Configuration… (also on the toolbar)
explicitly selects an existing YAML. “Project configuration: Unsaved” means no
project YAML is associated yet. A/B/C can be completed entirely in memory.

After Apply and A/B/C validation, Save opens Save As once for an unsaved project.
Subsequent Save atomically updates that same file; Save As is an explicit copy.
Unapplied site-group drafts and unvalidated changes cannot be saved through these
GUI actions. Recomputing preparation retains the associated YAML path.

Modern saves contain data_prep, atomic_model and qcbed inputs, with internal JSON
preparation references for exact replay. Those JSON/data caches are artifacts,
not additional project YAMLs. Arrays, mask-layer diagnostics, validation receipts,
compiled state and histories are not embedded in YAML. Compatibility-only fields
are accepted on load but omitted from modern saves. There is no separate raw/A/B/C
YAML save button. Explicit CLI configuration paths remain supported.

## Refinement drafts and committed starting values

Binding selects a target; Active identifies a complete committed selection.
Draft edits are never active until Apply succeeds. Missing D1/D2/PLAMDA fields
are named explicitly. Starting values come from the committed physical model;
an edited requested-depth draft is annotated until the atomic Apply succeeds.

D1 and D2 are dimensionless legacy scaling factors, not hard bounds:
u = (x − D1*x0) / ((D2−D1)*x0). PLAMDA displaces an initial simplex vertex in u,
corresponding to physical displacement PLAMDA*(D2−D1)*x0. These fields retain
their numerical semantics; the grouped ANNEAL path uses its adaptive proposals.
Lower/Upper columns specify separate physical hard bounds.

Grouped Run supports exactly requested depth plus one shared site-group
Debye–Waller target, both with positive finite ordered bounds containing their
starting values. Validate C records external file/output preflight as well as
configuration checks; no extra generic Validate action is required afterward.

Execution readiness details (below the Refinement table) reports visible drafts
separately from committed bindings, including starting values, D1/D2/PLAMDA and
inactive reasons. It lists prepared-data readability, recorded cache integrity,
CIF loading, output-directory safety and creatability separately. A Grouped Run
capability failure is an execution issue, never evidence that the prepared
experimental file failed to load. Raw CBED is needed to reprepare data, not to
execute a refinement using the prepared input.
