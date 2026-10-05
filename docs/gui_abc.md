# A/B/C configuration and D execution workspace

> Historical milestone/audit record. For current controls and the A/B/C/D workflow,
> use the [User Manual](user_manual.md); dated counts and replay inputs below retain
> their original provenance.

Both the source GUI and rebuilt macOS app read and save the compact version-2
schema, with external preparation assets resolved from their document bases.

Launch the source GUI with `python -m qcbed.gui` or the rebuilt `dist/QCBED.app`.
See [frozen-app validation](../packaging/macos/VALIDATION.md). The workflow tabs are:

- **A. Data Prep**: linked Pattern / Sigma / Weight views, experimental input file,
  optional full CBED overview and manual exclusion editing.
- **B. Atomic Model**: structure file, original/conventional zone/hkl inputs,
  Front/Top/Left views, physical site groups and capability information. Legacy
  documents retain their pre-sliced region/structure editor; no missing crystal
  structure or zone-axis geometry is fabricated for them.
- **C. Main QCBED Parameters**: dynamic refinement targets, requested depth,
  microscope/signal mode, scan, optimizer, execution and output settings.
- **D. Run & Results**: readiness dashboard and prominent Run QCBED, Resume
  QCBED, Cancel and Open Results buttons. Settings remain exclusively in C.

Fresh preparation projects start at requested depth 1000 Å; loaded values are
never replaced. Applying site groups commits pending B basis/orientation fields
and preserves independent C drafts. Partial B validation likewise retains C
drafts. Invalidating validation receipts does not reset these user inputs.
Configured physical hard bounds are visible on target generation; missing
canonical limits are explicitly identified rather than guessed. See the
[hard-bound audit](modern_refinement_catalogue.md#hard-bound-audit-and-immediate-target-display).

Existing Run/Resume/Cancel, Compatibility and Results controls remain available.
Results after Cancel and recorded latest/best/retained final states are unchanged.

## Data binding and manual exclusions

Data Prep reads the **configured experimental file using the production reader**.
Public arrays retain `(tilt, beam)` orientation and ordered entries, including
repeated beam indices. Within a beam `tilt = row * columns + column`. Pattern is
the already-prepared experimental signal. Sigma is displayed as the square root
of that file's variance; the original variance is never replaced by a square-root/
square round trip. No A4 provenance is claimed unless supplied externally.

Select a beam, click a pixel in any of the three views or drag a rectangle, review
the integer selection, then **Exclude selection**. The half-open bounds use x =
column, y = row. Rectangle selection includes pixel centers inside the drawn
rectangle; the visible integer bounds can be adjusted before applying. Only the
selected beam's weights become zero. Pattern and variance remain bitwise equal
when reread. No interpolation, denoising, normalization, signal conversion or
per-beam fitting occurs. Active-point count updates after each committed edit.

**Restore selection** restores that region's original weights; **Reset manual
exclusions** restores all original weights. Neither can activate pixels that the
source already excluded. Undo/Redo restores immutable configuration snapshots.
Legacy nonbinary weights are inspected unchanged; manual masking is disabled
rather than silently converting their meaning.

Each committed edit writes a new content-addressed headerless experimental body
under the document directory's `.qcbed-data/`, and updates the configuration's
experimental path. The original input is never overwritten. A `.dat.prep.json`
companion records the original path/hash, exclusions and optional overview
references. Those files are gitignored generated project assets, **not disposable
caches while a saved configuration references them**. Keep the data body for
execution and its companion/original source for Restore and overview provenance.

Save/Save As uses the existing canonical YAML serializer and absolute resolved
external paths. Save As does not copy external assets into its new directory.
Reload therefore preserves the exact mask. Runs consume the same ordinary
experimental body as CLI runs; no GUI-only numerical override exists. Existing
configuration/data fingerprint checks apply: Resume accepts the same mask and
rejects a changed mask. To refine with different exclusions, choose new run paths.
Active runs lock editor operations and retain their immutable input snapshot.

## Full CBED overview

Use **New Data Preparation** to load a full CBED image and generate the preparation
live, or **Load Previous Preparation…** to reopen a saved preparation. The saved
file restores its overview, beam regions and provenance together; no separate
attachment of image and metadata is needed. Existing saved formats remain supported.

The overview uses recorded ROI centers/radii as labelled bounding boxes and the
recorded `reciprocal_frame.zone_axis_origin_px` as its marker. It never infers an
origin from the central beam or a reciprocal index. Click a box to select that
ordered beam in all three linked views. For overlapping boxes, the closest ROI
center selects the beam. All image dimensions remain dynamic. No metadata means
an explicit unavailable message, not invented boxes or an invented origin.

## Atomic model and refinement

The three plots use the existing canonical beam-frame projection, without moving
atoms: Front = `(x_beam,y_beam)` (zone-axis view), Top = `(x_beam,depth)`, Left =
`(depth,y_beam)`, all in Å. Structure summary retains lattice lengths/angles, volume,
atom count, strict B5 compatibility and requested/effective depth separately.

Site-group membership, occupancy and absorption stay read-only; supported shared
B editing remains above the table. Refinement status and resolved active bounds
are displayed alongside group identity. Inactive specific bounds are labelled as
such. Use C → Refinement Targets for physical selection and bounds; search metadata is optional Expert / Legacy settings. The table remains
catalogue-driven, without a fixed parameter count. Modern execution accepts any nonempty supported subset, with depth optional,
using ANNEAL/AMOEBA/PARAMIN/SIMPARA. See the [parameter audit](modern_refinement_catalogue.md)
and [user guide](user_manual.md) for physical mappings and explicitly deferred families.

## Fresh runs and resume

**Fresh Run** and **Resume** remain separate. In D / Execution & Output, select or
type a directory and click **Validate Output Directory**. This checks only that it
can be created and written to; an existing nonempty folder is valid. Fresh Run owns
a new uniquely named subfolder and leaves existing files untouched. Resume separately
checks the newest managed checkpoint, with support for older direct checkpoints.
See [output workflow](user_manual.md#d-run--results).

## Validation record

Source validation for this interface update: five focused data/Qt tests passed;
94 GUI-focused tests passed; full suite **2844 passed**. Two Matplotlib
`constrained_layout` warnings occurred in offscreen GUI tests; no numerical test
failed. `git diff --check` was clean. The real-data inspection snapshot at
`tmp/abc_interface_validation/abc.yml` attaches the existing overview without manual
exclusions and retains 7281 active gradient points with active weights `[1.0]`.
The standalone app was not rebuilt.

## Sequential validation workflow

The configuration workspace now has three actual tabs. Each has an explicit
Validate A/B/C button. Tabs may be inspected in any order, but validation is
ordered: A → B → C → Fresh Run/Resume. Changing A clears all receipts; changing
atomic-model inputs clears B/C; other run-setting edits clear C. Opening a saved
document requires validation again. Receipts hold immutable configuration snapshots
and input identities; they are editor confirmations, not checkpoint identities.

A requires a full pattern with matching ROI metadata, finite signal/variance,
binary inclusion weights with positive active variance, and a selected or explicitly
confirmed zone-axis position. Choose the reference beam, enable **Select symmetry
center**, then click its local centrosymmetric feature on the full pattern. This
is a manual symmetry-center choice, not an automatic fitted result. It uses the
existing symmetry-center relation and A10 coordinate conversion; missing physical selection,
binning, central-center or basis metadata blocks this operation rather than guessing.
**Confirm recorded zone-axis position** accepts an existing symmetry-center result explicitly.
The Zone axis marker is distinct from the central/lattice origin. Labels show
experimental (n1,n2) and original/conventional hkl when available.

Inspect linked Pattern/Sigma/Weight, select a pixel/rectangle, and Exclude,
Restore selected or Reset manual exclusions. Only weights change; excluded signal
remains visible. Validate A captures those mask changes for subsequent stages.
B validates the existing crystallography/strict-B5 rules; C checks configuration,
execution capability and external files without evaluating an objective.

Modern Fresh Run uses unique `qcbed-run-<timestamp>-<id>` subfolders; it requires
no replacement confirmation. Resume retains canonical identity and lock checks.
Legacy GUI output replacement still uses explicit confirmation and sibling backups.
CLI overwrite behavior is unchanged. No refinement was rerun for GUI acceptance;
the frozen smoke path exercises editing and readiness only.

## Editable A/B/C YAML

Modern Save/Save As uses `data_prep`, `atomic_model`, and `qcbed`, mapping directly
to these tabs. The preparation descriptor preserves mask/origin/source references;
validation receipts themselves remain transient. See [modern schema](modern_yaml.md)
and [the compact example](../examples/fcc_al.yml). Old verbose modern YAML still
opens and is converted on Save. Legacy/pre-sliced saving is unchanged.

## Live experimental preparation

A. Data Prep opens with **New Data Preparation**. Import a full, single-channel
CBED image (numeric text, NPY, TIFF or PNG). An aperture image is optional.
Geometry precedence is explicit/manual override → optional disk-profile override
→ threshold geometry from the supplied aperture. With no image-derived geometry,
draw a diameter on the CBED and press Accept. A separate disk-profile image is
not required. No radius comes from array dimensions or the historical 38 px case.

With an aperture, existing A2/A3/A4 calibrate the detector PSF/OTF, correct the
whole pattern, and fit noise. Calibration and geometry remain separate backend
operations; their use of the same supplied aperture is explicit in provenance.
Without an aperture there is no correction or calibrated noise fit: the selected
Poisson/intensity approximation is sigma=sqrt(max(I,0)), variance=max(I,0), with
positive-intensity validity. Calibration failures never silently select this fallback.

Draw a rough centrosymmetric ROI, inspect its integer-grid preview, and press
Accept/Snap. Beam extraction dimensions are derived from the physical ellipse
diameters expanded upward to the next binning-compatible integer size; the measured
floating-point radii remain unchanged for arcs, ellipses, masks and provenance.
For binning 2, diameters 75.6 and 76.2 give array extents 76 and 78. Symmetry crop half sizes are half the accepted ROI width/height. Both appear
as read-only diagnostics. Binning is editable under Advanced / Data Processing
(default 2). Arbitrary drawn ROI bounds automatically preview an enclosing even-sized pixel-grid
rectangle. The center is preserved as closely as possible; near boundaries the
rectangle shifts inward minimally to keep all symmetry-search crops in frame.
Confirm the overlay with Accept. Only an impossible fit asks for a smaller region.
Finalized ROI and crop sizes are saved for exact replay. Physical radii are never rounded. Replay retains recorded sizes even for older explicit-size preparations.
Run Data Preparation invokes the existing A5–A10 routines in a spawned
process. Inspect the resulting beam boxes and refined Zone axis marker, explicitly
confirm the origin, inspect Pattern/Sigma/Weight, then Exclude/Restore bad pixels
and Validate A. Masking changes only weights. A may validate before B is complete;
an empty project starts with a clearly incomplete Atomic Model draft, never the
FCC Al reference model. Validate B and C remain required before Run.

**Load Previous Preparation** is the secondary replay path. Successful preparation
writes an exclusive `.qcbed-data/preparation-…/` cache containing ordinary
experimental blocks, full corrected image, uncertainty arrays, a data preparation metadata record,
selection/noise provenance and a compact prepared descriptor. Manual edits create
new immutable snapshots. Existing compact/verbose/legacy documents remain supported.

Raw-input drafts can be saved and reopened without prepared data:

```yaml
schema_version: 2
data_prep:
  pattern: full_cbed.dat
  aperture: aperture.tif  # optional
```

Existing atomic_model/qcbed sections may accompany the draft. After preparation,
normal modern Save retains raw references plus `data_prep.preparation` for replay;
selections belong to the cache, not an expanded YAML schema. An unprepared draft
cannot be compiled or run by the CLI. Paths are document-relative on input and
resolved on Save As.

Limitations are explicit: the current OTF backend requires matching array grids;
A4 retains its frame/fit requirements; A9 requires integer crop origins and does
not silently round refined beam centers. Fractional geometry or insufficient
lattice evidence produces an error. Live sparse A6 evidence (three to five arcs)
now reuses the existing two-validated-arc seed, requiring every detected candidate
to agree within the existing lattice residual tolerance. Deterministic selection
uses geometric confidence; detector thresholds and upstream radii are unchanged.
For `pattern.dat`, explicitly supplied 38 px radii give center (127,126), arc centers
(118,208), (136,44), (208,135), and seed (81,9), (-9,82). Sixteen predicted sites
produce nine A9 beams with full physical ROI coverage (76×76 selection, binning 2).
A7 insufficient/missing evidence remains nonfatal; seven sites fail A9 coverage.
No archived centers, catalogue, or origin are imported. The one-pixel seed difference
from the reviewed historical pair is retained, not tuned away. This is geometry/
extraction validation, not validation of the no-aperture Poisson noise approximation.


### Live A → B → C handoff

Validate A publishes an immutable, content-addressed `*.prepared.json` descriptor
in the document's `.qcbed-data` directory, using the same descriptor publisher as
Save. It links the confirmed symmetry-center origin, mask snapshot, beam identities/order,
Pattern/variance/Weight, ROIs, live lattice provenance and current A10 corners.
Sigma remains the square root of the recorded variance; full-frame uncertainty
products are retained by live preparation. B/C consume the validated A fields;
changes to beam order/labels, scan transform, scan geometry or source artifacts
invalidate the A receipt. Live A validation rejects mismatched A10 corners.
The descriptor is a preparation output/cache and can be opened with Load Previous
Preparation; it is not required to begin with raw CBED.

The real 38 px-radius live smoke uses g1=(81,9), g2=(-9,82), central/symmetry-center origin
(127,126), physical selection 76×76 and binning 2. Its identity-basis corners are
A=(-0.5008180871634687,-0.3962516733601071),
B=(0.40175516882344187,-0.4953145917001339),
C=(-0.40175516882344187,0.4953145917001339). These use the current live basis,
not the historical DM sign/basis convention or archived corners. The explicit
[001], (200)/(020) B controls map every live (n1,n2) to (2n1,2n2,0).
The smoke stops at stage dispatch with zero optimizer evaluations: its checkpoint
is marked failed by the intentional test exception, not a completed refinement.
No-aperture Poisson uncertainty remains labelled as an approximation.


## Final new-user acceptance scope

The real FCC Al acceptance starts with raw `reference_data/pattern.dat`, no prepared
metadata record. Automated Qt control/event checks draw a 76 px diameter (radius is derived,
not a detector constant), accept the local symmetry ROI, and launch the live preparation
process. Nine linked beam views, zone-axis marker, weight-only pixel exclusion,
validated descriptor replay, physical HK mapping, atomic projections and A/B/C
validation are checked. A separate initialization worker stops before optimizer entry
with zero evaluations. This does not constitute another refinement campaign.

The same check runs in the frozen application. It automates the Qt interactions;
literal human Finder/file-picker usability is not claimed by this test. No matched
real aperture is available: aperture calibration/correction/noise are covered by
existing component and orchestration tests, not a matched real-data acceptance.
The real acceptance therefore uses the explicit no-aperture Poisson approximation.
A7 global lattice fitting remains unavailable with sparse arcs; predicted sites
remain valid for extraction subject to unchanged ROI coverage requirements.

Final milestone acceptance (2026-10-04): **PASS** for the real raw-CBED
manual-diameter workflow in source and frozen QCBED.app. Focused tests: 128;
GUI tests: 102; full suite: 2880 passed. Ad-hoc deep/strict signature verified.
Run initialization performed zero optimizer evaluations; its deliberate pre-stage
stop is labelled failed in the test checkpoint, not a completed refinement.
Evidence: `tmp/final_new_user_acceptance/` (generated, not required runtime input).
The aperture and manual file-picker limitations above remain explicit.


### Subpixel beam centers in live preparation

Live preparation samples the nearest integer detector-grid window (NumPy ties-to-even),
without interpolation. The physical center, radii, reciprocal basis and scan corners
remain floating point and unchanged. Each beam's ideal origin, sampled bounds and
registration offset are recorded in preparation provenance; the mask ellipse retains
its measured position within the sampled window. This is an explicit nearest-pixel
registration approximation, with offsets at most half a pixel per axis, not a claim
of exact subpixel registration or legacy DigitalMicrograph interpolation equivalence.
The standalone legacy extraction API defaults to strict integer origins as before.


### Selected-beam inspection and masking

After live preparation, A opens the inspection page automatically. The whole-pattern
overview remains above the selected beam's vertically stacked Pattern, Sigma and
Weight panels. Clicking a beam updates all three. Click a pixel or drag a small
region, then use Exclude, Restore Selected or Reset Manual Mask. Exclusion changes
only final inclusion weight: the defect remains visible in Pattern and the original
variance/Sigma are untouched. Weight and active-point counts refresh immediately.
Validate A captures the overrides; Save/Save As and prepared replay retain them.
The compact replay header leaves room for the scrollable inspector; switching to
New Data Preparation retains the disk and centrosymmetry selection tools.


### Resizable whole-pattern workspace

A uses a vertical splitter: drag its separator to allocate space between live
preparation and selected-beam inspection. New preparation gives priority to the
whole-CBED image (minimum canvas 520×520), with compact controls alongside. The
canvas grows with the window instead of being capped at a small fixed height.
Both the interactive image and prepared overview have Home, Pan and Zoom tools.
Turn Pan/Zoom off before drawing a diameter, selecting a feature or picking a beam.
Replay prioritizes the large prepared overview and the vertical Pattern/Sigma/Weight
inspector; the splitter remains adjustable. Numerical behavior is unchanged.


### Modern data ownership

Modern configurations have no separate Experimental Data panel or low-level data
file picker. A owns raw CBED/aperture inputs, saved-preparation loading, derived scan
geometry, beam inspection, symmetry selection and manual masking. C owns accelerating
voltage and signal-mode settings. Changing the signal-mode setting does not transform
existing arrays: prepare or reopen matching data in A. The legacy configuration
editor and internal experimental-file/FOR007 loaders remain compatible.

### C sub-tabs

Validate A and B, then **Simulation Setup → Validate Simulation Setup → Refinement
Targets → Optimizer → Validate C**, followed by **D → Execution & Output**. Check
Refine and enter physical hard bounds; recommended search metadata is supplied
automatically. The normal target table hides D1/D2/PLAMDA. Expert overrides are in
the collapsed Expert / Legacy settings section. Compatible committed bindings
survive setup changes, but require re-validation. D provides separate writable-directory validation and Fresh Run / Resume status. See [the user guide](user_manual.md).

### Parameter visibility milestone

All expandable controls use one display-only section widget. A's disk-profile
override, binning and exact manual-mask rectangle coordinates are Advanced; live
selection and mask layers remain Basic. B's format override/site indices are
Advanced and complete compiler diagnostics are Expert. C's numerical controls,
optimizer budgets and execution policy are Advanced; full canonical state and
D1/D2/PLAMDA are Expert. Results still read existing artifacts only.

Slice boundaries are always visible in Simulation Setup. Following the
[backend audit](one_slice_simulation.md), `0, 1` is a valid one-slice period for
supported strict geometry; the earlier incomplete-draft rule is superseded.
Transverse mesh and cutoff remain Advanced. The
[control audit](gui_parameter_classification.md) records classifications.

C validation covers simulation, targets and optimizer configuration without requiring
output selection. D exclusively owns output validation and Advanced execution settings
(mode, batch size and run limits). A/B/C validity alone does not enable Run: D must
validate the writable output directory. Resume compatibility is checked separately.
