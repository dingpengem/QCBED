# Example study: FCC Al from raw CBED to QCBED results

[Quick Start](QUICK_START.md) · [User Manual](user_manual.md)

## Scientific goal and status

Prepare an experimental CBED dataset, attach a conventional FCC aluminium atomic
model, select requested sample depth and group Debye–Waller B for refinement, then
inspect recorded results. This is an **instructional workflow**, not a claim of a
scientifically converged aluminium refinement.

The runnable GUI template is [examples/fcc_al.yml](../examples/fcc_al.yml). It starts
from raw input, so a fresh checkout needs no prior `tmp/` campaigns. Its two-evaluation
safety limit is intentional. Release acceptance uses an even stricter stop before
the first optimizer evaluation; it does not manufacture a new convergence history.

## Exact reference inputs

| File | Role |
|---|---|
| [reference_data/pattern.dat](../reference_data/pattern.dat) | Raw full CBED image, numeric 256 × 256 detector array |
| [reference_data/Al_FCC.cif](../reference_data/Al_FCC.cif) | Explicit four-site conventional cubic FCC Al cell, a = 4.03893 Å |
| [examples/fcc_al.yml](../examples/fcc_al.yml) | Raw-project template and explicit study choices |

The CIF lists four explicit sites using a P1 symmetry listing; their coordinates
form the conventional FCC motif. Do not replace the supplied cell with a primitive
cell while retaining the same hkl labels. The study uses the **manual disk-diameter
path**, with **no aperture or detector PSF correction file**. `reference_data/sigma.dat`
and legacy output files are preserved reference material, not substitutes for the
newly prepared gradient uncertainty.

## Open the template

Launch QCBED and choose **Load Configuration… → examples/fcc_al.yml**. The raw
project opens A's interactive pattern view. From source the equivalent launch is:

```sh
python -m qcbed.gui examples/fcc_al.yml
```

The template stores physical choices but no derived matrices, optimizer history,
validation receipts, editable output filename or preparation arrays. Later choose
**Save As** to save your prepared project, leaving this template unchanged.

## A. Data Prep: reproduce the selections

1. The loaded **Full CBED** is `pattern.dat`. Select **Draw Disk Diameter**.
2. In detector x/y coordinates, drag from **(89, 126)** to **(165, 126)**. Click
   **Accept geometry / Snap ROI preview**. This selects a 76-pixel diameter, radius
   38 pixels; it is the study's observed disk geometry, not a global default.
3. Select **Select Centrosymmetric Feature**. Drag the rectangle with corners
   **(111, 105)** and **(155, 151)** and accept the snapped preview. The outward
   even-grid rectangle has the same bounds here. These are raw-detector coordinates,
   before binning. The symmetry search subsequently measures its origin.
4. Leave the established binning at **2** and signal mode at **gradient**; click
   **Run Data Preparation**. The 76 × 76 extraction becomes a **38 × 38** scan.
5. Inspect the detected/refined **nine-beam** catalogue. It is measured from the
   current input; the software does not hardcode nine beams for arbitrary images.
   Inspect the full-pattern reference beam and origin, then click
   **Confirm recorded zone-axis position**. That explicit human confirmation is
   part of A, not a fitted output parameter.
6. Inspect **Pattern**, **Sigma**, and **Final Weight** for each beam. Review
   **Geometric Mask**, **Processing Validity** and **Manual Mask** separately.
   For this baseline study add no manual exclusions. If you add exclusions for a
   scientific reason, record them and treat the resulting dataset as a different
   input state. Click **Validate A**.

The prepared image, uncertainty, ordered beam data, origin and mask provenance are
written into `.qcbed-data` beside the project. Numerical signal values and final
active-point counts are determined by this preparation, not copied from historical
B7 campaigns. A different drawn diameter/ROI can change the results. Use the
recorded descriptor if exact replay is needed; do not claim equivalence merely
because two previews look similar.

## B. Atomic Model

The template references **Al_FCC.cif** with:

- zone axis **[0, 0, 1]**;
- **g1 hkl [2, 0, 0]**, **g2 hkl [0, 2, 0]** in the conventional cell;
- physical site group **Al_fcc**, original sites **0, 1, 2, 3**;
- **Debye–Waller B = 0.86 Å²**, occupancy **1**, absorption **disabled**.

If building the project from a blank session instead, browse to that CIF, enter
these basis vectors and manually enter B = 0.86 Å². Loading the CIF alone does not
supply a universal thermal factor. Inspect the lattice/sites and group membership,
click **Apply Site-Group Parameters**, then **Validate B**. Confirm that g1/g2 stay
[2,0,0]/[0,2,0] after applying the group.

## C. Simulation Setup

The template and this tutorial use exactly:

| Control | Study setting |
|---|---|
| Accelerating voltage | 200 kV |
| Signal mode | Gradient, matching A |
| Requested depth | 1000 Å (the fresh GUI starting value) |
| Slice boundaries | `0, 1` |
| Transverse simulation mesh X/Y | 16 / 16 points |
| Reciprocal cutoff | 2.5 Å⁻¹ |

`0, 1` means **one stored slice** for the periodic cell. `0, 0.5, 1` would mean
**two stored slices** and is a separate numerical choice; do not switch silently
when comparing runs. Requested depth is expanded through the existing whole-slice
propagation rule. Strict transverse periodicity still applies in either case.

Mesh X/Y samples the zone-axis-recast transverse projected-potential/wavefield cell,
not detector pixels or the 38 × 38 experimental scan. The 16 × 16 starting mesh is
an inexpensive teaching choice already used by repository examples, **not a proven
converged mesh**. Leave the normal numerical settings for this workflow exercise;
a scientific study must check mesh, cutoff and slicing sensitivity independently.
Click **Validate Simulation Setup**.

## Refinement Targets and optimizer

Select **Refine** for both supported physical targets:

| Target | Starting value | Example study lower / upper bounds |
|---|---:|---:|
| Requested sample depth | 1000 Å | 600 / 1200 Å |
| Al_fcc Debye–Waller B | 0.86 Å² | 0.5 / 1.5 Å² |

These are **example study bounds**, drawn from the project's documented bring-up
settings; they are not canonical/universal aluminium recommendations. Starting
values come from B and Simulation Setup. Click **Validate Refinement Targets**.

Choose **AMOEBA** and retain **Use recommended optimizer settings**. The compact
YAML target entries `{}` request the existing target-specific metadata, avoiding
redundant D1/D2/PLAMDA in the teaching file. No optimizer mathematics or defaults
are changed. Click **Validate C**. A/B/C can now be valid before choosing output.

## D. Output, execution and Run

The template's output container is **runs/fcc_al**, relative to the repository root.
You may choose another writable folder in D; click **Validate Output Directory**.
It may be missing or nonempty. QCBED creates missing folders and each fresh GUI run
gets a new managed child, preserving existing files.

Check **Output directory: Valid**, **Run: Ready**, and, in a fresh folder,
**Resume: Not available — no previous run found**. Save your prepared project with
**Save As**, for example `runs/fcc_al_project.yml`. Keep its cache/source dependencies.

Optional: **Advanced execution settings** shows **batched / 64**. This is the normal
modern GUI setting; no user action is required. The template also explicitly sets
**maximum total evaluations = 2** and **maximum wall time = 60 seconds**. The normal
optimizer defaults are retained, but the safety cap stops this instructional run
before a meaningful convergence search. The cap is per invocation, including Resume.

Click **Run QCBED** only when you want that bounded illustrative evaluation. It is
normal to finish with **Interrupted** due to the safety budget. A wall-time limit
is checked at safe boundaries, not by killing a mid-propagation calculation.
**Cancel** also waits for a safe boundary. The release acceptance harness stops
before evaluation; it creates no new scientifically interpreted fit.

## Interpret Results without overstating them

**Open Results** follows the actual `qcbed-run-…` child. Inspect:

- **Overview:** status, active data, signal mode, evaluation count and budget reason.
- **Refined Parameters:** initial/latest/best/retained values and bounds. A one-point
  history cannot establish a parameter optimum or uncertainties.
- **χ² History** and **Physical / YNORM history:** agreement and normalization over
  recorded evaluations. Compare only like-for-like data, variance and masks.
- **Beam details / Pattern Matching:** experimental, calculated and residual signals
  for each recorded beam/evaluation. Excluded pixels do not enter the objective.
- **Report / Files:** provenance and recorded artifacts; not a command to rerun physics.

A zero-evaluation initialization has no simulated comparison. A cancelled/budget-
limited run may have latest/best values but no completed retained final result.
**Resume QCBED** requires the compatible checkpoint and unchanged identity; it adds
another invocation's budget. **Run QCBED** creates a separate fresh child instead.

Historical numerical results are documented in [B7 bring-up](b7_refinement_bringup.md)
and its referenced configurations. They used different depth starts, slice choices,
preparation/provenance and budgets; they are not outputs of this new tutorial and
are not evidence that this two-evaluation study converged. The user need not possess
those local campaign directories to complete preparation, save/reload or the
bounded teaching workflow here.
