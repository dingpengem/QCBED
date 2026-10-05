# Quick Start: first QCBED project

[Install first](INSTALLATION.md). Allow roughly 5–10 minutes to read this guide;
data preparation and runs take additional time. For repeatable selections and
values use the [FCC Al study](EXAMPLE_FCC_AL.md). No Advanced/Expert settings are
needed for the normal path.

## 1. Start or load a project

Open QCBED. A blank session starts at **A. Data Prep → New Data Preparation**.
There is no separate "New Project" button: starting a new preparation creates the
project. To use the teaching template instead, choose **Load Configuration…** and
open `examples/fcc_al.yml` from the source/reference archive.

Keep your raw input, structure and project cache available throughout the study.
Use **Save As** for your own copy of a shipped example.

## 2. A. Data Prep

1. Browse **Full CBED** to a single-channel pattern. For FCC Al choose
   `reference_data/pattern.dat`. The example uses no aperture calibration file.
2. Select **Draw Disk Diameter**, drag across a complete disk, then click
   **Accept geometry / Snap ROI preview**.
3. Select **Select Centrosymmetric Feature**, drag a rectangle around a suitable
   local symmetry feature, then accept its snapped preview. This measures the
   origin; it is different from estimating disk radius.
4. Choose **conventional** or **gradient** for the data you intend to refine. The
   FCC Al teaching project uses **gradient**. Click **Run Data Preparation**.
   A new blank project asks for a folder in which to store its preparation cache.
5. Inspect the full-pattern beam catalogue and selected symmetry position.
   Click **Confirm recorded zone-axis position** when correct. Inspect each beam's
   Pattern, Sigma and Final Weight. Exclude detector defects if appropriate;
   exclusion changes weight, not the measured signal.
6. Click **Validate A**. Fix preparation/selection errors before proceeding.

## 3. B. Atomic Model

1. Load the structure file. In the study it is `reference_data/Al_FCC.cif`.
2. Set zone axis **0, 0, 1**, **g1 hkl = 2, 0, 0**, and **g2 hkl = 0, 2, 0** for
   this example. These indices refer to the supplied conventional cell.
3. Inspect the atomic model and site groups. For the Al group, enter
   **Debye–Waller B = 0.86 Å²** only as the example's reference value; it is not
   inferred from the CIF or a universal aluminium constant.
4. Click **Apply Site-Group Parameters**, then **Validate B**. A different specimen
   requires its own lattice, basis, thermal factors, occupancy and absorption model.

## 4. C. Main QCBED Parameters

1. **Simulation Setup:** check voltage, signal mode, requested depth and slice
   boundaries. The teaching project uses **200 kV**, **gradient**, **1000 Å** and
   **0, 1** (one stored slice). Click **Validate Simulation Setup**.
2. **Refinement Targets:** select the physical quantities to refine and supply hard
   lower/upper bounds containing each starting value. The example selects depth
   with **600–1200 Å** and Al-group B with **0.5–1.5 Å²**. These are **example study
   bounds**, not universal recommendations. Click **Validate Refinement Targets**.
3. **Optimizer:** use **AMOEBA** for the teaching example and leave recommended
   settings selected. Click **Validate C**.

A/B/C can all be Valid without an output directory. C does not contain output or
execution controls. Changing scientific inputs invalidates the relevant later
validation steps; revalidate the changed stage before running.

## 5. D. Run & Results

1. Enter or browse to an **Output directory**, then click **Validate Output Directory**.
   QCBED can create a missing directory. A writable folder may already contain files.
2. Check the separate statuses: **Output directory: Valid**, **Run: Ready**, and,
   for a fresh folder, **Resume: Not available — no previous run found**.
3. Choose **Save As** to save your prepared project YAML, retaining its cache/input
   references. Keep the generated `.qcbed-data` folder and referenced source files.
4. Click **Run QCBED**. Each fresh GUI run uses its own subfolder and preserves
   existing files in the selected directory.

The shipped study has a deliberate **two-evaluation safety limit**. An interrupted
partial result is expected, not evidence of convergence. A new project outside
that teaching template uses the normal configured budget. Use **Cancel** to stop
at a safe boundary, and wait for the worker to finish.

## 6. Inspect and continue

Use **Open Results**. Inspect recorded parameters, χ² history and beam pattern
comparisons; excluded pixels do not contribute to the objective. Initial, latest,
best and retained values have different meanings. A cancelled or budget-limited
run can have useful partial results; a run stopped before its first evaluation
has no fitted pattern to display.

**Resume QCBED** is enabled only for compatible existing checkpoint/data identity.
It continues the recorded run. **Run QCBED** starts a new independent run. Neither a
short run nor a smooth history proves a physically correct or converged refinement.
Continue with the [User Manual](user_manual.md) and [FCC Al study](EXAMPLE_FCC_AL.md).
