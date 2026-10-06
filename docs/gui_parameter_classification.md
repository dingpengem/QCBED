# GUI parameter classification audit

> Historical milestone/audit record. For current controls and the A/B/C/D workflow,
> use the [User Manual](user_manual.md); dated counts and replay inputs below retain
> their original provenance.

Audit captured before the app-wide visibility changes (2026-10-05). Sources: `gui/document.py`, `preparation_view.py`, `data_view.py`, `modern_widgets.py`, `window.py`, `optimizer_panel.py`, `results_view.py`.

Original audit counts were Basic: 37, Expert: 13, Advanced: 21 control families (not dynamic row/widget instances). The later slicing correction moves Slice boundaries to Basic: updated counts are Basic: 38, Expert: 13, Advanced: 20. Repeated beam/site/target instances retain dynamic counts. Current values below are the loaded FCC Al example, **not universal defaults**. Compound rows identify coupled actions/display selectors rather than physical scalars. Backend-only controls not exposed today are listed separately below.

| Current label | Backend field | Meaning/section | Units | Current value/default | Safe default? | Before visibility | Classification | Rationale |
|---|---|---|---|---|---|---|---|---|
| Project title | title | Project | text | FCC Al grouped gradient baseline | Existing stored/default value; no reset | Visible | Basic | Physical/project choice; keep visible. |
| Energy (kV) | microscope.energy_kv | Microscope | kV | 200.0 | Existing stored/default value; no reset | Visible | Basic | Physical/project choice; keep visible. |
| Scan columns / NUMP | scan.scan_columns | Scan | tilt columns | 38 | Existing stored/default value; no reset | Visible | Expert | Prepared geometry; editing can invalidate A. Prefer recorded live geometry. |
| Scan rows / NUML | scan.scan_rows | Scan | tilt rows | 38 | Existing stored/default value; no reset | Visible | Expert | Prepared geometry; editing can invalidate A. Prefer recorded live geometry. |
| Corner A (2 components) | scan.corner_a | Scan | reciprocal-basis coefficients | -0.50623966320854, 0.4005412719891745 | Existing stored/default value; no reset | Visible | Expert | Prepared geometry; editing can invalidate A. Prefer recorded live geometry. |
| Corner B (2 components) | scan.corner_b | Scan | reciprocal-basis coefficients | 0.39497819876710266, 0.5006765899864682 | Existing stored/default value; no reset | Visible | Expert | Prepared geometry; editing can invalidate A. Prefer recorded live geometry. |
| Corner C (2 components) | scan.corner_c | Scan | reciprocal-basis coefficients | -0.39497819876710266, -0.5006765899864682 | Existing stored/default value; no reset | Visible | Expert | Prepared geometry; editing can invalidate A. Prefer recorded live geometry. |
| Execution mode | execution.execution_backend | D / Advanced execution settings | enum | batched | Existing stored/default value; no reset | Visible | Advanced | Specialist numerical/search or file-format setting; preserve configured value. |
| Batch size (blank for serial) | execution.batch_size | D / Advanced execution settings | tilts per batch | 64 | Existing stored/default value; no reset | Visible | Advanced | Specialist numerical/search or file-format setting; preserve configured value. |
| Output directory | run.output_directory | Outputs | filesystem path | tmp/modern_grouped_gui_run | Existing stored/default value; no reset | Visible | Basic | Physical/project choice; keep visible. |
| Structure file | crystallography.structure.path | Structure | filesystem path | examples/fcc_al_conventional.extxyz | No universal default | Visible | Basic | Physical/project choice; keep visible. |
| Structure format (blank: infer from file) | crystallography.structure.format | Structure | ASE format name | extxyz | Existing stored/default value; no reset | Visible | Advanced | Specialist numerical/search or file-format setting; preserve configured value. |
| Zone axis [u,v,w] | crystallography.zone_axis | Structure | integer [u,v,w] | 0, 0, 1 | No universal default | Visible | Basic | Physical/project choice; keep visible. |
| g1 hkl (original/conventional cell) | crystallography.cbed_basis.g1_hkl | Structure | integer (h,k,l) | 2, 0, 0 | No universal default | Visible | Basic | Physical/project choice; keep visible. |
| g2 hkl (original/conventional cell) | crystallography.cbed_basis.g2_hkl | Structure | integer (h,k,l) | 0, 2, 0 | No universal default | Visible | Basic | Physical/project choice; keep visible. |
| Transverse simulation mesh X (points) | crystallography.mesh.mesh_x | Numerical / Simulation cell | points | 16 | Existing stored/default value; no reset | Visible | Advanced | Specialist numerical/search or file-format setting; preserve configured value. |
| Transverse simulation mesh Y (points) | crystallography.mesh.mesh_y | Numerical / Simulation cell | points | 16 | Existing stored/default value; no reset | Visible | Advanced | Specialist numerical/search or file-format setting; preserve configured value. |
| Reciprocal cutoff (Å⁻¹) | crystallography.dsap_inverse_angstrom | Numerical / Simulation cell | Å⁻¹ | 2.5 | Existing stored/default value; no reset | Visible | Advanced | Specialist numerical/search or file-format setting; preserve configured value. |
| Slice boundaries | crystallography.slice_boundaries_fractional | Simulation Setup | fraction of propagation period | Fresh: 0, 1; loaded values retained | Valid one-slice choice; convergence requires review | Always visible | Basic | Required definition of the slicing model; strict geometry proof retained. |
| Requested depth (Å) | crystallography.depth_angstrom | Numerical / Simulation cell | Å | 840.0 | No universal default | Visible | Basic | Physical/project choice; keep visible. |
| Signal mode | experiment.gradient_option | Microscope | conventional/gradient | 1 | Existing stored/default value; no reset | Visible | Basic | Physical/project choice; keep visible. |
| Workflow | refinement.optimizer.method | Optimizer / Run | workflow enum | ANSI | Existing stored/default value; no reset | Visible | Basic | Physical/project choice; keep visible. |
| Loops | refinement.optimizer.loops | Optimizer / Run | stage cycles | 1 | Existing stored/default value; no reset | Visible | Advanced | Specialist numerical/search or file-format setting; preserve configured value. |
| Stage iteration / soft evaluation budget | refinement.optimizer.maxfsimp | Optimizer / Run | iterations / soft evaluations | 59 | Existing stored/default value; no reset | Visible | Advanced | Specialist numerical/search or file-format setting; preserve configured value. |
| Coordinate-search budget | refinement.optimizer.maxfpara | Optimizer / Run | coordinate-search cycles | 3 | Existing stored/default value; no reset | Visible | Advanced | Specialist numerical/search or file-format setting; preserve configured value. |
| Random seed | refinement.optimizer.seed | Optimizer / Run | integer RNG seed | 2026 | Existing stored/default value; no reset | Visible | Advanced | Specialist numerical/search or file-format setting; preserve configured value. |
| ANNEAL maximum temperature | refinement.optimizer.maxtemp | Optimizer / Run | optimizer temperature | 0.05 | Existing stored/default value; no reset | Visible | Advanced | Specialist numerical/search or file-format setting; preserve configured value. |
| ANNEAL minimum temperature | refinement.optimizer.mintemp | Optimizer / Run | optimizer temperature | 0.0001 | Existing stored/default value; no reset | Visible | Advanced | Specialist numerical/search or file-format setting; preserve configured value. |
| Cooling rate | refinement.optimizer.coolrate | Optimizer / Run | dimensionless ratio | 0.9 | Existing stored/default value; no reset | Visible | Advanced | Specialist numerical/search or file-format setting; preserve configured value. |
| Objective-spread tolerance | refinement.optimizer.ftol | Optimizer / Run | relative objective tolerance | 1e-06 | Existing stored/default value; no reset | Visible | Advanced | Specialist numerical/search or file-format setting; preserve configured value. |
| Minimum objective | refinement.optimizer.chimin | Optimizer / Run | objective value | 1e-06 | Existing stored/default value; no reset | Visible | Advanced | Specialist numerical/search or file-format setting; preserve configured value. |
| Full CBED | NewPreparationView.inputs.pattern | Input image | path | None | User input | Visible | Basic | Separate optional override; never infer unrelated image scales. |
| Aperture (optional) | NewPreparationView.inputs.calibration_aperture | Input image | path | None | User input | Visible | Basic | Separate optional override; never infer unrelated image scales. |
| Disk-profile override (optional) | NewPreparationView.inputs.disk_profile | Input image | path | None | User input | Visible | Advanced | Separate optional override; never infer unrelated image scales. |
| New Data Preparation / Load Previous Preparation | preparation.mode | Fresh preparation or replay | choice | New | No numerical substitution | Visible | Basic | Physical choice, operation, or read-only result. |
| Draw Disk Diameter / Select Centrosymmetric Feature | preparation.tool | Measured geometry and local symmetry | px | User drawn | No numerical substitution | Visible | Basic | Physical choice, operation, or read-only result. |
| Accept geometry / Snap ROI preview | preparation.accept | Confirm finalized geometry | px | None | No numerical substitution | Visible | Basic | Physical choice, operation, or read-only result. |
| Signal mode | preparation.signal | Prepare conventional or gradient signal | choice | conventional | No numerical substitution | Visible | Basic | Physical choice, operation, or read-only result. |
| Run Data Preparation | preparation.start | Recompute live arrays | action | None | No numerical substitution | Visible | Basic | Physical choice, operation, or read-only result. |
| Selected beam | data_panel.beams | Select beam for linked images | identity | First | No numerical substitution | Visible | Basic | Physical choice, operation, or read-only result. |
| Final Weight / Geometric Mask / Processing Validity / Manual Mask | data_panel.mask_layer | Explain binary inclusion | 0/1 | Final Weight | No numerical substitution | Visible | Basic | Physical choice, operation, or read-only result. |
| Exclude / Restore Selected / Reset Manual Mask | data_panel.mask | Edit only manual weight overrides | pixels | No overrides | No numerical substitution | Visible | Basic | Physical choice, operation, or read-only result. |
| Select symmetry center / Confirm recorded zone-axis position | data_panel.select_zone | Explicit origin confirmation | px | Unconfirmed | No numerical substitution | Visible | Basic | Physical choice, operation, or read-only result. |
| Group ID | site_groups.id | Stable physical group name | text | Structure generated | No numerical substitution | Visible | Basic | Physical choice, operation, or read-only result. |
| Element | site_groups.element | Chemical element | symbol | CIF | No numerical substitution | Visible | Basic | Physical choice, operation, or read-only result. |
| Debye–Waller B (Å²) | site_groups.debye_waller_b | Thermal parameter | Å² | Blank unless supplied | No numerical substitution | Visible | Basic | Physical choice, operation, or read-only result. |
| Occupancy | site_groups.occupancy | Site occupation | fraction | CIF, otherwise existing 1 | No numerical substitution | Visible | Basic | Physical choice, operation, or read-only result. |
| Absorption (a,b,c; blank=disabled) | site_groups.absorption | Physical absorption coefficients | legacy coefficient units | Blank/loaded | No numerical substitution | Visible | Basic | Physical choice, operation, or read-only result. |
| Apply Site-Group Parameters | apply_site_groups | Commit reviewed physical metadata | action | None | No numerical substitution | Visible | Basic | Physical choice, operation, or read-only result. |
| Front / Top / Left / 3D | atomic_axes/atomic_3d | Same beam-frame structural views | Å | Current structure | No numerical substitution | Visible | Basic | Physical choice, operation, or read-only result. |
| Refine | refinement.parameters | Physical target selection | boolean | Loaded/unchecked | No numerical substitution | Visible | Basic | Physical choice, operation, or read-only result. |
| Physical parameter / Starting value / Status | available_refinement_targets | Committed physical target | target units | Derived | No numerical substitution | Visible | Basic | Physical choice, operation, or read-only result. |
| Hard lower bound / Hard upper bound | refinement.bounds.specific | Allowed physical domain | target units | Explicit only | No numerical substitution | Visible | Basic | Physical choice, operation, or read-only result. |
| Run / Resume / Cancel | RunController | Execute immutable snapshot safely | action | Disabled until ready | No numerical substitution | Visible | Basic | Physical choice, operation, or read-only result. |
| Inspect evaluation / physical target / beam | ResultsView | Read recorded results | identity | Latest | No numerical substitution | Visible | Basic | Physical choice, operation, or read-only result. |
| Common Experimental/Calculated scale | ResultsView.common_scale | Display-only scale | boolean | True | No numerical substitution | Visible | Basic | Physical choice, operation, or read-only result. |
| Open output / report / result / checkpoint / snapshot | ResultsView.file_buttons | Open artifacts without modification | action | Unavailable until recorded | No numerical substitution | Visible | Basic | Physical choice, operation, or read-only result. |
| Load / Save / Save As / Undo / Redo / Revert | ConfigWindow.document | One project and immutable editing history | action | Unsaved | No numerical substitution | Visible | Basic | Physical choice, operation, or read-only result. |
| Validate A / B / Simulation / Refinement / C | Workflow | Stage receipts and readiness | action | Pending | No numerical substitution | Visible | Basic | Physical choice, operation, or read-only result. |
| Binning | LiveSelection.final_binning | Block processing size | pixels | 2 in GUI | Existing value | Visible | Advanced | Specialist precision/display control; basic drawing remains available. |
| x/y start/stop | data_panel.box | Exact manual-mask rectangle | pixels | (0,0,1,1) | Existing value | Visible | Advanced | Specialist precision/display control; basic drawing remains available. |
| Original sites | site_groups.site_indices | Review/edit group membership | indices | Structure generated | Existing value | Visible | Advanced | Specialist precision/display control; basic drawing remains available. |
| Zoom / pan / save figure | matplotlib toolbar | Display-only navigation | display | Current view | Existing value | Visible | Advanced | Specialist precision/display control; basic drawing remains available. |
| Search-range factor D1 | ParameterBinding.range_factors[0] | Scaling; PARAMIN TOP | dimensionless | Category profile | No replacement of explicit input | Collapsed Advanced | Expert | Internal/legacy representation or debugging; explicit opt-in. |
| Search-range factor D2 | ParameterBinding.range_factors[1] | Scaling; PARAMIN BOT | dimensionless | Category profile | No replacement of explicit input | Collapsed Advanced | Expert | Internal/legacy representation or debugging; explicit opt-in. |
| Perturbation parameter PLAMDA | ParameterBinding.plamda | Optimizer-specific perturbation | scaled or physical | Category profile | No replacement of explicit input | Collapsed Advanced | Expert | Internal/legacy representation or debugging; explicit opt-in. |
| All fields (read-only) | QcbedConfig recursive view | Complete canonical/internal state | mixed | Derived/loaded | No replacement of explicit input | Visible navigation | Expert | Internal/legacy representation or debugging; explicit opt-in. |
| Model diagnostics | bundle.diagnostics | Recast/transformation/indexing diagnostics | mixed | Derived | No replacement of explicit input | Visible B summary | Expert | Internal/legacy representation or debugging; explicit opt-in. |
| Capabilities | CAPABILITIES | Implemented/deferred mapping audit | text | Backend defined | No replacement of explicit input | Visible navigation | Expert | Internal/legacy representation or debugging; explicit opt-in. |
| Compatibility import/export/report | legacy menu | Legacy migration and representability | files | Explicit only | No replacement of explicit input | Separate menu | Expert | Internal/legacy representation or debugging; explicit opt-in. |
| Log / traceback / raw event details | run_log | Execution diagnostics | text | Recorded | No replacement of explicit input | Visible log | Expert | Internal/legacy representation or debugging; explicit opt-in. |

## Not currently editable GUI parameters

Gradient support, beam-radius factor, lattice thresholds, PSF/noise tuning, worker count and checkpoint frequency have no modern edit widgets. Their existing backend/default/recorded values remain unchanged; this milestone does not invent new scientific controls. They classify as Advanced if exposed later. Internal HK, transform matrices, recast coordinates and parameter indices are derived/read-only Expert information. Unsupported vacancy/original-site/original-HKL mappings stay disabled with reasons; no new physics is enabled.

## Default limitation discovered

At audit time a new preparation draft used the existing 16×16 mesh and 2.5 Å⁻¹ cutoff, but **empty slice boundaries** (superseded by the explicit user decision recorded below). No structure-independent slice partition exists. A loaded, explicitly validated numerical profile can be retained as recommended settings. An incomplete fresh profile must expose the numerical panel and ask for explicit boundaries; the FCC Al reference uses 0, 0.5, 1, but that must not be injected for an arbitrary CIF. Mesh/cutoff convergence is not guaranteed by a default.

## Frozen pre-milestone check

The four-tab optimizer build succeeded, but its first smoke input (`examples/fcc_al.yml`) did not contain a replayable full-pattern overview. Validate A correctly rejected it. This is an acceptance-input deficiency, not justification to bypass A validation; final packaged acceptance will use a live prepared case.

## Implemented outcome

The audited families retain the classifications above. Shared collapsed sections
now contain Data Processing, exact mask coordinates, structure format/member
indices, numerical settings, optimizer budgets, execution settings, compiler
internals/canonical state and detailed logs. Geometry interaction and result
inspection remain accessible. Native plot zoom/pan stays directly available as a
display tool for accurate interaction rather than becoming a scientific input.

The Simulation Setup correction supersedes the earlier incomplete `0, 1` draft
policy. Slice boundaries are Basic and always visible; `(0,1)` is valid one-slice
support after the [backend audit](one_slice_simulation.md). The exact transverse
periodicity proof remains unchanged. Mesh X/Y and cutoff remain Advanced, with
compiled transverse P/Q dimensions and sampling steps when available.

Fresh modern drafts use the existing batched default (64) and rejecting bounds;
existing project settings are never silently replaced. Historical recovery records
in [macOS validation](../packaging/macos/VALIDATION.md) describe the policy at that
time; the later one-slice entry records the current behavior.

Execution mode, batch size, run safety limits and output-directory validation are
owned exclusively by D. C contains Simulation Setup, Refinement Targets and
Optimizer. Modern grouped checkpoint/result cadence is fixed at every completed
evaluation; legacy cadence metadata remains preserved in project YAML.
