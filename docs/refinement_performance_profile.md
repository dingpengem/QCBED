# Refinement performance profile

Current batched depth-plus-scan measurements on Apple Silicon, including exact
reuse opportunities, are recorded in [the Stage 1 baseline](apple_silicon_refinement_baseline.md).
No cache has been implemented.

> Follow-up: the wave signs have now been corrected. Numerical findings below
> describe the pre-correction implementation unless explicitly stated otherwise.
> See [correction validation](fft_sign_correction.md). Old profiling timings
> remain measurements, not corrected-physics equivalence claims. Thickness
> acceleration remains deferred.

## Thickness-only refinement

This is an audit and instrumented profile, not an optimization. No production
code, Fortran source, reference data, optimizer, cache or process pool changed.
The existing `evaluate_refinement_trial(..., objective_mode="masked")` executes
all three trials. Profiling instrumentation is confined to
`benchmarks/profile_thickness_only.py`.

### Exact parameter and source semantics

The real FOR007 case has 18 active parameters, one region and one stored slice
type. `depth`, target `(0,)`, is original packed Python index **0**, corresponding
to Fortran **P(1)/PT(1) -> DEPTH(1)** because voltage refinement is disabled.
Starting region depth is **840 angstrom**. The stored slice thickness DELTAZ(1)
is **2.0245 angstrom**, distinct from region depth. Mesh is `(16,16)` and the
scan comprises 1521 tilts and nine measured beams.

An immutable copy switches off every non-depth refinement mask: absorption
A/B/C (including INELB/INELC), DW, scan A/B/C, reflection multipliers, vacancies,
dilation, atom x/y and voltage. Their physical values remain fixed, not reset.
The depth's original D1/D2 range pair is retained. Exactly one active parameter
is asserted, with all other mapped physical fields identical when it changes.
The baseline input file and parsed baseline remain unchanged.

Evidence, in priority order:

* `private predecessor implementation`, ARRANGEPARS: voltage precedes thickness; TMASK selects
  DEPTH entries. `private predecessor implementation`, UPDPARS: assigns `DEPTH(I)=DABS(PT(COUNT+1))`.
* `private predecessor implementation`, ONEITPARA: three parameter/sequence passes, resetting a
  region's THICK to zero and cycling its stored types with MOD. Each occurrence
  adds `DELTAZ(L)*APPLIEDDILATION(L)` while THICK < DEPTH. It sets RSLICES/SLICES
  and total MS. `private predecessor implementation` has the parallel ONEIT construction.
* `private predecessor implementation`: phase gratings rebuilt by stored type, not occurrence.
  Lines 112-127 build PROPG with **DELTAZ(K)**, not DEPTH. `private predecessor implementation`
  uses `DZ*APPLIEDDILATION(IS)` in its propagator phase.
* `private predecessor implementation`: initial grating is MSEQ(1); execution samples beams
  after the final occurrence's reciprocal propagation, when J == MS. Source
  statements following that sample perform an unused transform/next-grating
  access. The existing Python path already stops at the sampled reciprocal
  wave and does not reproduce that unused source tail; this task changes none
  of those established conventions.
* `private predecessor implementation`: general legacy vacancy, dilation and
  surface-DW equations can depend on depth/RSLICES. They prevent a blanket
  thickness-independence claim for all legacy models. The current trial path
  explicitly rejects vacancy/free-surface hazards and G-atom randomization.
* Python `refinement._entries/apply_refinement_parameters`,
  `refinement_trial._prepare_trial/evaluate_refinement_trial`,
  `simulation.prepare_multislice_sample/build_tilt_propagators`, and
  `multislice.simulate_slice_sequence` provide the active execution dataflow.
  No private predecessor implementation-only declaration is used as authoritative evidence.

For this supported subset, applied dilation is one and the three source
sequence passes have the same result. The initial sequence is 415 occurrences
of slice type 0. Increasing depth by **2.0245 angstrom** to **842.0245 angstrom**
adds one occurrence (416). It does not modify per-slice thickness, nor create
a partial slice. A structural-only check at 840.1 angstrom retains 415
occurrences; that fourth depth is not an additional timed forward evaluation.
Use the literal cumulative while-loop when constructing sequence keys: a new
rounding/ceil rule could differ at floating-point boundaries.

Changing depth changes occurrence count, total executed MSEQ length and final
beam-sampling position. A multi-type region cycles the same fixed grating types
and can change the terminal type; here the sole grating type remains 0. No
independent scalar "region repetition" multiplier or final partial propagation
exists. Increasing requested thickness within the same sequence interval leaves
the supported deterministic forward calculation unchanged. The thickness
objective is therefore piecewise constant under these conventions, not a
smooth scalar-rescaling model.

### Dependency table

All answers below apply only to this fixed S-atom, zero-vacancy, no-surface case.

| Quantity | Thickness dependent? | Current behavior / reason |
| --- | --- | --- |
| Electron wavelength | No | Fixed voltage; recomputed in preparation and called through other physics helpers. |
| Reciprocal lattice and grids | No | Fixed direct cell, mesh, ordering and applied dilation; reconstructed per trial. |
| Elastic structure factors | No | Fixed atoms, occupancies, DW, geometry and SFFAC; recalculated per trial. |
| Absorptive structure factors | No | Fixed absorption coefficients, voltage, DW and geometry; recalculated per trial. |
| Phase gratings | No | Fixed projected structure factors; no region-depth scalar is an argument; rebuilt per trial. |
| Slice sequence / occurrence count | Yes, discretely | Literal cumulative depth loop; 415 -> 416 here. |
| Stored grating selection | Types unchanged; occurrence schedule can change | Single type here; multi-type last selection depends on sequence. |
| Per-slice propagator thickness | No | Stored 2.0245 angstrom, not requested region depth. |
| Propagators | No | Fixed slice geometry and tilts; nevertheless rebuilt for all 1521 tilts every trial. |
| Scan tilts | No | A/B/C and bases fixed; generated again per trial. |
| Beam mapping | No | H/K list, mesh and ordering fixed; validated/mapped again during preparation and per tilt. |
| Complex waves / intensities | Yes when sequence changes | Different number of sequential grating/propagator operations. |
| GRADIENT preprocessing | Yes through calculated intensity | Recomputed over simulated maps after all tilts; not applied to experimental arrays. |
| Experimental mask | No | Immutable weights; masked-objective selection rebuilt each trial. |
| YNORM | Yes through simulation | Fixed IREF=-1 selects weighted equal-count normalization over the same valid set. |
| Objective | Yes through simulation/YNORM | Weighted normalized root error over valid points only. |
| Surface correction | Not active | FSTYPE=0; unsupported general surface/vacancy paths remain rejected. |

### Profiling method

Use exactly:

```powershell
Set-Location C:\data\QCBED
$env:OPENBLAS_NUM_THREADS = "1"
$env:OMP_NUM_THREADS = "1"
& "C:\Users\peng863\AppData\Local\miniconda3\envs\qcbed-python\python.exe" benchmarks/profile_thickness_only.py
```

The script calls the unchanged production evaluator three times: initial,
identical repeat, and changed depth. There is no optimizer, cache, reduced mesh,
reduced scan, persistent worker or process startup/serialization. Experimental
parsing occurs once, outside the timed evaluations. First-call wall time is
not full Python startup time: imports and parsing have already completed.

`perf_counter` wrappers time original functions and count calls; temporary
patches are restored after each case. FFT timing wraps NumPy fft2/ifft2, not a
replacement transform. These are **instrumented wall times**, including
wrapper/hashing overhead, not calibrated uninstrumented speed benchmarks.
Nested timings must not be added to their parents. Sample-preparation children
include scattering and grating work; NumPy FFT times include the preparation
FFT as well as the multislice FFTs. The sequence timer includes its existing
applied-parameter/absorption-coefficient metadata setup. All-tilt forward work
is kernel construction plus multislice execution, with outer array assignments
and observer overhead in the unattributed wall remainder.

There is one wall-time sample per case, with no statistical confidence interval.
The full pytest validation overlapped part of the repeated evaluation. Therefore
cross-case wall-time differences are not a controlled speed comparison, nor
evidence of caching or of a cheaper changed-depth calculation. Cache-cost
estimates below use the within-evaluation breakdown of the first case instead.

SHA-256 fingerprints check phase gratings, reciprocal STL grids, scan arrays
and every tilt's complex propagator; wavelength and beam mapping are compared
directly. Both amplitudes and intensities of the repeated trial must match
exactly. The changed trial is not asserted to match the first wave. No forward
result is reused by this profiling script.

### Measured results (2026-09-16)

Interpreter: `C:\Users\peng863\AppData\Local\miniconda3\envs\qcbed-python\python.exe`;
Python 3.11.16, NumPy 2.4.6, Windows build 26200. The complete machine-readable
measurements, counts and hashes are in `tmp/thickness_profile/profile.json`.

| Case | Requested depth (angstrom) | Occurrences | Realized depth (angstrom, rounded) | Wall seconds | Masked objective | YNORM |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Initial | 840 | 415 | 840.1675 | 173.727997 | 1.4154022826662773 | 4723684.377058765 |
| Exact repeat | 840 | 415 | 840.1675 | 166.193672 | 1.4154022826662773 | 4723684.377058765 |
| Changed | 842.0245 | 416 | 842.1920 | 168.032103 | 1.4043594728622661 | 4712936.058252941 |

The repeat matched complex amplitudes and intensities bit-for-bit (maximum
amplitude difference 0). The changed amplitude maximum absolute difference was
0.0176695736949594, a physical thickness change, not an equivalence error.
The lower changed objective is one prescribed probe, **not an optimization**.
All three cases had identical grating, reciprocal-STL-grid, scan and all-tilt
propagator fingerprints, wavelength and beam positions. Baseline and parameter
immutability checks passed.

| Timed component (seconds) | Initial | Repeat | Changed |
| --- | ---: | ---: | ---: |
| Trial mapping | 0.000387 | 0.000845 | 0.000371 |
| Sequence + applied-state metadata | 0.000340 | 0.000610 | 0.000202 |
| Sample preparation | 0.018948 | 0.014695 | 0.011227 |
| Scan generation | 0.000198 | 0.000183 | 0.000100 |
| All propagator construction | 0.215973 | 0.204265 | 0.207568 |
| All-tilt multislice execution | 173.423262 | 165.916422 | 167.752564 |
| GRADIENT preprocessing | 0.001311 | 0.000604 | 0.000634 |
| Masked YNORM + objective | 0.019045 | 0.010881 | 0.013792 |
| Unattributed outer work / instrumentation | 0.048533 | 0.045167 | 0.045645 |
| **All-tilt kernels + multislice subtotal** | **173.639236** | **166.120686** | **167.960131** |
| Nested NumPy fft2 + ifft2 calls, including preparation FFT | 54.830176 | 52.364214 | 52.886708 |
| Nested elastic structure factors | 0.009569 | 0.006283 | 0.005854 |
| Nested absorptive structure factors | 0.007389 | 0.006776 | 0.004436 |
| Nested phase-grating construction | 0.000796 | 0.000531 | 0.000407 |
| Nested wavelength calls | 0.000027 | 0.000031 | 0.000017 |
| Nested beam mapping (preparation + per tilt) | 0.132304 | 0.127280 | 0.128521 |

The first nine rows partition wall time. The subtotal and nested rows must not
be added again. Approximately 118.59 s of first-case multislice time was outside
NumPy FFT calls: this includes products, centered-bin reorderings, allocation,
copies, input checks, sequence validation, beam extraction and Python control
flow. These contributors were **not individually timed**, so this profile does
not establish that any single one dominates that remainder.

### Repeated-work audit and measured speedup limits

Every trial, including the exact repeat, made these calls:

* 1 mapping, 1 sequence/applied-state construction, 1 sample preparation;
* 1 elastic SF, 1 absorptive SF, 1 grating construction, 3 wavelength calls;
* 1 scan generation, 1521 propagator-build calls, 1521 multislice calls;
* 3044 beam mappings (two during preparation, two per tilt);
* 1 GRADIENT and 1 combined masked YNORM/objective computation.

The initial/repeated trial performs 631215 ifft2 calls and 629695 fft2 calls.
One fft2 is grating preparation; multislice therefore performs
1521*(415+414) = **1,260,909 FFTs**. The changed trial performs 632736 ifft2 and
631216 fft2 calls, or **1,263,951 multislice FFTs**: precisely 3042 extra FFTs,
two per tilt for the added occurrence. No final unused inverse FFT is added.

There is no numerical-result cache: unchanged trials rebuild everything and
repropagate the full sequence. Reciprocal grids are also rebuilt inside sample
preparation (not separately timed). The current loop repeats ordering checks,
finite-array validation, sequence validation and beam-index resolution. None
has been bypassed. The experimental file itself is not re-read or transformed
per trial, but masked-objective selection/validation is recomputed.

Measured first-case preparation + propagators + scan generation totals only
**0.235119 s** of 173.727997 s (about **0.1353%**). Even eliminating that whole
cost for free gives an Amdahl upper estimate of only **1.00136x**, not a material
speedup. Wavelength/SF/grating/grid costs are already included in preparation;
adding them again would exaggerate savings. Beam mapping adds roughly 0.13 s
of potential overhead reduction, with a small overlapping preparation part.

The main opportunity is avoiding repeated propagation, which occupied **99.82%**
of the first trial. Exact-sequence forward-result reuse could avoid that work
on repeats or same-bin thickness requests; a resulting end-to-end speedup is
not measured because no cache was implemented. For a one-occurrence extension,
a correct prefix-wave cache could reduce propagated occurrences from 416 to 1;
that is an operation-count opportunity, **not a claimed 416x runtime gain**.
Cache lookup/copying, remaining validation, transformations, mask/YNORM/objective
and memory costs still need measurement in a separately authorized implementation.

Validation: setup check passed; complete suite **1559 passed in 68.60 s**.
`git diff --check` returned 1 only for pre-existing whitespace in protected
`private predecessor implementation` at lines 4, 101 and 105. Its SHA-256 remained
`7905B6D659D2B7EDA135FFBD737C042C5FF6154C2A34DB099102045A9EFD3B51`.
The three new task files pass separate whitespace/newline checks. No files were
staged or committed; existing untracked `reference_data/bestin.dm4` and
`tatus --short` were left untouched.

### Safest prospective reuse architecture (not implemented)

1. Define an immutable fixed-physics preparation keyed by voltage, lateral
   geometry, stored slice thickness/shifts, atoms/occupancies/species, applied
   DW/vacancy/dilation, absorption, reflection multipliers, mesh/order/boundary
   and phase-grating method. Reuse gratings, reciprocal grids and wavelength.
   The existing `PreparedMultisliceSample` also embeds a sequence in its
   parameters; a future separation must not leave that sequence stale.
2. Key scan tilts, beam mapping and per-type/per-tilt propagators by their actual
   inputs. Region depth is not a propagator input in this restricted case.
   No need to rebuild kernels merely because an occurrence is added. Preserve
   source indexing, FFT signs/normalization and immutable complex128 arrays.
3. Rebuild the cheap exact sequence from mapped depth on each trial; use its
   actual tuple (not rounded depth) to select the forward result. Repeated
   identical sequences may reuse amplitudes/intensities if every fixed-physics
   key matches. Always return the new trial's mapped physical state, even when
   its diffraction result is reused. Experimental data/IREF/mode must also be
   keyed if reusing GRADIENT, YNORM or objectives.
4. Start with prepared-state reuse and exact-sequence result reuse. For genuinely
   different sequences, recompute propagation initially. Only after equivalence
   tests should a separate prefix-wave cache be considered: retain the complex
   reciprocal wave for **every tilt and every grid point**, not just the nine
   beam intensities. Resume a longer sequence with the next inverse transform,
   grating and propagator, preserving terminal sampling. Shortening requires a
   cached earlier prefix or recomputation; do not attempt inverse propagation.

A full endpoint wave for this case occupies 1521*16*16*16 = **6,230,016 bytes**
(5.94 MiB), excluding Python metadata. A table of all fixed one-type propagators
has the same raw size. These are storage estimates, not an implemented cache or
a measured speedup. Retaining all 415 prefixes would multiply endpoint storage
by 415; a bounded checkpoint policy would need separate design and tests.

Do not apply this plan to the general legacy vacancy/surface/G-atom paths.
Thickness may alter their derived applied state; the current explicit rejection
must remain. No approximation, JIT, GPU, optimizer change or process pool is
part of this profile. The correctness test for any later fast path must compare
intermediate complex waves and beam amplitudes, not only objective values.
