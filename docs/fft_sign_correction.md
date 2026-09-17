# Source-backed wave FFT correction

The 2026-09-16 correction changes only physical wave-transform roles. It does
not tune physics parameters, change the objective, optimize execution, or
resume thickness acceleration. Potential construction is unchanged.

## Source evidence and operation table

`private predecessor implementation` dispatches ISGN >= 0 to CFTFSUB and negative ISGN to
CFTBSUB. The executable four-point butterflies at 2838-2876 give positive
and negative exponential sums respectively. `private predecessor implementation`
applies the same sign to both axes without normalization. Their hashes are
recorded in [the reconciliation audit](fft_replacement_audit.md).

`private predecessor implementation` packs real/imaginary pairs into DATA and unpacks them
without conjugation, scaling or transposition. The normalized Python grid is
`(mesh_y,mesh_x)`, axis 0 = y/K and axis 1 = x/H, with C-order flattening.
Fortran operates in place; Python retains immutable input ownership and
returns independent arrays. This does not change the mathematical transform.

| Operation | Exponent | Normalization | NumPy operation |
| --- | --- | --- | --- |
| Potential coefficients to V | negative | none | fft2, norm="backward" |
| Real wave to reciprocal wave | negative | 1/(mesh_x*mesh_y) | fft2, norm="forward" |
| Reciprocal wave to real wave | positive | none | ifft2, norm="forward" |

These operations differ because their callers differ. SFMAP2 -> OSHIFT ->
FFT(-1) -> PGMAP (`private predecessor implementation`, `private predecessor implementation`) has no RMT.
The wave loop uses FFT(-1), then RMT*PMAP*PRARRAY, samples beams, then
FFT(+1) if propagation continues (`private predecessor implementation`,
`private predecessor implementation`). RMT=1/DBLE(MT) at `private predecessor implementation`.
Python incorporates this factor in the first wave transform. Its floating-point
operation order need not be bitwise identical to external Fortran scaling.

## Unchanged reciprocal and physics interpretation

Beam positions remain `(K mod mesh_y,H mod mesh_x)` from
`private predecessor implementation`, without H/K negation. PROPG at `private predecessor implementation`
uses `A1=2*TILT(1)+H`, `B1=2*TILT(2)+K`, negative SG, and
`TEMP=-TWOPI*(SG*DZ*APPLIEDDILATION+H*X0+K*Y0)`.
Thus tilt coefficients, propagation phase, negative fractional-shift phase,
OSHIFT ordering, and excluded negative-Nyquist row/column remain unchanged.
There is no compensating beam reversal or visual flip.

Odd/non-power-of-two meshes are mathematical Python extensions tested with
scalar DFTs. The restored FFT library itself requires power-of-two dimensions;
these extension tests do not establish odd-grid Fortran support.

## Implementation and independent validation

`fft_backend.py` adds explicit single-grid and batched wave helpers. Existing
generic potential-sum/inverse APIs retain their contracts. `wave_propagation.py`
delegates its existing public functions to the new wave helpers;
`tilt_series.py` uses the batched wave helpers. The process path already calls
the batched engine and needs no scheduling change. `phase_grating.py` has only
docstring corrections; its numerical expressions are unchanged.

The literal source CFTB040 example `[0,1,0,0]` maps to
`[.25,-.25j,-.25,.25j]` after RMT. An asymmetric complex input `[1,1j,0,0]`
has labeled intensities `[.125,.25,.125,0]`; this distinguishes the signs
without relying on a round trip or a symmetric diffraction pattern.

Existing nested scalar DFT expectations now use the source signs. Explicit
positive spatial plane waves sample positive H/K bins. Independent tests use
odd/even rectangular meshes, varying gratings, non-even kernels, nonzero signed
tilts, shifts, repeated/multiple slice types, complete complex intermediate
waves and beam sampling. Printed-coefficient small models remain independent
of production scattering/grating/propagator construction. Real Al FOR007
geometry, atoms, absorption and 16x16 grid are also tested for three occurrences
at zero and positive/negative tilts against scalar nested DFT propagation.

Measured maximum absolute residuals on asymmetric small cases:

| Comparison | Residual |
| --- | ---: |
| Real -> reciprocal versus scalar DFT | 1.052e-15 |
| Reciprocal -> real versus scalar DFT | 3.198e-14 |
| Single-step reciprocal wave | 8.192e-16 |
| Single-step returned real wave | 4.232e-15 |
| Single-step intensity | 2.221e-16 |
| Four-slice reciprocal wave | 9.355e-16 |
| Four-slice intensity | 8.882e-16 |
| Serial versus batched amplitudes/intensities, measured seven-tilt case | 0 / 0 |
| Serial versus spawned-process amplitudes/intensities, same case | 0 / 0 |

The broader multi-tilt matrix passes amplitude atol=3e-15 and intensity
atol=4e-15, rtol=0. It verifies order, duplicate beams/tilts and no cross-talk.
Real-Al short-sequence bounds are 8e-15 complex and 1e-14 intensity.
Measured real-Al maxima were 1.814e-15 for beam amplitudes, 2.443e-15 for
intensities and 3.597e-15 over the complete complex reciprocal wave.
These checks establish agreement with the source-derived mathematical oracle,
not compiled historical CDFT execution.

## Full reference and visualization

`examples/validate_fft_correction.py` runs the full 1521-tilt, nine-beam
reference at 840 angstrom and one deterministic depth trial at 842.0245
angstrom, using the existing serial trial evaluator. It saves complex
amplitudes, raw/GRADIENT/scaled arrays, source-format output, per-beam PNGs and
summary JSON under `tmp/fft_correction/`. No optimizer or reference-file edits
are involved. `--compare-saved` additionally compares the retained prior
`tmp/end_to_end_validation/comparison.npz` without recomputing physics.

| Baseline metric | Before correction | After correction |
| --- | ---: | ---: |
| firstout maximum absolute difference | 0.33890555339166895 | 0.3389055519364774 |
| firstout relative L2 difference | 4.120887163251928e-7 | 4.120887167430503e-7 |
| YNORM | 4723684.377058765 | 4723684.377058742 |
| Masked objective | 1.4154022826662773 | 1.4154022826663093 |

Historical agreement is **comparable**, not materially improved or worsened.
No historical executable-equivalence tolerance has been established. The old
and new saved baseline raw intensities differ by at most 9.565e-14; scaled
maps differ by at most 1.953e-7 (relative L2 1.766e-13).
This centrosymmetric Al case can obscure sign errors; asymmetric complex,
shifted-atom and signed-tilt tests are essential and are not replaced by it.

Identity-oriented firstout comparison has relative L2 4.121e-7. Flipping rows
or columns gives about .957657; inversion gives 1.050755; transposition gives
1.068583. No mirroring, inversion or H/K compensation is justified. All nine
beam comparisons are regenerated using original scan row/column mapping.
Before/after difference panels have their own explicitly labeled color scale,
so tiny numerical differences are not presented as large pattern changes.

The supplied 117 zero variances still make the unchanged legacy objective
undefined (NaN). Legacy YNORM remains finite and agrees with masked YNORM in
this IREF=-1 case. Masked mode retains its explicitly non-legacy policy and
uses the same valid points for YNORM and objective. Neither variances nor
weights were modified.

## Validation status and remaining limits

Interpreter: `C:\Users\peng863\AppData\Local\miniconda3\envs\qcbed-python\python.exe`.
Focused suite: **570 passed in 9.76 s**.
Complete `-m pytest -q`: **1567 passed, 1 failed in 118.25 s**.
The sole failure remains the explicit scattering-table provenance test:
the replaced private predecessor implementation no longer contains ATOMICSCATFACT or its 98 table rows.
Coefficients are unchanged and were previously checked exactly against Git
source `de5b2b5cf08c4a34532b2ce636d15e42caa3907d`. The test was not weakened,
removed, skipped or made to claim current-source verification.

During validation, the user also replaced private predecessor implementation and private predecessor implementation. Their
line-normalized texts match private predecessor implementation and private predecessor implementation respectively, so they
introduce no conflicting wave/GRADIENT equations. They were not edited here;
the duplicated source definitions require future explicit build selection.

The original compiler settings, literal rounding, mutable INT(IP) workspace
interface, duplicate definitions, scattering-source provider and historical
linkage remain unresolved. Original executable complex intermediates and a
historical final refined out.dat are unavailable. Long optimization and all
old fitted histories need separate reruns; a short test suite does not certify
their final fitted parameters. Old performance measurements describe the old
implementation; thickness-response acceleration remains deferred.
