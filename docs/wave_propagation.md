# One multislice wave step

## Source audit: category B

The operation order is established by the only readable ONEIT
(`private predecessor implementation`) and ONEITPARA (`private predecessor implementation`) loops. Their
propagation bodies are identical. Historical CDFT sign/normalization and
linkage equivalence remain unverified. private predecessor implementation is only a related-version
declaration reference; no capacities or include-only constants are used.

The exact loop, including its unused tail, is:

```text
JPG = MSEQ(1)
PMAP = PGARRAY[JPG]                 # unit real-space incident wave times first grating
J = 1
for N = 1 .. MS:
    PMAP = FFT(-1, PMAP)            # source physical role: real -> reciprocal
    PMAP = RMT * PMAP * PRARRAY[JPG]
    if J == MS:
        UTEST(J)                   # sum |PMAP|^2; does not alter PMAP
        BEAM[:,T] = |PMAP[IB[:]]|^2
    PMAP = FFT(+1, PMAP)            # reciprocal -> real
    J = J + 1
    JPG = MSEQ(J)
    PMAP = PMAP * PGARRAY[JPG]      # next grating, including unused invalid final lookup
```

The initial PMAP is already a real-space wave after the first phase grating;
no explicit incident plane-wave array is constructed there. Tilt is supplied
through PRARRAY. Subsequent gratings are multiplied after the preceding return
transform, equivalently before their own reciprocal propagation. The new
single-step API accepts the wave **before** its grating. Regions affect MSEQ
construction, with no reset or special region-boundary branch in this loop.

Final-slice PRARRAY multiplication is used by beam sampling. Despite the
comment "EXIT HERE", no early exit exists: the final return transform,
MSEQ(MS+1) lookup and extra PGARRAY multiplication still execute, but their
PMAP result is unused. MSEQ is zero-initialized and only 1:MS populated;
the tail can reference PGARRAY(0,...) or exceed sequence capacity. This API
returns one valid step's real wave and does not perform that extra grating or
sequence lookup. Final-slice transform elision belongs to future sequencing.

## Fourier mapping and normalization

FFT packs PMAP(I+(J-1)*MESHX) into complex pairs without a transpose, shift,
conjugation or scaling (`private predecessor implementation`). CDFT2D and CDFT2D_SUB pass ISIGN
to both axis transforms (`private predecessor implementation`), with no wrapper scale.
The restored CDFT dispatch (`private predecessor implementation`) and butterflies
(`private predecessor implementation`) establish negative exponent for -1 and positive for +1,
both unnormalized. The 2D executable statements pass that sign along both
axes. Python now follows those wave-domain signs without changing beam labels.
See [correction validation](fft_sign_correction.md) for downstream results.

The source-backed wave convention is:

```text
u = real_wave * phase_grating
c[H,K] = (1/M) sum_r,c u[r,c] exp(-2*pi*i*(H*c/Nx + K*r/Ny))
s[H,K] = c[H,K] * propagator[H,K]
v[r,c] = sum_H,K s[H,K] exp(+2*pi*i*(H*c/Nx + K*r/Ny))
M = Nx*Ny
```

Real -> reciprocal uses np.fft.fft2(norm="forward"), and reciprocal ->
real uses np.fft.ifft2(norm="forward"). A positive-exponent spatial plane
wave maps to the positive signed H/K bin. The round trip has unit scaling.
Potential construction separately remains a negative-exponent unnormalized
sum; the generic potential/inverse helpers have not changed.

The driver sets MT=MESHX*MESHY and RMT=1/DBLE(MT)
(`private predecessor implementation`). RMT compensates the unnormalized pair described in
the historical specification. Python's first transform already divides by M;
no explicit RMT is required or exposed. Adding it would incorrectly attenuate
every wave by 1/M per step. This establishes Python normalization, not the
original executable's sign or scaling behavior.

## Domains, ordering and ownership

Every array has nonempty shape (mesh_y,mesh_x), complex128, with y/K on rows
and x/H on columns. Real-space samples start at fractional (0,0); phase
gratings are never shifted. PRARRAY is unshifted after PROPG's OSHIFT. The
driver converts signed beam indices to unshifted bins and then one-based
flat IB (`private predecessor implementation`). Thus unshifted is the natural reciprocal
working order, while centered is supported explicitly through backend helpers.

All wave transform calls require ordering. Reciprocal multiplication requires
both wave_ordering and propagator_ordering; unequal labels are rejected.
Convert explicitly using existing centered_to_unshifted/unshifted_to_centered
helpers. Shape equality alone cannot prove ordering; callers must supply
truthful metadata. The high-level ordering selects both its reciprocal
intermediate and the expected propagator order. No ordering is inferred from
values. Generic odd meshes are supported; exact legacy boundary restrictions
remain in propagation_kernel, not the wave transforms.

The source reuses PMAP in place. Python operations preserve all inputs and
return independent arrays. WaveStepResult is frozen and stores defensive
copies backed by immutable bytes, so its arrays cannot be made writable.
Standalone operation outputs are new caller-owned arrays. Inputs must be
finite NumPy complex128 arrays; no dtype promotion, broadcasting, nonfinite
replacement, or overflow repair is performed.

## Norm, absorption and boundaries

The tested real-space norm is N(u)=sum(abs(u)**2). Parseval for this pair gives
N(u)=M*sum(abs(c)**2). Unit-modulus gratings and propagators preserve N through
a step (up to roundoff), without additional normalization. Reciprocal norm
alone is smaller by M; it must not be equated with N. UTEST sums reciprocal
intensities at the sampling point; its commented-out normalization is inactive.

PGMAP gives |T|=exp(SIGMA*Im(V)). Negative imaginary potential attenuates;
positive imaginary potential can amplify as written. The wave API neither
clips nor assumes that every complex grating is attenuating. It multiplies
the grating once, without inserting thickness, shifts or another potential.

The propagator may already contain exact negative-Nyquist zeros. Multiplying
by it removes those reciprocal components and can reduce norm even without
absorption. No new mask or aperture is applied by this module. With boundary
disabled and unit-modulus multipliers, odd meshes conserve norm too.

## Public API

```python
apply_phase_grating(real_wave, phase_grating)
real_to_reciprocal_wave(real_wave, *, ordering)
apply_propagator(reciprocal_wave, propagator, *, wave_ordering, propagator_ordering)
reciprocal_to_real_wave(reciprocal_wave, *, ordering)
propagate_one_step(real_wave, phase_grating, propagator, *, ordering, propagator_ordering)
# -> WaveStepResult(real_wave, reciprocal_wave, ordering)
```

The result's reciprocal_wave is **after** propagator multiplication and the
Python inverse normalization, immediately before the return transform. No
beam extraction is implemented. The result's real_wave is after that return
transform and before any next grating.

## Validation limits

Tests use fixed complex arrays, explicit roots of unity, independently nested
scalar DFTs, and existing grating/propagator constructors. They compare complex
intermediates and real outputs, not just intensities. Complete direct/FFT
comparisons cover (Ny,Nx)=(4,6),(3,5),(4,4), in both reciprocal orderings.
Absolute acceptance tolerances are 2e-14 for the real wave, 8e-15 for the
reciprocal intermediate, and 3e-14 for real-space norm conservation. These
allow accumulated scalar trigonometric/summation roundoff on these small
meshes; tests do not normalize the outputs to obtain agreement.

Measured with the required qcbed-python interpreter:

| Quantity | Maximum absolute residual | Maximum pointwise relative residual |
| --- | --- | --- |
| Real wave: FFT vs nested DFT | 2.443e-15 | 1.201e-14 |
| Sampling intermediate: FFT vs nested DFT | 5.794e-16 | 1.001e-14 |
| Real-space norm: output vs input | 1.777e-15 | 1.840e-16 |

Worst absolute real-wave residual: shape (4,6), unshifted, row=0,column=5.
Worst absolute reciprocal residual: shape (4,6), unshifted, row=3,column=3
(K=-1,H=-3). Worst real-wave relative residual is shape (4,4), unshifted,
row=1,column=0. Worst norm residual is shape (3,5), centered. Relative
residuals divide by nonzero reference magnitudes; absolute tolerances govern
near-zero cases. Norm is sum(abs(real_wave)**2), not its per-pixel average.

No saved intermediate
PMAP wave grid was identified in the readable WRITE paths or prior reference
audit. Saved final intensities require sequence/tilt orchestration and cannot
validate this isolated step. Historical executable equivalence is not claimed.
MSEQ iteration, multiple slices/regions, beam extraction, incident setup, tilt
loops, multiprocessing, GPU support and refinement remain deferred.
