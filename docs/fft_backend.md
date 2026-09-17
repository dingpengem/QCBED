# NumPy CPU Fourier backend

## Specification and algebraic mapping

The numerical specification is the validated Python direct sum:

```text
V[r,c] = sum_HK C[K,H] exp(-2*pi*i*(H*c/Nx + K*r/Ny))
```

There is no normalization or output shift; x=c/Nx, y=r/Ny. Shapes are
(Ny,Nx)=(mesh_y,mesh_x), with row/K/y on axis 0 and column/H/x on axis 1.

The installed NumPy 2.4.6 `numpy.fft` documentation defines the DFT with a
negative exponent and its inverse with a positive exponent and 1/N scale.
Applying the two 1D definitions on axes 0 and 1 yields:

```text
fft2(U)[r,c] = sum_q=0..Ny-1 sum_p=0..Nx-1
               U[q,p] exp(-2*pi*i*(p*c/Nx + q*r/Ny))
```

Set p=H mod Nx, q=K mod Ny, and U[q,p]=C[K,H]. For integer pixel indices,
replacing H by H+m*Nx (or K by K+n*Ny) multiplies the exponential by
exp(-2*pi*i*integer)=1. Thus every summand matches the direct specification.
This determines the transform algebraically, without choosing a routine by
trial against tests:

```python
V = np.fft.fft2(U, axes=(0, 1), norm="backward")
U = np.fft.ifft2(V, axes=(0, 1), norm="backward")
```

NumPy's explicit `norm="backward"` means fft2 is unscaled and ifft2 divides
by Nx*Ny. The inverse follows from the orthogonality sum over spatial pixels:
the product of axis sums is Nx*Ny for matching bins, zero for distinct bins.
The inverse therefore recovers the coefficients exactly in mathematics,
within rounding in floating-point arithmetic. This inverse is defined for
the Python operation; it is **not** yet a validated legacy multislice wave
transform. No physical area, wavelength, thickness or propagation factor
belongs in either backend transform.

## Ordering, odd sizes and Nyquist

The existing `grids.reciprocal_index_grid` defines the signed frequencies.
The backend obtains centered H/K from that helper and scatters each value
to its modulo bin; the reverse function gathers from those same bins.
There is no independent centered-frequency definition.

- Size 4: centered [-2,-1,0,1], unshifted [0,1,-2,-1].
- Size 5: centered [-2,-1,0,1,2], unshifted [0,1,2,-2,-1].
- Each axis uses its own size; rectangular and mixed-parity meshes work
  without transpose. Length-one axes are also supported.

Tests establish that this explicit permutation equals NumPy ifftshift on
centered input and fftshift in reverse, including odd dimensions where the
two permutations differ by one position. Even meshes also agree with the
existing OSHIFT helper. The backend does not call OSHIFT on odd grids or
extend its legacy half-swap contract.

For even N the project labels the Nyquist bin as -N/2. Modulo mapping sends
it to N/2; the coefficient is preserved, including its imaginary part. No
Hermitian symmetry, real-valuedness or zero boundary is assumed. Reordering
and Fourier transformation apply no physical masks. The separate source
Nyquist mask in phase_grating.py suppresses the designated reciprocal row
and column only when invoked through phase-grating composition.

## Public API and integration

All functions below are available from qcbed and in their respective modules:

```python
centered_to_unshifted(coefficients) -> complex128 array
unshifted_to_centered(coefficients) -> complex128 array
reciprocal_to_real_sum(coefficients, *, ordering="unshifted") -> complex128 array
real_to_reciprocal_coefficients(real_space, *, ordering="unshifted") -> complex128 array
projected_potential_fft(coefficients, reciprocal_grid) -> complex128 array
```

The two reorder functions do not transform or normalize. In
reciprocal_to_real_sum, ordering describes the reciprocal **input**. In
real_to_reciprocal_coefficients it describes the reciprocal **output**.
Real-space arrays always start at fractional (0,0); no spatial centering is
accepted. Ordering must be exactly "centered" or "unshifted".

Inputs must be finite nonempty 2D NumPy complex128 arrays. Float, complex64,
lists, wrong dimensionality, nonfinite values and invalid ordering are
rejected instead of silently converted. Outputs are new C-contiguous
complex128 arrays. Inputs may be read-only or noncontiguous and are not
mutated. Nonfinite FFT results raise ValueError; no NaN replacement occurs.

projected_potential_fft reads the ordering from the existing reciprocal grid
and requires matching shape. Like projected_potential_direct it applies no
masks and accepts even, odd and mixed-parity meshes. The existing high-level
function now accepts:

```python
phase_grating_from_structure_factors(
    elastic, absorptive, reciprocal_grid, *, voltage_kv,
    dsap_inverse_angstrom, method="direct"
)
```

The default remains the direct oracle. "fft" explicitly selects the backend;
unknown methods fail. Both paths share the exact same combination of elastic
and absorptive coefficients, inclusive STL<=DSAP mask, negative-Nyquist mask,
interaction constant and phase-grating exponentiation. This high-level path
still rejects odd dimensions because its source-backed boundary mask requires
even dimensions; that restriction does not apply to the mathematical backend.

## Validation and timing

Tests cover DC normalization, signed H and K sinusoids against explicit roots
of unity, Friedel cosines, modulo-bin placement, even OSHIFT agreement, odd
permutations, Nyquist preservation versus masking, complete gratings, inverse
normalization, input immutability and invalid inputs. Deterministic complex
data use seed 20260916 on meshes 4x4, 5x5, 6x4, 5x4, 4x5, 5x3, 1x1 and
1x4, in both orderings. A separate nested scalar sum and the existing direct
array sum are both compared. No timing assertion belongs to pytest.

For the small random meshes, transform comparisons use atol=8e-14 and
rtol=2e-14, allowing accumulation of direct trigonometric rounding. Round
trips use atol=rtol=2e-15. Complete phase gratings use atol=8e-15,
rtol=2e-14. Near-zero reference entries are controlled by absolute tolerance.

Measured with the required qcbed-python interpreter and NumPy 2.4.6:

| Comparison | Maximum absolute residual | Maximum pointwise relative residual |
| --- | --- | --- |
| FFT vs direct, small random meshes | 1.313e-14 | 8.035e-15 |
| Reciprocal -> real -> reciprocal | 8.951e-16 | not used for acceptance alone |
| Complete direct vs FFT gratings | 1.111e-16 | comparable magnitude (gratings near unit modulus) |

The worst FFT/direct case is centered 6x4; the worst round trip is unshifted
5x3. Pointwise relative residual divides by the magnitude of the nonzero
reference; these finite random cases contain no exact-zero denominator.

`benchmarks/benchmark_fft.py` is a standalone informational timing helper.
Run it with the same interpreter as pytest. It times the public reconstruction
functions, including validation and centered-bin reordering, on centered
32x32 complex128 random input. One warmup is followed by five calls per
function, reporting the median. Measured in this environment:

- Direct: 0.022336200 seconds (22.3362 ms).
- NumPy FFT: 0.000076900 seconds (0.0769 ms).
- Maximum absolute discrepancy: 4.975e-13; relative L2 discrepancy: 3.195e-15.

These timings are informational, include Python overhead and depend on the
machine/load. They are not a performance guarantee or evidence of historical
Fortran equivalence. No multiprocessing, thread management, Numba, external
FFT dependency or GPU support was added.

## Historical Fortran boundary

The readable `private predecessor implementation` OSHIFT, `private predecessor implementation` FFT wrapper,
`private predecessor implementation` CDFT2D and `private predecessor implementation` CDFT2D_SUB were inspected.
The wrappers preserve H-fastest packing, pass the same sign through both
axes and contain no normalization. Their mathematical specification agrees
with the direct sum, and the restored 1D CDFT now corroborates its sign. Original linkage remains
unverified. private predecessor implementation is only a related-version declaration
reference. This work validates the NumPy implementation against the **Python
mathematical specification**, not the historical CDFT executable. The source
audit's category B limitation remains unchanged.

No saved complex phase-grating pixels were identified in the preceding audit;
there is no new historical FFT comparison. Existing saved structure-factor
and rounded SIGMA checks remain in the complete suite. Wave propagation,
PROPG, multislice, probes, detectors, tilt scanning and refinement are deferred.

## Source-backed wave transforms (2026-09-16)

Dedicated wave helpers use negative exponent / N and positive exponent with
no normalization, for both single grids and batches. Generic potential sums
and their inverse retain their original contracts. See
[FFT sign correction](fft_sign_correction.md) for source evidence and validation.
