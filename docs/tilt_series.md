# Independent tilt series

## Source audit and confidence

**Category B:** source-backed independent-tilt control flow, validated in the
established Python Fourier convention. Category A executable equivalence is
not established; this is not a category C unresolved active-definition case.
The absent readable historical 1D CDFT and original compiler/literal behavior
still prevent claiming agreement with the original executable. private predecessor implementation
is only a related-version declaration reference. No include-only capacity
or COMMON layout is needed by this module.

Primary evidence is `private predecessor implementation`, `private predecessor implementation:ONEIT` (lines 70-185),
`private predecessor implementation:ONEITPARA` (lines 60-173), and the existing PROPG equations.
The latter two routines use the same explicit serial tilt/slice loops.
ONEITPARA is a parameter/objective variant; its name does not establish
parallel execution. Optimizer calls rerun the objective, including applied
state and sequence reconstruction, rather than scheduling independent tilts.

The driver reads scan basis components AY,BEE,CEE,DEE, then corner pairs
A,B,C. Experimental input supplies NUMP,NUML,GRADIENTOPT. NORMFACTOR is
LENGTHY/LENGTHX for the supplied scan basis; INVNORMFACTOR is its inverse.
The objective first multiplies the second coordinate of A,B,C by NORMFACTOR,
then divides both coordinates by each stored type's APPLIEDDILATION.
Writing the resulting corners as Ad,Bd,Cd, the tilt-list construction is:

```
outer scan row K = 1..NUML
  inner scan column N = 1..NUMP
    M = (K-1)*NUMP + N
    TILTARRAY(type,M,:) = Ad
        + (K-1)/(NUML-1) * (Cd-Ad)
        + (N-1)/(NUMP-1) * (Bd-Ad)
```

Both endpoints are included when each dimension exceeds one. The source
divides by zero for singleton scan dimensions; Python does not invent a
grid rule for that case. Zero tilt is not guaranteed. This is row-major
scan-pixel order, not a guarantee that either lattice H or K stays constant
along a scan row. Both physical components may vary along either scan axis.

Immediately before PROPG, the second coordinate is multiplied by
INVNORMFACTOR and the two coordinates u,v are transformed to the cell basis:

```
tilt_h = AY*u + CEE*v
tilt_k = BEE*u + DEE*v
```

These are dimensionless fractional reciprocal-cell coefficients, not degrees,
radians, milliradians or inverse angstroms. Dilation division commutes with
the linear basis transformation. The Python input supplies the transformed
**undilated** coefficients h0,k0. For stored type i, PROPG receives
`h0/d_i, k0/d_i`, exactly once, as in `simulate_one_tilt`. The reciprocal
grid retains its own established dilation scaling. See [one tilt](single_tilt_simulation.md)
and [propagator equations](propagator.md) for that distinction.

The source builds PGARRAY once per stored type outside the tilt loop. For
each T, it rebuilds PRARRAY for every stored type, starts PMAP at the first
phase grating (unit illumination), and executes MSEQ sequentially. At the
last reciprocal sampling point it sets
`BEAM(I,T)=real(PMAP(IB(I)))**2+imag(PMAP(IB(I)))**2` for every requested beam.
Every beam uses the same T list. There is no coherent sum across T.

BEAM is indexed `(beam,tilt)`. THEORINTS flattens it with beam outermost and
T innermost. Experimental XPINTS, VARIANCE and WEIGHTS are read as
`NUML*MBOUT` records of NUMP pixels: beam blocks, scan rows, then scan columns.
The Python result instead explicitly exposes `(tilt,beam)`; it is not silently
reshaped or flattened into experimental file order. The driver maps signed
beam H/K to one-based unshifted PMAP IB; the existing validated Python mapping
is reused, including request order and duplicates.

The objective changes shared applied parameters, MSEQ and PGARRAY before the
tilt loop. Within it PMAP is overwritten, PRARRAY is rebuilt, BEAM is filled,
FFTCOUNT/IFFTCOUNT advance, and UTEST updates TOTALINT. The local Python
execution owns its working arrays and never mutates prepared state. No COMMON
storage, random generator, counters or global cache is introduced.

## Public API and execution

```python
serial = simulate_tilt_series_serial(prepared_sample, tilts)
batched = simulate_tilt_series_batched(prepared_sample, tilts, batch_size=8)
```

`tilts` must be a finite float64 NumPy array of shape `(n_tilts,2)`, with at
least one row. Column 0 is h0, column 1 is k0 as defined above. Lists, bool,
integer/float32 arrays and malformed shapes are rejected rather than coerced.
`batch_size` must be a positive integer excluding bool. One, larger than the
tilt count, and final incomplete chunks are supported.

`TiltSeriesResult` holds immutable bytes-backed snapshots of `tilts`,
complex128 `amplitudes` and float64 `intensities`, together with requested
`beam_indices_hk`, `execution_mode`, and the requested maximum `batch_size`
(None for serial). Amplitudes and intensities have shape `(n_tilts,n_beams)`.
Duplicate tilts and duplicate beams retain their positions. Intensities are
`real(amplitude)**2+imag(amplitude)**2` independently per experiment.

The serial correctness reference calls `simulate_one_tilt` once per row.
The batched path keeps separate waves of shape `(B,mesh_y,mesh_x)` and uses
dedicated backend transforms with axes **(-2,-1)**. Reciprocal coefficients
are the positive-exponent DFT divided by mesh_x*mesh_y; real-space sums use
the unnormalized negative-exponent transform. There is no batch-size scaling,
extra RMT factor, intensity fitting or normalization to unit total intensity.
Only nonfinal occurrences transform back to real space.

Both centered and unshifted prepared grids are supported. The batched engine
uses unshifted bins internally: one cached permutation converts each native
kernel, the shared mask and beam positions. Real-space gratings are unchanged.
No batch-axis shift or transform occurs. Kernels still use the existing
propagation routine, including slice shifts, with its mask disabled so the
shared mask can be applied without rebuilding it for each tilt.

Preparation caches reciprocal grids, phase gratings, wavelength, stored-type
geometry, applied parameters, slice sequence and native beam positions. Each
series call computes its boundary mask and internal ordering/beam mapping
once. Each chunk builds one propagator per stored type per tilt, then executes
the sequential occurrences using broadcast gratings. A function scope releases
chunk waves and kernels before starting the next chunk. With S stored types
and M mesh points, principal working storage is O(B*M*(S+constant)); returned
tilts and beam results necessarily use O(n_tilts*n_beams+n_tilts). No full
all-tilt wave history, all-tilt propagator cache or automatic memory sizing is
created. Propagator construction remains a small Python loop, deliberately
reusing the audited physical calculation.

This represents independent plane-wave experiments. It does not represent a
spherical/convergent coherent probe. No disk-grid helper is exposed: the
legacy scan and experimental sections remain unparsed, and callers must
provide explicit cell-basis tilts and retain any external grid metadata.

## Validation and benchmark

`tests/test_tilt_series.py` covers ordered signed and duplicate tilts/beams,
A-B-A-B with two stored types, direct/FFT gratings, absorption on/off, shifts,
nonunit dilation, square/rectangular meshes, centered/unshifted storage, chunk
sizes 1/2/greater than count, immutable results, invalid inputs and empty beam
requests. A change to the middle tilt must leave both neighboring complex
amplitudes and intensities exactly unchanged. Instrumentation checks one mask
per call, propagator counts, bounded chunks, release of previous chunk arrays
and absence of the final reciprocal-to-real transform.

An independent 2x4 calculation uses printed Al scattering values
281.953595001 at STL=0 and 275.340445007875 at STL=.025 inverse angstrom,
independently checked wavelength/sigma, scalar atomic sums, Cartesian
reciprocal vectors and nested DFTs. No production scattering, grating,
propagator or FFT routine generates its expected values. Amplitude tolerance
is 3e-15 absolute; intensity tolerance is 4e-15. The uniform elastic DC-only
case independently checks intensity one and rules out another RMT factor.
Intermediate transmitted real-space and post-propagator reciprocal complex
waves are also compared at every used transform boundary against nested DFTs.

The real Al FOR007 input uses five explicit tilts, nine recorded beams,
16x16 mesh and 415 occurrences at its explicit initial applied parameters.
Serial/batched comparisons use 3e-14 absolute tolerance. This exercises real
input, not executable-output equivalence: fitted YNORM outputs cannot serve
as raw wave/intensity references, and historical CDFT is still unresolved.

Run `benchmarks/benchmark_tilt_series.py` with the required qcbed-python
interpreter. It reports preparation and median-of-three execution times after
warmup for 32 tilts, 16x16 mesh, two types, 20 A-B occurrences and five beams,
with absorption, shifts and unequal dilation. Comparisons include repeated
one-tilt calls, the serial wrapper, and batch sizes 1,2,4,8,16,64. There are no
pytest timing assertions or unmeasured scaling claims.

Measured on 2026-09-16 with
`C:\Users\peng863\AppData\Local\miniconda3\envs\qcbed-python\python.exe`
and NumPy 2.4.6: preparation took 0.016501 s. These are local informational
measurements, not portable performance guarantees.

| Execution | Total seconds | Milliseconds per tilt |
| --- | ---: | ---: |
| Repeated one-tilt calls | 0.094489 | 2.952791 |
| Serial wrapper | 0.097158 | 3.036184 |
| Batch 1 | 0.039083 | 1.221350 |
| Batch 2 | 0.028589 | 0.893394 |
| Batch 4 | 0.019261 | 0.601909 |
| Batch 8 | 0.013086 | 0.408941 |
| Batch 16 | 0.010854 | 0.339187 |
| Batch 64 (one 32-tilt chunk) | 0.010235 | 0.319837 |

Measured maximum absolute serial/batched amplitude and intensity residuals
were both **0** across the test matrix, Al subset, and benchmark. Against the
independent tiny-model DFT, maxima were **4.478e-16** for complex amplitudes
and **4.441e-16** for intensities. Neighboring experiments in the no-cross-talk
test were exactly unchanged. These statements describe the measured cases;
tests retain nonzero tolerances for numerical portability.

Validation with that interpreter: focused tilt-series tests **87 passed in
1.47 s**; complete pytest suite **1135 passed in 17.37 s**. The complete run
includes existing reference-data checks; no new executable-Fortran comparison
was possible. `git diff --check` passed.

Deferred: scan-grid construction, multiprocessing/threads, coherent probes,
GPU, experimental-intensity parsing, YNORM, chi-square, refinement and CLI.
