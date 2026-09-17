# One-tilt multislice sequence and beam extraction

## Source audit and decision

**Category B.** Sequence selection and beam sampling are established by the
readable ONEIT (`private predecessor implementation`) and ONEITPARA (`private predecessor implementation`). The
Python Fourier convention remains as documented in wave_propagation.md;
historical CDFT equivalence is not claimed. private predecessor implementation is a related-version
declaration reference, not the verified original include.

MSEQ is not read from FOR007. Each objective function declares a local
MSEQ(MSL), initializes it to zero, and constructs it in three UPDPARS/sequence
update passes before the tilt loop (`private predecessor implementation`, `private predecessor implementation`):

```text
K = 0; SLN = 0
for region J = 1..NREGIONS:
    THICK = 0; I = 0
    while THICK < DEPTH(J):
        K += 1; I += 1
        L = MOD(I, RNSL(J)); if L == 0: L = RNSL(J)
        L += SLN
        MSEQ(K) = L
        THICK += DELTAZ(L) * APPLIEDDILATION(L)
    SLN += RNSL(J)
    RSLICES(J) = I; SLICES(J) = DBLE(I)
MS = K
```

MSEQ entries are **one-based stored slice-type IDs**: the same JPG selects
PGARRAY(JPG,:) and PRARRAY(JPG,:). They are not region IDs or occurrence
numbers. NSLICE counts stored slice types, MS counts executed occurrences,
RNSL counts stored types assigned to each region, and RSLICES counts executed
occurrences within a region. No readable NSLICES variable was found. J is the
one-based occurrence cursor in the propagation loop; other routines reuse J
as a local loop variable with no shared sequence meaning.

Region types occupy consecutive blocks in the stored arrays. A region with
two stored types cycles A,B,A,B,...; the next region starts at its own first
type, offset by SLN. For example regions with stored types (A,B) and (C), and
occurrence counts four and two, produce legacy (1,2,1,2,3,3). The wave is not
reset at region boundaries. This engine accepts that explicit sequence; it
does not derive counts from depth, vacancy, dilation or parsed inputs.

The driver reads NREGIONS, DEPTH/TMASK and RNSL (`private predecessor implementation`) and
performs diagnostic/pre-dilation counting (`174-188`, `313-329`) but does not
construct MSEQ there. The second preliminary loop uses THICK<=DEPTH, unlike
the effective objective loop's strict <. It must not be substituted for the
objective's sequence construction. Refinement callers pass ONEIT/ONEITPARA
as objective functions; UPDPARS changes depth/applied parameters before each
rebuild. MSEQ is fixed across tilts within that evaluation. The parser
currently validates then discards region counts/depths and beam requests;
this task adds no parser fields or configuration-to-array construction.

PROPG is called for each stored type K at each tilt with DELTAZ(K),
XSHIFT(K),YSHIFT(K) (`private predecessor implementation`, `private predecessor implementation`). Its own
geometry, applied dilation, internal tilt and fractional shifts determine
PRARRAY(K,:). No preceding/following slice is passed. The phase grating and
propagator therefore both belong to the **current** selected slice.

Only MSEQ(1:MS) is active. MSL is storage capacity, not MS+1 semantics, and
must not be copied into Python. No intentional sentinel at MS+1 is assigned.
The source still increments J and reads MSEQ(MS+1) after final sampling,
normally obtaining the zero-initialized inactive entry and using PGARRAY(0,:);
at capacity the sequence lookup itself can be out of bounds. The reason for
this tail is the unconditional loop body, not an extra physical slice. Its
resulting PMAP is not used by reported beams. Runtime failure from an invalid
access is possible; equivalence of that invalid execution is not asserted.

## Effective loop and incident wave

Both sources initialize PMAP=PGARRAY(MSEQ(1),:). This is exactly a unit
real-space wave times the first transmission grating. No focused probe,
nonuniform incident array, multi-beam illumination or explicit incident
reciprocal delta is constructed in either readable path. A unit DC delta
is mathematically equivalent under the Python transform pair but is not a
separate source initialization branch. Tilt enters the precomputed propagators.

The effective Python loop is:

```text
real_wave = incident_wave                 # unit real-space plane wave by default
for occurrence, slice_id in sequence_zero_based:
    transmitted = real_wave * phase_gratings[slice_id]
    reciprocal = real_to_reciprocal_wave(transmitted)
    reciprocal *= propagators[slice_id]
    if occurrence is final:
        extract requested complex beams and intensities
        return final reciprocal wave
    real_wave = reciprocal_to_real_wave(reciprocal)
```

There is no final return transform, MS+1 lookup, dummy sequence entry or next
grating. To avoid an unnecessary all-ones allocation/multiplication, the
default path passes the first selected grating directly to the nonmutating
transform. An optional custom incident wave is a documented Python extension:
it must be finite complex128, have the same shape, and is not mutated.

## Beam mapping and intensities

`private predecessor implementation` reads MBOUT and all (HIN,KIN) pairs in one list-directed
READ. Lines 342-355 copy to HOUT/KOUT, add MESHX/MESHY to negative indices,
add one to each coordinate, and calculate
IB=(KOUT-1)*MESHX+HOUT. IB is a **one-based** PMAP index. Source sampling is
unshifted after PROPG's OSHIFT; wave transforms apply no further centering.

For Nx columns and Ny rows:

| Quantity | Mapping |
| --- | --- |
| Dimensionless beam | (H,K), distinct from physical reciprocal coordinates |
| Centered row,column | (K+Ny//2,H+Nx//2) |
| Unshifted row,column | (K mod Ny,H mod Nx) |
| Zero-based C-order linear offset | row*Nx+column |
| Legacy IB | unshifted linear offset + 1 |

The implementation reuses validate_mesh_dimensions and row_col_to_pmap_index;
it does not duplicate PMAP flattening arithmetic. Axes are independent on
rectangular meshes; no transpose is used. Accepted signed ranges follow the
existing grid convention: -N//2 through (N-1)//2 (where the lower bound means
-(N//2)). Even positive Nyquist +N/2 is rejected, with -N/2 the canonical bin.
Out-of-grid values are rejected before modulo mapping, including otherwise
aliased values. The source performs only one negative wrap and has no bounds
check; rejection is explicit Python validation, not a recovered Fortran check.

Both loops sample after PRARRAY and RMT, before returning to real space
(`private predecessor implementation`, `private predecessor implementation`). Each amplitude is complex128;
intensity is float64 real(amplitude)**2+imag(amplitude)**2. Input beam order
and duplicates are preserved; no special scaling is applied to the direct
beam. An empty beam list is allowed (MBOUT=0), retaining the final wave.

UTEST (`private predecessor implementation`) only updates TOTALINT=sum |PMAP|^2. Its WRITE is
commented out, J=J does nothing, and it neither normalizes nor changes PMAP.
The surrounding PMAP/TOTALINT normalization is also commented out.

MT=Nx*Ny and RMT=1/MT (`private predecessor implementation`). The Python real-to-reciprocal
fft2(norm="forward") divides by MT; reciprocal-to-real ifft2(norm="forward") is unnormalized. No
extra RMT or fitted intensity scale is applied. For lossless multipliers,
sum |final reciprocal|^2 = sum |incident real wave|^2 / MT. Thus unit incident
illumination has total reciprocal intensity one. Absorptive gratings and
propagator boundary zeros can each reduce this norm; neither is renormalized.
Duplicate selected beams should not be summed as a total-intensity diagnostic.

## API, validation and memory

```python
legacy_slice_sequence_to_zero_based(sequence_one_based, slice_type_count)
unit_incident_wave(mesh_x, mesh_y)
beam_indices_to_row_col(beam_indices_hk, *, mesh_x, mesh_y, ordering)
beam_indices_to_legacy_ib(beam_indices_hk, *, mesh_x, mesh_y)
extract_beam_amplitudes(reciprocal_wave, beam_indices_hk, *, ordering)
beam_intensities(amplitudes)
simulate_slice_sequence(
    phase_gratings, propagators, slice_sequence_zero_based, *,
    beam_indices_hk, reciprocal_ordering, propagator_ordering,
    incident_wave=None,
) -> MultisliceResult
```

The conversion helper validates one-based entries before subtracting one.
The engine accepts only nonempty zero-based sequences within collection
bounds, never infers the base. Integers include NumPy integers but exclude
bool, floats, strings and nonfinite values. Collections must be nonempty,
equal length and consist of finite complex128 arrays with identical nonempty
(mesh_y,mesh_x) shapes. Every stored array is validated, including unused
types. Broadcasting and automatic dtype conversion are not allowed.
reciprocal_ordering and propagator_ordering are required and must agree;
use existing FFT reorder helpers explicitly for a different stored order.
Real-space gratings are never shifted. Odd meshes are allowed when supplied
arrays permit them; boundary policy remains in the propagator constructor.

The frozen result contains final_reciprocal_wave, amplitudes, intensities,
beam_indices_hk, slice_sequence_zero_based, slice_type_count and
reciprocal_ordering. Array fields are defensive immutable snapshots, with
complex128 waves/amplitudes and float64 intensities. No per-step wave history
or unused final real wave is retained. Stored arrays are referenced without
stacking/copying the entire collection; working memory is O(mesh pixels),
plus final snapshots, O(sequence length) integer metadata and O(beam count)
extraction data. No legacy fixed capacity limits are introduced.

## Numerical validation

Independent two-point arithmetic checks two distinct slices, A-B-A-B reuse,
and two later occurrences of an absorptive third type. These hard-coded
complex results agree exactly in the tested binary-representable cases.
Separate nested scalar DFTs check four-step sequences on (Ny,Nx)=(2,3),
(3,5),(4,6), in both reciprocal orderings, with nonuniform incident waves
and spatially varying unit-modulus gratings and propagators. Expectations
never call the production sequence or FFT backend.

Measured with the required qcbed-python interpreter:

| Comparison | Maximum absolute residual | Maximum relative residual |
| --- | --- | --- |
| Final complex reciprocal wave vs nested DFT | 1.078e-15 | 1.609e-14 |
| Total reciprocal intensity vs incident real norm / MT | 1.111e-15 | 9.800e-16 |

The worst absolute complex residual occurs on shape (4,6), centered engine
output, at unshifted row=3,column=0 (K=-1,H=0; centered row=1,column=3).
The worst pointwise relative residual is the same mesh/order at unshifted
row=1,column=4 (K=1,H=-2). Relative residual divides by the nonzero reference
amplitude magnitude. Norm residuals are worst on the same mesh/order.
Acceptance uses absolute tolerance 8e-15 for the final wave and 4e-15 for
the normalized norm, allowing accumulated scalar DFT rounding. It does not
rescale the waves to force conservation. Tests also cover nonzero source
PROPG phase, source PGMAP attenuation, exact boundary zeros, explicit/default
illumination, invalid data, immutable snapshots and absence of the final
return transform. No benchmark or historical executable comparison is implied.

## Reference boundary

The current FOR007 requests nine beams on a 16x16 mesh: (0,0), (1,1),
(1,-1), (-1,-1), (-1,1), (2,0), (0,-2), (-2,0), (0,2). Tests verify these
requests against independent one-based IB values. This is a reference-input
mapping check, not a numerical propagation comparison.

Saved out.dat and bestout.dat are empty. firstout.dat contains text output,
but the driver writes YNORM*THEORINTS (`private predecessor implementation`), and objective
output likewise writes that product (`private predecessor implementation`, `private predecessor implementation`).
YNORM is selected/fitted in CHISQCAL, separately from raw beam extraction.
No saved precomputed PGARRAY/PRARRAY or final complex PMAP fixture tied to
that output was identified. Reproducing it requires configuration, refinement
state and tilt construction excluded from this task. Consequently no honest
per-beam Python/reference intensity error can be reported for this isolated
engine. List-directed text precision is compiler-dependent (the observed
firstout records have about 15 significant digits), not a full-precision
reference. No intensity values or missing precomputed inputs are guessed.

Historical CDFT signs/normalization, linkage, original literal kinds and
executable equivalence remain unresolved. Sequence construction from region
depths, parsed-input numerical setup, tilt grids, experimental data parsing,
fitted normalization/refinement, parallelism, GPU support and CLI remain
deferred.
