# Prepared sample and one incident tilt

## Decision and source setup

**Category B:** the supported numerical path is source-backed, but the parser
does not retain all preparation state. Applied values are explicit, not
guessed from masks or reconstructed by a refinement engine. Historical CDFT
equivalence remains unresolved; the established Python transform pair is
unchanged. private predecessor implementation is a related-version declaration reference only.

The readable setup path is:

```text
PROGRAM MULTIS: read KV/zone, one global mesh, per-slice geometry, species
    -> CELL2R from each CELL2D
    -> DSAP, regions/depths/RNSL, vacancy/dilation/thermal/absorption input
    -> SFMAP reads atoms and initializes applied state for each stored type
    -> read beams, scan basis/corners, reflection/optimizer parameters
    -> ONEIT or ONEITPARA
        -> three UPDPARS + MSEQ rebuild passes
        -> SFMAP2, OSHIFT, FFT(-1), PGMAP for each stored type
        -> for each tilt: PROPG for each stored type
        -> execute MSEQ, sample final reciprocal beams
```

Primary source locations: `private predecessor implementation` (voltage, mesh, geometry,
species, DSAP), `155-329` (regional input and preliminary counting),
`332`/`1739-1919` (SFMAP), `341-355` (beams), `398-462` (scan and reflection
input), `private predecessor implementation` and `private predecessor implementation` (objective setup/loops),
`private predecessor implementation` (SFMAP2/PGMAP), and `682-723` (PROPG).

Each `LegacyInput.slices[i]` and `atoms_by_slice[i]` defines one stored type.
Species are global; atom species indices remain legacy one-based and are
validated through the existing conversion helper. SFMAP randomizes G sites;
these unresolved records remain unsupported. Fixed S coordinates and starting
occupancies are used with the explicit applied vacancy correction once.

MSEQ contains one-based stored type IDs, selecting both PGARRAY and PRARRAY.
The effective sequence repeats each region's contiguous RNSL block until
its depth is reached using DELTAZ*dilation, then advances to the next block.
The source assigns disjoint stored-type blocks to regions, not shared IDs
between regions. Identical material may therefore appear in distinct stored
types. Repeated occurrences within a region reuse one type. This API accepts
an explicit zero-based sequence without reconstructing regions or assuming
that arbitrary ID reuse describes the original region input.

## Parameter ownership and compatibility

| Scope | Source state / Python boundary |
| --- | --- |
| Global | One MESHX/MESHY, applied voltage and wavelength, DSAP, species list, raw ABSORB, INELA/B/C per species, requested beams. |
| Stored slice type | Direct/reciprocal cell, thickness, fractional shifts, atoms, applied dilation, vacancy and APPLIEDDW per species. |
| Region | Input DW per species, depth/repetition, vacancy/surface models and reflection refinement selection; supply their applied per-slice results here. |
| Tilt | Propagation phase/kernel and all propagated waves/intensities. |

No per-slice mesh is read by the source. This integration requires one shared
even mesh. The source stores a separate lateral cell for every type, averages
some quantities for scan setup, and contains no explicit lateral-equality
check or wave-resampling operation. To avoid inventing an interface between
different undilated periodic cells, the **supported high-level subset requires
exactly matching undilated a,b,gamma** across types. This is an explicit Python
restriction, not a recovered Fortran check. Per-slice applied dilation remains
allowed: the source explicitly uses it in area, reciprocal metric, distance
and tilt construction while propagating on the same fractional mesh. No
physical Cartesian wave resampling is added for those dilation changes.

Input geometry, mesh and atom association are taken from LegacyInput.
Preparation accepts one immutable SamplePreparationParameters with:

- `applied_voltage_kv`: explicit current KV (use model.energy_kv for an initial
  snapshot). Wavelength and interaction constants come from existing helpers.
- `slices`: one SliceAppliedParameters for every stored type, containing
  `debye_waller_angstrom_squared` in global species order, `vacancy_fraction`
  and `dilation`. These are APPLIEDDW/VAC/APPLIEDDILATION, not bulk DW,
  vacancy-profile coefficients or the raw input DILATION coefficient.
- `dsap_inverse_angstrom`: global inclusive STL cutoff, not |g| cutoff or DSMAX.
- `absorption_coefficients`: global per-species applied INELA/B/C values, or
  explicit None when disabled. Enabling follows model.absorption_raw[:2]=='YE'
  exactly, without trimming, uppercasing or accepting quoted YES.
- `slice_sequence_zero_based` and `beam_indices_hk`: explicit validated
  occurrence IDs and signed dimensionless beam indices. The named legacy
  sequence converter can convert MSEQ before preparation.
- `reciprocal_ordering` and `phase_grating_method`: explicitly centered or
  unshifted, and direct or fft. Gratings remain unshifted real-space samples.
- `boundary_policy='legacy'`: this preparation path uses the existing even-mesh
  phase-grating and PROPG Nyquist masks; unsupported policies are rejected.
- `reflection_multipliers='unity'`: explicitly selects the implemented raw
  atomic sums. Nonunit SFFAC/ABSFAC and their symmetry/refinement rules are
  unsupported, not silently assumed absent because the parser stopped early.

All fields are required. Disabled absorption skips the absorptive sum and uses
an exact zero coefficient grid; no made-up absorption coefficients or positive
B are inserted. B>=0 is valid for elastic-only preparation; absorption requires
B>0 because the existing source FAATM formula is singular at zero. Signed
finite absorption weights retain the lower-level contract. Vacancy is in [0,1],
occupancy is nonnegative, dilation and thickness are positive, shifts and atom
coordinates are finite fractional values without wrapping.

Preflight checks every stored type, all species/atoms, mesh/order, common cell,
sequence bounds, beam representability, parameter counts, physical domains and
derived wavelength before any structure-factor sum. Errors identify a slice
type using a zero-based index where applicable. The current raw scattering
helpers evaluate the full grid before masks; accordingly the full grid must
remain in their atomic scattering domains, even where DSAP would exclude
points. No hidden extrapolation or fake STL substitution is performed.

## Preparation and immutable cache boundary

```python
prepared = prepare_multislice_sample(model, parameters)
```

For each stored type exactly once:

```text
SliceGeometry -> DirectLattice2D -> ReciprocalLattice2D/ReciprocalVectorGrid
    -> elastic_structure_factors
    -> optional absorptive_structure_factors
    -> phase_grating_from_structure_factors(method=explicit method)
       [existing combination, DSAP/Nyquist masks, direct/FFT sum, PGMAP]
```

No equation from those modules is reimplemented. Thickness and fractional
slice shifts are not inserted into the structure factors or phase grating.
The result holds the applied configuration, mesh, raw title/absorption metadata,
wavelength, beam positions and one PreparedSliceType per stored ID. Each type
holds immutable geometry, reciprocal grid and a defensive immutable grating
snapshot. Derived grids use existing immutable byte storage. No original mutable
atom/species list or numerical work array is retained. Preparing a changed
structural/applied state produces a new sample; there is no global cache.

Grids, raw sums, masks, beam mappings and gratings are tilt-independent at fixed
applied state. Raw sums and mask temporaries are discarded after producing
the gratings; they need not be cached as additional arrays. Beam positions are
available as metadata; the existing sequence executor independently validates
beam mapping when extracting its result. Memory scales with stored types and
mesh size, not repeated MSEQ occurrences or a retained wave history.

## One-tilt convention and execution

```python
kernels = build_tilt_propagators(
    prepared, tilt_h_undilated_fractional=..., tilt_k_undilated_fractional=...,
)
result = simulate_one_tilt(
    prepared, tilt_h_undilated_fractional=..., tilt_k_undilated_fractional=...,
)
```

The pair is dimensionless in the common undilated reciprocal-cell basis,
after any user scan-basis transformation but **before** division by each
slice's applied dilation. These are not radians, degrees, mrad, tangents or
physical inverse-angstrom components. The names distinguish this boundary
from PropagationParameters.tilt_h/k_fractional, which arrive at PROPG after
the caller's dilation division.

In `private predecessor implementation` and `private predecessor implementation`, scan corners and hence linear
interpolated TILTARRAY are divided by each slice's d. INVNORM and the AY/BEE/
CEE/DEE transform then act linearly (`private predecessor implementation`). Therefore a common
already-transformed undilated pair (h0,k0) reaches PROPG as (h0/d,k0/d).
This setup division is applied once. It does not replace PROPG's separate
reciprocal metric /d^2 and propagation distance DZ*d. At d=1, the supplied
pair equals the lower-level PROPG pair. No scan-grid or angle conversion is
implemented. The same fractional H/K array labels are preserved for oblique
cells; existing reciprocal geometry handles the cross term.

Each call constructs one kernel per stored type with its own thickness,
shifts and dilation, using the prepared wavelength and order. The public
builder returns immutable kernels for inspection. simulate_one_tilt builds
them and calls simulate_slice_sequence with unit real-space illumination;
this high-level API does not expose the lower-level custom-probe extension.
The existing MultisliceResult returns the final complex reciprocal wave and
requested amplitudes/intensities, preserving duplicates and input order.
No final return transform, MS+1 read, extra RMT, YNORM or normalization is added.

The sequence is fixed across tilts of a prepared sample. Structural/applied
changes (including voltage, atom coordinates, vacancy, thermal factors,
dilation, cutoff or reflection policy) require preparation again. In readable
refinement, mesh/species/region assignment/beam requests and cutoff are fixed
by input, whereas depth, applied material state, scan parameters and selected
reflection factors may change after UPDPARS. The original energy-refinement
path may leave VOLTS stale while updating KV/WAVEL; this stateless API uses
one explicit voltage consistently and does not emulate that historical issue.

## Parser and validation limits

The parser consumes but does not retain DSAP, region depths/RNSL, vacancy/
dilation coefficients, absorption weights, regional DW/surface parameters,
or output beam pairs. Later scan basis/corners, reflection SFFAC/ABSFAC/SYMOPT,
optimizer/refinement sections and experimental data are not parsed. Their
effects cannot be inferred from LegacyInput alone. This module executes an
explicit supported snapshot; it is not a complete FOR007/refinement runner.

The reference-input integration uses the actual Al input and explicitly
corroborated initial state: 200 kV, B=.86, vacancy=0, applied d=1, DSAP=2.5,
absorption weights (0,1,1), unity reflection factors. Depth=840 and
DELTAZ=2.0245 give 415 occurrences under the effective strict-< depth loop,
all stored ID 0 in Python. The nine recorded beam requests are preserved.
This initial zero-tilt calculation checks finite gratings, propagators and
nonnegative intensities; it does not reproduce a fitted historical scan.

Saved firstout.dat is scaled by YNORM and does not identify matching prepared
arrays and a zero-tilt wave; out.dat/bestout.dat are empty. No historical
numerical equivalence is claimed. Independent direct-DFT and lower-level
composition comparisons provide the current integration validation.

## Numerical validation

The independent complete tiny model has a 2x4 mesh, common a=20,b=10,
gamma=90-degree cell, Al and H sites with non-origin fractional coordinates,
fractional occupancies, per-species B=.86/1.2, vacancy=.2, and distinct slice
thicknesses/shifts. The legacy mask leaves K=0,H=-1,0,1; independent source
table values at STL=0,.025 therefore suffice for all retained coefficients.
The test constructs the atomic sums, potential, gratings and propagators
independently and uses nested scalar DFTs for every sequence occurrence.
It uses independently checked wavelength/sigma literals, never the production
high-level API to generate expectations. Cases include one occurrence, two
distinct types, A-B-A-B, both grating methods and four signed tilt choices.

Separate manual composition checks use the existing lower-level routines on
square/rectangular meshes, 60-degree cells, enabled absorption, different
per-slice dilation, both reciprocal orders, multiple beams and duplicates.
These check wiring and are distinct from the independent numerical oracle.
Instrumentation verifies exactly one grating per stored ID during preparation
and one propagator per stored ID per tilt, irrespective of occurrence count.

Measured maximum absolute residuals with the required qcbed-python interpreter:

| Comparison | Residual |
| --- | --- |
| Prepared grating vs independent tiny reference | 2.169e-19 |
| Final complex wave vs independent tiny reference | 4.463e-16 |
| Total reciprocal intensity vs independent reference | 6.662e-16 |
| Direct vs FFT prepared gratings | 3.904e-18 |
| Final waves using direct vs FFT preparation | 2.221e-16 |

The worst tiny-model wave residual is the FFT A-B-A-B case at tilt (.2,0),
unshifted row=0,column=3 (H=-1,K=0). Maximum pointwise relative residual is
3.259e-13, excluding reference amplitudes <=1e-15 from division; absolute
tolerance handles those small/zero entries. Direct/FFT maxima occur on the
6x4 orthogonal mesh with centered storage. Acceptance tolerances are 5e-16
for tiny gratings, 3e-15 for tiny final waves, 4e-15 for beam intensities and
2e-15 for final-wave composition/direct-vs-FFT checks. Norm comparison here
uses sum(abs(final reciprocal wave)**2); it compares independent calculations,
not an assumption of unit norm after boundary suppression. The uniform
origin-only case separately verifies unit intensity and absence of extra RMT.

Focused tests also check early validation before structure-factor sums,
immutable storage and metadata, invalid cross-slice compatibility, raw flag
semantics, invalid species/sequence references and nonfinite parameters.
No timings or performance claims are made.

Deferred: nonunit reflection refinement, gas randomization, regional applied
state construction, incompatible lateral-cell interfaces, scan/tilt grids,
multi-tilt batching, multiprocessing, probes, GPU, experimental-intensity
parsing, YNORM fitting, chi-square, parameter refinement and CLI.
