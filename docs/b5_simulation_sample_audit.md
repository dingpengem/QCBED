# B5 simulation-sample audit — strict policy

This original audit precedes implementation. The minimum-slice policy was later
revised by the [one-slice audit](one_slice_simulation.md); strict cell geometry is
unchanged. The original audit does not change propagation, slicing,
configuration, refinement, experimental data or legacy/reference files.

## Existing engine contract

- `LegacyInput` contains global species, `SliceGeometry` records and ordered
  `AtomSite` groups per stored slice. An atom has species index, fractional x/y,
  starting occupancy and record type; it has no z coordinate.
- `SamplePreparationParameters` supplies applied voltage, per-slice species B,
  vacancy, dilation, absorption coefficients, requested HK beams and an explicit
  stored-slice sequence. `prepare_multislice_sample()` creates immutable
  `PreparedMultisliceSample` gratings. It is not a 3-D slicing algorithm.
- All stored types share one lateral a,b,gamma and an even mesh. `DirectLattice2D`
  supports oblique gamma; the FFT grid is a fractional lattice-coordinate grid,
  not a requirement that physical a and b be rectangular.
- The 2-D Cartesian convention puts a on +x and b at positive y. The reciprocal
  vectors are its dual, with no 2π. This is not necessarily B2's viewing frame.
- Structure factors use +i*6.283185307*(H*x+K*y), scattering at STL=|g|/2, B in Å²
  via exp(-B*STL²), and starting occupancy*(1-vacancy). B is applied per species
  per stored slice, not independently per atom. Atom summation retains source order.
- `ReciprocalVectorGrid` uses H*a*+K*b*. The same 2-D metric enters propagator SG,
  with tilt expressed in H/K coefficients. No fixed Miller pair is hard-coded,
  but consistency of the supplied slice metric, atomic fractional coordinates,
  and beam HK is the caller's responsibility.
- Depth is in Å. Slice `thickness_angstrom` gives propagation distance. Region
  sequence construction cycles ordered stored slices while accumulated thickness
  is strictly less than requested depth, including the final whole slice that can
  overshoot. It does not truncate a slice or infer atom z.
- Existing slice shifts enter PROPG as a transverse phase. They do not implement
  depth-dependent transverse boundary conditions for an oblique 3-D cell.
- Existing atom-coordinate refinement acts on stored-slice fractional x/y indexed
  by slice and atom. It is not original 3-D fractional-coordinate refinement.
  Site-to-slice ownership must be explicit in any new adapter; existing bindings
  must not be reinterpreted.
- Low-level applied dilation scales reciprocal geometry and propagation distance.
  The supported refinement trial path explicitly uses APPLIEDDILATION=1 for its
  zero-vacancy/no-free-surface subset. B5 cannot reinterpret it as a general
  original-cell lattice refinement without changing existing behavior.

Primary implementation references: `simulation.py`, `structure_factors.py`,
`reciprocal.py`, `propagator.py`, `refinement_trial.py`, `refinement.py`.
The existing source audits `structure_factors.md`, `reciprocal_geometry.md`,
`multislice_sequence.md` and `wave_propagation.md` remain authoritative for legacy
arithmetic and boundary conventions. No kernel change is proposed here.

## Geometry that can be derived without a slicing decision

Use B3 to obtain recast A,B,C and unwrapped atom coordinates; use B2 to obtain
z=C/|C|. Define transverse projections:

```text
P = A - (A dot z)*z
Q = B - (B dot z)*z
```

P,Q span the beam plane, with a positive oriented area because det(A,B,C)>0.
Their 2-D dual embedded in the beam plane is precisely B3's A*,B*. Thus a lateral
metric |P|, |Q|, angle(P,Q) preserves the intended internal HK reciprocal vectors.
A proper simulation Cartesian frame can put P on +x, z along C and y=z cross x;
its rotation relative to B2 must be recorded, not absorbed into HK sign changes.
Atoms retain x/y fractional coefficients relative to P,Q. Their physical depth is
`f_A*(A dot z)+f_B*(B dot z)+f_C*|C|`, not merely f_C*|C|.

These statements establish coordinate and metric consistency. They do not prove
that a depth-resolved slice is periodic under P,Q.

## Blocking boundary-condition distinction

In general A and B include longitudinal components. Translating by A shifts both
transverse position by P and depth by A dot z. Therefore a slab at fixed physical
z is not necessarily invariant under transverse translation P at unchanged z.
The same applies to B/Q. Repeating each independently binned slab with the engine's
ordinary 2-D FFT periodicity can create a different crystal.

A pure transverse lattice translation with the same projected period P requires
A-m*C with integer m and (A dot C)/(C dot C)=m (or an independently established
smaller crystal translation from the atomic motif). A larger commensurate cell
may exist for rational longitudinal offsets, but changes transverse reciprocal
orders and requires an explicit beam-index embedding. For a fully general metric
there is no guarantee of a finite exact beam-perpendicular crystallographic
supercell. A float tolerance or angle approximation is not such a construction.

Integrating a complete C period yields a periodic column-projected potential on
P,Q, but using that as one multislice grating discards intra-period depth ordering.
It is an explicitly chosen projection approximation, not a generic equivalent
replacement for depth-resolved multislice. B1/B3 Cartesian atom-position invariance
alone cannot establish propagation equivalence.

## Selected policy: strict depth-resolved geometry

The canonical B5 policy rejects unsupported transverse periodicity. Full-C-period
projection is not a fallback for incompatible geometry. The later one-slice audit
permits an explicit `(0,1)` partition when the strict simulation-cell proof passes.

`ZoneAxisCell` remains the B3 indexing/recast crystallographic basis.
`MultisliceSimulationCell` is a distinct type: the currently supported constructor
uses P=A-mC, Q=B-nC for exactly integer longitudinal offsets m,n. Its transform
is unimodular and preserves the internal HK reciprocal basis. A future constructor
may choose different integer lattice translations, requiring explicit reciprocal
index embedding and possibly atom replication. This implementation does not search
motif translations or larger commensurate supercells.

Compatibility is proved with rational arithmetic on the supplied floating-point
original lattice and exact integer B3 transform. No near-integer snapping or
metric tolerance is used. Consequently nominally orthogonal cells containing
roundoff longitudinal components can be rejected conservatively. Such rejection
must not be repaired by silently rounding the scientific input.

`build_depth_resolved_sample()` requires explicit strictly increasing fractional
C boundaries from 0 to 1 with at least two boundary values (one or more slices;
the original two-slab minimum is superseded). It partitions point sites by
physical depth in half-open intervals; it does not approximate a continuous atomic
potential or select an optimal slice thickness. Integer C translations place sites
in one depth period and are recorded by original site index. Transverse coordinates
are unwrapped. Original structure, site indices, occupancy and B metadata survive.
Missing B remains missing; no scattering/absorption defaults are invented.

`build_zone_axis_simulation_sample()` bridges these slabs to the unchanged
`LegacyInput` / `SamplePreparationParameters` / `prepare_multislice_sample()`
path and retains both provenance and the prepared snapshot. Species follow first
appearance, atom order is preserved within each slab, duplicate beams are retained,
and each supplied B4 beam HK is checked against its original hkl and the attached
B3 cell. The simulation frame rows are x=P/|P|, z=C/|C|, y=z cross x;
Cartesian row vectors multiply frame.T. This frame may differ from B2's viewing
frame; neither frame changes HK identities.

The adapter requires explicit mesh, energy, DSAP, absorption coefficients (or
None), depth, grating method and reciprocal ordering. It supports explicit
identical B across all sites of a chemical species; missing or differing site B
is rejected instead of silently regrouping species. Vacancy=0, dilation=1 and
unity reflection multipliers are the declared subset. Slice shifts are zero.
Depth in Å cycles the explicit slabs with the existing whole-slice overshoot
convention; requested and effective depth are reported separately.

B4b production configuration integration remains deferred: this is an opt-in
backend bridge, with no remapping of existing refinement bindings. Original site
indices identify slab ownership, but legacy atom x/y refinement still acts on
stored-slice fractional coordinates. No arbitrary original-structure refinement
is introduced. Full-period projection is not available through this adapter.

Focused tests establish exact complex-wave, amplitude and intensity agreement
with an independently hand-constructed compatible legacy sample. They do not
prove a particular slice discretization converged to a continuous 3-D potential.
The legacy/reference-input path remains unchanged. No real-data objective or
optimization is executed. Run diagnostics with
`python examples/b5_strict_geometry.py`.

### Measured coordinate checks

For the cubic [001], orthorhombic [001], hexagonal [001] and integer-C-shear
examples, original versus simulation-cell reciprocal vectors and B2 versus
simulation beam direction had maximum measured differences of zero. The cubic
site reconstruction also had zero measured Cartesian error. These are discrete
fixtures, not universal error bounds.

The engine's unchanged `DirectLattice2D` uses the source pi literal
`3.14159265358` to convert gamma degrees. Reconstructing its reciprocal basis in
original Cartesian coordinates therefore differs slightly from the crystallographic
matrix inverse: maximum measured basis-vector error across these examples was
`2.9018731861896185e-12 Å^-1`. B5 does not compensate gamma or modify this legacy
arithmetic. Exact wave agreement is tested against the same legacy metric
convention, not claimed against a different full-precision trigonometric kernel.
