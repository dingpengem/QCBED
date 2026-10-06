# Patch B1–B4: crystallography and explicit experimental-beam indexing

B1 is additive and independent of the existing two-dimensional multislice lattice
implementation. It does not create slices, assign experimental reflections, or
modify canonical run configuration, propagation, refinement or Patch A.

## Architecture and conventions

| Module | Responsibility |
|---|---|
| `qcbed.crystallography.lattice` | Immutable direct and reciprocal 3-D lattices |
| `qcbed.crystallography.structure` | Ordered immutable fractional atomic sites |
| `qcbed.crystallography.ase_adapter` | Copy ASE inputs into canonical models |
| `qcbed.crystallography.visualization` | Optional Matplotlib 3-D diagnostics |

`DirectLattice.vectors` has shape `(3,3)` with **rows** a,b,c, in Å.
`volume` is positive in Å³. Left-handed, singular, nonfinite and incorrectly shaped
cells are rejected, never repaired. `metric_tensor = A @ A.T` has units Å²;
`lengths` returns (a,b,c); `angles` returns (alpha,beta,gamma) in **degrees**,
where alpha is between b,c; beta between a,c; gamma between a,b.

`ReciprocalLattice.vectors` stores rows a*,b*,c* in Å⁻¹. The crystallographic
convention is **no 2π**:

```python
B = np.linalg.inv(A).T
B @ A.T == np.eye(3)  # within floating-point roundoff
cartesian_position = fractional_position @ A
cartesian_g = np.array([h, k, l]) @ B
d_hkl = 1 / np.linalg.norm(cartesian_g)
```

`g_vector(hkl)` requires three integer indices; its zero vector is allowed.
`d_spacing((0,0,0))` is rejected. A geometric spacing is not a statement about
systematic absences or reflection intensity. Extremely ill-conditioned cells can
lose floating-point accuracy; B1 does not perform cell reduction or introduce an
arbitrary condition-number cutoff. Representable positive finite volume and
finite inverse are required.

All model arrays are owned float64 snapshots with immutable byte storage. Returned
metric, length, angle and g-vector arrays are also read-only. Sites preserve atom
order and duplicates. Fractional positions outside [0,1) remain unchanged; no
wrapping, merging, symmetry expansion or coordinate transform is hidden.
`AtomicSite.occupancy` satisfies 0 < occupancy <= 1; optional `debye_waller_b` is
finite and nonnegative in Å² (isotropic B). It is metadata only, with no B/U
conversion or coupling to existing thermal-factor/refinement code.

## ASE frontend and dependencies

Install the optional structure-input extra (the development conda environment also
includes ASE):

```bash
python -m pip install -e '.[crystallography]'
```

ASE requirement: `ase>=3.23,<4`; tested with ASE 3.29.0 and Python 3.11.
ASE itself depends on SciPy and Matplotlib. Importing the crystallography backend
does not import ASE or Matplotlib; frontend imports are deferred until use.
No additional crystallography libraries are introduced.

```python
from ase.build import bulk
from qcbed.crystallography import atomic_structure_from_ase, plot_atomic_structure

atoms = bulk('Al', 'fcc', a=4.05, cubic=True)
structure = atomic_structure_from_ase(atoms)
reciprocal = structure.direct_lattice.reciprocal()
figure, axes = plot_atomic_structure(structure, repeats=(2, 2, 2))
```

For files, callers may use `ase.io.read(path)` explicitly and then the same adapter.
B1 has no format autodetection beyond ASE's own frontend. Conversion copies
`atoms.cell.array`, `atoms.get_scaled_positions(wrap=False)` and chemical symbols.
A complete right-handed cell is mandatory, including for nonperiodic ASE objects.
Fractional round trips through ASE Cartesian storage are subject to floating-point
roundoff. The canonical model retains no `ase.Atoms` object.

The basic adapter defaults occupancies to 1 and B to None. It does not interpret
arbitrary ASE arrays/info, PBC flags, calculators, velocities or disordered-CIF
occupancy conventions. Build explicit canonical `AtomicSite` records when those
scientific attributes matter; B1 does not promise a lossless import of all ASE
metadata or CIF disorder. Chemical symbols use the 118 element symbols; dummy
atoms and isotope aliases are not accepted.

## Diagnostic visualization

`plot_atomic_structure(structure, repeats=(1,1,1), ax=None)` returns `(figure,axes)`.
It draws each unit-cell frame and Cartesian atoms, distinguishes species by
color/marker with a legend, and labels lengths in Å. Replication is translation
by integer direct lattice vectors, preserving original fractional coordinates.
It does not add duplicate boundary atoms, infer bonds, encode occupancy as atom
counts, orient a zone axis, call `show()`, or modify the model. The caller owns
saving/displaying/closing the figure. Scientific classes contain no visual style.

Run the offline deterministic example:

```bash
python -m examples.b1_crystal_structure --output tmp/b1_crystal_structure.png
```

The conventional FCC Al cell has four atoms, a=4.05 Å, volume=66.430125 Å³ and
reciprocal axes of magnitude 0.2469135802469136 Å⁻¹. The 2×2×2 plot contains 32
atoms. Geometric d100=4.05 Å, d110=2.86378246381 Å, d111=2.33826859022 Å and
d200=2.025 Å. These are not fitted to experimental or DigitalMicrograph values.

## Ownership boundary and later stages

Patch A owns stable experimental beam IDs, detector geometry, empirical (n1,n2),
signal, uncertainty and support. **(n1,n2) is not automatically physical (h,k).**
B1 owns only the physical 3-D atomic structure and direct/reciprocal geometry.
No integration with the refinement engine or new GUI controls is introduced.

B2 implements arbitrary integer zone axes and Cartesian beam-frame projection.
B3 now implements exact integer recast cells and original-hkl conversion (below).
Physical reflection enumeration is not implemented. B4 explicit experimental-beam annotation is documented below; automatic indexing
and production configuration migration remain outside these stages.

## Validation

Focused tests cover cubic and orthorhombic analytic results; a triclinic inverse
relation and independent cross-product reciprocal oracle; wrong shape, nonfinite,
singular and left-handed cells; invalid sites and spacings; ASE order/unwrapped
positions and copy isolation; and headless plots with replication. No Fortran
routine is translated or replaced in B1; the existing full suite supplies
regression coverage for the unchanged simulation/refinement backend.


## B2: zone axes and beam-frame projection

New backend modules are `orientation.py` and `projection.py`. `replication.py`
provides the shared direct-cell translation logic used by B1 and B2 viewers.
All plotting remains in `visualization.py`, with deferred Matplotlib imports.
No new dependency is required.

### Axis and orientation convention

`ZoneAxis(uvw)` accepts exactly three integers (not floats, strings or booleans),
rejects zero, divides by their greatest common divisor, then makes the first
nonzero component positive. Thus [-2,-2,0] becomes [1,1,0], and [0,-2,2] becomes
[0,1,-1]. This is an **unoriented crystallographic axis**: opposite incidence
signs are intentionally identified. A future directed beam API must be separate.

`zone_axis_cartesian(lattice, zone_axis)` returns a read-only Cartesian unit vector:

```text
r = [u,v,w] @ A = u*a + v*b + w*c
z = r / |r|
```

It does not normalize the integer triple as though it were Cartesian. Unequal
lengths and nonorthogonal cells are handled by the direct lattice, in Å.

`CrystalOrientation(direct_lattice, zone_axis)` computes immutable `beam_basis`
and `beam_direction_cartesian`; callers cannot inject inconsistent matrices.
Its exact frame rule is:

1. Find the Cartesian reference axis least parallel to z (smallest |z component|).
2. Break ties with x, then y, then z priority. Values within 32 float64 machine
   eps of the minimum count as ties, preventing roundoff-dependent choices at
   symmetric directions. This dimensionless tie band does not alter z.
3. `x_beam = normalize(cross(reference,z))`.
4. `y_beam = normalize(cross(z,x_beam))`.
5. Store `R = [x_beam; y_beam; z_beam]` as **rows**.

The chosen reference is recorded as `reference_axis_index`. The reference is
always well separated from parallel; no tiny cross-product normalization is used.
This is deterministic and piecewise smooth, but not globally continuous across
reference-axis selection boundaries. No global smoothness is promised.

For row vectors (or arrays with final dimension 3):

```text
v_beam = v_cart @ R.T
v_cart = v_beam @ R
R @ R.T ≈ I; det(R) ≈ +1; cross(x_beam,y_beam) ≈ z_beam
```

The methods `cartesian_to_beam()` and `beam_to_cartesian()` preserve units and
return immutable arrays. Positions transform in Å, reciprocal vectors in Å⁻¹.
For cubic [001], the fixed convention is x_beam=(0,-1,0), y_beam=(1,0,0),
z_beam=(0,0,1), so Cartesian (x,y,z) becomes (-y,x,z). No translation is applied;
positions remain relative to the crystallographic origin.

### Atomic projection and replication

`project_atomic_structure(structure, orientation, repeats=(1,1,1))` returns a tuple
of immutable `ProjectedAtomicSite` records with `source_site_index`, `cell_offset`,
`source_site`, `symbol`, `x`, `y`, and `depth`. x/y/depth are in Å. Keeping the
original immutable source site also retains occupancy, B, and unwrapped fractional
coordinates. The orientation and structure must have identical direct lattices.

Replication accepts positive integers only. Cell offsets are lexicographic with
c fastest; within each offset, original site order is retained:

```text
position_cartesian = (fractional_position + cell_offset) @ A
position_beam = position_cartesian @ R.T
```

No sites are wrapped, merged, sorted by depth, or removed. Overlapping projected
columns retain all individual atoms and depths. The replicated direct cell is not
an automatically constructed commensurate zone-axis supercell. Optional column
grouping is deferred; there is no hidden spatial tolerance.

### Supplied reciprocal vectors and zone law

`project_hkl(reciprocal_lattice, orientation, hkl)` returns an immutable
`ProjectedReciprocalVector` with `hkl`, `g_cartesian`, `g_beam`, `g_projected_xy`
and `beam_parallel_component`. It handles one supplied HKL, without enumeration,
selection or claiming allowed-reflection intensity. Reciprocal convention remains
no 2π. Supply the dual of the orientation's direct lattice when testing the zone law.

`zone_law_value(hkl, zone_axis)` evaluates integer `h*u+k*v+l*w` exactly. For dual
lattices, the quantitative consistency relation is:

```text
gz = (h*u + k*v + l*w) / |[u,v,w] @ A|
```

Hence an exact zone reflection has negligible gz at floating-point precision,
including for nonorthogonal cells. Zero and nonzero cases are tested. This helper
does not implement zone filtering or reflection enumeration.

**The orthonormal beam frame is a Cartesian viewing frame, not a 2-D integer
crystallographic reciprocal basis.** Projected a*,b*,c* can be inspected, but none
is selected as an indexing basis because of visual alignment. Patch A (n1,n2)
remains unrelated to physical HK until a later explicit mapping stage.

### Viewer and example

`plot_projected_atomic_structure(structure, orientation, repeats=(1,1,1), ax=None)`
returns a Matplotlib figure/2-D axes, with equal x/y aspect, element colors/markers,
Å labels and a `Zone axis [u v w]` title. It does not call `show()`, group atoms,
change orientation or mutate inputs. Coincident symbols may visually overlap;
underlying projected sites are retained. No periodic 2-D cell boundary is inferred.

```bash
python -m examples.b2_zone_axis_projection --output tmp/b2_zone_axis_projection
```

This creates FCC Al [001], [110], [111] images `zone_001.png`, `zone_110.png`,
`zone_111.png`, plus `geometry_diagnostics.json`. The output directory must be new.
It prints each beam direction/frame and supplied in-zone reciprocal vectors.
The geometry examples do not assert allowed diffraction intensities.

Measured example checks over cubic, orthorhombic and triclinic cells:

- Maximum orthonormality error: 2.220446049250313e-16.
- Beam-frame determinant range: 0.9999999999999998 to 1.0000000000000004.
- Maximum Cartesian round-trip error: 8.881784197001252e-16 Å for the test vectors.
- Maximum in-zone reciprocal parallel component: 4.7007413780846676e-17 Å⁻¹.

B1+B2 focused tests cover 105 cases, including [100], [1 -1 0], [2 1 1], invalid
axis types, sign/gcd canonicalization, symmetric tie handling, replica provenance,
immutability and headless plotting. B2 changes no Patch A, canonical runner,
multislice, objective or refinement behavior. B3 integer machinery is documented below; B4 is not started.


## B3: original Miller identity and internal recast cell

**Users input true original-cell (h,k,l). Users do not input recast (H,K).**
B3 is crystallographic algebra only. It does not modify configuration schemas,
PreparedExperimentalDataset, experimental beam records, refinement beam tables,
or the existing numerical engine. No real objective is run.

Four concepts must remain distinct:

| Concept | Meaning |
|---|---|
| Original Miller (h,k,l) | Canonical user-facing physical reflection identity |
| User CBED basis hkl_1,hkl_2 | Two ordered original-cell reciprocal reflections |
| Experimental (n1,n2) | Patch A's detector-lattice coordinates |
| Internal (H,K,0) | Recast reciprocal coordinates, shortened to internal (H,K) |

Neither experimental (n1,n2) nor the first two original Miller components are
implicitly internal (H,K). Changing hkl to HKL is a basis change, not a change in
the physical reflection. No new user workflow requires recast indices.

### Deterministic unimodular completion

`build_zone_axis_cell(lattice, zone_axis)` returns immutable `ZoneAxisCell`, keeping
`original_lattice`, canonical `zone_axis`, integer `transform`, integer
`inverse_transform`, recast `direct_lattice` and `reciprocal_lattice`.
The constructor derives all transformed fields; callers cannot inject inconsistent
matrices. The third transform row is exactly the B2-canonical primitive [u,v,w].

The implementation uses signed extended Euclid, via ordinary Euclidean division
on absolute values, followed by sign restoration of Bézout coefficients. No
floating-point search, symbolic dependency or Cartesian perpendicular basis is
used. Specifically:

```text
g = gcd(u,v)
s*u + t*v = g
alpha*g + beta*w = 1

p = (t,-s,0)
q = (beta*u/g, beta*v/g, -alpha)
z = (u,v,w)
T = [p; q; z]
```

All divisions in this construction are exact integer divisions. If g=0, canonical
primitive z is [001], and T is the identity. Otherwise:

```text
cross(p,q) = (s*alpha,t*alpha,beta)
det(T) = cross(p,q) dot z = alpha*g + beta*w = +1
```

Thus every primitive axis admits this completion; no supercell determinant is
accepted. The exact inverse is the adjugate with columns `q×z`, `z×p`, `p×q`.
Algebra uses Python integers. Published transform arrays use immutable int64 byte
storage; out-of-range coefficients raise an explicit OverflowError before any
truncation. Integer reflection transforms remain Python integer arithmetic and
support indices beyond int64. Physical arrays remain float64, with its usual
conditioning and representability limits. Large axes can yield poorly conditioned
cells; B3 does not silently reduce or optimize this deterministic basis.

### Row-vector transformations and derivation

Let original direct rows be a,b,c, recast rows A,B,C. Then:

```text
A_zone = T @ A_original
C = u*a + v*b + w*c
B_original = inv(A_original).T
B_zone = inv(A_zone).T = inv(T).T @ B_original
```

All direct units are Å and reciprocal units Å⁻¹, with **no 2π**. det(T)=+1
preserves handedness and cell volume. Direct A and B need not be perpendicular to
C; reciprocal A* and B* are perpendicular to C by duality. The recast cell is an
integer crystal basis, distinct from B2's orthonormal Cartesian viewing frame.

Equating physical reciprocal vectors gives
`hkl @ B_original = HKL @ inv(T).T @ B_original`, hence:

```text
HKL_zone = hkl_original @ T.T
hkl_original = HKL_zone @ inv(T).T
```

`original_hkl_to_zone_HKL()` and `zone_HKL_to_original_hkl()` use exact integer
products with the exact integer inverse. In particular the third component is
**L=h*u+k*v+l*w**. `is_zone_reflection_original()` and `is_zolz_zone_index()`
therefore agree exactly, independent of lattice metric. Zero reflection counts as
in-zone algebraically; this does not define a finite d-spacing or an intensity.

`ZoneReflectionHK(H,K)` explicitly means internal (H,K,0), with `as_HKL()`.
`original_hkl_to_internal_HK()` rejects out-of-zone reflections; the inverse
`internal_HK_to_original_hkl()` restores the user-facing original identity.
Raw two-tuples are not accepted in place of the explicit internal type.

For fractional row positions, equating Cartesian positions gives:

```text
f_original @ A_original = f_zone @ T @ A_original
f_zone = f_original @ inv(T)
f_original = f_zone @ T
```

`original_fractional_to_zone()` and `zone_fractional_to_original()` apply these
maps to (...,3) coordinates without wrapping. `recast_structure_to_zone_cell()`
keeps site count, order, duplicates, species, occupancy and B metadata, transforming
only fractional coordinates and the basis. It checks original-lattice provenance.
No atoms are replicated or merged. Float64 Cartesian coordinates are invariant
within roundoff; integer Miller round trips are exact.

### User original-hkl basis primitives

`OriginalReflectionBasis(hkl_1,hkl_2,zone_axis)` requires two nonzero, integer,
exactly independent original reflections satisfying the selected zone law.
Independence is checked by their exact integer cross product, rejecting all
collinear pairs, including rational multiples. Nonprimitive reflections such as
002 are deliberately preserved; a supplied basis may span a sublattice and need
not span every zone reflection. User basis ordering/sign is not canonicalized.

`reflection_from_basis_coordinates(basis,n1,n2)` returns exactly:

```text
hkl = n1*hkl_1 + n2*hkl_2
```

`original_reflection_basis_vectors(original_reciprocal,basis)` returns two
Cartesian reciprocal rows `g1=hkl_1@B_original`, `g2=hkl_2@B_original`. Lengths
and angle are preserved, not orthonormalized. B2 can project these with
`orientation.cartesian_to_beam(vectors)` when needed.

This is the intended future B4 interpretation of Patch A's (n1,n2), but no Patch A
object is imported, connected, relabelled or changed in B3. Existing Patch A
basis fields retain their meaning; a future configuration migration must be explicit.

For [110], the implemented completion is:

```text
T = [[1,0,0], [0,0,-1], [1,1,0]]
(H,K,L) = (h,-l,h+k)
```

For user basis `(1,-1,0)` and `(0,0,2)`:

| Integer basis coordinates | Original hkl | Internal HK0 |
|---|---|---|
| (1,0) | (1,-1,0) | (1,0,0) |
| (0,1) | (0,0,2) | (0,-2,0) |
| (1,1) | (1,-1,2) | (1,-2,0) |
| (-1,1) | (-1,1,2) | (-1,-2,0) |

These internal labels depend on the documented deterministic completion; physical
g vectors do not change. They are not user-specified or inferred detector indices.

### Example and validation

```bash
python -m examples.b3_zone_axis_cell --output tmp/b3_zone_axis_cell
```

The directory must be new. The example prints original/recast lattices, T,
det(T), original→zone→original indices for [001]/[110]/[111], internal HK for
in-zone reflections, and the [110] basis algebra above. It writes
`zone_cell_diagnostics.json`. No remote input or ASE dependency is needed.

B1–B3 focused coverage: 169 tests. Across seven axes and three metrics, every
reflection in [-3,3]^3 is checked (7,203 exact round-trip and zone-law/L=0 checks).
An additional signed-axis grid checks completion, including gcd reductions and
zero components. Measured example maxima:

| Check | Maximum absolute error |
|---|---:|
| Reciprocal-vector invariance | 4.996003610813204e-16 Å⁻¹ |
| Cartesian atom-position invariance | 1.7763568394002505e-15 Å |
| B2 z versus normalized recast C | 2.220446049250313e-16 |
| In-zone reciprocal beam-parallel component | 1.022317817829952e-16 Å⁻¹ |
| Recast reciprocal duality | 9.992007221626409e-16 |
| Reciprocal transformation relation | 2.220446049250313e-16 Å⁻¹ |

B3 does not enumerate physical reflections, apply extinction rules, associate
experimental beams, construct multislice input, run objectives or modify refinement.
The deterministic cell is not a shortest/reduced or Cartesian-orthogonal cell.
B4 is not started.

## B4: explicit experimental-axis to physical-reflection bridge

B4 implements deterministic annotation and refinement-input adaptation. It does
not perform indexing search, fit detector calibration, run an objective, or change
Patch A detection/extraction. The exact correspondence supplied by the caller is:

```text
experimental pixel g1 axis ↔ user original-cell hkl_1
experimental pixel g2 axis ↔ user original-cell hkl_2
(n1,n2) -> n1*hkl_1+n2*hkl_2 -> original hkl -> internal (H,K,0)
```

No basis swap, sign flip, 90-degree rotation or symmetry-equivalent replacement
is inferred. Reversing or negating a user basis changes the output correspondingly.
Original hkl remains the physical/user-facing identity; internal HK is derived.
Pixel vectors such as (81,9),(-10,81) are in detector pixels, not Å⁻¹. Semantic
axis correspondence supplies no reciprocal-space detector scale calibration.

### Public interfaces and immutable metadata

`beam_indexing.py` supplies:

- `IndexedBeamReflection(stable_id, experimental_index, original_hkl, internal_hk)`.
- `index_experimental_beam(...)`, with keyword-only identity, basis and zone cell.
- `index_experimental_beams(beams, reflection_basis, zone_cell)`, accepting an
  ordered iterable of `(stable_id,(n1,n2))` pairs.

These reuse B3's `reflection_from_basis_coordinates()` and
`original_hkl_to_internal_HK()` without duplicating the crystallographic algebra.
Basis and zone-cell canonical axes must agree. Stable IDs must be nonempty/unique
in batches; duplicate experimental coordinates under different IDs are preserved.
Order is never sorted. Empty batches are supported. The direct beam maps through
exact zero-index algebra to original 000 and internal 00, without d-spacing.
B3 validation rejects noninteger, zero, collinear or out-of-zone basis reflections.

`prepared_indexing.py` is the narrow Patch A boundary:

- `index_beam_catalogue()` reads only site `beam_id` and `integer_indices`, in
  catalogue order. It does not filter statuses or promote unconfirmed geometry.
- `index_prepared_dataset()` returns immutable `IndexedPreparedDataset`, containing
  the unchanged A1 dataset, authoritative reflection basis, original-lattice/zone
  provenance, and one derived `indexed_beams` tuple.
- `internal_scan_basis()` derives internal scan coefficients using an explicit
  A10 basis transform.
- `adapt_indexed_prepared_dataset()` passes internal HK assignments and derived
  scan basis to the existing `adapt_prepared_dataset()` API.

**PreparedExperimentalDataset is not extended or relabelled.** Signal, variance,
weights, ROIs, center provenance, scan corners and ordering remain untouched.
The annotation wrapper references the immutable source, rather than copying large
arrays or storing competing assignments inside unrelated records. It requires
explicit `lattice_coordinates`; it never treats `BeamObservation.indices` as
experimental (n1,n2). Indexing works when scientific variance is absent. The
existing variance readiness error still blocks adaptation in that case.

The pure indexing backend does not import Patch A. Imports of A1/engine adapter
code are deferred until the integration functions are used. No plotting or ASE
dependency is introduced to indexing.

### Scan-coordinate and simulation-frame boundary

Beam identities and A10 scan coefficients are distinct. A10 can have used an
explicit integer matrix M with `scan columns = experimental basis columns @ M`.
B4 therefore requires `scan_basis_transform` explicitly; it does not infer M from
corners, pixel lengths or prose provenance. Let Q contain the internal HK of
hkl_1/hkl_2 as columns. The engine scan basis columns are `Q @ M`.
Only the explicit scan-coordinate conversion applies M; beam indexing always
uses the authoritative original hkl pair in its supplied order.

For the existing A10 convention M=diag(1,-1), the [110] demonstration derives scan
basis rows `(1,0,0)` and `(0,2,0)`. A10 corners are unchanged. M must be a 2×2
integer unimodular matrix; internal scan coefficients must be exactly representable
when converted to the existing float64 scan model. Omitting M is an error, not an
implicit identity choice.

```python
from qcbed.crystallography import (
    ZoneAxis, OriginalReflectionBasis, build_zone_axis_cell,
    index_prepared_dataset, adapt_indexed_prepared_dataset,
)

zone = ZoneAxis((1,1,0))
basis = OriginalReflectionBasis((1,-1,0), (0,0,2), zone)
zone_cell = build_zone_axis_cell(original_structure.direct_lattice, zone)
indexed = index_prepared_dataset(dataset, basis, zone_cell)
# Requires validated variance, and the explicit transform actually used by A10:
inputs = adapt_indexed_prepared_dataset(
    indexed, scan_basis_transform=conversion.canonical_to_scan_basis,
    objective_mode='masked',
)
```

The historical engine interprets reflection indices as Fourier orders in its
stored two-dimensional slice reciprocal cell; `generate_scan_tilts()` returns HK
in that same basis. Therefore its simulation sample MUST already represent the
same zone-recast reciprocal frame. B4 does not build or verify a multislice sample
from the B1 3-D structure. In particular oblique recast A/B vectors cannot simply
be copied as perpendicular 2-D slice lengths. A physically consistent slice/metric
adapter is separate work. A successful experimental adapter call is not proof that
an arbitrary existing sample is compatible, and is not permission to run a real
objective. No structure factors, coupling, propagation, objective, YNORM, packing
or optimizer code is changed.

### Configuration audit: B4a delivered, B4b deferred explicitly

No existing QcbedConfig/YAML fields are reinterpreted or added in B4a:

- `ScanConfig.basis_1/basis_2` compile directly to legacy HK coefficient vectors.
- `MicroscopeConfig.zone_axis` is retained existing microscope/legacy metadata.
- `SampleConfig` owns stored 2-D slices and x/y atoms, not B1's original 3-D cell.
- The compiler has no B1-structure-to-recast-slice adapter.

Adding only a `crystallography` YAML block would leave scientifically incomplete
runtime behavior and also affect canonical serialization/checkpoint identities.
B4b configuration integration is therefore deferred until original-structure
ownership and the recast multislice boundary are specified. Current original-hkl
inputs are the typed B3 `OriginalReflectionBasis`, validated on construction;
there is no workflow requiring manually calculated internal HK. The future YAML
concept `crystallography.cbed_basis.g1_hkl/g2_hkl` is NOT accepted by today's strict
loader. Existing Python/YAML/legacy workflows remain unchanged.

### Nine-beam example and numerical verification

```bash
python -m examples.b4_beam_indexing --output tmp/b4_beam_indexing
```

The directory must be new. It prints and writes `indexed_beams.csv` and
`indexed_beams.json`, including all three index identities, original g, basis-sum
g and difference norm. This is a **demonstration** of a supplied [110] basis
`(1,-1,0),(0,0,2)`, not an experimentally established assignment for pattern.dat.

| Stable ID | Experimental (n1,n2) | Original hkl | Internal HK |
|---|---|---|---|
| lattice_0_0 | (0,0) | (0,0,0) | (0,0) |
| lattice_-1_0 | (-1,0) | (-1,1,0) | (-1,0) |
| lattice_1_0 | (1,0) | (1,-1,0) | (1,0) |
| lattice_0_-1 | (0,-1) | (0,0,-2) | (0,2) |
| lattice_0_1 | (0,1) | (0,0,2) | (0,-2) |
| lattice_-1_-1 | (-1,-1) | (-1,1,-2) | (-1,2) |
| lattice_1_1 | (1,1) | (1,-1,2) | (1,-2) |
| lattice_-1_1 | (-1,1) | (-1,1,2) | (-1,-2) |
| lattice_1_-1 | (1,-1) | (1,-1,-2) | (1,2) |

Tests cover negative/central indices, duplicates, ordering, invalid bases, mismatched
zones, explicit sign/swap changes, unchanged A9 arrays and A10 corners, gradient
metadata, missing variance and native adapter HK delivery. Cubic, orthorhombic,
hexagonal (a=b, gamma=120°, c≠a) and triclinic metrics are included. Hexagonal input
uses original-cell **three-index** notation; Miller–Bravais hkil/uvtw is a future
explicit adapter, never silently accepted as three indices.

For 49 integer-coordinate pairs in each of the four example metrics, measured
maximum `|g(original_hkl) - (n1*g1+n2*g2)|` is 2.220446049250313e-16 Å⁻¹.
The focused suite (crystallography, prepared bridge and existing prepared integration)
contains 207 passing tests. B4 performs no objective evaluation or optimization.

### Remaining real-data steps

This task explicitly confirms `tests/fixtures/legacy_reference/sigma.dat` as a legacy per-pixel
sigma map, so a future explicit integration uses variance=sigma². This supersedes
the earlier Patch A closure's unresolved sigma-versus-variance status; historical
reports remain unmodified. B4 does not load/square the file, fabricate A4 masks,
or imply that a noise-characterization provenance chain has been established.

Before a real objective: explicitly supply/confirm the actual crystal, zone and
ordered original-hkl basis for the experiment; integrate the confirmed sigma map
through the existing uncertainty/extraction path; retain the actual A10 transform;
and establish a matching recast simulation sample/metric. Configuration integration
is B4b follow-up. Automatic indexing/calibration, four-index notation, full real
objective and optimization remain outside this patch.

### B5 strict simulation-cell geometry

B5 uses strict transverse-periodicity checks before explicit depth slicing. B3
`ZoneAxisCell` is an indexing/recast basis, distinct from
`MultisliceSimulationCell`. The restricted constructor supports exact integer C
shears; unsupported geometries are rejected, without full-period projection.
`build_depth_resolved_sample()` retains original site provenance and requires
explicit depth boundaries. `build_zone_axis_simulation_sample()` prepares the
supported subset through the existing propagation backend with explicit scattering
settings. Production configuration/refinement integration remains deferred. See [the B5 audit](b5_simulation_sample_audit.md)
for conservative compatibility rules and remaining limitations.

## B4b production configuration integration

`examples/crystallography_config.yml` is a deterministic synthetic demonstration,
not an assignment of a zone axis or hkl basis to `reference_data/pattern.dat`.
Users supply an original structure, original direct-lattice zone direction, and
ordered original-cell `g1_hkl/g2_hkl`. Users do not enter recast HK.

The two input modes are mutually exclusive:

- Existing `sample` species/regions/stored slices, without `crystallography`.
- `sample: null` plus `crystallography`. This mode uses `CrystalScanConfig`
  (dimensions/corners only) and `CrystalExperimentalDataConfig` (file metadata
  only). Legacy `basis_1/basis_2` and `beam_indices_hk` inputs are rejected here.

The new section contains:

```yaml
crystallography:
  structure: {path: b4b_demo.extxyz, format: extxyz}
  zone_axis: [0, 0, 1]             # original direct-lattice direction
  cbed_basis:
    g1_hkl: [1, 0, 0]             # original-cell Miller indices
    g2_hkl: [0, 1, 0]
  mesh: {mesh_x: 8, mesh_y: 8}
  depth_angstrom: 5.0
  slice_boundaries_fractional: [0.0, 0.5, 1.0]
  debye_waller:
    - {species: Si, value_angstrom_squared: 0.3}
  dsap_inverse_angstrom: 0.5
  beams:
    - {id: direct, lattice_coordinates: [0, 0]}
    - {id: g1, lattice_coordinates: [1, 0]}
  scan_basis_transform: [[1, 0], [0, -1]]
  absorption: {enabled: false}
```

`microscope.zone_axis` must describe the same original direction. CBED basis
signs/order and beam list order are never corrected. Duplicate beam coordinates
are retained, but IDs must be unique. Only 3-index input is supported; 4-index
Miller–Bravais input is rejected. The explicit integer unimodular scan transform
has exactly the B4a/A10 meaning; corners are retained without transformation.

### File boundary and pure compiler

```python
from qcbed import load_yaml_config, compile_config_from_sources
config = load_yaml_config("examples/crystallography_config.yml")
compiled = compile_config_from_sources(config)  # reads local structure; no objective
print(compiled.crystallography.diagnostics)
```

The YAML loader validates the file-backed configuration and resolves paths relative
to the document. `resolve_config_paths` itself only resolves paths for this mode.
`compile_config_from_sources` reads the resolved absolute structure path, then
calls `compile_config(config, structure=canonical_structure)`. The latter remains
pure: no filesystem access or ASE dependency. Direct Python frontends can supply
an existing `AtomicStructureModel` through this argument. No `ase.Atoms` survives
loading. `load_config_structure(config)` exposes the input boundary separately.
Default structured validation also remains filesystem-independent: supply
`structure=` for compile diagnostics, or explicitly request the FILES phase to
load sources. Crystallography errors carry section/field locations.

ASE is optional (`pip install 'qcbed[crystallography]'`). File loading uses one
local frame, optionally an explicit ASE format. CIF, POSCAR/CONTCAR and extended
XYZ are accepted when ASE can read a valid full right-handed lattice. Plain XYZ
without a lattice, multiple frames, unreadable files and partial/disordered CIF
occupancies are rejected. The existing basic ASE adapter imports full occupancy
and does not infer B from arbitrary metadata. Explicit B per chemical symbol is
therefore mandatory. Direct canonical input retains occupancy, and conflicting
explicit site B is rejected. No online lookup occurs.

### Compiled integration

Compilation reuses B1–B5 to derive orientation, zone cell, strict simulation cell,
reflection basis, indexed beams, depth slices and a prepared sample. It then feeds
an equivalent pre-sliced configuration into the existing compiler. The immutable
`CompiledConfig.crystallography` retains original structure, source path,
orientation, basis, zone cell, B5 sample, scan transform and scan configuration.
`CompiledConfig.configuration` is the resolved pre-sliced execution configuration;
`source_configuration` retains the user's original declarations. Legacy-mode
compilation retains its previous configuration identity.

Generated IDs are deterministic: region `crystal`, species by chemical symbol,
slices `crystal-slice-0`, etc., atoms `site-<original index>`. Bindings use the
existing compiler/layout, not original 3-D coordinate refinement. The B5 supported
subset remains zero vacancy, unit applied dilation, unity reflection multipliers,
one region, uniform B per species and explicit slices. Linked absorption C is
currently rejected in this bridge; existing pre-sliced linked-C behavior is
unchanged. Strict incompatible geometry fails without a projection fallback.

`compiled.crystallography.index_dataset(dataset)` provides the B4a immutable
annotation wrapper. `adapt_dataset(dataset, objective_mode=...)` calls the existing
refinement adapter with derived HK, checks configured beam identities/order and
scan dimensions/corners, and preserves the scientific dataset. Variance readiness
is still enforced by that adapter. Compilation does not consume experimental
arrays or evaluate YNORM/objective.

The canonical run snapshot loader uses the same source-compilation frontend;
there is no separate crystallographic runner or changed checkpoint algorithm.
This patch tests snapshot preparation, not objective execution. YAML Save/Save As
preserves all declarations and writes resolved paths; it cannot overwrite the
structure source. Structure files remain external and should be archived with
YAML for reproducibility.

GUI editing is deferred to B4c: opening this mode in the current pre-sliced editor
fails clearly without replacing the current document. Legacy imports retain the
old mode and do not invent original hkl. Crystallography-derived FOR007 export is
explicitly unsupported, since original structure/hkl provenance is not represented.
Existing legacy import/export remains unchanged.

Real `pattern.dat` zone/basis confirmation and separate variance integration still
precede any real-objective work. No real variance integration, objective or
optimization was performed in B4b.

B4b verification: 91 focused configuration/indexing tests and 2678 full-suite
tests passed. Synthetic nine-beam mapping preserved IDs/order and signs; direct
B5 versus compiled preparation matched parameters and complex phase gratings
exactly. YAML Save/Save As/reload preserved canonical declarations. These checks
include snapshot/data loading, without new objective or optimizer execution.

## B6a preflight gate

The production compiled bundle can now be checked against prepared experimental
input with `qcbed.input_preflight.preflight_inputs`. It requires exact attached
ZoneAxisCell provenance, validates integer indexing and floating reciprocal/axis
invariants, and verifies the unchanged B4a scan/array adapter before propagation.
Synthetic PASS is distinct from real-input confirmation. The real assignment is
now confirmed as conventional FCC Al [001], (200)/(020), and real B6a passes;
see [the B6a contract and diagnostics](b6a_input_preflight.md).
No crystallographic transformations or sample-builder physics changed.

## Canonical cell/indexing policy and FCC Al closure

Use the standard user-facing crystallographic cell (conventional where practical).
User-facing zone axes and original hkl refer to that cell. This does not require
all materials to use a nonprimitive cell. Primitive, recast and legacy cells are
internal representations, never substitutes for an explicitly supplied indexing
convention. B3 `ZoneAxisCell` is an indexing/recast basis; B5's distinct
`MultisliceSimulationCell` supplies strict transverse-periodic depth slicing.

The confirmed real case uses a conventional cubic four-site FCC Al cell,
a0=4.049 Å, [001], original g1=(200), g2=(020). Experimental Patch A (n1,n2)
therefore maps to (2*n1,2*n2,0), not to the legacy (110)/(1-10) directions.
B3 derives internal refinement HK; for this case T is identity and HK=(h,k).
The experimental pixel basis remains (81,9)/(-10,81); the explicit A10 transform
is diag(1,-1). No pixel geometry or scan corners changed.

Compatibility-only APIs in `qcbed.crystallography.legacy_for007_basis`:

- `legacy_for007_001_basis_to_conventional_hkl(p,q)` returns (p+q,p-q,0).
- `conventional_hkl_to_legacy_for007_001_basis(hkl)` requires l=0 and same-parity
  h,k, returning ((h+k)//2,(h-k)//2). Both reject noninteger values and booleans.
- `fcc_conventional_reflection_diagnostic(hkl)` returns direct/allowed/forbidden.
  Allowed means all even or all odd. This reporting helper never filters beams,
  gratings or dynamical channels.

Legacy FOR007 internal reciprocal coordinates (1,1)/(1,-1) correspond to canonical
(200)/(020); (1,0)/(0,1) instead correspond to (110)/(1-10). The legacy lateral
length a0/sqrt(2) is not the conventional cubic a0. See the provenance and measured
preflight in [B6a](b6a_input_preflight.md). No objective or propagation is authorized
by a successful input gate.

## B7a physical site groups

Opt-in `crystallography.site_groups` owns original physical sites and their B,
occupancy and independent absorption metadata. The existing target catalogue and
compiler resolve group B aliases into the existing dynamic layout; no simulation
atom index is user input. Recast/slab copies retain original index and group ID.
See [B7 bring-up](b7_refinement_bringup.md) for the implemented schema, one-group-per-
element limitation, explicit bounds, exact baseline tests and version2 identity.
Existing non-group semantics and checkpoint identities remain unchanged.
Group optimizer/resume execution and GUI editing are deferred.

## Explicit one-slice simulation periods

The later [one-slice simulation audit](one_slice_simulation.md) supersedes the
original adapter's two-slab minimum. A proven `MultisliceSimulationCell` accepts
`[0,1]` as one stored slice spanning |C| Å, with the same atom ownership and
cyclic depth sequence as finer partitions. It does not bypass the integer-C-shear
transverse-periodicity proof. The audit also documents that Mesh X/Y sample the
recast P/Q lattice (columns/rows), not the original [001] cell in general.
